#ifndef CNV_H
#define CNV_H

#include <boost/filesystem.hpp>
#include <boost/multi_array.hpp>
#include <boost/progress.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/sam.h>


namespace torali
{

  struct CNV {
    int32_t chr;
    int32_t start;
    int32_t end;
    int32_t ciposlow;
    int32_t ciposhigh;
    int32_t ciendlow;
    int32_t ciendhigh;
    uint32_t srleft;
    uint32_t srright;
    uint32_t nsnps;
    float cn;
    float rdsupport;
    float penalty;
    float mappable;
    double maf;
    
  CNV(int32_t const c, int32_t const s, int32_t const e, int32_t const cil, int32_t const cih, int32_t const cel, int32_t ceh, float const estcn, float const sp, float const pty, float const mp) : chr(c), start(s), end(e), ciposlow(cil), ciposhigh(cih), ciendlow(cel), ciendhigh(ceh), srleft(0), srright(0), nsnps(0), cn(estcn), rdsupport(sp), penalty(pty), mappable(mp), maf(0) {}
  };

    template<typename TCNV>
    struct SortCNVs : public std::binary_function<TCNV, TCNV, bool>
      {
	inline bool operator()(TCNV const& sv1, TCNV const& sv2) {
	  return ((sv1.chr<sv2.chr) || ((sv1.chr==sv2.chr) && (sv1.start<sv2.start)) || ((sv1.chr==sv2.chr) && (sv1.start==sv2.start) && (sv1.end<sv2.end)) || ((sv1.chr==sv2.chr) && (sv1.start==sv2.start) && (sv1.end==sv2.end) && (sv1.penalty > sv2.penalty)));
	}
      };


  struct BpCNV {
    int32_t start;
    int32_t end;
    double zscore;

    BpCNV(int32_t const s, int32_t const e, double const z) : start(s), end(e), zscore(z) {}
  };
  
  template<typename TConfig, typename TGcBias, typename TCoverage>
  inline void
  callCNVs(TConfig const& c, std::pair<uint32_t, uint32_t> const& gcbound, std::vector<uint16_t> const& gcContent, std::vector<uint16_t> const& uniqContent, TGcBias const& gcbias, TCoverage const& cov, bam_hdr_t const* hdr, int32_t const refIndex) {

    // Parameters
    int32_t smallestWin = 2500;
    int32_t biggestWin = 15000;
    uint32_t chain = 10;
    float sdundo = 2;
    
    // Scanning window sizes
    std::vector<int32_t> winsize;
    int32_t wsize = smallestWin;
    while (wsize < biggestWin) {
      winsize.push_back(wsize);
      wsize *= 2;
    }

    // Iterate window sizes, largest to smallest
    typedef int32_t TCnVal;
    typedef std::vector<TCnVal> TCN;
    typedef std::vector<int32_t> TChrPos;
    std::vector<BpCNV> bpvec;
    for(uint32_t idx = 0; idx < winsize.size(); ++idx) {
      TCN cnvec;
      TChrPos wpos;
      uint32_t wstart = 0;
      while(wstart < hdr->target_len[refIndex]) {
	double covsum = 0;
	double expcov = 0;
	int32_t winlen = 0;
	uint32_t pos = wstart;
	while ((winlen < winsize[idx]) && (pos < hdr->target_len[refIndex])) {
	  if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	    covsum += cov[pos];
	    expcov += gcbias[gcContent[pos]].coverage;
	    ++winlen;
	  }
	  ++pos;
	}
	if (winlen == winsize[idx]) {
	  // Full window
	  cnvec.push_back((int32_t) boost::math::round(c.ploidy * covsum / expcov * 100.0));
	  wpos.push_back(wstart);
	}
	wstart = pos;
      }

      // Identify breakpoints
      TCN pre(chain, -1);
      TCN suc(chain, -1);
      TChrPos prep(chain, 0);
      TChrPos sucp(chain, 0);
      uint32_t idxbp = 0;
      for(uint32_t k = 0; k < cnvec.size(); ++k) {
	if (k < chain) {
	  pre[k % chain] = cnvec[k];
	  prep[k % chain] = wpos[k];
	  if (k + 1 < cnvec.size()) {
	    if (idx == 0 ) bpvec.push_back(BpCNV(wpos[k], wpos[k+1], 0));
	    else idxbp += 2*idx;
	  }
	} else if (k < 2 * chain) {
	  suc[k % chain] = cnvec[k];
	  sucp[k % chain] = wpos[k];
	} else {
	  // Midpoint
	  TCnVal val = suc[k%chain];	  
	  int32_t pos = sucp[k%chain];
	  int32_t posNext = sucp[(k+1)%chain];
	  suc[k%chain] = cnvec[k];
	  sucp[k%chain] = wpos[k];

	  // Debug
	  //for(uint32_t m = 0; m < pre.size(); ++m) std::cerr << prep[m] << '\t' << pre[m] << std::endl;
	  //std::cerr << "M:" << pos << '\t' << val << std::endl;
	  //for(uint32_t m = 0; m < suc.size(); ++m) std::cerr << sucp[m] << '\t' << suc[m] << std::endl;

	  // Any shift in CN?
	  boost::accumulators::accumulator_set<TCnVal, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> accpre;
	  boost::accumulators::accumulator_set<TCnVal, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> accsuc;
	  for(uint32_t m = 0; m < pre.size(); ++m) accpre(pre[m]);
	  for(uint32_t m = 0; m < suc.size(); ++m) accsuc(suc[m]);
	  double diff = std::abs(boost::accumulators::mean(accsuc) - boost::accumulators::mean(accpre));
	  // Breakpoint candidate
	  double zscore = 0;
	  if ((diff > sdundo *sqrt(boost::accumulators::variance(accpre))) && (diff > sdundo * sqrt(boost::accumulators::variance(accsuc)))) {
	    zscore = diff / std::max(sqrt(boost::accumulators::variance(accpre)), sqrt(boost::accumulators::variance(accsuc)));
	  }
	  if (idx == 0) bpvec.push_back(BpCNV(pos, posNext, zscore));
	  else {
	    for(uint32_t sub = idxbp; sub < idxbp + 2 * idx; ++sub) bpvec[sub].zscore += zscore;
	    idxbp += 2*idx;
	  }
	  pre[k%chain] = val;
	  prep[k%chain] = pos;
	}
      }

      /*
      // Select local maxima
      std::vector<BpCNV> bpmax;
      if (bpvec.size()) {
	int32_t pos = bpvec[0].start;
	int32_t posNext = bpvec[0].end;
	double bestDiff = bpvec[0].diff;
	for(uint32_t n = 1; n < bpvec.size(); ++n) {
	  if (bpvec[n].start != bpvec[n-1].end) {
	    bpmax.push_back(BpCNV(pos, posNext, bestDiff));
	    pos = bpvec[n].start;
	    posNext = bpvec[n].end;
	    bestDiff = bpvec[n].diff;
	  } else {
	    if (bpvec[n].diff > bestDiff) {
	      pos = bpvec[n].start;
	      posNext = bpvec[n].end;
	      bestDiff = bpvec[n].diff;
	    }
	  }
	}
	bpmax.push_back(BpCNV(pos, posNext, bestDiff));
      }
      */
    }
    
    // Merge local maxima
    for(uint32_t n = 0; n < bpvec.size(); ++n) {
      std::cerr << bpvec[n].start << '-' << bpvec[n].end << ':' << bpvec[n].zscore << std::endl;
    }
     

  }
  

}

#endif
