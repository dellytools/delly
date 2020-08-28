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
      
  
  template<typename TConfig, typename TGcBias, typename TCoverage>
  inline void
  callCNVs(TConfig const& c, std::pair<uint32_t, uint32_t> const& gcbound, std::vector<uint16_t> const& gcContent, std::vector<uint16_t> const& uniqContent, TGcBias const& gcbias, TCoverage const& cov, bam_hdr_t const* hdr, int32_t const refIndex) {
    
    // cnv tiling windows, window_offset smallest window, window_size largest window
    uint32_t numTiling = 0;
    {
      uint32_t winbound = c.window_offset;
      while (winbound < c.window_size) {
	++numTiling;
	winbound *= 2;
      }
    }

    typedef std::vector<int16_t> TCNVec;
    typedef std::vector<TCNVec> TTilingCN;
    TTilingCN cnvec(numTiling, TCNVec());
    for(uint32_t i = 0; i < numTiling; ++i) cnvec[i].resize(hdr->target_len[refIndex]/50, -1);
    for(uint32_t start = 0; start < hdr->target_len[refIndex]; start = start + 50) {
      if (start + c.window_offset < hdr->target_len[refIndex]) {
	double covsum = 0;
	double expcov = 0;
	uint32_t winlen = 0;
	uint32_t winbound = c.window_offset;
	uint32_t tilingPos = 0;
	for(uint32_t pos = start; ((pos < start + c.window_size) && (pos < hdr->target_len[refIndex])); ++pos) {
	  if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	    covsum += cov[pos];
	    expcov += gcbias[gcContent[pos]].coverage;
	    ++winlen;
	  }
	  // Multiple of window size?
	  if ((pos - start) == winbound) {
	    if (winlen >= c.fracWindow * (pos - start)) {
	      cnvec[tilingPos][start/hdr->target_len[refIndex]] = (int32_t) boost::math::round(c.ploidy * covsum / expcov * 100.0);
	    }
	    winbound *= 2;
	    ++tilingPos;
	  }
	}
      }
    }

    // Iterate CNVs
    for(int32_t i = numTiling - 1; i>=0; --i) {
      for(uint32_t k = 0; k < cnvec[i].size(); ++k) {
	std::cout << cnvec[i][k] << ',';
      }
      std::cout << std::endl;
    }
  }
  

}

#endif
