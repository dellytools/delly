#ifndef CNV_H
#define CNV_H

#include <boost/filesystem.hpp>
#include <boost/multi_array.hpp>
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
#include "util.h"


namespace torali
{

  struct SVBreakpoint {
    int32_t pos;
    int32_t cilow;
    int32_t cihigh;
    int32_t qual;

    explicit SVBreakpoint(int32_t const p) : pos(p), cilow(0), cihigh(0), qual(0) {}
    SVBreakpoint(int32_t const p, int32_t const cil, int32_t const cih, int32_t q) : pos(p), cilow(cil), cihigh(cih), qual(q) {}
  };


  template<typename TSVBp>
  struct SortSVBreakpoint : public std::binary_function<TSVBp, TSVBp, bool>
  {
    inline bool operator()(TSVBp const& sv1, TSVBp const& sv2) {
      return ((sv1.pos<sv2.pos) || ((sv1.pos==sv2.pos) && (sv1.qual<sv2.qual)));
    }
  };

  
  struct BpCNV {
    int32_t start;
    int32_t end;
    double zscore;

    BpCNV(int32_t const s, int32_t const e, double const z) : start(s), end(e), zscore(z) {}
  };


  template<typename TConfig>
  inline void
  mergeCNVs(TConfig const& c, std::vector<CNV>& chrcnv, std::vector<CNV>& cnvs) {
    // Merge neighboring segments if too similar
    bool merged = true;
    std::vector<CNV> newcnv;
    while(merged) {
      int32_t k = -1;
      for(int32_t i = 0; i < (int32_t) chrcnv.size(); ++i) {
	if (i <= k) continue;
	k = i;
	for(int32_t j = i + 1; j < (int32_t) chrcnv.size(); ++j) {
	  bool allValid = true;
	  for(int32_t pre = i; pre < j; ++pre) {
	    double diff = std::abs(chrcnv[pre].cn - chrcnv[j].cn);
	    if (diff >= c.cn_offset) {
	      allValid = false;
	      break;
	    }
	  }
	  if (allValid) k = j;
	  else break;
	}
	if (k > i) {
	  // Merge
	  double cn = (chrcnv[i].cn + chrcnv[k].cn) / 2.0;
	  double mp = (chrcnv[i].mappable + chrcnv[k].mappable) / 2.0;
	  newcnv.push_back(CNV(chrcnv[i].chr, chrcnv[i].start, chrcnv[k].end, chrcnv[i].ciposlow, chrcnv[i].ciposhigh, chrcnv[k].ciendlow, chrcnv[k].ciendhigh, cn, mp));	  
	} else {
	  newcnv.push_back(chrcnv[i]);
	}
      }
      if (newcnv.size() == chrcnv.size()) merged = false;
      else {
	chrcnv = newcnv;
	newcnv.clear();
      }
    }

    // Insert into global CNV vector
    for(uint32_t i = 0; i < chrcnv.size(); ++i) {
      cnvs.push_back(chrcnv[i]);
      //std::cerr << chrcnv[i].chr << '\t' << chrcnv[i].start << '\t' << chrcnv[i].end << "\tMerged" << std::endl;
    }
  }


  template<typename TConfig, typename TGcBias, typename TCoverage, typename TGenomicBreakpoints>
  inline void
  breakpointRefinement(TConfig const& c, std::pair<uint32_t, uint32_t> const& gcbound, std::vector<uint16_t> const& gcContent, std::vector<uint16_t> const& uniqContent, TGcBias const& gcbias, TCoverage const& cov, bam_hdr_t const* hdr, int32_t const refIndex, TGenomicBreakpoints const& svbp, std::vector<CNV>& cnvs) {
    typedef typename TGenomicBreakpoints::value_type TSVs;
    
    // Estimate CN shift
    for(uint32_t n = 1; n < cnvs.size(); ++n) {
      if ((cnvs[n-1].chr != refIndex) || (cnvs[n].chr != refIndex)) continue;
      double precovsum = 0;
      double preexpcov = 0;
      double succovsum = 0;
      double sucexpcov = 0;
      int32_t pos = cnvs[n-1].start;
      while((pos < cnvs[n].end) && (pos < (int32_t) hdr->target_len[refIndex])) {
	if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	  if (pos < cnvs[n-1].end) {
	    precovsum += cov[pos];
	    preexpcov += gcbias[gcContent[pos]].coverage;
	  } else {
	    succovsum += cov[pos];
	    sucexpcov += gcbias[gcContent[pos]].coverage;
	  }
	}
	++pos;
      }
      double precndiff = std::abs((c.ploidy * precovsum / preexpcov) - (c.ploidy * succovsum / sucexpcov));

      // Intersect with delly SVs
      typename TSVs::const_iterator itbest = svbp[refIndex].end();
      int32_t searchStart = std::max(0, std::min(cnvs[n-1].ciendlow, cnvs[n-1].end - 1000));
      int32_t searchEnd = std::max(cnvs[n].ciposhigh, cnvs[n].start + 1000);
      int32_t midpoint = (int32_t) ((cnvs[n-1].start + cnvs[n-1].end) / 2);
      if (searchStart < midpoint) searchStart = midpoint;
      midpoint = (int32_t) ((cnvs[n].start + cnvs[n].end) / 2);
      if (searchEnd > midpoint) searchEnd = midpoint;
      // Current CNV start for this breakpoint
      typename TSVs::const_iterator itsv = std::lower_bound(svbp[refIndex].begin(), svbp[refIndex].end(), SVBreakpoint(searchStart), SortSVBreakpoint<SVBreakpoint>());
      for(; itsv != svbp[refIndex].end(); ++itsv) {
	if (itsv->pos > searchEnd) break;
	if ((itbest == svbp[refIndex].end()) || (itsv->qual > itbest->qual)) itbest = itsv;
      }
      if ((itbest != svbp[refIndex].end()) && (itbest->qual >= 50)) {
	// Check refined CNV
	precovsum = 0;
	preexpcov = 0;
	succovsum = 0;
	sucexpcov = 0;
	pos = cnvs[n-1].start;
	while((pos < cnvs[n].end) && (pos < (int32_t) hdr->target_len[refIndex])) {
	  if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	    if (pos < itbest->pos) {
	      precovsum += cov[pos];
	      preexpcov += gcbias[gcContent[pos]].coverage;
	    } else {
	      succovsum += cov[pos];
	      sucexpcov += gcbias[gcContent[pos]].coverage;
	    }
	  }
	  ++pos;
	}
	double postcndiff = std::abs((c.ploidy * precovsum / preexpcov) - (c.ploidy * succovsum / sucexpcov));
	//std::cerr << cnvs[n-1].end << ',' << itbest->pos << ',' << precndiff << ',' << postcndiff << std::endl;
	if ((precndiff < postcndiff + c.cn_offset) && (std::abs(cnvs[n].start - itbest->pos) < 50000)) {
	  // Accept new breakpoint
	  cnvs[n-1].end = itbest->pos;
	  cnvs[n].start = itbest->pos;
	  cnvs[n-1].ciendlow = itbest->pos + itbest->cilow;
	  cnvs[n-1].ciendhigh = itbest->pos + itbest->cihigh;
	  cnvs[n].ciposlow = itbest->pos + itbest->cilow;
	  cnvs[n].ciposhigh = itbest->pos + itbest->cihigh;
	}
      }
    }
  }
  

  template<typename TConfig, typename TGcBias, typename TCoverage>
  inline void
  breakpointRefinement2(TConfig const& c, std::pair<uint32_t, uint32_t> const& gcbound, std::vector<uint16_t> const& gcContent, std::vector<uint16_t> const& uniqContent, TGcBias const& gcbias, TCoverage const& cov, bam_hdr_t const* hdr, int32_t const refIndex, std::vector<CNV>& cnvs) {

    int32_t maxbpshift = 10000;
	
    // Breakpoint refinement
    for(uint32_t n = 1; n < cnvs.size(); ++n) {
      int32_t prehalf = (cnvs[n-1].start + cnvs[n-1].end) / 2;
      prehalf = std::max(cnvs[n-1].end - maxbpshift, prehalf);
      int32_t suchalf = (cnvs[n].start + cnvs[n].end) / 2;
      suchalf = std::min(cnvs[n].start + maxbpshift, suchalf);
      double precovsum = 0;
      double preexpcov = 0;
      double succovsum = 0;
      double sucexpcov = 0;
      int32_t pos = cnvs[n-1].start;
      std::vector<int32_t> validpos;
      while((pos < cnvs[n].end) && (pos < (int32_t) hdr->target_len[refIndex])) {
	if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	  if (pos < prehalf) {
	    precovsum += cov[pos];
	    preexpcov += gcbias[gcContent[pos]].coverage;
	  } else {
	    if (pos <= suchalf) validpos.push_back(pos);
	    succovsum += cov[pos];
	    sucexpcov += gcbias[gcContent[pos]].coverage;
	  }
	}
	++pos;
      }
      double precn = c.ploidy * precovsum / preexpcov;
      double succn = c.ploidy * succovsum / sucexpcov;
      // Shift Bp
      std::vector<double> diffcn(validpos.size(), 0);
      for(uint32_t idx = 0; idx < validpos.size(); ++idx) {
	if ((preexpcov > 0) && (sucexpcov > 0)) {
	  precn = c.ploidy * precovsum / preexpcov;
	  succn = c.ploidy * succovsum / sucexpcov;
	  diffcn[idx] = std::abs(precn - succn);
	  //if (validpos[idx] == cnvs[n-1].end) std::cerr << "-->";
	  //std::cerr << validpos[idx] << ',' << precn << ',' << succn << ',' << diffcn[idx] << std::endl;
	}
	// Add to pre, remove from suc
	precovsum += cov[validpos[idx]];
	preexpcov += gcbias[gcContent[validpos[idx]]].coverage;
	succovsum -= cov[validpos[idx]];
	sucexpcov -= gcbias[gcContent[validpos[idx]]].coverage;
      }
      // Find best
      int32_t bestIdx = -1;
      for(uint32_t idx = 0; idx < validpos.size(); ++idx) {
	if ((bestIdx == -1) || (diffcn[idx] > diffcn[bestIdx])) bestIdx = idx;
      }
      if (bestIdx != -1) {
	// Update breakpoint
	cnvs[n-1].end = validpos[bestIdx];
	cnvs[n].start = validpos[bestIdx];
      }
    }
    //for(uint32_t n = 0; n < cnvs.size(); ++n) std::cerr << hdr->target_name[cnvs[n].chr] << '\t' << cnvs[n].start << '\t' << cnvs[n].end << "\tRefinement" << std::endl;
  }
  

  template<typename TConfig, typename TGcBias, typename TCoverage>
  inline void
  genotypeCNVs(TConfig const& c, std::pair<uint32_t, uint32_t> const& gcbound, std::vector<uint16_t> const& gcContent, std::vector<uint16_t> const& uniqContent, TGcBias const& gcbias, TCoverage const& cov, bam_hdr_t const* hdr, int32_t const refIndex, std::vector<CNV>& cnvs) {
    for(uint32_t n = 0; n < cnvs.size(); ++n) {
      if (cnvs[n].chr != refIndex) continue;
      double covsum = 0;
      double expcov = 0;
      int32_t winlen = 0;
      int32_t pos = cnvs[n].start;
      while((pos < cnvs[n].end) && (pos < (int32_t) hdr->target_len[refIndex])) {
	if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	  covsum += cov[pos];
	  expcov += gcbias[gcContent[pos]].coverage;
	  ++winlen;
	}
	++pos;
      }
      double cn = c.ploidy;
      if (expcov > 0) cn = c.ploidy * covsum / expcov;
      double mp = (double) winlen / (double) (cnvs[n].end - cnvs[n].start);
      cnvs[n].cn = cn;
      cnvs[n].mappable = mp;

      // Estimate SD
      boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance> > acc;
      uint32_t wsz = winlen / 10;
      if (wsz > 1) {
	covsum = 0;
	expcov = 0;
	winlen = 0;
	pos = cnvs[n].start;
	while((pos < cnvs[n].end) && (pos < (int32_t) hdr->target_len[refIndex])) {
	  if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	    covsum += cov[pos];
	    expcov += gcbias[gcContent[pos]].coverage;
	    ++winlen;
	    if (winlen % wsz == 0) {
	      double cn = c.ploidy;
	      if (expcov > 0) cn = c.ploidy * covsum / expcov;
	      acc(cn);
	      covsum = 0;
	      expcov = 0;
	    }
	  }
	  ++pos;
	}
	cnvs[n].sd = sqrt(boost::accumulators::variance(acc));
	if (cnvs[n].sd < 0.025) cnvs[n].sd = 0.025;
      } else {
	// Invalid
	cnvs[n].cn = -1;
	cnvs[n].sd = 0.025;
      }
    }
  }
  
  template<typename TConfig, typename TGcBias, typename TCoverage>
  inline void
  callCNVs(TConfig const& c, std::pair<uint32_t, uint32_t> const& gcbound, std::vector<uint16_t> const& gcContent, std::vector<uint16_t> const& uniqContent, TGcBias const& gcbias, TCoverage const& cov, bam_hdr_t const* hdr, int32_t const refIndex, std::vector<CNV>& cnvs) {

    // Parameters
    int32_t smallestWin = c.minCnvSize / 10;
    int32_t biggestWin = smallestWin * 200;
    uint32_t chain = 10;

    // Find breakpoints
    std::vector<BpCNV> bpmax;
    if (bpmax.empty()) {
      // Scanning window sizes
      std::vector<int32_t> winsize;
      int32_t wsize = smallestWin;
      while (wsize < biggestWin) {
	winsize.push_back(wsize);
	wsize *= 2;
      }

      // Iterate window sizes
      typedef int32_t TCnVal;
      typedef std::vector<TCnVal> TCN;
      typedef std::vector<int32_t> TChrPos;
      std::vector<BpCNV> bpvec;
      for(uint32_t idx = 0; idx < winsize.size(); ++idx) {
	uint32_t idxOffset = winsize[idx] / winsize[0];
	//std::cerr << idx << ',' << winsize[idx] << ',' << idxOffset << ',' << bpvec.size() << ',' << hdr->target_len[refIndex] << std::endl;
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
	    if (expcov > 0) cnvec.push_back((int32_t) boost::math::round(c.ploidy * covsum / expcov * 100.0));
	    else cnvec.push_back((int32_t) boost::math::round(c.ploidy * 100.0));
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
	      else idxbp += idxOffset;
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
	    boost::accumulators::accumulator_set<TCnVal, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance> > accpre;
	    boost::accumulators::accumulator_set<TCnVal, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance> > accsuc;
	    for(uint32_t m = 0; m < pre.size(); ++m) accpre(pre[m]);
	    for(uint32_t m = 0; m < suc.size(); ++m) accsuc(suc[m]);
	    double diff = std::abs(boost::accumulators::mean(accsuc) - boost::accumulators::mean(accpre));
	    // Breakpoint candidate
	    double zscore = 0;
	    if ((diff > c.stringency * sqrt(boost::accumulators::variance(accpre))) && (diff > c.stringency * sqrt(boost::accumulators::variance(accsuc)))) {
	      zscore = diff / std::max(sqrt(boost::accumulators::variance(accpre)), sqrt(boost::accumulators::variance(accsuc)));
	    }
	    if (idx == 0) bpvec.push_back(BpCNV(pos, posNext, zscore));
	    else {
	      for(uint32_t sub = idxbp; sub < idxbp + idxOffset; ++sub) bpvec[sub].zscore += zscore;
	      idxbp += idxOffset;
	    }
	    pre[k%chain] = val;
	    prep[k%chain] = pos;
	  }
	}
      }
    
      // Local maxima
      if (bpvec.size()) {
	int32_t pos = bpvec[0].start;
	int32_t posNext = bpvec[0].end;
	double bestDiff = bpvec[0].zscore;
	for(uint32_t n = 1; n < bpvec.size(); ++n) {
	  //std::cerr << "B:" << bpvec[n].start << '-' << bpvec[n].end << ':' << bpvec[n].zscore << std::endl;
	  if (bpvec[n].zscore == 0) {
	    if (bestDiff != 0) {
	      //std::cerr << "M:" << pos << '-' << posNext << ':' << bestDiff << std::endl;
	      bpmax.push_back(BpCNV(pos, posNext, bestDiff));
	      pos = bpvec[n].start;
	      posNext = bpvec[n].end;
	      bestDiff = bpvec[n].zscore;
	    }
	  } else {
	    if (bpvec[n].zscore > bestDiff) {
	      // Replace local max
	      pos = bpvec[n].start;
	      posNext = bpvec[n].end;
	      bestDiff = bpvec[n].zscore;
	    } else if (bpvec[n].zscore == bestDiff) {
	      // Extend local max
	      posNext = bpvec[n].end;
	    }
	  }
	}
      }
    }

    // Breakpoints
    for(uint32_t n = 0; n <= bpmax.size(); ++n) {
      int32_t cil = 0;
      int32_t cih = 0;
      if (n > 0) {
	cil = bpmax[n-1].start;
	cih = bpmax[n-1].end;
      }
      int32_t cel = hdr->target_len[refIndex] - 1;
      int32_t ceh = hdr->target_len[refIndex] - 1;
      if (n < bpmax.size()) {
	cel = bpmax[n].start;
	ceh = bpmax[n].end;
      }
      //std::cerr << (cih - cil) << ';' << (ceh - cel) << std::endl;
      int32_t cnvstart = (int32_t) ((cil + cih)/2);
      int32_t cnvend = (int32_t) ((cel + ceh)/2);
      int32_t estcnvstart = -1;
      int32_t estcnvend = -1;
      double covsum = 0;
      double expcov = 0;
      int32_t winlen = 0;
      int32_t pos = cnvstart;
      while((pos < cnvend) && (pos < (int32_t) hdr->target_len[refIndex])) {
	if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	  if (estcnvstart == -1) estcnvstart = pos;
	  estcnvend = pos;
	  covsum += cov[pos];
	  expcov += gcbias[gcContent[pos]].coverage;
	  ++winlen;
	}
	++pos;
      }
      if ((estcnvstart != -1) && (estcnvend != -1) && (estcnvend - estcnvstart > 0)) {
	double cn = c.ploidy;
	if (expcov > 0) cn = c.ploidy * covsum / expcov;
	double mp = (double) winlen / (double) (estcnvend - estcnvstart);
	cnvs.push_back(CNV(refIndex, estcnvstart, estcnvend, cil, cih, cel, ceh, cn, mp));
	//std::cerr << hdr->target_name[refIndex] << '\t' << estcnvstart << '\t' << estcnvend << '\t' << '(' << cil << ',' << cih << ')' << '\t' << '(' << cel << ',' << ceh << ')' << '\t' << cn << '\t' << mp << std::endl;
      }
    }
  }


  // Parse Delly CNV VCF file
  template<typename TConfig>
  inline void
  parseVcfCNV(TConfig const& c, bam_hdr_t* hd, std::vector<CNV>& cnvs) {
    // Load bcf file
    htsFile* ifile = bcf_open(c.genofile.string().c_str(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);
    bcf1_t* rec = bcf_init();

    // Parse bcf
    int32_t nsvend = 0;
    int32_t* svend = NULL;
    int32_t ncipos = 0;
    int32_t* cipos = NULL;
    int32_t nmp = 0;
    float* mp = NULL;
    int32_t nsvt = 0;
    char* svt = NULL;
    int32_t nmethod = 0;
    char* method = NULL;
    uint16_t wimethod = 0; 
    while (bcf_read(ifile, hdr, rec) == 0) {
      bcf_unpack(rec, BCF_UN_INFO);

      // Delly BCF file?
      if (!wimethod) {
	wimethod = 2;
	if (bcf_get_info_string(hdr, rec, "SVMETHOD", &method, &nmethod) > 0) {
	  std::string mstr = std::string(method);
	  if ((mstr.size() >= 10) && (mstr.substr(0, 10)  == "EMBL.DELLY")) wimethod = 1;
	}
      }

      // Delly
      if (wimethod == 1) {
	// Fill SV record
	CNV cnv;
	std::string chrName = bcf_hdr_id2name(hdr, rec->rid);
	int32_t tid = bam_name2id(hd, chrName.c_str());
	cnv.chr = tid;
	cnv.start = rec->pos - 1;
	cnv.qval = rec->qual;

	// Parse CNV type
	if (bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt) > 0) {
	  if (std::string(svt) != "CNV") continue;
	} else continue;
      
	// Parse INFO
	if (bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend) > 0) cnv.end = *svend;
	else continue;
	if (bcf_get_info_int32(hdr, rec, "CIPOS", &cipos, &ncipos) > 0) {
	  cnv.ciposlow = cnv.start + cipos[0];
	  cnv.ciposhigh = cnv.start + cipos[1];
	} else {
	  cnv.ciposlow = cnv.start - 50;
	  cnv.ciposhigh = cnv.start + 50;
	}
	if (bcf_get_info_int32(hdr, rec, "CIEND", &cipos, &ncipos) > 0) {
	  cnv.ciendlow = cnv.end + cipos[0];
	  cnv.ciendhigh = cnv.end + cipos[1];
	} else {
	  cnv.ciendlow = cnv.end - 50;
	  cnv.ciendhigh = cnv.end + 50;
	}
	if (bcf_get_info_float(hdr, rec, "MP", &mp, &nmp) > 0) cnv.mappable = (double) *mp;
	else cnv.mappable = 0;

	cnvs.push_back(cnv);
      }
    }
    // Clean-up
    free(svend);
    free(svt);
    free(method);
    free(cipos);
    free(mp);
    
    // Close VCF
    bcf_hdr_destroy(hdr);
    bcf_close(ifile);
    bcf_destroy(rec);
  }

  
  template<typename TConfig>
  inline void
  cnvVCF(TConfig const& c, std::vector<CNV> const& cnvs) {
    // Open one bam file header
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* bamhd = sam_hdr_read(samfile);

    // Output all copy-number variants
    std::string fmtout = "wb";
    if (c.outfile.string() == "-") fmtout = "w";
    htsFile *fp = hts_open(c.outfile.string().c_str(), fmtout.c_str());
    bcf_hdr_t *hdr = bcf_hdr_init("w");

    // Print vcf header
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    boost::gregorian::date today = now.date();
    std::string datestr("##fileDate=");
    datestr += boost::gregorian::to_iso_string(today);
    bcf_hdr_append(hdr, datestr.c_str());
    bcf_hdr_append(hdr, "##ALT=<ID=CNV,Description=\"copy-number variants\">");
    bcf_hdr_append(hdr, "##FILTER=<ID=LowQual,Description=\"Poor quality copy-number variant\">");
    bcf_hdr_append(hdr, "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END\">");
    bcf_hdr_append(hdr, "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS\">");
    bcf_hdr_append(hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the copy-number variant\">");
    bcf_hdr_append(hdr, "##INFO=<ID=MP,Number=1,Type=Float,Description=\"Mappable fraction of CNV\">");
    bcf_hdr_append(hdr, "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise copy-number variant\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect CNV\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Integer copy-number\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=CNL,Number=.,Type=Float,Description=\"Log10-scaled copy-number likelihoods\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Per-sample genotype filter\">");

    bcf_hdr_append(hdr, "##FORMAT=<ID=RDCN,Number=1,Type=Float,Description=\"Read-depth based copy-number estimate\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=RDSD,Number=1,Type=Float,Description=\"Read-depth standard deviation\">");

    // Add reference
    std::string refloc("##reference=");
    refloc += c.genome.string();
    bcf_hdr_append(hdr, refloc.c_str());
    for (int i = 0; i<bamhd->n_targets; ++i) {
      std::string refname("##contig=<ID=");
      refname += std::string(bamhd->target_name[i]) + ",length=" + boost::lexical_cast<std::string>(bamhd->target_len[i]) + ">";
      bcf_hdr_append(hdr, refname.c_str());
    }
    // Add samples
    bcf_hdr_add_sample(hdr, c.sampleName.c_str());
    bcf_hdr_add_sample(hdr, NULL);
    if (bcf_hdr_write(fp, hdr) != 0) std::cerr << "Error: Failed to write BCF header!" << std::endl;

    uint32_t cnvid = 0;
    if (!cnvs.empty()) {
      // Genotype arrays
      int32_t *gts = (int*) malloc(bcf_hdr_nsamples(hdr) * 2 * sizeof(int));
      int32_t *gqval = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
      int32_t *cnval = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
      float *cnrdval = (float*) malloc(bcf_hdr_nsamples(hdr) * sizeof(float));
      float *cnsdval = (float*) malloc(bcf_hdr_nsamples(hdr) * sizeof(float));
      float *cnl = (float*) malloc(bcf_hdr_nsamples(hdr) * MAX_CN * sizeof(float));

      std::vector<std::string> ftarr;
      ftarr.resize(bcf_hdr_nsamples(hdr));
    
      // Iterate all structural variants
      now = boost::posix_time::second_clock::local_time();
      std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Genotyping" << std::endl;
      bcf1_t *rec = bcf_init();
      for(uint32_t i = 0; i < cnvs.size(); ++i) {
	// Invalid CNV?
	if ((!c.hasGenoFile) && (cnvs[i].cn == -1)) continue;

	// Integer copy-number
	int32_t absCN = (int32_t) boost::math::round(cnvs[i].cn);
	if ((!c.segmentation) && (absCN == c.ploidy)) continue;
      
	// Output main vcf fields
	rec->rid = bcf_hdr_name2id(hdr, bamhd->target_name[cnvs[i].chr]);
	int32_t svStartPos = cnvs[i].start + 1;
	int32_t svEndPos = cnvs[i].end;
	if (svEndPos >= (int32_t) bamhd->target_len[cnvs[i].chr]) svEndPos = bamhd->target_len[cnvs[i].chr] - 1;
	rec->pos = svStartPos;
	std::string id("CNV");
	std::string padNumber = boost::lexical_cast<std::string>(++cnvid);
	padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
	id += padNumber;
	bcf_update_id(hdr, rec, id.c_str());
	std::string svtype = "CNV";
	std::string alleles = "N,<" + svtype + ">";
	bcf_update_alleles_str(hdr, rec, alleles.c_str());
      
	// Add INFO fields
	bcf_update_info_flag(hdr, rec, "IMPRECISE", NULL, 1);

	bcf_update_info_string(hdr, rec, "SVTYPE", svtype.c_str());
	std::string dellyVersion("EMBL.DELLYv");
	dellyVersion += dellyVersionNumber;
	bcf_update_info_string(hdr,rec, "SVMETHOD", dellyVersion.c_str());
	int32_t tmpi = svEndPos;
	bcf_update_info_int32(hdr, rec, "END", &tmpi, 1);
	int32_t ciend[2];
	ciend[0] = cnvs[i].ciendlow - cnvs[i].end;
	ciend[1] = cnvs[i].ciendhigh - cnvs[i].end;
	int32_t cipos[2];
	cipos[0] = cnvs[i].ciposlow - cnvs[i].start;
	cipos[1] = cnvs[i].ciposhigh - cnvs[i].start;
	bcf_update_info_int32(hdr, rec, "CIPOS", cipos, 2);
	bcf_update_info_int32(hdr, rec, "CIEND", ciend, 2);
	float tmpf = cnvs[i].mappable;
	bcf_update_info_float(hdr, rec, "MP", &tmpf, 1);

	// Genotyping
	cnval[0] = absCN;
	cnrdval[0] = cnvs[i].cn;
	cnsdval[0] = cnvs[i].sd;
	gts[0] = bcf_gt_missing;
	gts[1] = bcf_gt_missing;
	int32_t qval = _computeCNLs(c, cnvs[i].cn, cnvs[i].sd, cnl, gqval);
	if (c.hasGenoFile) rec->qual = cnvs[i].qval;  // Leave site quality in genotyping mode
	else rec->qual = qval;
	tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
	if (rec->qual < 15) tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "LowQual");
	bcf_update_filter(hdr, rec, &tmpi, 1);
	
	if (gqval[0] < 15) ftarr[0] = "LowQual";
	else ftarr[0] = "PASS";
	std::vector<const char*> strp(bcf_hdr_nsamples(hdr));
	std::transform(ftarr.begin(), ftarr.end(), strp.begin(), cstyle_str());	
	bcf_update_genotypes(hdr, rec, gts, bcf_hdr_nsamples(hdr) * 2);
	bcf_update_format_int32(hdr, rec, "CN", cnval, bcf_hdr_nsamples(hdr));
	bcf_update_format_float(hdr, rec, "CNL",  cnl, bcf_hdr_nsamples(hdr) * MAX_CN);
	bcf_update_format_int32(hdr, rec, "GQ", gqval, bcf_hdr_nsamples(hdr));
	bcf_update_format_string(hdr, rec, "FT", &strp[0], bcf_hdr_nsamples(hdr));
	bcf_update_format_float(hdr, rec, "RDCN",  cnrdval, bcf_hdr_nsamples(hdr));
	bcf_update_format_float(hdr, rec, "RDSD",  cnsdval, bcf_hdr_nsamples(hdr));
	bcf_write1(fp, hdr, rec);
	bcf_clear1(rec);
      }
      bcf_destroy1(rec);
    
      // Clean-up
      free(gts);
      free(gqval);
      free(cnval);
      free(cnrdval);
      free(cnsdval);
      free(cnl);
    }

    // Close BAM file
    bam_hdr_destroy(bamhd);
    sam_close(samfile);
    
    // Close VCF file
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    
    // Build index
    if (c.outfile.string() != "-") bcf_index_build(c.outfile.string().c_str(), 14);
  }
 

}

#endif
