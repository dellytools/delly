#ifndef CNV_H
#define CNV_H

#include <limits>
#include <algorithm>
#include <fstream>

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
    int32_t support;

    explicit SVBreakpoint(int32_t const p) : pos(p), cilow(0), cihigh(0), qual(0), support(0) {}
    SVBreakpoint(int32_t const p, int32_t const cil, int32_t const cih, int32_t q, int32_t const sup) : pos(p), cilow(cil), cihigh(cih), qual(q), support(sup) {}

    bool operator<(const SVBreakpoint& sv2) const {
      return ((pos<sv2.pos) || ((pos==sv2.pos) && (qual<sv2.qual)));
    }
  };

  // Segment boundary
  struct CnvBoundary {
    int32_t w;
    int32_t bp;
    int32_t sr;

    CnvBoundary(int32_t const w_, int32_t const bp_, int32_t const sr_) : w(w_), bp(bp_), sr(sr_) {}

    bool operator<(const CnvBoundary& o) const {
      // keep better-supported junction for same window
      return ((w < o.w) || ((w == o.w) && (sr > o.sr)));
    }
  };

  // Collect candidate CNV boundaries
  template<typename TConfig, typename TGcBias, typename TCoverage>
  inline void
  collectBreakpoints(TConfig const& c, std::pair<uint32_t, uint32_t> const& gcbound, std::vector<uint16_t> const& gcContent, std::vector<uint16_t> const& uniqContent, TGcBias const& gcbias, TCoverage const& cov, bam_hdr_t const* hdr, int32_t const refIndex, std::vector<int32_t>& clips, std::vector<SVBreakpoint>& chrbp) {
    if (clips.empty()) return;
    std::sort(clips.begin(), clips.end());
    int32_t bpTol = (int32_t) (2 * c.minClip);
    double flankExpTarget = (c.targetExpCov > 0) ? c.targetExpCov : 1000.0;
    int32_t maxFlank = 1000000;
    uint32_t i = 0;
    while (i < clips.size()) {
      // Cluster split-reads
      uint32_t j = i;
      int64_t possum = clips[i];
      uint32_t support = 1;
      while ((j + 1 < clips.size()) && (clips[j+1] - clips[j] <= bpTol)) {
	++j;
	possum += clips[j];
	++support;
      }
      if (support >= c.minBpSupport) {
	int32_t bppos = (int32_t) (possum / support);
	// Read-depth shift left/right of breakpoint
	double lcov = 0;
	double lexp = 0;
	int32_t lspan = 0;
	int32_t pos = bppos - 1;
	while ((pos >= 0) && (lexp < flankExpTarget) && (lspan < maxFlank)) {
	  if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	    lcov += cov[pos];
	    lexp += gcbias[gcContent[pos]].coverage;
	  }
	  --pos;
	  ++lspan;
	}
	double rcov = 0;
	double rexp = 0;
	int32_t rspan = 0;
	pos = bppos;
	while ((pos < (int32_t) hdr->target_len[refIndex]) && (rexp < flankExpTarget) && (rspan < maxFlank)) {
	  if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	    rcov += cov[pos];
	    rexp += gcbias[gcContent[pos]].coverage;
	  }
	  ++pos;
	  ++rspan;
	}
	if ((lexp >= 0.5 * flankExpTarget) && (rexp >= 0.5 * flankExpTarget)) {
	  double cnL = c.ploidy * lcov / lexp;
	  double cnR = c.ploidy * rcov / rexp;
	  if (std::abs(cnL - cnR) >= c.minCnShift) {
	    int32_t qual = 50 + (int32_t) std::min(support, (uint32_t) 40);
	    chrbp.push_back(SVBreakpoint(bppos, -bpTol, bpTol, qual, (int32_t) support));
	  }
	}
      }
      i = j + 1;
    }
    std::sort(chrbp.begin(), chrbp.end());
  }


  template<typename TConfig, typename TGcBias, typename TCoverage>
  inline void
  genotypeCNVs(TConfig const& c, std::pair<uint32_t, uint32_t> const& gcbound, std::vector<uint16_t> const& gcContent, std::vector<uint16_t> const& uniqContent, TGcBias const& gcbias, TCoverage const& cov, TCoverage const& covUniq, TCoverage const& covMap, bam_hdr_t const* hdr, int32_t const refIndex, std::vector<CNV>& cnvs) {
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

      // Uniquely-mappable
      double ufrac = -1;
      if (c.hasMapFile) {
	// Mappability map
	double umapsum = 0;
	int32_t uspan = 0;
	for(int32_t p = cnvs[n].start; (p < cnvs[n].end) && (p < (int32_t) hdr->target_len[refIndex]); ++p) {
	  umapsum += (double) uniqContent[p];
	  ++uspan;
	}
	if ((uspan > 0) && (c.meanisize > 0)) ufrac = std::min(1.0, umapsum / ((double) uspan * (double) c.meanisize));
      } else {
	// Map free
	double ucov = 0;
	double tcov = 0;
	for(int32_t p = cnvs[n].start; (p < cnvs[n].end) && (p < (int32_t) hdr->target_len[refIndex]); ++p) {
	  ucov += covUniq[p];
	  tcov += covMap[p];
	}
	if (tcov > 0) ufrac = ucov / tcov;
      }
      cnvs[n].uniqfrac = ufrac;

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
  
  // Merge segments
  inline void
  mergeAdjacentSameCN(std::vector<CNV>& cnvs, double const mergeTol) {
    if (cnvs.empty()) return;
    std::vector<CNV> out;
    out.push_back(cnvs[0]);
    for(uint32_t i = 1; i < cnvs.size(); ++i) {
      CNV& prev = out.back();
      CNV const& cur = cnvs[i];
      // Same CN?
      bool sameCN = false;
      if ((prev.cn >= 0) && (cur.cn >= 0)) {
	double zl = std::log2(std::max(prev.cn, 0.03));
	double zr = std::log2(std::max(cur.cn, 0.03));
	sameCN = (std::abs(zl - zr) < mergeTol);
      }
      if (sameCN && (prev.chr == cur.chr)) {
	double w1 = (double) (prev.end - prev.start);
	double w2 = (double) (cur.end - cur.start);
	double wsum = (w1 + w2 > 0) ? (w1 + w2) : 1.0;
	prev.cn = (prev.cn * w1 + cur.cn * w2) / wsum;
	prev.mappable = (prev.mappable * w1 + cur.mappable * w2) / wsum;
	prev.sd = (prev.sd * w1 + cur.sd * w2) / wsum;
	prev.ciendlow = cur.ciendlow;
	prev.ciendhigh = cur.ciendhigh;
	prev.end = cur.end;
	prev.srright = cur.srright;
      } else out.push_back(cur);
    }
    cnvs.swap(out);
  }

  // piecewise segmentation
  inline void
  cnvSegment(std::vector<double> const& y, double const beta, int32_t const kmin, std::vector<int32_t>& bnd) {
    int32_t N = (int32_t) y.size();
    bnd.clear();
    if (N < 2 * kmin) return; // one segment
    std::vector<double> S1(N + 1, 0.0);
    std::vector<double> S2(N + 1, 0.0);
    for(int32_t i = 0; i < N; ++i) {
      S1[i+1] = S1[i] + y[i];
      S2[i+1] = S2[i] + y[i] * y[i];
    }
    std::vector<double> F(N + 1, 0.0);
    std::vector<int32_t> prev(N + 1, 0);
    F[0] = -beta;
    std::vector<int32_t> R;
    for(int32_t t = kmin; t <= N; ++t) {
      if ((t == kmin) || (t >= 2 * kmin)) R.push_back(t - kmin);
      double best = std::numeric_limits<double>::max();
      int32_t bestS = 0;
      for(uint32_t ri = 0; ri < R.size(); ++ri) {
	int32_t s = R[ri];
	double n = (double) (t - s);
	double sm = S1[t] - S1[s];
	double val = F[s] + ((S2[t] - S2[s]) - sm * sm / n) + beta;
	if (val < best) {
	  best = val;
	  bestS = s;
	}
      }
      F[t] = best;
      prev[t] = bestS;
      
      // Prune
      std::vector<int32_t> Rn;
      Rn.reserve(R.size());
      for(uint32_t ri = 0; ri < R.size(); ++ri) {
	int32_t s = R[ri];
	double n = (double) (t - s);
	double sm = S1[t] - S1[s];
	if (F[s] + ((S2[t] - S2[s]) - sm * sm / n) <= F[t]) Rn.push_back(s);
      }
      R = Rn;
    }
    
    // Backtrack internal boundaries
    std::vector<int32_t> rev;
    int32_t t = N;
    while (t > 0) {
      int32_t s = prev[t];
      if (s > 0) rev.push_back(s);
      if (s >= t) break;
      t = s;
    }
    for(int32_t i = (int32_t) rev.size() - 1; i >= 0; --i) bnd.push_back(rev[i]);
  }

  // Segment read-depth
  template<typename TConfig, typename TGcBias, typename TCoverage>
  inline void
  segmentRD(TConfig const& c, std::pair<uint32_t, uint32_t> const& gcbound, std::vector<uint16_t> const& gcContent, std::vector<uint16_t> const& uniqContent, TGcBias const& gcbias, TCoverage const& cov, bam_hdr_t const* hdr, int32_t const refIndex, std::vector<SVBreakpoint> const& chrbp, std::vector<CNV>& cnvs) {
    int32_t reflen = (int32_t) hdr->target_len[refIndex];
    int32_t kmin = 4;
    int32_t bpTol = (int32_t) (2 * c.minClip);
    
    // Coverage-adaptive segmentation
    double pcfTargetExp = (c.targetExpCov > 0) ? (c.targetExpCov / 8.0) : 0.0;
    int32_t pcfWinBases = std::max(1, (int32_t) (c.minCnvSize / 10));

    // Per-window profile using log2 scale
    double const rFloor = 1.0 / 64.0;
    std::vector<double> z;
    std::vector<double> wcov;
    std::vector<double> wexp;
    std::vector<int32_t> ws;
    std::vector<int32_t> we;
    {
      double covsum = 0;
      double expcov = 0;
      int32_t winlen = 0;
      int32_t start = -1;
      int32_t last = -1;
      for(int32_t pos = 0; pos < reflen; ++pos) {
	if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	  if (start < 0) start = pos;
	  covsum += cov[pos];
	  expcov += gcbias[gcContent[pos]].coverage;
	  last = pos;
	  ++winlen;
	  bool close = (pcfTargetExp > 0) ? (expcov >= pcfTargetExp) : (winlen >= pcfWinBases);
	  if (close) {
	    double r = (expcov > 0) ? (covsum / expcov) : 1.0;
	    z.push_back(std::log2(std::max(r, rFloor)));
	    wcov.push_back(covsum);
	    wexp.push_back(expcov);
	    ws.push_back(start);
	    we.push_back(pos + 1);
	    covsum = 0;
	    expcov = 0;
	    winlen = 0;
	    start = -1;
	  }
	}
      }
      if ((winlen > 0) && (start >= 0)) {
	double r = (expcov > 0) ? (covsum / expcov) : 1.0;
	z.push_back(std::log2(std::max(r, rFloor)));
	wcov.push_back(covsum);
	wexp.push_back(expcov);
	ws.push_back(start);
	we.push_back(last + 1);
      }
    }
    int32_t N = (int32_t) z.size();
    if (N < 1) return;

    // per-window noise
    double sigma = 0.02;
    if (N > 1) {
      std::vector<double> diff;
      diff.reserve(N - 1);
      for(int32_t i = 1; i < N; ++i) diff.push_back(std::abs(z[i] - z[i-1]));
      std::sort(diff.begin(), diff.end());
      int32_t m = std::max(1, (2 * (int32_t) diff.size()) / 3);
      double s = 0;
      for(int32_t i = 0; i < m; ++i) s += diff[i];
      sigma = (s / m) / 1.128379;
    }
    if (sigma < 0.02) sigma = 0.02;

    // Segmentation
    double beta = c.penalty * sigma * sigma * std::log((double) std::max(N, 2));
    std::vector<int32_t> pcfbnd;
    cnvSegment(z, beta, kmin, pcfbnd);

    // Boundary set
    std::vector<CnvBoundary> B;
    B.push_back(CnvBoundary(0, -1, 0));
    for(uint32_t i = 0; i < pcfbnd.size(); ++i) B.push_back(CnvBoundary(pcfbnd[i], -1, 0));
    B.push_back(CnvBoundary(N, -1, 0));

    // Fuse split-read breakpoints
    for(uint32_t k = 0; k < chrbp.size(); ++k) {
      int32_t bppos = chrbp[k].pos;
      int32_t lo = 0, hi = N;
      while (lo < hi) {
	int32_t mid = (lo + hi) / 2;
	if (ws[mid] < bppos) lo = mid + 1;
	else hi = mid;
      }
      int32_t wi = lo;
      if ((wi <= 0) || (wi >= N)) continue;
      int32_t bi = 0;
      for(uint32_t b = 1; b + 1 < B.size(); ++b) {
	if (std::abs(B[b].w - wi) < std::abs(B[bi].w - wi)) bi = b;
      }
      if ((bi > 0) && (std::abs(B[bi].w - wi) <= 1)) {
	B[bi].w = wi;
	B[bi].bp = bppos;
	B[bi].sr = chrbp[k].support;
      } else {
	B.push_back(CnvBoundary(wi, bppos, chrbp[k].support));
      }
    }
    std::sort(B.begin(), B.end());
    B.erase(std::unique(B.begin(), B.end(), [](CnvBoundary const& a, CnvBoundary const& b){ return a.w == b.w; }), B.end());

    // Segment sums
    uint32_t ns = B.size() - 1;
    std::vector<double> segcov(ns, 0);
    std::vector<double> segexp(ns, 0);
    std::vector<int32_t> segnw(ns, 0);
    for(uint32_t s = 0; s < ns; ++s) {
      for(int32_t w = B[s].w; w < B[s+1].w; ++w) {
	segcov[s] += wcov[w];
	segexp[s] += wexp[w];
	++segnw[s];
      }
    }

    // Merge neighbours
    double const zK = 3.0;
    double const zFloor = c.cnMergeTol;
    while (ns > 1) {
      int32_t best = -1;
      double bestDz = 0;
      for(uint32_t s = 0; s + 1 < ns; ++s) {
	if (B[s+1].bp >= 0) continue;  // keep split-read boundaries
	double cnL = (segexp[s] > 0) ? (c.ploidy * segcov[s] / segexp[s]) : (double) c.ploidy;
	double cnR = (segexp[s+1] > 0) ? (c.ploidy * segcov[s+1] / segexp[s+1]) : (double) c.ploidy;
	double dz = std::abs(std::log2(std::max(cnL / c.ploidy, rFloor)) - std::log2(std::max(cnR / c.ploidy, rFloor)));
	double se = sigma * std::sqrt(1.0 / (double) std::max(segnw[s], 1) + 1.0 / (double) std::max(segnw[s+1], 1));
	double tol = std::max(zFloor, zK * se);
	if ((dz < tol) && ((best < 0) || (dz < bestDz))) { best = (int32_t) s; bestDz = dz; }
      }
      if (best < 0) break;
      uint32_t s = (uint32_t) best;
      segcov[s] += segcov[s+1];
      segexp[s] += segexp[s+1];
      segnw[s] += segnw[s+1];
      B.erase(B.begin() + s + 1);
      segcov.erase(segcov.begin() + s + 1);
      segexp.erase(segexp.begin() + s + 1);
      segnw.erase(segnw.begin() + s + 1);
      --ns;
    }

    // Output segments
    for(uint32_t s = 0; s < ns; ++s) {
      int32_t wa = B[s].w;
      int32_t wb = B[s+1].w;
      if (wb <= wa) continue;
      int32_t start = (B[s].bp >= 0) ? B[s].bp : ws[wa];
      int32_t end = (B[s+1].bp >= 0) ? B[s+1].bp : we[wb-1];
      int32_t cil = (B[s].bp >= 0) ? (start - bpTol) : ws[wa];
      int32_t cih = (B[s].bp >= 0) ? (start + bpTol) : (we[wa] - 1);
      int32_t cel = (B[s+1].bp >= 0) ? (end - bpTol) : ws[wb-1];
      int32_t ceh = (B[s+1].bp >= 0) ? (end + bpTol) : (we[wb-1]);
      double cn = (segexp[s] > 0) ? (c.ploidy * segcov[s] / segexp[s]) : (double) c.ploidy;
      CNV cnvRec(refIndex, start, end, cil, cih, cel, ceh, cn, 1.0);
      cnvRec.srleft = B[s].sr;
      cnvRec.srright = B[s+1].sr;
      cnvs.push_back(cnvRec);
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
    int32_t nsrl = 0;
    int32_t* srl = NULL;
    int32_t nsrr = 0;
    int32_t* srr = NULL;
    int32_t nsvt = 0;
    char* svt = NULL;
    int32_t nmethod = 0;
    char* method = NULL;
    uint16_t wimethod = 0; 
    while (bcf_read(ifile, hdr, rec) == 0) {
      bcf_unpack(rec, BCF_UN_SHR);

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
	cnv.start = rec->pos;
	cnv.qval = rec->qual;
	if ((rec->d.id != NULL) && (std::string(rec->d.id) != ".")) cnv.id = std::string(rec->d.id);

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
	if (bcf_get_info_int32(hdr, rec, "SRL", &srl, &nsrl) > 0) cnv.srleft = *srl;
	else cnv.srleft = 0;
	if (bcf_get_info_int32(hdr, rec, "SRR", &srr, &nsrr) > 0) cnv.srright = *srr;
	else cnv.srright = 0;

	cnvs.push_back(cnv);
      }
    }
    // Clean-up
    free(svend);
    free(svt);
    free(method);
    free(cipos);
    free(mp);
    free(srl);
    free(srr);
    
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
    bcf_hdr_append(hdr, "##INFO=<ID=MP,Number=1,Type=Float,Description=\"Callable (GC- and mappability-passing) fraction of the CNV span\">");
    bcf_hdr_append(hdr, "##INFO=<ID=UNIQ,Number=1,Type=Float,Description=\"Uniquely-mappable fraction\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SRL,Number=1,Type=Integer,Description=\"Split-read support at the left breakpoint\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SRR,Number=1,Type=Integer,Description=\"Split-read support at the right breakpoint\">");
    bcf_hdr_append(hdr, "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise copy-number variant\">");
    bcf_hdr_append(hdr, "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise copy-number variant\">");
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

    // Segmentation BED output
    std::ofstream segOut;
    if (c.hasSegFile) segOut.open(c.segfile.string().c_str());

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

	// Segmentation
	if (c.hasSegFile) segOut << bamhd->target_name[cnvs[i].chr] << '\t' << cnvs[i].start << '\t' << cnvs[i].end << "\tSEG" << (i + 1) << '\t' << cnvs[i].cn << '\n';

	// Only true CNVs, unless in genotyping mode
	if ((!c.hasGenoFile) && (absCN == c.ploidy)) continue;

	// Output main vcf fields
	rec->rid = bcf_hdr_name2id(hdr, bamhd->target_name[cnvs[i].chr]);
	int32_t svStartPos = cnvs[i].start;
	int32_t svEndPos = cnvs[i].end;
	if (svEndPos >= (int32_t) bamhd->target_len[cnvs[i].chr]) svEndPos = bamhd->target_len[cnvs[i].chr] - 1;
	rec->pos = svStartPos;
	std::string id;
	if ((c.hasGenoFile) && (!cnvs[i].id.empty())) id = cnvs[i].id;
	else {
	  id = "CNV";
	  std::string padNumber = boost::lexical_cast<std::string>(++cnvid);
	  padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
	  id += padNumber;
	}
	bcf_update_id(hdr, rec, id.c_str());
	std::string svtype = "CNV";
	std::string alleles = "N,<" + svtype + ">";
	bcf_update_alleles_str(hdr, rec, alleles.c_str());
      
	// Add INFO fields
	if ((cnvs[i].srleft > 0) && (cnvs[i].srright > 0)) bcf_update_info_flag(hdr, rec, "PRECISE", NULL, 1);
	else bcf_update_info_flag(hdr, rec, "IMPRECISE", NULL, 1);

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
	float uniqf = cnvs[i].uniqfrac;
	bcf_update_info_float(hdr, rec, "UNIQ", &uniqf, 1);
	int32_t srltmp = cnvs[i].srleft;
	bcf_update_info_int32(hdr, rec, "SRL", &srltmp, 1);
	int32_t srrtmp = cnvs[i].srright;
	bcf_update_info_int32(hdr, rec, "SRR", &srrtmp, 1);

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
