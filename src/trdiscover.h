#ifndef TRDISCOVER_H
#define TRDISCOVER_H

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>

#include <htslib/sam.h>

#include "tags.h"
#include "util.h"
#include "split.h"
#include "gotoh.h"
#include "tandem.h"

namespace torali
{

  // Observed TR
  struct TrObs {
    int32_t net;
    int32_t mapq;
    int32_t qa;
    int32_t qb;
    uint32_t file;
    std::size_t id;
    std::string seq;
  };

  // Tandem repeat allele
  struct TrAllele {
    double net;
    uint32_t support;
    int32_t status;   // 0: not tried, 1: no consensus, 2: no anchor, 3: no gap, 4: invalid split, 5: low identity, 8: gap length, 9: emitted
    int32_t gap;
    int32_t consNet;
    int32_t consLen;
    double ident;
    std::vector<uint32_t> obs;

    TrAllele(double const n, uint32_t const s) : net(n), support(s), status(0), gap(0), consNet(0), consLen(0), ident(0) {}
  };

  // Tandem repeat locus
  struct TrLocus {
    int32_t ts;
    int32_t te;
    int32_t ws;
    int32_t we;
    uint32_t nspan;
    bool capped;
    std::vector<int32_t> members;
    std::vector<int32_t> nets;
    std::vector<TrObs> obs;
    std::vector<TrAllele> alleles;

    TrLocus(int32_t const s, int32_t const e) : ts(s), te(e), ws(s), we(e), nspan(0), capped(false) {}
  };

  // Query offset
  inline int32_t
  _findSeqBpStrict(bam1_t* rec, int32_t const pos, int32_t const minOp, int32_t const guard) {
    int32_t rp = rec->core.pos;
    int32_t sp = 0;
    int32_t hit = -1;
    const uint32_t* cigar = bam_get_cigar(rec);
    for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
      int32_t op = bam_cigar_op(cigar[i]);
      int32_t len = bam_cigar_oplen(cigar[i]);
      if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF)) {
	if ((hit < 0) && (rp <= pos) && (pos < rp + len)) hit = sp + (pos - rp);
	rp += len;
	sp += len;
      } else if ((op == BAM_CDEL) || (op == BAM_CREF_SKIP)) {
	if ((rp <= pos) && (pos < rp + len)) return -1;
	if ((len >= minOp) && ((std::abs(pos - rp) < guard) || (std::abs(pos - (rp + len)) < guard))) return -1;
	rp += len;
      } else if (op == BAM_CINS) {
	if ((len >= minOp) && (std::abs(pos - rp) < guard)) return -1;
	sp += len;
      } else if (op == BAM_CSOFT_CLIP) {
	sp += len;
      }
    }
    return hit;
  }

  // Tandem repeat period
  template<typename TConfig>
  inline int32_t
  _trCandidatePeriod(TConfig const& c, char const* seq, int32_t const seqlen, StructuralVariantRecord const& sv) {
    int32_t const a = sv.svStart + 1;
    if ((sv.svt == 2) && (sv.svEnd - sv.svStart <= 2000)) {
      int32_t p = _trWindowPeriod(seq, seqlen, sv.svStart, sv.svEnd, c.trMinFrac);
      if (p > 0) return p;
    }
    static const int32_t win[3] = {50, 100, 200};
    for(int32_t k = 0; k < 3; ++k) {
      int32_t p = _trWindowPeriod(seq, seqlen, std::max(0, a - win[k]), a, c.trMinFrac);
      if (p > 0) return p;
      p = _trWindowPeriod(seq, seqlen, a, std::min(seqlen, a + win[k]), c.trMinFrac);
      if (p > 0) return p;
    }
    return 0;
  }

  // Tandem repeat loci
  template<typename TConfig>
  inline void
  _trBuildLoci(TConfig const& c, char const* seq, int32_t const seqlen, int32_t const refIndex, std::vector<StructuralVariantRecord> const& svs, std::vector<TrLocus>& loci) {
    typedef std::pair<int32_t, int32_t> TTract;
    typedef std::pair<TTract, int32_t> TTractSv;
    std::vector<TTractSv> tr;
    for(uint32_t i = 0; i < svs.size(); ++i) {
      if ((svs[i].chr != refIndex) || ((svs[i].svt != 2) && (svs[i].svt != 4))) continue;
      if ((svs[i].svt == 2) && (svs[i].svEnd - svs[i].svStart > c.indelsize)) continue;
      if ((svs[i].svt == 4) && (svs[i].insLen > c.indelsize)) continue;
      int32_t a = svs[i].svStart + 1;
      int32_t p = _trCandidatePeriod(c, seq, seqlen, svs[i]);
      if (p <= 0) continue;
      int32_t ts = -1;
      int32_t te = -1;
      if ((!_trTract(seq, seqlen, a, p, ts, te)) && (!_trTract(seq, seqlen, a - p, p, ts, te))) continue;
      int32_t svE = (svs[i].svt == 4) ? svs[i].svStart : svs[i].svEnd;
      if ((svs[i].svStart < ts - c.trFlank) || (svE > te + c.trFlank)) continue;
      tr.push_back(std::make_pair(std::make_pair(ts, te), (int32_t) i));
    }
    std::sort(tr.begin(), tr.end());

    // Merge overlapping repeats
    for(uint32_t k = 0; k < tr.size(); ++k) {
      int32_t ts = tr[k].first.first;
      int32_t te = tr[k].first.second;
      if ((!loci.empty()) && (ts <= loci.back().te + c.trFlank) && (std::max(te, loci.back().te) - loci.back().ts <= 10000)) loci.back().te = std::max(loci.back().te, te);
      else loci.push_back(TrLocus(ts, te));
      loci.back().members.push_back(tr[k].second);
    }

    for(uint32_t k = 0; k < loci.size(); ++k) {
      loci[k].ws = loci[k].ts;
      loci[k].we = loci[k].te;
      for(uint32_t m = 0; m < loci[k].members.size(); ++m) {
	StructuralVariantRecord const& sv = svs[loci[k].members[m]];
	int32_t svE = (sv.svt == 4) ? sv.svStart : sv.svEnd;
	loci[k].ws = std::min(loci[k].ws, sv.svStart);
	loci[k].we = std::max(loci[k].we, svE);
      }
      loci[k].ws = std::max(0, loci[k].ws - c.trFlank);
      loci[k].we = std::min(seqlen, loci[k].we + c.trFlank);
    }
  }

  // Repeat length
  template<typename TConfig>
  inline void
  _trScanRead(TConfig const& c, bam1_t* rec, uint32_t const file_c, std::vector<TrLocus>& loci, uint32_t& li) {
    if (rec->core.qual < c.minMapQual) return;
    int32_t const rpos = rec->core.pos;
    int32_t const rend = bam_endpos(rec);
    while ((li < loci.size()) && (loci[li].ws < rpos)) ++li;
    std::string sequence;
    std::vector<std::pair<int32_t, uint32_t> > bps;
    for(uint32_t k = li; ((k < loci.size()) && (loci[k].ws <= rend)); ++k) {
      TrLocus& l = loci[k];
      if (l.we > rend) continue;
      if (l.nspan >= c.maxGenoReadCount) {
	l.capped = true;
	continue;
      }
      bps.push_back(std::make_pair(l.ws, 2 * k));
      bps.push_back(std::make_pair(l.we, 2 * k + 1));
    }
    if (bps.empty()) return;
    std::sort(bps.begin(), bps.end());
    std::vector<int32_t> qoff(bps.size(), -1);
    {
      int32_t const minOp = c.minRefSep;
      int32_t const guard = 10;
      int32_t rp = rec->core.pos;
      int32_t sp = 0;
      uint32_t bi = 0;
      const uint32_t* cigar = bam_get_cigar(rec);
      auto invalidate = [&](int32_t lo, int32_t hi) {
	std::vector<std::pair<int32_t, uint32_t> >::iterator it = std::lower_bound(bps.begin(), bps.end(), std::make_pair(lo, (uint32_t) 0));
	for(; (it != bps.end()) && (it->first <= hi); ++it) qoff[it - bps.begin()] = -2;
      };
      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	int32_t op = bam_cigar_op(cigar[i]);
	int32_t len = bam_cigar_oplen(cigar[i]);
	if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF)) {
	  while ((bi < bps.size()) && (bps[bi].first < rp + len)) {
	    if ((bps[bi].first >= rp) && (qoff[bi] == -1)) qoff[bi] = sp + (bps[bi].first - rp);
	    ++bi;
	  }
	  rp += len;
	  sp += len;
	} else if ((op == BAM_CDEL) || (op == BAM_CREF_SKIP)) {
	  while ((bi < bps.size()) && (bps[bi].first < rp + len)) {
	    if (bps[bi].first >= rp) qoff[bi] = -2;   // inside a deletion
	    ++bi;
	  }
	  if (len >= minOp) { invalidate(rp - guard + 1, rp + guard - 1); invalidate(rp + len - guard + 1, rp + len + guard - 1); }
	  rp += len;
	} else if (op == BAM_CINS) {
	  if (len >= minOp) invalidate(rp - guard + 1, rp + guard - 1);
	  sp += len;
	} else if (op == BAM_CSOFT_CLIP) {
	  sp += len;
	}
      }
    }
    std::vector<int32_t> qa_of(loci.size() - li, -1);
    std::vector<int32_t> qb_of(loci.size() - li, -1);
    for(uint32_t b = 0; b < bps.size(); ++b) {
      uint32_t k = bps[b].second / 2;
      if (bps[b].second % 2 == 0) qa_of[k - li] = (qoff[b] >= 0) ? qoff[b] : -1;
      else qb_of[k - li] = (qoff[b] >= 0) ? qoff[b] : -1;
    }
    for(uint32_t k = li; ((k < loci.size()) && (loci[k].ws <= rend)); ++k) {
      TrLocus& l = loci[k];
      if ((l.we > rend) || (l.capped)) continue;
      int32_t qa = qa_of[k - li];
      if (qa < 0) continue;
      int32_t qb = qb_of[k - li];
      if (qb < 0) continue;
      ++l.nspan;
      int32_t net = (qb - qa) - (l.we - l.ws);
      l.nets.push_back(net);
      if ((std::abs(net) >= (int32_t) c.minRefSep) && (l.obs.size() < 4 * c.maxReadPerSV)) {
	if (sequence.empty()) {
	  sequence.resize(rec->core.l_qseq);
	  const uint8_t* seqptr = bam_get_seq(rec);
	  for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	}
	int32_t sPos = std::max(0, qa - c.minConsWindow);
	int32_t ePos = std::min((int32_t) sequence.size(), qb + c.minConsWindow);
	TrObs o;
	o.net = net;
	o.mapq = rec->core.qual;
	o.qa = qa - sPos;
	o.qb = qb - sPos;
	o.file = file_c;
	o.id = hash_lr(rec);
	o.seq = sequence.substr(sPos, ePos - sPos);
	l.obs.push_back(std::move(o));
      }
    }
  }

  // Read noise
  inline double
  _trNoiseScale(int32_t const tractLen) {
    if (tractLen < 300) return 8.0;
    if (tractLen < 1000) return 10.0;
    return 14.0;
  }

  inline double
  _trMedian(std::vector<int32_t>& v) {
    std::sort(v.begin(), v.end());
    std::size_t n = v.size();
    if (n % 2) return v[n / 2];
    return ((double) v[n / 2 - 1] + (double) v[n / 2]) / 2.0;
  }

  // Repeat alleles of a locus
  template<typename TConfig>
  inline void
  _trCallAlleles(TConfig const& c, TrLocus& l) {
    l.alleles.clear();
    if (l.nets.size() < c.minCliqueSize) return;
    double const b = _trNoiseScale(l.te - l.ts);
    double const minSupp = c.trAlleleFrac * (double) l.nspan;
    int32_t lo = *std::min_element(l.nets.begin(), l.nets.end()) - 3 * (int32_t) b;
    int32_t hi = *std::max_element(l.nets.begin(), l.nets.end()) + 3 * (int32_t) b;

    // Kernel density
    std::vector<double> dens(hi - lo + 1, 0);
    for(int32_t x = lo; x <= hi; ++x) {
      for(uint32_t r = 0; r < l.nets.size(); ++r) dens[x - lo] += std::exp(-std::abs((double) x - (double) l.nets[r]) / b);
    }

    std::vector<int32_t> modes;
    for(int32_t x = lo + 1; x < hi; ++x) {
      if ((dens[x - lo] >= dens[x - lo - 1]) && (dens[x - lo] > dens[x - lo + 1])) {
	uint32_t mass = 0;
	for(uint32_t r = 0; r < l.nets.size(); ++r) {
	  if (std::abs((double) l.nets[r] - (double) x) <= 2 * b) ++mass;
	}
	if ((mass >= c.minCliqueSize) && ((double) mass >= minSupp)) modes.push_back(x);
      }
    }

    typedef std::pair<double, uint32_t> TAllele;
    std::vector<TAllele> alle;
    for(uint32_t m = 0; m < modes.size(); ++m) {
      std::vector<int32_t> mem;
      for(uint32_t r = 0; r < l.nets.size(); ++r) {
	int32_t d = std::abs(l.nets[r] - modes[m]);
	if ((double) d > 5 * b) continue;
	bool nearest = true;
	for(uint32_t k = 0; k < modes.size(); ++k) {
	  if (std::abs(l.nets[r] - modes[k]) < d) {
	    nearest = false;
	    break;
	  }
	}
	if (nearest) mem.push_back(l.nets[r]);
      }
      if ((mem.size() >= c.minCliqueSize) && ((double) mem.size() >= minSupp)) alle.push_back(std::make_pair(_trMedian(mem), (uint32_t) mem.size()));
    }
    std::sort(alle.begin(), alle.end());

    // Merge alleles
    std::vector<TAllele> merged;
    for(uint32_t a = 0; a < alle.size(); ++a) {
      if (!merged.empty()) {
	uint32_t nmin = std::min(alle[a].second, merged.back().second);
	double resol = std::max(10.0, std::max(0.03 * std::abs((double) (l.te - l.ts) + alle[a].first), 2.83 * b / std::sqrt(std::max(1.0, (double) nmin / 2.0))));
	if (std::abs(alle[a].first - merged.back().first) <= resol) {
	  uint32_t tot = merged.back().second + alle[a].second;
	  merged.back().first = (merged.back().first * merged.back().second + alle[a].first * alle[a].second) / (double) tot;
	  merged.back().second = tot;
	  continue;
	}
      }
      merged.push_back(alle[a]);
    }
    for(uint32_t a = 0; a < merged.size(); ++a) {
      if (std::abs(merged[a].first) >= (double) c.minRefSep) l.alleles.push_back(TrAllele(merged[a].first, merged[a].second));
    }

    // Read support per allele
    for(uint32_t o = 0; o < l.obs.size(); ++o) {
      int32_t best = -1;
      double bestD = 5 * b;
      for(uint32_t a = 0; a < l.alleles.size(); ++a) {
	double d = std::abs((double) l.obs[o].net - l.alleles[a].net);
	if (d <= bestD) {
	  if ((best < 0) || (d < bestD)) {
	    best = a;
	    bestD = d;
	  }
	}
      }
      if (best >= 0) l.alleles[best].obs.push_back(o);
    }
  }

  // Repeat Allele Discovery
  template<typename TConfig>
  inline void
  _trCallLoci(TConfig const& c, std::vector<TrLocus>& loci) {
    std::vector<int32_t> span;
    for(uint32_t k = 0; k < loci.size(); ++k) span.push_back(loci[k].nspan);
    double spanMed = (span.empty()) ? 0 : _trMedian(span);
    for(uint32_t k = 0; k < loci.size(); ++k) {
      TrLocus& l = loci[k];
      if ((l.capped) || ((double) l.nspan < 0.5 * spanMed) || ((double) l.nspan > 2.5 * spanMed)) continue;
      _trCallAlleles(c, l);
    }
  }

  // Consensus
  template<typename TConfig>
  inline bool
  _trAlleleConsensus(TConfig const& c, char const* seq, int32_t const seqlen, TrLocus const& l, TrAllele const& a, std::string& cons) {
    if (a.obs.size() < 2) return false;
    std::vector<uint32_t> idx(a.obs);
    std::sort(idx.begin(), idx.end(), [&](uint32_t x, uint32_t y) {
      double dx = std::abs((double) l.obs[x].net - a.net);
      double dy = std::abs((double) l.obs[y].net - a.net);
      if (dx != dy) return dx < dy;
      int32_t bx = std::min(l.obs[x].qa, (int32_t) l.obs[x].seq.size() - l.obs[x].qb);
      int32_t by = std::min(l.obs[y].qa, (int32_t) l.obs[y].seq.size() - l.obs[y].qb);
      if (bx != by) return bx > by;
      return x < y; });
    std::vector<std::string> sps;
    for(uint32_t i = 0; (i < idx.size()) && ((int32_t) i < (int32_t) c.maxReadPerSV); ++i) sps.push_back(l.obs[idx[i]].seq);
    std::string prefix = boost::to_upper_copy(std::string(seq + std::max(l.ws - c.minConsWindow, 0), seq + l.ws));
    std::string suffix = boost::to_upper_copy(std::string(seq + l.we, seq + std::min(seqlen, l.we + c.minConsWindow)));
    cons.clear();
    msaWfa(c, sps, cons, prefix, suffix);
    return (!cons.empty());
  }

  // Anchor reference flank
  inline bool
  _trAnchorFlank(std::string const& flank, std::string const& cons, int32_t& cStart, int32_t& cEnd, double& ident) {
    if ((flank.size() < 50) || (cons.size() < flank.size())) return false;
    EdlibAlignResult cig = edlibAlign(flank.c_str(), flank.size(), cons.c_str(), cons.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    bool ok = ((cig.status == EDLIB_STATUS_OK) && (cig.numLocations > 0));
    if (ok) {
      cEnd = infixEnd(cig) + 1;
      cStart = infixStart(cig);
      ident = 1.0 - (double) cig.editDistance / (double) flank.size();
    }
    edlibFreeAlignResult(cig);
    return ok;
  }

  // Align consensus
  template<typename TConfig>
  inline int32_t
  _trAlignAllele(TConfig const& c, char const* seq, int32_t const seqlen, int32_t const refIndex, TrLocus const& l, TrAllele& a, StructuralVariantRecord& sv, int32_t& gap) {
    std::string const& cons = sv.consensus;
    sv.chr = refIndex;
    sv.chr2 = refIndex;
    sv.svt = (a.net < 0) ? 2 : 4;

    // Flank anchors
    int32_t const anchorLen = std::max(100, c.minConsWindow / 2);
    std::string prefix = boost::to_upper_copy(std::string(seq + std::max(l.ws - anchorLen, 0), seq + l.ws));
    std::string suffix = boost::to_upper_copy(std::string(seq + l.we, seq + std::min(seqlen, l.we + anchorLen)));
    int32_t pS = 0, pE = 0, sS = 0, sE = 0;
    double pid = 0, sid = 0;
    a.consLen = (int32_t) cons.size();
    if (!_trAnchorFlank(prefix, cons, pS, pE, pid)) return 2;
    if (!_trAnchorFlank(suffix, cons, sS, sE, sid)) return 2;
    a.ident = std::min(pid, sid);
    if (sS <= pE) return 4;
    sv.srAlignQuality = std::min(pid, sid);
    if (sv.srAlignQuality < c.flankQuality) return 5;

    // Window alignment
    std::string refWin = boost::to_upper_copy(std::string(seq + l.ws, seq + l.we));
    std::string altWin = cons.substr(pE, sS - pE);
    double tol = std::max(2.0 * _trNoiseScale(l.te - l.ts), 0.05 * std::abs(a.net));
    int32_t net = (int32_t) altWin.size() - (int32_t) refWin.size();
    a.consNet = net;
    if (std::abs((double) net - a.net) > tol) return 8;
    typedef boost::multi_array<char, 2> TAlign;
    TAlign align;
    AlignConfig<false, false> global;
    DnaScore<int> sc(5, -4, -5000, -1);
    gotoh(altWin, refWin, align, global, sc);

    // Longest gap (consensus vs reference)
    AlignDescriptor ad;
    int32_t refIdx = 0;
    int32_t varIdx = 0;
    int32_t gapStartRefIdx = 0;
    int32_t gapStartVarIdx = 0;
    bool inGap = false;
    for(int32_t j = 0; j < (int32_t) align.shape()[1]; ++j) {
      if (align[0][j] != '-') ++varIdx;
      if (align[1][j] != '-') ++refIdx;
      if (((align[0][j] == '-') || (align[1][j] == '-')) && (refIdx > 0) && (varIdx > 0)) {
	if (!inGap) {
	  gapStartVarIdx = (align[0][j] != '-') ? (varIdx - 1) : varIdx;
	  gapStartRefIdx = (align[1][j] != '-') ? (refIdx - 1) : refIdx;
	  inGap = true;
	}
      } else {
	if ((inGap) && (_checkSVGap((refIdx - gapStartRefIdx), (ad.rEnd - ad.rStart), (varIdx - gapStartVarIdx), (ad.cEnd - ad.cStart), sv.svt))) {
	  ad.rStart = gapStartRefIdx;
	  ad.rEnd = refIdx;
	  ad.cStart = gapStartVarIdx;
	  ad.cEnd = varIdx;
	}
	inGap = false;
      }
    }
    if (ad.rEnd <= ad.rStart) return 3;
    if (!_validSRAlignment(ad.cStart, ad.cEnd, ad.rStart, ad.rEnd, sv.svt)) return 4;
    gap = (sv.svt == 2) ? (ad.rEnd - ad.rStart) : (ad.cEnd - ad.cStart - 1);
    if (std::abs((double) gap - std::abs(a.net)) > tol) return 8;

    // Exact alleles
    std::string refVCF;
    std::string altVCF;
    int32_t cpos = 0;
    bool inSV = false;
    for(uint32_t j = 0; j < align.shape()[1]; ++j) {
      if (align[0][j] != '-') {
	++cpos;
	if (cpos == ad.cStart) inSV = true;
	else if (cpos == ad.cEnd) inSV = false;
      }
      if (inSV) {
	if (align[0][j] != '-') altVCF += align[0][j];
	if (align[1][j] != '-') refVCF += align[1][j];
      }
    }

    // Left-alignment
    {
      int32_t shift = 0;
      int32_t rS = ad.rStart, rE = ad.rEnd, cS = ad.cStart, cE = ad.cEnd;
      int32_t const refN = (int32_t) refWin.size();
      int32_t const altN = (int32_t) altWin.size();
      if (sv.svt == 4) {
	// Insertion
	while ((rS > 1) && (cS > 1) && (cE - 1 <= altN) && (refWin[rS - 2] == altWin[cS - 2]) && (altWin[cE - 2] == refWin[rS - 1])) {
	  --rS;
	  --rE;
	  --cS;
	  --cE;
	  ++shift;
	}
      } else {
	// Deletion
	while ((rS > 1) && (cS > 1) && (rE <= refN) && (refWin[rS - 2] == altWin[cS - 2]) && (refWin[rE - 2] == refWin[rS - 1])) {
	  --rS;
	  --rE;
	  --cS;
	  --cE;
	  ++shift;
	}
      }
      if (shift > 0) {
	ad.rStart = rS; ad.rEnd = rE; ad.cStart = cS; ad.cEnd = cE;
	if (sv.svt == 4) {
	  refVCF = refWin.substr(rS - 1, 1);
	  altVCF = altWin.substr(cS - 1, cE - cS + 1);
	} else {
	  refVCF = refWin.substr(rS - 1, rE - rS);
	  altVCF = altWin.substr(cS - 1, 1);
	}
      }
    }
    sv.alleles = _addAlleles(refVCF, altVCF);
    sv.precise = true;
    sv.svStart = l.ws + ad.rStart;
    sv.svEnd = l.ws + ad.rEnd;
    sv.insLen = ad.cEnd - ad.cStart - 1;
    sv.consBp = pE + ad.cStart;

    // Homology
    std::string svRefStr = prefix + refWin + suffix;
    ad.cStart += pE;
    ad.cEnd += pE;
    ad.rStart += (int32_t) prefix.size();
    ad.rEnd += (int32_t) prefix.size();
    _findHomology(cons, svRefStr, ad, sv.svt);
    sv.homLen = std::max(0, ad.homLeft + ad.homRight - 2);
    int32_t ci_wiggle = std::max(ad.homLeft, ad.homRight);
    sv.ciposlow = -ci_wiggle;
    sv.ciposhigh = ci_wiggle;
    sv.ciendlow = -ci_wiggle;
    sv.ciendhigh = ci_wiggle;
    return 0;
  }

  // Output TR alleles
  template<typename TConfig, typename TSRStore>
  inline void
  _trEmitLocus(TConfig const& c, char const* seq, int32_t const seqlen, int32_t const refIndex, TrLocus& l, TSRStore const& srStoreChr, std::vector<StructuralVariantRecord>& svs, std::vector<int32_t> const& chrSv, std::vector<StructuralVariantRecord>& out, int32_t& alleleCounter, uint32_t& nCalled, uint32_t& nEmitted, uint32_t& nReplaced) {
    if (l.alleles.empty()) return;
    if (l.we - l.ws > c.indelsize) return;
    std::vector<StructuralVariantRecord> rec;
    std::set<std::size_t> alleleReads;
    for(uint32_t a = 0; a < l.alleles.size(); ++a) {
      ++nCalled;
      StructuralVariantRecord sv;
      if (!_trAlleleConsensus(c, seq, seqlen, l, l.alleles[a], sv.consensus)) {
	l.alleles[a].status = 1;
	continue;
      }
      int32_t reason = _trAlignAllele(c, seq, seqlen, refIndex, l, l.alleles[a], sv, l.alleles[a].gap);
      if (reason) {
	l.alleles[a].status = reason;
	continue;
      }
      l.alleles[a].status = 9;
      double mq = 0;
      for(uint32_t o = 0; o < l.alleles[a].obs.size(); ++o) {
	mq += l.obs[l.alleles[a].obs[o]].mapq;
	alleleReads.insert(l.obs[l.alleles[a].obs[o]].id);
      }
      mq /= (double) l.alleles[a].obs.size();
      sv.srSupport = l.alleles[a].support;
      sv.peSupport = 0;
      sv.srMapQuality = (int32_t) mq;
      sv.mapq = (int32_t) (mq * l.alleles[a].support);
      sv.precise = true;
      sv.trStart = std::min(l.ts, sv.svStart);
      sv.trEnd = std::max(l.te, (sv.svt == 2) ? sv.svEnd : sv.svStart);
      sv.alleleid = alleleCounter;
      rec.push_back(sv);
    }
    if (rec.empty()) return;
    ++alleleCounter;
    for(uint32_t r = 0; r < rec.size(); ++r) {
      rec[r].nallele = rec.size();
      out.push_back(rec[r]);
      ++nEmitted;
    }

    // Candidates inside the repeat tract
    for(uint32_t i = 0; i < chrSv.size(); ++i) {
      StructuralVariantRecord& cand = svs[chrSv[i]];
      if (cand.svStart < l.ts - c.trFlank) continue;
      if (cand.svStart > l.te + c.trFlank) break;
      if (cand.srSupport == 0) continue;
      int32_t candE = (cand.svt == 4) ? cand.svStart : cand.svEnd;
      if (candE > l.te + c.trFlank) continue;
      int32_t candLen = (cand.svt == 4) ? cand.insLen : (cand.svEnd - cand.svStart);
      bool explained = false;
      for(uint32_t r = 0; (r < rec.size()) && (!explained); ++r) {
	if (rec[r].svt != cand.svt) continue;
	int32_t recLen = (rec[r].svt == 4) ? rec[r].insLen : (rec[r].svEnd - rec[r].svStart);
	int32_t lengthvar = std::max(15, (int32_t) (0.1 * std::min(candLen, recLen)));
	if (std::abs(candLen - recLen) < lengthvar) explained = true;
      }
      if (explained) {
	cand.srSupport = 0;
	cand.consensus.clear();
	++nReplaced;
      }
    }

    std::map<int32_t, std::set<std::size_t> > memberReads;
    for(uint32_t m = 0; m < l.members.size(); ++m) memberReads[l.members[m]] = std::set<std::size_t>();
    for(typename TSRStore::const_iterator it = srStoreChr.begin(); it != srStoreChr.end(); ++it) {
      for(uint32_t k = 0; k < it->second.size(); ++k) {
	typename std::map<int32_t, std::set<std::size_t> >::iterator mr = memberReads.find(it->second[k].svid);
	if (mr != memberReads.end()) mr->second.insert(it->first.second);
      }
    }
    for(typename std::map<int32_t, std::set<std::size_t> >::const_iterator mr = memberReads.begin(); mr != memberReads.end(); ++mr) {
      if (mr->second.empty()) continue;
      uint32_t hit = 0;
      for(std::set<std::size_t>::const_iterator ri = mr->second.begin(); ri != mr->second.end(); ++ri) {
	if (alleleReads.find(*ri) != alleleReads.end()) ++hit;
      }
      if ((2 * hit >= mr->second.size()) && (svs[mr->first].srSupport > 0)) {
	svs[mr->first].srSupport = 0;
	svs[mr->first].consensus.clear();
	++nReplaced;
      }
    }
  }

}

#endif
