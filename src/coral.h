#ifndef CORAL_H
#define CORAL_H

#include <limits>
#include <iomanip>
#include <cmath>

#include <boost/icl/split_interval_map.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/math/special_functions/round.hpp>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "bed.h"
#include "scan.h"
#include "gcbias.h"
#include "cnv.h"
#include "version.h"

namespace torali
{

  struct CountDNAConfig {
    bool basecov;
    bool somatic;
    bool adaptive;
    bool hasStatsFile;
    bool hasScanFile;
    bool noScanWindowSelection;
    bool regionalGc;
    bool hasSegFile;
    bool hasGenoFile;
    bool hasExcludeFile;
    bool hasSnvSites;
    uint32_t nchr;
    uint32_t minSnvDepth;
    uint8_t minBaseQual;
    uint32_t minClip;
    uint32_t minRefSep;
    uint32_t minBpSupport;
    float cnMergeTol;
    float penalty;
    uint32_t meanisize;
    uint32_t window_size;
    uint32_t window_offset;
    uint32_t scanWindow;
    uint32_t minChrLen;
    uint32_t minCnvSize;
    uint32_t targetReads;
    uint32_t minSegWin;
    double targetExpCov;
    uint16_t minQual;
    uint16_t mapqUniq;
    uint16_t mad;
    bool ploidyAuto;
    float ploidy;
    float ctrlPloidy;
    float expectedCN;
    float purity;
    float exclgc;
    float uniqueToTotalCovRatio;
    float fracWindow;
    float fragmentUnique;
    int32_t cnvMinQual;
    float cnvDelConfirm;
    float cnvDelRatio;
    std::string sampleName;
    boost::filesystem::path segfile;
    boost::filesystem::path genofile;
    boost::filesystem::path outfile;
    boost::filesystem::path covfile;
    boost::filesystem::path genome;
    boost::filesystem::path statsFile;
    boost::filesystem::path bamFile;
    boost::filesystem::path scanFile;
    boost::filesystem::path exclude;
    boost::filesystem::path snvSites;
    std::set<int32_t> refIdx;
  };

  struct CovWin {
    uint32_t start;
    uint32_t end;
    uint32_t winlen;
    double covsum;
    double expcov;
    double ucov;
    double tcov;
    double aall;
    double eall;
    bool valid;

    CovWin(uint32_t const s, uint32_t const e, uint32_t const w, double const cs, double const ec, double const uc, double const tc, double const aa, double const ea, bool vld) : start(s), end(e), winlen(w), covsum(cs), expcov(ec), ucov(uc), tcov(tc), aall(aa), eall(ea), valid(vld) {}
  };

  struct CountDNAConfigLib {
    uint16_t madCutoff;
    uint16_t madNormalCutoff;
    boost::filesystem::path genome;
    std::vector<boost::filesystem::path> files;
  };
  
  template<typename TConfig>
  inline int32_t
  bamCount(TConfig const& c, LibraryInfo const& li, std::vector<GcBias> const& gcbias, std::pair<uint32_t, uint32_t> const& gcbound, std::vector<double> const& regcorr, uint32_t const regWin) {
    // Load bam file
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Parse BAM file
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Count fragments" << std::endl;

    // Open output files
    boost::iostreams::filtering_ostream dataOut;
    if (!c.covfile.empty()) {
      dataOut.push(boost::iostreams::gzip_compressor());
      dataOut.push(boost::iostreams::file_sink(c.covfile.string(), std::ios_base::out | std::ios_base::binary));
      dataOut << "chr\tstart\tend\t" << c.sampleName << "_uniqfrac\t" << c.sampleName << "_logR\t" << c.sampleName << "_CN" << std::endl;
    }

    // CNVs
    std::vector<CNV> cnvs;
    if (c.hasGenoFile) parseVcfCNV(c, hdr, cnvs);

    // Iterate chromosomes
    faidx_t* faiRef = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      if ((!c.hasGenoFile) && (chrNoData(c, refIndex, idx))) continue;

      // Haploid chromosome?
      float chrCtrlPloidy = c.ctrlPloidy;
      float chrPloidy = c.ploidy;
      if (c.refIdx.find(refIndex) != c.refIdx.end()) {
	chrCtrlPloidy = c.ctrlPloidy - 1;
	chrPloidy = c.ploidy - 1;
      }
      
      // Reference sequence
      std::string tname(hdr->target_name[refIndex]);
      int32_t seqlen = faidx_seq_len(faiRef, tname.c_str());
      if (seqlen == - 1) continue;
      else seqlen = -1;
      char* ref = faidx_fetch_seq(faiRef, tname.c_str(), 0, faidx_seq_len(faiRef, tname.c_str()), &seqlen);
      if (ref == NULL) continue;

      // Get GC content
      std::vector<uint16_t> uniqContent(hdr->target_len[refIndex], 0);
      std::vector<uint16_t> gcContent(hdr->target_len[refIndex], 0);
      {
	// GC map
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet gcref(hdr->target_len[refIndex], false);
	for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	  if ((ref[i] == 'c') || (ref[i] == 'C') || (ref[i] == 'g') || (ref[i] == 'G')) gcref[i] = 1;
	}

	// Sum across fragment
	int32_t halfwin = (int32_t) (c.meanisize / 2);
	int32_t gcsum = 0;
	for(int32_t pos = halfwin; pos < (int32_t) hdr->target_len[refIndex] - halfwin; ++pos) {
	  if (pos == halfwin) {
	    for(int32_t i = pos - halfwin; i<=pos+halfwin; ++i) gcsum += gcref[i];
	  } else {
	    gcsum -= gcref[pos - halfwin - 1];
	    gcsum += gcref[pos + halfwin];
	  }
	  gcContent[pos] = gcsum;
	}
      }

      // Broad window for replication wave
      std::vector<float> tileFac;
      if ((!regcorr.empty()) && (regWin > 0)) {
	uint32_t ntile = hdr->target_len[refIndex] / regWin + 1;
	tileFac.resize(ntile, 1.0);
	for(uint32_t t = 0; t < ntile; ++t) {
	  uint32_t s = t * regWin;
	  uint32_t e = std::min((uint32_t) hdr->target_len[refIndex], s + regWin);
	  double gcnum = 0;
	  uint32_t winlen = 0;
	  for(uint32_t pos = s; pos < e; ++pos) {
	    if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second)) { gcnum += gcContent[pos]; ++winlen; }
	  }
	  if (winlen > 0) tileFac[t] = (float) regCorrFactor(regcorr, (gcnum / (double) winlen) / (double) c.meanisize);
	}
      }

      // Coverage track
      typedef uint16_t TCount;
      uint32_t maxCoverage = std::numeric_limits<TCount>::max();
      typedef std::vector<TCount> TCoverage;
      TCoverage cov(hdr->target_len[refIndex], 0);
      TCoverage covUniq(hdr->target_len[refIndex], 0);
      TCoverage covAll(hdr->target_len[refIndex], 0);
      TCoverage covTot;
      if (!c.basecov) covTot.resize(hdr->target_len[refIndex], 0);
      TCoverage& covMap = (!c.basecov) ? covTot : cov;

      // Split-read breakpoints
      std::vector<int32_t> clips;
      {
	// Mate map
	typedef boost::unordered_map<std::size_t, bool> TMateMap;
	TMateMap mateMap;

	// Count reads
	hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	int32_t lastAlignedPos = 0;
	std::set<std::size_t> lastAlignedPosReads;
	while (sam_itr_next(samfile, iter, rec) >= 0) {
	  if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
	  if ((rec->core.flag & BAM_FPAIRED) && ((rec->core.flag & BAM_FMUNMAP) || (rec->core.tid != rec->core.mtid))) continue;

	  // Base coverage
	  addBaseCoverage3(rec, covAll, covMap, covUniq, c.minQual, c.mapqUniq, hdr->target_len[refIndex], maxCoverage);

	  // Filter by quality
	  if (rec->core.qual < c.minQual) continue;

	  // Collect split-read breakpoints
	  if (rec->core.qual >= c.mapqUniq) addSplitReadBreakpoints(rec, c.minClip, c.minRefSep, hdr->target_len[refIndex], clips);

	  // Base-level counting done
	  if (c.basecov) continue;

	  // Fragment coverage
	  int32_t midPoint = rec->core.pos + halfAlignmentLength(rec);
	  if (rec->core.flag & BAM_FPAIRED) {
	    std::size_t seed = hash_sr(rec);
	    
	    // Clean-up the read store for identical alignment positions
	    if (rec->core.pos > lastAlignedPos) {
	      lastAlignedPosReads.clear();
	      lastAlignedPos = rec->core.pos;
	    }
	    
	    if (_firstPairObs(rec, lastAlignedPosReads)) {
	      // First read
	      lastAlignedPosReads.insert(seed);
	      std::size_t hv = hash_pair(rec);
	      mateMap[hv] = true;
	      continue;
	    } else {
	      // Second read
	      std::size_t hv = hash_pair_mate(rec);
	      auto itMM = mateMap.find(hv);
	      if ((itMM == mateMap.end()) || (!itMM->second)) continue; // Mate discarded
	      mateMap.erase(itMM);
	    }
	    
	    // update midpoint
	    int32_t isize = (rec->core.pos + alignmentLength(rec)) - rec->core.mpos;
	    if ((li.minNormalISize < isize) && (isize < li.maxNormalISize)) midPoint = rec->core.mpos + (int32_t) (isize/2);
	  }
	  
	  // Count fragment
	  if ((midPoint >= 0) && (midPoint < (int32_t) hdr->target_len[refIndex]) && (cov[midPoint] < maxCoverage - 1)) ++cov[midPoint];
	}
	// Clean-up
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }

      // Callable positions: unique if mostly high-MAPQ reads
      for(uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
	bool u;
	if (covMap[pos] == 0) u = ((ref[pos] != 'N') && (ref[pos] != 'n'));
	else u = (2 * (uint32_t) covUniq[pos] >= (uint32_t) covMap[pos]);
	uniqContent[pos] = (u ? (uint16_t) c.meanisize : 0);
      }
      // hom-del or unmappable?
      uint32_t maxHomDel = 1000000;
      uint32_t rstart = 0;
      while (rstart < hdr->target_len[refIndex]) {
	if (covMap[rstart] == 0) {
	  uint32_t rend = rstart;
	  while ((rend < hdr->target_len[refIndex]) && (covMap[rend] == 0)) ++rend;
	  bool leftOK = (rstart > 0) && (uniqContent[rstart - 1] > 0);
	  bool rightOK = (rend < hdr->target_len[refIndex]) && (uniqContent[rend] > 0);
	  if ((!leftOK) || (!rightOK) || (rend - rstart > maxHomDel)) {
	    for(uint32_t k = rstart; k < rend; ++k) uniqContent[k] = 0;
	  }
	  rstart = rend;
	} else ++rstart;
      }
      // Genome-wide read-depth windows
      DepthTrack const dt = uniqueTrack();
      std::vector<CovWin> wins;
      if (c.adaptive) {
	double covsum = 0;
	double expraw = 0;
	double expcor = 0;
	double ucov = 0;
	double tcov = 0;
	double aall = 0;
	double eall = 0;
	uint32_t winlen = 0;
	uint32_t start = 0;
	for(uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
	  ucov += covUniq[pos];
	  tcov += covMap[pos];
	  if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second)) {
	    aall += covAll[pos];
	    eall += gcbias[gcContent[pos]].coverageTotal;
	  }
	  if (_posUsed(dt, c, gcContent, uniqContent, gcbound, pos)) {
	    double e1 = _expCov(dt, gcbias, gcContent[pos]);
	    covsum += cov[pos];
	    expraw += e1;
	    expcor += e1 * (tileFac.empty() ? 1.0 : (double) tileFac[pos / regWin]);
	    ++winlen;
	    if (expraw >= c.targetExpCov) {
	      wins.push_back(CovWin(start, pos + 1, winlen, covsum, expcor, ucov, tcov, aall, eall, true));
	      covsum = 0;
	      expraw = 0;
	      expcor = 0;
	      ucov = 0;
	      tcov = 0;
	      aall = 0;
	      eall = 0;
	      winlen = 0;
	      start = pos + 1;
	    }
	  }
	}
      } else {
	// Fixed windows (non-overlapping)
	for(uint32_t start = 0; start < hdr->target_len[refIndex]; start = start + c.window_offset) {
	  if (start + c.window_size < hdr->target_len[refIndex]) {
	    double covsum = 0;
	    double expcov = 0;
	    double ucov = 0;
	    double tcov = 0;
	    double aall = 0;
	    double eall = 0;
	    uint32_t winlen = 0;
	    for(uint32_t pos = start; pos < start + c.window_size; ++pos) {
	      ucov += covUniq[pos];
	      tcov += covMap[pos];
	      if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second)) {
		aall += covAll[pos];
		eall += gcbias[gcContent[pos]].coverageTotal;
	      }
	      if (_posUsed(dt, c, gcContent, uniqContent, gcbound, pos)) {
		covsum += cov[pos];
		expcov += _expCov(dt, gcbias, gcContent[pos]) * (tileFac.empty() ? 1.0 : (double) tileFac[pos / regWin]);
		++winlen;
	      }
	    }
	    bool valid = (winlen >= c.fracWindow * c.window_size);
	    wins.push_back(CovWin(start, start + c.window_size, winlen, covsum, expcov, ucov, tcov, aall, eall, valid));
	  }
	}
      }

      // Separate true hom. dels from unmappable
      uint32_t nw = wins.size();
      std::vector<bool> naFlag(nw, false);
      std::vector<bool> suspect(nw, false);
      std::vector<bool> strong(nw, false);
      double lowFrac = 0.1;            // close to CN0
      double flankFrac = 0.5;          // clear no CN0
      uint32_t maxHomDelWin = 1000000; // max size for hom. DEL
      for(uint32_t i = 0; i < nw; ++i) {
	if ((!wins[i].valid) || (wins[i].expcov <= 0)) {
	  naFlag[i] = true;
	  continue;
	}
	double r = wins[i].covsum / wins[i].expcov;
	suspect[i] = (r < lowFrac);
	strong[i] = (r >= flankFrac);
      }
      for(uint32_t i = 0; i < nw; ) {
	if (naFlag[i] || (!suspect[i])) {
	  ++i;
	  continue;
	}
	uint32_t a = i;
	uint32_t b = i;
	while ((b + 1 < nw) && (!naFlag[b+1]) && suspect[b+1]) ++b;
	uint32_t runBp = wins[b].end - wins[a].start;
	bool leftStrong = (a > 0) && (!naFlag[a-1]) && strong[a-1];
	bool rightStrong = (b + 1 < nw) && (!naFlag[b+1]) && strong[b+1];
	bool keepDel = leftStrong && rightStrong && (runBp <= maxHomDelWin);
	if (!keepDel) {
	  for(uint32_t k = a; k <= b; ++k) naFlag[k] = true;
	}
	i = b + 1;
      }

      // Total-depth deletion signal
      double delRatio = c.cnvDelRatio;

      // Flag non-unique windows
      bool uniqGate = c.basecov;
      if (uniqGate) {
	for(uint32_t i = 0; i < nw; ++i) {
	  if (naFlag[i]) continue;
	  double rTot = (wins[i].eall > 0) ? (wins[i].aall / wins[i].eall) : 1.0;
	  if (rTot < delRatio) continue;
	  if ((wins[i].tcov > 0) && (wins[i].ucov <= c.uniqueToTotalCovRatio * wins[i].tcov)) naFlag[i] = true;
	}
      }

      // Flag low callable windows
      if ((c.adaptive) && (nw > 4)) {
	std::vector<int32_t> spans(nw);
	for(uint32_t i = 0; i < nw; ++i) spans[i] = (int32_t) (wins[i].end - wins[i].start);
	std::vector<int32_t> tmp(spans);
	std::sort(tmp.begin(), tmp.end());
	int32_t medSpan = tmp[tmp.size() / 2];
	for(uint32_t i = 0; i < nw; ++i) tmp[i] = std::abs(spans[i] - medSpan);
	std::sort(tmp.begin(), tmp.end());
	int32_t madSpan = tmp[tmp.size() / 2];
	int32_t maxSpan = std::max(medSpan + (int32_t) c.mad * madSpan, 2 * medSpan);
	for(uint32_t i = 0; i < nw; ++i) {
	  if (spans[i] > maxSpan) {
	    double rTot = (wins[i].eall > 0) ? (wins[i].aall / wins[i].eall) : 1.0;
	    if (rTot >= delRatio) naFlag[i] = true;
	  }
	}
      }

      // Exclude NA windows
      std::vector<std::pair<int32_t, int32_t> > naiv;
      for(uint32_t i = 0; i < nw; ++i) {
	if (naFlag[i]) {
	  for(uint32_t k = wins[i].start; k < wins[i].end; ++k) uniqContent[k] = 0;
	  int32_t nas = (int32_t) wins[i].start;
	  int32_t nae = (int32_t) wins[i].end;
	  if ((!naiv.empty()) && (naiv.back().second >= nas)) naiv.back().second = std::max(naiv.back().second, nae);
	  else naiv.push_back(std::make_pair(nas, nae));
	}
      }

      // Split-read breakpoints
      std::vector<SVBreakpoint> chrbp;
      collectBreakpoints(c, gcbound, gcContent, uniqContent, gcbias, cov, hdr, refIndex, clips, chrbp);

      // CNV discovery
      if (!c.hasGenoFile) segmentRD(c, gcbound, gcContent, uniqContent, gcbias, tileFac, regWin, cov, hdr, refIndex, chrbp, uniqueTrack(), naiv, cnvs);

      // CNV genotyping
      genotypeCNVs(c, gcbound, gcContent, uniqContent, gcbias, tileFac, regWin, cov, covUniq, covMap, covAll, ref, hdr, refIndex, cnvs);
      if (ref != NULL) free(ref);

      // Write windows
      if (!c.covfile.empty()) {
	std::string chrn(hdr->target_name[refIndex]);
	for(uint32_t i = 0; i < nw; ++i) {
	  double uniqFrac;
	  if (uniqGate) uniqFrac = (wins[i].tcov > 0) ? (wins[i].ucov / wins[i].tcov) : -1.0;
	  else uniqFrac = (wins[i].end > wins[i].start) ? ((double) wins[i].winlen / (double) (wins[i].end - wins[i].start)) : -1.0;
	  if (naFlag[i]) {
	    dataOut << chrn << "\t" << wins[i].start << "\t" << wins[i].end << "\t" << uniqFrac << "\tNA\tNA" << std::endl;
	  } else {
	    double cn = chrPloidy;
	    double logR = 0;
	    if (wins[i].expcov > 0) {
	      cn = (c.expectedCN * wins[i].covsum / wins[i].expcov - chrCtrlPloidy * (1 - c.purity)) / c.purity;
	      logR = std::log2((wins[i].covsum + 1.0) / (wins[i].expcov + 1.0));
	    }
	    dataOut << chrn << "\t" << wins[i].start << "\t" << wins[i].end << "\t" << uniqFrac << "\t" << logR << "\t" << cn << std::endl;
	  }
	}
      }
    }

    // Sort CNVs
    sort(cnvs.begin(), cnvs.end());

    // Merge CNVs
    if (!c.hasGenoFile) mergeAdjacentSameCN(cnvs, c.cnMergeTol);

    // Exclude regions
    if (c.hasExcludeFile) {
      typedef boost::icl::interval_set<uint32_t> TChrIntervals;
      typedef TChrIntervals::interval_type TIVal;
      std::vector<TChrIntervals> validRegions;
      if (_parseExcludeIntervals(c, hdr, validRegions)) {
	std::vector<CNV> keep;
	keep.reserve(cnvs.size());
	for(uint32_t i = 0; i < cnvs.size(); ++i) {
	  int32_t s = cnvs[i].start;
	  int32_t e = cnvs[i].end;
	  if ((e <= s) || (cnvs[i].chr < 0) || (cnvs[i].chr >= (int32_t) validRegions.size())) { keep.push_back(cnvs[i]); continue; }
	  TChrIntervals ci;
	  ci.insert(TIVal::right_open((uint32_t) s, (uint32_t) e));
	  TChrIntervals inter = ci & validRegions[cnvs[i].chr];
	  uint32_t validbp = 0;
	  for(TChrIntervals::const_iterator it = inter.begin(); it != inter.end(); ++it) validbp += (it->upper() - it->lower());
	  if ((double) validbp >= 0.5 * (double) (e - s)) keep.push_back(cnvs[i]);
	}
	cnvs.swap(keep);
      }
    }

    // Genotype CNVs
    cnvVCF(c, cnvs);

    // clean-up
    fai_destroy(faiRef);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
    if (!c.covfile.empty()) {
      dataOut.pop();
      dataOut.pop();
    }
    
    return 0;
  }

  
  // Purity/ploidy fit
  struct FitBin {
    int32_t refIndex;
    double ratio;
    double baf;
    double hetFrac;

    FitBin(int32_t const ri, double const r) : refIndex(ri), ratio(r), baf(-1), hetFrac(-1) {}
  };

  // Noise estimation
  inline double
  _binNoise(std::vector<double> const& diffs, double const minNoise) {
    if (diffs.size() < 10) return minNoise;
    std::vector<double> d(diffs);
    for(uint32_t i = 0; i < d.size(); ++i) d[i] = std::abs(d[i]);
    std::sort(d.begin(), d.end());
    double const mad = d[d.size() / 2] * 1.4826 / std::sqrt(2.0);
    return std::max(minNoise, mad);
  }

  // BAF
  inline double
  _bafObserved(double const b, double const bal) {
    return 0.5 - std::sqrt((0.5 - b) * (0.5 - b) + (0.5 - bal) * (0.5 - bal));
  }

  // Nearest integer CN
  inline double
  _binChi2(FitBin const& bin, double const purity, double const expCN, double const ctrl, double const bal, double const sigR, double const sigB, double const minFrac, int32_t& n, int32_t& nB) {
    double const cnEst = (expCN * bin.ratio - ctrl * (1.0 - purity)) / purity;
    int32_t const n0 = std::max(1, (int32_t) std::floor(cnEst));
    double best = std::numeric_limits<double>::max();
    n = n0;
    nB = -1;
    for(int32_t cn = n0; cn <= n0 + 1; ++cn) {
      double const total = purity * cn + (1.0 - purity) * ctrl;
      double const rExp = total / expCN;
      double const chiR = (bin.ratio - rExp) / sigR;
      if (bin.baf < 0) {
	if (chiR * chiR < best) { best = chiR * chiR; n = cn; nB = -1; }
	continue;
      }
      // CN
      int32_t b0 = 0;
      int32_t b1 = 0;
      if (bin.baf > 0) {
	double const dev2 = (0.5 - bin.baf) * (0.5 - bin.baf) - (0.5 - bal) * (0.5 - bal);
	double const bafTrue = (dev2 > 0) ? (0.5 - std::sqrt(dev2)) : 0.5;
	double const bEst = (bafTrue * expCN * bin.ratio - 0.5 * ctrl * (1.0 - purity)) / purity;
	b0 = std::max(0, std::min(cn / 2, (int32_t) std::floor(bEst)));
	b1 = std::max(0, std::min(cn / 2, b0 + 1));
      }
      for(int32_t b = b0; b <= b1; ++b) {
	double const bExp = _bafObserved((total > 0) ? ((purity * b + 0.5 * ctrl * (1.0 - purity)) / total) : 0.5, bal);
	double chiB = 0;
	if (bin.baf > 0) chiB = (bin.baf - bExp) / sigB;
	else if (bExp > minFrac) chiB = (bExp - minFrac) / sigB;
	double const chi2 = chiR * chiR + chiB * chiB;
	if (chi2 < best) { best = chi2; n = cn; nB = b; }
      }
    }
    return best;
  }

  // Fit
  struct FitNorm {
    double gauss1;
    double gauss2;
    double outlier1;
    double outlier2;

    FitNorm(double const sigR, double const sigB, double const maxRatio, double const bafRange) {
      double const pi = std::acos(-1.0);
      gauss1 = std::log(std::sqrt(2.0 * pi) * sigR);
      gauss2 = std::log(2.0 * pi * sigR * sigB);
      outlier1 = std::log(maxRatio);
      outlier2 = std::log(maxRatio * bafRange);
    }
  };

  // Negative log-likelihood of a CN state
  inline double
  _binNLL(FitBin const& bin, double const purity, double const expCN, double const ctrl, double const bal, double const sigR, double const sigB, double const minFrac, FitNorm const& norm, int32_t& n, int32_t& nB) {
    double const chi2 = _binChi2(bin, purity, expCN, ctrl, bal, sigR, sigB, minFrac, n, nB);
    double const nll = 0.5 * chi2 + ((bin.baf < 0) ? norm.gauss1 : norm.gauss2);
    double const out = (bin.baf < 0) ? norm.outlier1 : norm.outlier2;
    if (nll <= out) return nll;
    n = -1;
    nB = -1;
    return out;
  }

  // Score fit
  inline double
  _fitScore(std::vector<FitBin> const& bins, double const purity, double const ploidy, double const ctrl, double const bal, double const sigR, double const sigB, double const minFrac, FitNorm const& norm) {
    double const expCN = purity * ploidy + (1.0 - purity) * ctrl;
    typedef std::pair<int32_t, int32_t> TState;
    std::vector<TState> st(bins.size());
    std::vector<double> nll(bins.size());
    std::map<TState, uint32_t> occ;
    for(uint32_t i = 0; i < bins.size(); ++i) {
      nll[i] = _binNLL(bins[i], purity, expCN, ctrl, bal, sigR, sigB, minFrac, norm, st[i].first, st[i].second);
      ++occ[st[i]];
    }
    double sum = 0;
    for(uint32_t i = 0; i < bins.size(); ++i) sum += nll[i] - std::log((double) occ[st[i]] / (double) bins.size());
    return sum / (double) bins.size();
  }

  // Joint purity and ploidy estimation
  template<typename TConfig>
  inline void
  fitPurityPloidy(TConfig const& c, bam_hdr_t const* hdr, std::vector<RegionRatio> const& regRatio, uint32_t const regWin, std::vector< std::vector<AlleleCount> > const& snv, bool const fixedPurity, float& purity, float& ploidy) {
    uint32_t const binSize = std::max((uint32_t) 1000000, 5 * regWin);
    uint32_t const minBins = 100;
    uint32_t const minSites = 100;
    uint32_t const minHet = 20;
    double const minFrac = 0.1;
    double const minBal = 0.4;
    double const balStep = 0.02;
    double const lohFrac = 0.25;
    double const knownHet = 0.5;
    double const minNoise = 0.01;
    double const minSep = 4.0;
    double const minRange = 1.5;
    double const tolFrac = 0.03;
    double const family = 0.5;
    double const gainFrac = 0.07;
    double const minState = 0.05;
    double const ctrl = c.ctrlPloidy;
    double const defPurity = purity;
    double const defPloidy = ploidy;

    // Autosomal bins
    std::vector<FitBin> bins;
    std::vector<int32_t> binStart;
    std::vector<double> cur;
    int32_t curRef = -1;
    int64_t curBin = -1;
    for(uint32_t i = 0; i <= regRatio.size(); ++i) {
      bool flush = (i == regRatio.size());
      if (!flush) {
	std::string tname(hdr->target_name[regRatio[i].refIndex]);
	if ((hdr->target_len[regRatio[i].refIndex] < c.minChrLen) || (c.refIdx.find(regRatio[i].refIndex) != c.refIdx.end())) continue;
	if ((tname == "chrX") || (tname == "chrY") || (tname == "X") || (tname == "Y") || (tname == "chrM") || (tname == "MT")) continue;
	flush = ((regRatio[i].refIndex != curRef) || ((int64_t) (regRatio[i].start / binSize) != curBin));
      }
      if (flush) {
	if (cur.size() >= 3) {
	  std::sort(cur.begin(), cur.end());
	  bins.push_back(FitBin(curRef, cur[cur.size() / 2]));
	  binStart.push_back(curBin * binSize);
	}
	cur.clear();
	if (i == regRatio.size()) break;
	curRef = regRatio[i].refIndex;
	curBin = (int64_t) (regRatio[i].start / binSize);
      }
      cur.push_back(regRatio[i].ratio);
    }
    if (bins.size() < minBins) {
      std::cerr << "Not enough bins to estimate purity and ploidy!" << std::endl;
      return;
    }

    // BAF and heterozygosity per bin
    HetBand const hb(minFrac);
    for(uint32_t i = 0; i < bins.size(); ++i) {
      BinBaf bb = binBaf(snv[bins[i].refIndex], binStart[i], binStart[i] + binSize, minHet, hb);
      if (bb.nsites >= minSites) bins[i].hetFrac = (double) bb.nhet / (double) bb.nsites;
      bins[i].baf = bb.baf;
    }

    // LOH
    {
      std::vector<double> hf;
      for(uint32_t i = 0; i < bins.size(); ++i) {
	if (bins[i].hetFrac >= 0) hf.push_back(bins[i].hetFrac);
      }
      if (!hf.empty()) {
	std::sort(hf.begin(), hf.end());
	double const typical = hf[(uint32_t) (0.9 * (double) (hf.size() - 1))];
	if (typical >= knownHet) {
	  for(uint32_t i = 0; i < bins.size(); ++i) {
	    if ((bins[i].hetFrac >= 0) && (bins[i].hetFrac < lohFrac * typical)) bins[i].baf = 0;
	  }
	}
      }
    }

    // Noise estimation
    double sigR = minNoise;
    double sigB = minNoise;
    {
      std::vector<double> dr;
      std::vector<double> db;
      for(uint32_t i = 1; i < bins.size(); ++i) {
	if (bins[i].refIndex != bins[i-1].refIndex) continue;
	dr.push_back(bins[i].ratio - bins[i-1].ratio);
	if ((bins[i].baf > 0) && (bins[i-1].baf > 0)) db.push_back(bins[i].baf - bins[i-1].baf);
      }
      sigR = _binNoise(dr, minNoise);
      sigB = _binNoise(db, minNoise);
    }

    // Outliers
    double maxRatio = minRange;
    {
      std::vector<double> r;
      for(uint32_t i = 0; i < bins.size(); ++i) r.push_back(bins[i].ratio);
      std::sort(r.begin(), r.end());
      maxRatio = std::max(maxRatio, r[(uint32_t) (0.99 * (double) (r.size() - 1))]);
    }
    FitNorm const norm(sigR, sigB, maxRatio, 0.5 - minFrac);
    double const gap = norm.outlier2 - norm.gauss2;   // log-likelihood difference of an outlier bin and a bin at its state
    double const tol = tolFrac * gap;
    double const gain = gainFrac * gap;

    // Grid search over purity and ploidy
    struct GridPoint {
      double score;
      double purity;
      double ploidy;
      double bal;

      GridPoint(double const s, double const pu, double const pl, double const b) : score(s), purity(pu), ploidy(pl), bal(b) {}
    };
    double const pMin = (fixedPurity) ? defPurity : 0.1;
    double const pMax = (fixedPurity) ? defPurity : 1.0;
    double best = std::numeric_limits<double>::max();
    std::vector<GridPoint> grid;
    for(double pl = 1.0; pl <= 6.0 + 1e-9; pl += 0.01) {
      for(double pu = pMin; pu <= pMax + 1e-9; pu += 0.01) {
	double const expCN = pu * pl + (1.0 - pu) * ctrl;
	if (pu / expCN < minSep * sigR) continue;
	for(double bal = minBal; bal <= 0.5 + 1e-9; bal += balStep) {
	  double const sc = _fitScore(bins, pu, pl, ctrl, bal, sigR, sigB, minFrac, norm);
	  grid.push_back(GridPoint(sc, pu, pl, bal));
	  if (sc < best) best = sc;
	}
      }
    }

    // Pick smallest ploidy
    double fitPurity = defPurity;
    double fitPloidy = defPloidy;
    double fitBal = 0.5;
    double fitScore = std::numeric_limits<double>::max();
    bool found = (!grid.empty());
    double const maxScore = best + tol;
    double minPloidy = 0;
    for(uint32_t k = 0; k < grid.size(); ++k) {
      if (grid[k].score <= maxScore) {
	minPloidy = grid[k].ploidy;
	break;
      }
    }
    for(uint32_t k = 0; k < grid.size(); ++k) {
      if ((grid[k].ploidy <= minPloidy + family) && (grid[k].score <= fitScore)) { 
	fitScore = grid[k].score;
	fitPurity = grid[k].purity;
	fitPloidy = grid[k].ploidy;
	fitBal = grid[k].bal;
      }
    }
    double defScore = std::numeric_limits<double>::max();
    double defBal = 0.5;
    for(double bal = minBal; bal <= 0.5 + 1e-9; bal += balStep) {
      double const sc = _fitScore(bins, defPurity, defPloidy, ctrl, bal, sigR, sigB, minFrac, norm);
      if (sc < defScore) {
	defScore = sc;
	defBal = bal;
      }
    }
    if (!found) {
      fitScore = defScore;
      fitBal = defBal;
    }

    // At least 2 CN
    uint32_t states = 0;
    {
      double const expCN = fitPurity * fitPloidy + (1.0 - fitPurity) * ctrl;
      std::map<std::pair<int32_t, int32_t>, uint32_t> occ;
      for(uint32_t i = 0; i < bins.size(); ++i) {
	int32_t n = 0;
	int32_t nB = 0;
	_binNLL(bins[i], fitPurity, expCN, ctrl, fitBal, sigR, sigB, minFrac, norm, n, nB);
	if (n >= 0) ++occ[std::make_pair(n, nB)];
      }
      for(std::map<std::pair<int32_t, int32_t>, uint32_t>::const_iterator it = occ.begin(); it != occ.end(); ++it) {
	if ((double) it->second >= minState * (double) bins.size()) ++states;
      }
    }
    if ((found) && (states >= 2) && (defScore - fitScore >= gain)) {
      purity = fitPurity;
      ploidy = fitPloidy;
    } else fitBal = defBal;
    uint32_t nbaf = 0;
    for(uint32_t i = 0; i < bins.size(); ++i) {
      if (bins[i].baf >= 0) ++nbaf;
    }
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Purity/ploidy fit: purity " << std::fixed << std::setprecision(2) << purity << ", ploidy " << ploidy << std::endl;
  }

  int coral(int argc, char **argv) {
    CountDNAConfig c;
    std::string haploidChr;
    std::string mode;
    std::string ploidyStr;
    uint16_t minBaseQual;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome file")
      ("exclude,x", boost::program_options::value<boost::filesystem::path>(&c.exclude), "file with regions to exclude")
      ("quality,q", boost::program_options::value<uint16_t>(&c.minQual)->default_value(10), "min. mapping quality")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "BCF output file")
      ("covfile,c", boost::program_options::value<boost::filesystem::path>(&c.covfile), "gzipped coverage file")
      ("segmentation,u", boost::program_options::value<boost::filesystem::path>(&c.segfile), "segmentation BED output file")
      ;

    boost::program_options::options_description cnv("CNV calling");
    cnv.add_options()
      ("cnv-size,z", boost::program_options::value<uint32_t>(&c.minCnvSize)->default_value(1000), "min. CNV size")
      ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.genofile), "input CNV BCF file for re-genotyping")
      ("minclip", boost::program_options::value<uint32_t>(&c.minClip)->default_value(25), "min. clipping length")
      ("minrefsep", boost::program_options::value<uint32_t>(&c.minRefSep)->default_value(30), "min. reference separation")
      ("min-bp-support", boost::program_options::value<uint32_t>(&c.minBpSupport)->default_value(3), "min. split-read support")
      ("penalty", boost::program_options::value<float>(&c.penalty)->default_value(3), "segmentation penalty")
      ("cnv-merge", boost::program_options::value<float>(&c.cnMergeTol)->default_value(0.25), "min. log2 ratio to separate CNVs")
      ("cnv-qual", boost::program_options::value<int32_t>(&c.cnvMinQual)->default_value(5), "min. quality for PASS")
      ;

    boost::program_options::options_description cancer("Ploidy/purity correction");
    cancer.add_options()
      ("mode,m", boost::program_options::value<std::string>(&mode)->default_value("germline"), "CNV calling mode [germline|somatic]")
      ("ploidy,y", boost::program_options::value<std::string>(&ploidyStr)->default_value("2"), "sample ploidy (number or auto)")
      ("purity,p", boost::program_options::value<float>(&c.purity)->default_value(1), "sample purity [0.1, 1]")
      ("snv-sites,b", boost::program_options::value<boost::filesystem::path>(&c.snvSites), "common SNV sites in VCF/BCF format (req. for -y auto)")
      ("ctrl-ploidy", boost::program_options::value<float>(&c.ctrlPloidy)->default_value(2), "control ploidy")
      ("haploid-chr", boost::program_options::value<std::string>(&haploidChr), "haploid chromosomes, e.g. chrX,chrY")
      ;
    
    boost::program_options::options_description window("Read-depth windows");
    window.add_options()
      ("window,w", boost::program_options::value<uint32_t>(&c.window_size)->default_value(0), "window size in bp (0: automatic)")
      ("fraction-unique", boost::program_options::value<float>(&c.uniqueToTotalCovRatio)->default_value(0.8), "uniqueness filter [0,1]")
      ("cnv-del-confirm", boost::program_options::value<float>(&c.cnvDelConfirm)->default_value(1.5), "max. total-depth CN for DEL")
      ("cnv-del-ratio", boost::program_options::value<float>(&c.cnvDelRatio)->default_value(0.75), "total-depth ratio")
      ("basecov", "force base-level counting")
      ("fragmentcov", "force fragment-level counting")
      ("no-regional-gc", "disable broad GC correction")
      ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamFile), "input BAM/CRAM file")
      ("fragment", boost::program_options::value<float>(&c.fragmentUnique)->default_value(0.97), "min. fragment uniqueness [0,1]")
      ("statsfile", boost::program_options::value<boost::filesystem::path>(&c.statsFile), "gzipped stats output file (optional)")
      ("window-offset", boost::program_options::value<uint32_t>(&c.window_offset)->default_value(0), "window offset (0: window size)")
      ("fraction-window", boost::program_options::value<float>(&c.fracWindow)->default_value(0.25), "min. callable window fraction [0,1]")
      ("mapq-uniq", boost::program_options::value<uint16_t>(&c.mapqUniq)->default_value(20), "min. MAPQ for a uniquely-placed read")
      ("target-reads", boost::program_options::value<uint32_t>(&c.targetReads)->default_value(150), "target reads/window")
      ("min-windows", boost::program_options::value<uint32_t>(&c.minSegWin)->default_value(2), "min. windows per segment")
      ("scan-window", boost::program_options::value<uint32_t>(&c.scanWindow)->default_value(10000), "GC scanning window size")
      ("scan-regions", boost::program_options::value<boost::filesystem::path>(&c.scanFile), "GC scanning regions in BED format")
      ("mad-cutoff", boost::program_options::value<uint16_t>(&c.mad)->default_value(3), "median + 3 * mad count cutoff")
      ("snv-depth", boost::program_options::value<uint32_t>(&c.minSnvDepth)->default_value(10), "min. depth at SNV sites")
      ("snv-baseq", boost::program_options::value<uint16_t>(&minBaseQual)->default_value(20), "min. base quality at SNV sites")
      ("percentile", boost::program_options::value<float>(&c.exclgc)->default_value(0.0005), "excl. extreme GC fraction")
      ("no-window-selection", "no scan window selection")
      ;

    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);

    // Set the visibility
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(cnv).add(cancer).add(window).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(cnv).add(cancer).add(window);

    // Parse command-line
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
      std::cerr << std::endl;
      std::cerr << "Usage: delly " << argv[0] << " [OPTIONS] -g <genome.fa> <aligned.bam>" << std::endl;
      std::cerr << visible_options << "\n";
      return 1;
    }

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cerr << "delly ";
    for(int i=0; i<argc; ++i) { std::cerr << argv[i] << ' '; }
    std::cerr << std::endl;

    // Stats file
    if (vm.count("statsfile")) c.hasStatsFile = true;
    else c.hasStatsFile = false;

    // Scan regions
    if (vm.count("scan-regions")) c.hasScanFile = true;
    else c.hasScanFile = false;

    // Exclude regions
    if (vm.count("exclude")) c.hasExcludeFile = true;
    else c.hasExcludeFile = false;

    // Scan window selection
    if (vm.count("no-window-selection")) c.noScanWindowSelection = true;
    else c.noScanWindowSelection = false;
    c.regionalGc = (!vm.count("no-regional-gc"));

    // Adaptive windows
    c.adaptive = false;
    c.targetExpCov = 0;
    if (c.window_size == 0) c.adaptive = true;
    if (c.targetReads == 0) c.targetReads = 150;

    // Calling mode
    c.somatic = (mode == "somatic");
    if (c.somatic) {
      // More sensitive segmentation
      if (vm["penalty"].defaulted()) c.penalty = 1.25;
      if (vm["cnv-merge"].defaulted()) c.cnMergeTol = 0;
    }

    // Ploidy
    c.ploidyAuto = (ploidyStr == "auto");
    if (c.ploidyAuto) c.ploidy = 2;
    else {
      try {
	c.ploidy = boost::lexical_cast<float>(ploidyStr);
      } catch (boost::bad_lexical_cast const&) {
	std::cerr << "Sample ploidy needs to be a number or auto: " << ploidyStr << std::endl;
	return 1;
      }
    }

    // SNV sites
    c.hasSnvSites = vm.count("snv-sites");
    c.minBaseQual = (uint8_t) std::min((uint16_t) 255, minBaseQual);
    if (c.ploidyAuto != c.hasSnvSites) {
      std::cerr << "Please provide common SNV sites (-b) for an automatic purity and ploidy fit." << std::endl;
      return 1;
    }
    if (c.hasSnvSites) {
      SnvSitesFile sf;
      bool ok = sf.open(c.snvSites.string());
      sf.close();
      if (!ok) {
	bcf_index_build(c.snvSites.string().c_str(), 14);
	SnvSitesFile sf2;
	bool ok2 = sf2.open(c.snvSites.string());
	sf2.close();
	if (!ok2) {
	  std::cerr << "Fail to open indexed BCF file with SNV sites: " << c.snvSites.string() << std::endl;
	  return 1;
	}
      }
    }

    // Purity
    if (c.purity > 1) c.purity = 1;
    if (c.purity < 0.1) c.purity = 0.1;

    // Expected ploidy
    c.expectedCN = c.purity * c.ploidy + (1.0 - c.purity) * c.ctrlPloidy;
    
    // Segmentation
    if (vm.count("segmentation")) c.hasSegFile = true;
    else c.hasSegFile = false;

    // Window offset
    if ((c.window_offset == 0) || (c.window_offset > c.window_size)) c.window_offset = c.window_size;

    // Check input VCF file (CNV genotyping)
    if (vm.count("vcffile")) {
      if (!(boost::filesystem::exists(c.genofile) && boost::filesystem::is_regular_file(c.genofile) && boost::filesystem::file_size(c.genofile))) {
	std::cerr << "Input VCF/BCF file is missing: " << c.genofile.string() << std::endl;
	return 1;
      }
      htsFile* ifile = bcf_open(c.genofile.string().c_str(), "r");
      if (ifile == NULL) {
	std::cerr << "Fail to open file " << c.genofile.string() << std::endl;
	return 1;
      }
      bcf_hdr_t* hdr = bcf_hdr_read(ifile);
      if (hdr == NULL) {
	std::cerr << "Fail to open index file " << c.genofile.string() << std::endl;
	return 1;
      }
      bcf_hdr_destroy(hdr);
      bcf_close(ifile);
      c.hasGenoFile = true;
    } else c.hasGenoFile = false;

    // Check outfile
    if (!vm.count("outfile")) c.outfile = "-";
    else {
      if (c.outfile.string() != "-") {
	if (!_outfileValid(c.outfile)) return 1;
      }
    }
    
    // Check bam file
    LibraryInfo li;
    bool pairedLib = false;
    if (!(boost::filesystem::exists(c.bamFile) && boost::filesystem::is_regular_file(c.bamFile) && boost::filesystem::file_size(c.bamFile))) {
      std::cerr << "Alignment file is missing: " << c.bamFile.string() << std::endl;
      return 1;
    } else {
      // Get scan regions
      typedef boost::icl::interval_set<uint32_t> TChrIntervals;
      typedef typename TChrIntervals::interval_type TIVal;
      typedef std::vector<TChrIntervals> TRegionsGenome;
      TRegionsGenome scanRegions;

      // Open BAM file
      samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
      if (samfile == NULL) {
	std::cerr << "Fail to open file " << c.bamFile.string() << std::endl;
	return 1;
      }
      hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
      if (idx == NULL) {
	if (bam_index_build(c.bamFile.string().c_str(), 0) != 0) {
	  std::cerr << "Fail to open index for " << c.bamFile.string() << std::endl;
	  return 1;
	}
      }
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.bamFile.string() << std::endl;
	return 1;
      }
      c.nchr = hdr->n_targets;
      c.minChrLen = setMinChrLen(hdr, 0.95);
      std::string sampleName = "unknown";
      getSMTag(std::string(hdr->text), c.bamFile.stem().string(), sampleName);
      c.sampleName = sampleName;

      // Check special chromosomes
      if (haploidChr.size()) _selectedRefIndices(c, hdr, haploidChr);

      // Check matching chromosome names
      faidx_t* faiRef = fai_load(c.genome.string().c_str());
      uint32_t refFound = 0;
      for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
	std::string tname(hdr->target_name[refIndex]);
	if (faidx_has_seq(faiRef, tname.c_str())) ++refFound;
	else {
	  std::cerr << "Warning: BAM chromosome " << tname << " not present in reference genome!" << std::endl;
	}
      }
      fai_destroy(faiRef);
      if (!refFound) {
	std::cerr << "Reference genome chromosome naming disagrees with BAM file!" << std::endl;
	return 1;
      }

      // Estimate library params
      if (c.hasScanFile) {
	if (!_parseBedIntervals(c.scanFile.string(), c.hasScanFile, hdr, scanRegions)) {
	  std::cerr << "Warning: Couldn't parse BED intervals. Do the chromosome names match?" << std::endl;
	  return 1;
	}
      } else {
	scanRegions.resize(hdr->n_targets);
	for (int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
	  scanRegions[refIndex].insert(TIVal::right_open(0, hdr->target_len[refIndex]));
	}
      }
      typedef std::vector<LibraryInfo> TSampleLibrary;
      TSampleLibrary sampleLib(1, LibraryInfo());
      CountDNAConfigLib dellyConf;
      dellyConf.genome = c.genome;
      dellyConf.files.push_back(c.bamFile);
      dellyConf.madCutoff = 9;
      dellyConf.madNormalCutoff = c.mad;
      getLibraryParams(dellyConf, scanRegions, sampleLib);
      li = sampleLib[0];
      pairedLib = (li.median > 0);
      if (!li.median) {
	li.median = 250;
	li.mad = 15;
	li.minNormalISize = 0;
	li.maxNormalISize = 400;
      }
      c.meanisize = ((int32_t) (li.median / 2)) * 2 + 1;

      // Coverage-aware GC scan window
      bool autoScanWin = (!vm.count("scan-window")) || (vm["scan-window"].defaulted());
      if ((autoScanWin) && (idx != NULL)) {
	uint64_t totalMapped = 0;
	uint64_t genomeLen = 0;
	for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
	  uint64_t mapped = 0;
	  uint64_t unmapped = 0;
	  hts_idx_get_stat(idx, refIndex, &mapped, &unmapped);
	  if (mapped > 0) {
	    totalMapped += mapped;
	    genomeLen += hdr->target_len[refIndex];
	  }
	}
	if (pairedLib) totalMapped /= 2;   // fragments
	if ((totalMapped > 0) && (genomeLen > 0)) {
	  double fragPerBp = (double) totalMapped / (double) genomeLen;
	  uint32_t targetScanReads = 30; 
	  uint32_t autoScan = (uint32_t) (targetScanReads / fragPerBp);
	  if (autoScan < c.scanWindow) autoScan = c.scanWindow;
	  uint32_t maxScan = 1000000;
	  if (autoScan > maxScan) autoScan = maxScan;
	  if (autoScan > c.scanWindow) c.scanWindow = autoScan;
	}
      }

      // Clean-up
      bam_hdr_destroy(hdr);
      hts_idx_destroy(idx);
      sam_close(samfile);
    }

    // Counting model
    if (vm.count("basecov")) c.basecov = true;
    else if (vm.count("fragmentcov")) c.basecov = false;
    else c.basecov = ((!pairedLib) && (li.rs >= 500));
    if ((c.basecov) && (!vm.count("cnv-del-confirm"))) c.cnvDelConfirm = std::numeric_limits<float>::max();

    // GC bias estimation
    typedef std::pair<uint32_t, uint32_t> TGCBound;
    TGCBound gcbound;
    std::vector<GcBias> gcbias(c.meanisize + 1, GcBias());
    typedef std::vector<ScanWindow> TWindowCounts;
    typedef std::vector<TWindowCounts> TGenomicWindowCounts;
    TGenomicWindowCounts scanCounts(c.nchr, TWindowCounts());
    std::vector< std::vector<AlleleCount> > snv(c.nchr, std::vector<AlleleCount>());
    {
      // Scan genomic windows
      scan(c, li, scanCounts);

      // Check coverage
      {
	std::vector<uint32_t> sampleScanVec;
	for(uint32_t i = 0; i < scanCounts.size(); ++i) {
	  for(uint32_t j = 0; j < scanCounts[i].size(); ++j) {
	    sampleScanVec.push_back(scanCounts[i][j].cov);
	  }
	  if (sampleScanVec.size() > 1000000) break;
	}
	std::sort(sampleScanVec.begin(), sampleScanVec.end());
	if (sampleScanVec.empty()) {
	  std::cerr << "Not enough windows!" << std::endl;
	  return 1;
	}
	if (sampleScanVec[sampleScanVec.size()/2] < 5) {
	  std::cerr << "Coverage in the GC scan window is too low." << std::endl;
	  return 1;
	}
      }
    
      // Select stable windows
      selectWindows(c, scanCounts);

      // Estimate GC bias and collect allele counts
      gcBias(c, scanCounts, li, gcbias, gcbound, snv);

      // Statistics output
      if (c.hasStatsFile) {
	// Open stats file
	boost::iostreams::filtering_ostream statsOut;
	statsOut.push(boost::iostreams::gzip_compressor());
	statsOut.push(boost::iostreams::file_sink(c.statsFile.string(), std::ios_base::out | std::ios_base::binary));
	
	// Library Info
	statsOut << "LP\t" << li.rs << ',' << li.median << ',' << li.mad << ',' << li.minNormalISize << ',' << li.maxNormalISize << std::endl;
	
	// Scan window summry
	samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
	bam_hdr_t* hdr = sam_hdr_read(samfile);
	statsOut << "SW\tchrom\tstart\tend\tselected\tcoverage\tuniqcov" <<  std::endl;
	for(uint32_t refIndex = 0; refIndex < (uint32_t) hdr->n_targets; ++refIndex) {
	  for(uint32_t i = 0; i < scanCounts[refIndex].size(); ++i) {
	    statsOut << "SW\t" <<  hdr->target_name[refIndex] << '\t' << scanCounts[refIndex][i].start << '\t' << scanCounts[refIndex][i].end << '\t' << scanCounts[refIndex][i].select << '\t' << scanCounts[refIndex][i].cov << '\t' << scanCounts[refIndex][i].uniqcov << std::endl;
	  }
	}
	bam_hdr_destroy(hdr);
	sam_close(samfile);
	
	// GC bias summary
	statsOut << "GC\tgcsum\tsample\treference\tpercentileSample\tpercentileReference\tfractionSample\tfractionReference\tobsexp\tmeancoverage\tmeancoverageTotal\treferenceTotal" << std::endl;
	for(uint32_t i = 0; i < gcbias.size(); ++i) statsOut << "GC\t" << i << "\t" << gcbias[i].sample << "\t" << gcbias[i].reference << "\t" << gcbias[i].percentileSample << "\t" << gcbias[i].percentileReference << "\t" << gcbias[i].fractionSample << "\t" << gcbias[i].fractionReference << "\t" << gcbias[i].obsexp << "\t" << gcbias[i].coverage << "\t" << gcbias[i].coverageTotal << "\t" << gcbias[i].referenceTotal << std::endl;
	statsOut << "BoundsGC\t" << gcbound.first << "," << gcbound.second << std::endl;
	statsOut.pop();
	statsOut.pop();
      }
    }

    // Coverage-aware window size
    uint32_t effWin = (c.window_size > 0) ? c.window_size : 50000;
    if (c.adaptive) {
      // Mean expected coverage at CN2 over the callable GC range
      double covMean = 0;
      uint64_t refCnt = 0;
      for(uint32_t i = gcbound.first + 1; i < gcbound.second; ++i) {
	covMean += gcbias[i].coverage * (double) gcbias[i].reference;
	refCnt += gcbias[i].reference;
      }
      if (refCnt) covMean /= (double) refCnt;
      if (covMean <= 0) {
	// Fixed windows
	c.adaptive = false;
	c.window_size = 10000;
	c.window_offset = c.window_size;
      } else {
	double readLen = (li.rs > 0) ? (double) li.rs : (double) c.meanisize;
	double molPerBp = c.basecov ? (covMean / readLen) : covMean;
	if (molPerBp <= 0) molPerBp = 1e-9;
	double winBp = (double) c.targetReads / molPerBp;
	double minWin = std::max(100.0, 4.0 * readLen);
	double maxWin = 2000000.0;
	if (winBp < minWin) winBp = minWin;
	if (winBp > maxWin) winBp = maxWin;
	c.targetExpCov = covMean * winBp;
	effWin = (uint32_t) winBp;
	double effReads = molPerBp * winBp;
	double covDepth = c.basecov ? covMean : (covMean * readLen);
	now = boost::posix_time::second_clock::local_time();
	std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Auto window size: " << (uint32_t) winBp << " bp, " << (uint32_t) effReads << " reads/window (" << (c.basecov ? "base-level" : "fragment") << ", coverage " << std::fixed << std::setprecision(2) << covDepth << std::defaultfloat << "x)" << std::endl;
      }
    }

    // Regional GC correction
    std::vector<double> regcorr;
    uint32_t regWin = std::max((uint32_t) 50000, effWin);
    if ((c.regionalGc) || (c.ploidyAuto)) {
      std::vector<RegionRatio> regRatio;
      estimateRegionalGc(c, gcbound, gcbias, scanCounts, regWin, regcorr, regRatio);
      if (!c.regionalGc) regcorr.clear();

      // Sample purity and ploidy
      if (c.ploidyAuto) {
	samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
	bam_hdr_t* hdr = sam_hdr_read(samfile);
	fitPurityPloidy(c, hdr, regRatio, regWin, snv, !vm["purity"].defaulted(), c.purity, c.ploidy);
	c.expectedCN = c.purity * c.ploidy + (1.0 - c.purity) * c.ctrlPloidy;
	bam_hdr_destroy(hdr);
	sam_close(samfile);
      }
    }

    // Count reads
    if (bamCount(c, li, gcbias, gcbound, regcorr, regWin)) {
      std::cerr << "Read counting error!" << std::endl;
      return 1;
    }

    // Done
    now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;

    return 0;
  }

  
}

#endif
