#ifndef METHYL_H
#define METHYL_H

#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <math.h>
#include <stdio.h>

#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "edlib.h"
#include "tags.h"
#include "util.h"

namespace torali {

  // Half-window size (bp) around each SV breakpoint for methylation collection
  static const int32_t METHYL_WINDOW = 2000;
  // Minimum ML probability to call a base as methylated (128/256 ≈ 50%)
  static const uint8_t METHYL_PROB_THRESHOLD = 128;

  struct MethylInfo {
    float altSvStartL;
    float altSvStartR;
    float altSvRightL;
    float altSvRightR;
    float refSvStartL;
    float refSvStartR;
    float refSvRightL;
    float refSvRightR;

    MethylInfo() : altSvStartL(-1.0f), altSvStartR(-1.0f), altSvRightL(-1.0f), altSvRightR(-1.0f),
                   refSvStartL(-1.0f), refSvStartR(-1.0f), refSvRightL(-1.0f), refSvRightR(-1.0f) {}
  };

  // Per-SV read-level methylation accumulator (counts before fractions are computed)
  struct MethylAccum {
    uint32_t altStartL_m, altStartL_t;
    uint32_t altStartR_m, altStartR_t;
    uint32_t altRightL_m, altRightL_t;
    uint32_t altRightR_m, altRightR_t;
    uint32_t refStartL_m, refStartL_t;
    uint32_t refStartR_m, refStartR_t;
    uint32_t refRightL_m, refRightL_t;
    uint32_t refRightR_m, refRightR_t;

    MethylAccum() : altStartL_m(0), altStartL_t(0), altStartR_m(0), altStartR_t(0),
                    altRightL_m(0), altRightL_t(0), altRightR_m(0), altRightR_t(0),
                    refStartL_m(0), refStartL_t(0), refStartR_m(0), refStartR_t(0),
                    refRightL_m(0), refRightL_t(0), refRightR_m(0), refRightR_t(0) {}
  };


  struct ModHit {
    int32_t pos;
    char code;
    uint8_t prob;
    bool rev;
    char base;

    ModHit(int32_t const p, char const c, uint8_t const r, bool const s, char const b) : pos(p), code(c), prob(r), rev(s), base(b) {}
  };

    struct BedEntry {
    uint32_t start;
    uint32_t end;
    std::string id;
    
    BedEntry(uint32_t const s, uint32_t const e, std::string const& i) : start(s), end(e), id(i) {}
  };
    
  inline char
  complement_base(char b) {
    switch (std::toupper(static_cast<unsigned char>(b))) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    default: return b;
    }
  }

  
  // Parse BED file
  template<typename TConfig, typename TRegionsGenome>
  inline int32_t
  _parseBedIntervals(TConfig const& c, TRegionsGenome& bedRegions) {
    // Open file handles
    samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Typedef
    typedef typename TRegionsGenome::value_type TChrIntervals;
    
    int32_t intervals = 0;
    bedRegions.resize(hdr->n_targets, TChrIntervals());
    std::ifstream chrFile(c.bedfile.string().c_str(), std::ifstream::in);
    if (chrFile.is_open()) {
      while (chrFile.good()) {
	std::string chrFromFile;
	getline(chrFile, chrFromFile);
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t,;");
	Tokenizer tokens(chrFromFile, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter!=tokens.end()) {
	  std::string chrName = *tokIter++;
	  int32_t tid = bam_name2id(hdr, chrName.c_str());
	  if (tid >= 0) {
	    if (tokIter!=tokens.end()) {
	      int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
	      int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
	      if (start >= end) {
		std::cerr << "Warning: BED interval start >= end! 0-based BED intervals are interpreted as [start, end(" << std::endl;
		continue;
	      }
	      bedRegions[tid].push_back(BedEntry(start, end, *tokIter));
	      ++intervals;
	    }
	  } else {
	    std::cerr << "Warning: Chromosome not present in BAM! " << chrName << std::endl;
	  }
	}
      }
      chrFile.close();
    }

    // Clean-up
    bam_hdr_destroy(hdr);
    sam_close(samfile);

    return intervals;
  }



  inline double _medianVectorD(std::vector<double> myVector) {
    std::vector<double> myVectorCopy = myVector;
    const auto middleItr = myVectorCopy.begin() + myVectorCopy.size() / 2;
    std::nth_element(myVectorCopy.begin(), middleItr, myVectorCopy.end());
    if (myVectorCopy.size() % 2 == 0) {
      const auto leftMiddleItr = std::max_element(myVectorCopy.begin(), middleItr);
      return (*leftMiddleItr + *middleItr) / 2.0;
    } else {
      return *middleItr;
    }
  }
  
  // Parse MM+ML tags for a single read and fill methCall[fwdPos]:
  //   1  = methylated 5mC (prob >= probModTh)
  //   0  = unmodified C (explicit in MM or implicit)
  //  -1  = position is not a C in fwdseq (no call)
  // fwdseq is the 5'->3' molecule sequence: same as BAM SEQ for forward reads,
  // reverseComplement(BAM SEQ) for reverse reads.
  // Returns true if the MM tag was present and parsed.
  inline bool
  buildMethylCalls(bam1_t* rec, uint8_t probModTh, std::vector<int8_t>& methCall) {
    int32_t l = rec->core.l_qseq;
    methCall.assign(l, -1);

    uint8_t* mm_aux = bam_aux_get(rec, "MM");
    if (!mm_aux || *mm_aux != 'Z') return false;

    bool readRev = (rec->core.flag & BAM_FREVERSE);

    // Decode BAM sequence
    const uint8_t* seqptr = bam_get_seq(rec);
    std::string sequence(l, 'N');
    for (int32_t i = 0; i < l; ++i)
      sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

    // Build forward-oriented molecule sequence (5'->3')
    std::string fwdseq = sequence;
    if (readRev) reverseComplement(fwdseq);

    // Build per-base position index from fwdseq
    std::unordered_map<char, std::vector<int32_t>> basepos;
    for (int32_t i = 0; i < l; ++i)
      basepos[std::toupper(static_cast<unsigned char>(fwdseq[i]))].push_back(i);

    // Parse MM tag (semicolon-separated modification strings).
    // Track whether any 5mC entry uses the '?' status (unlisted positions unknown).
    // With '?' the pre-marking of all C positions as unmodified must be skipped so
    // that only explicitly listed positions contribute to the methylation fraction.
    const char* mmstr = reinterpret_cast<const char*>(mm_aux + 1);
    std::string mmString(mmstr);
    std::vector<ModHit> modhits;
    std::vector<std::string> tokens;
    bool mCHasSkipStatus = false;
    boost::split(tokens, mmString, boost::is_any_of(";"));
    for (const auto& tok : tokens) {
      if (tok.empty()) continue;
      std::size_t idx = 0;
      if (idx >= tok.size()) continue;
      char base = tok[idx++];
      if (idx >= tok.size()) continue;
      char strand = tok[idx++];
      bool revMod = (strand == '-');
      std::string mod_codes;
      bool hasSkip = false;
      while (idx < tok.size() && tok[idx] != ',') {
        char ch = tok[idx++];
        if (ch == '?' || ch == '.') hasSkip = true;
        else if (std::isalpha(static_cast<unsigned char>(ch))) mod_codes.push_back(ch);
      }
      if (hasSkip)
        for (char c : mod_codes)
          if (c == 'm' || c == 'M') mCHasSkipStatus = true;
      if (idx < tok.size() && tok[idx] == ',') {
        std::string pos_str = tok.substr(idx + 1);
        if (!pos_str.empty()) {
          std::vector<std::string> pos_tokens;
          boost::split(pos_tokens, pos_str, boost::is_any_of(","));
          int32_t current = -1;
          for (const auto& pt : pos_tokens) {
            if (pt.empty()) continue;
            current += std::stoi(pt) + 1;
            for (char mc : mod_codes)
              modhits.emplace_back(current, mc, static_cast<uint8_t>(255), revMod, base);
          }
        }
      }
    }

    // Pre-mark all C positions as unmodified only when '?' is absent (implicit
    // unmodified assumption per SAM spec).  When '?' is present, unlisted positions
    // are unknown and must remain -1 so they are excluded from the denominator.
    if (!mCHasSkipStatus) {
      for (int32_t i = 0; i < l; ++i)
        if (std::toupper(static_cast<unsigned char>(fwdseq[i])) == 'C')
          methCall[i] = 0;
    }

    // Parse ML tag for modification probabilities
    uint8_t* ml_aux = bam_aux_get(rec, "ML");
    if (ml_aux && *ml_aux == 'B') {
      char subtype = *reinterpret_cast<char*>(ml_aux + 1);
      if (subtype == 'C') {
        const uint8_t* p = ml_aux + 2;
        int32_t n = (int32_t)(p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24));
        const uint8_t* data = p + 4;
        int32_t assign = std::min<int32_t>(n, (int32_t)modhits.size());
        for (int32_t i = 0; i < assign; ++i) modhits[i].prob = data[i];
      }
    }

    // Map 5mC modhits to fwdseq positions and record methylation state
    for (const auto& mh : modhits) {
      if (mh.code != 'm' && mh.code != 'M') continue; // only 5mC
      char ub = std::toupper(static_cast<unsigned char>(mh.base));
      char target_base = mh.rev ? complement_base(ub) : ub;
      auto it = basepos.find(target_base);
      if (it == basepos.end()) continue;
      const auto& occs = it->second;
      if (mh.pos < 0 || static_cast<std::size_t>(mh.pos) >= occs.size()) continue;
      int32_t fwdPos = occs[mh.pos];
      methCall[fwdPos] = (mh.prob >= probModTh) ? 1 : 0;
    }

    return true;
  }

  // Collect methylation counts from a set of reference-coordinate windows in a
  // single CIGAR walk.  wins must be sorted by window start.  For each window,
  // methCounts[i] / totCounts[i] gives the methylated / total C count.
  inline void
  collectMethylFromWindows(bam1_t* rec,
                           const std::vector<int8_t>& methCall,
                           const std::vector<std::pair<int32_t,int32_t>>& wins,
                           std::vector<uint32_t>& methCounts,
                           std::vector<uint32_t>& totCounts) {
    methCounts.assign(wins.size(), 0);
    totCounts.assign(wins.size(), 0);
    if (wins.empty() || methCall.empty()) return;

    bool readRev = (rec->core.flag & BAM_FREVERSE);
    int32_t l = rec->core.l_qseq;
    // Maximum reference end across all windows for early exit
    int32_t maxEnd = 0;
    for (const auto& w : wins) maxEnd = std::max(maxEnd, w.second);

    int32_t rp = rec->core.pos;
    int32_t sp = 0;
    const uint32_t* cigar = bam_get_cigar(rec);
    for (uint32_t ci = 0; ci < rec->core.n_cigar; ++ci) {
      uint32_t op = bam_cigar_op(cigar[ci]);
      int32_t oplen = (int32_t)bam_cigar_oplen(cigar[ci]);
      if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
        for (int32_t k = 0; k < oplen; ++k, ++rp, ++sp) {
          if (rp >= maxEnd) return;
          int32_t fwdPos = readRev ? (l - sp - 1) : sp;
          if (fwdPos < 0 || fwdPos >= l) continue;
          int8_t call = methCall[fwdPos];
          if (call < 0) continue; // not a C in fwdseq
          for (std::size_t wi = 0; wi < wins.size(); ++wi) {
            if (rp >= wins[wi].first && rp < wins[wi].second) {
              ++totCounts[wi];
              if (call == 1) ++methCounts[wi];
            }
          }
        }
      } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
        rp += oplen;
        if (rp >= maxEnd) return;
      } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
        sp += oplen;
      }
      // BAM_CHARD_CLIP: bases not in SEQ, do not advance sp
    }
  }

  // Collect methylation counts from the insertion sequence of a read at a given
  // reference position.  useFirst=true takes the first windowLen bases of the
  // insertion; useFirst=false takes the last windowLen bases.
  inline void
  collectMethylFromInsertion(bam1_t* rec,
                              const std::vector<int8_t>& methCall,
                              int32_t insRefPos,
                              int32_t insLen,
                              bool useFirst,
                              int32_t windowLen,
                              uint32_t& methCount,
                              uint32_t& totCount) {
    methCount = 0;
    totCount = 0;
    if (methCall.empty() || insLen <= 0 || windowLen <= 0) return;

    bool readRev = (rec->core.flag & BAM_FREVERSE);
    int32_t l = rec->core.l_qseq;
    int32_t rp = rec->core.pos;
    int32_t sp = 0;
    const uint32_t* cigar = bam_get_cigar(rec);
    for (uint32_t ci = 0; ci < rec->core.n_cigar; ++ci) {
      uint32_t op = bam_cigar_op(cigar[ci]);
      int32_t oplen = (int32_t)bam_cigar_oplen(cigar[ci]);
      if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
        for (int32_t k = 0; k < oplen; ++k, ++rp, ++sp) {
          if (rp > insRefPos) return; // passed insertion site
        }
      } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
        rp += oplen;
        if (rp > insRefPos) return;
      } else if (op == BAM_CINS) {
        if (rp == insRefPos) {
          // Found the insertion at the right reference position
          int32_t wlen = std::min(windowLen, oplen);
          int32_t spBeg = useFirst ? sp : (sp + oplen - wlen);
          int32_t spEnd = useFirst ? (sp + wlen) : (sp + oplen);
          for (int32_t si = spBeg; si < spEnd && si < l; ++si) {
            int32_t fwdPos = readRev ? (l - si - 1) : si;
            if (fwdPos < 0 || fwdPos >= l) continue;
            int8_t call = methCall[fwdPos];
            if (call < 0) continue;
            ++totCount;
            if (call == 1) ++methCount;
          }
          return;
        }
        sp += oplen;
      } else if (op == BAM_CSOFT_CLIP) {
        sp += oplen;
      }
      // BAM_CHARD_CLIP: bases not in SEQ, do not advance sp
    }
  }

  // Compute final methylation fractions from accumulated counts.
  // Fields with zero total coverage are set to -1 (no data).
  inline void
  finalizeMethylInfo(const MethylAccum& ma, MethylInfo& mi) {
    auto frac = [](uint32_t m, uint32_t t) -> float {
      return t > 0 ? float(m) / float(t) : -1.0f;
    };
    mi.altSvStartL = frac(ma.altStartL_m, ma.altStartL_t);
    mi.altSvStartR = frac(ma.altStartR_m, ma.altStartR_t);
    mi.altSvRightL = frac(ma.altRightL_m, ma.altRightL_t);
    mi.altSvRightR = frac(ma.altRightR_m, ma.altRightR_t);
    mi.refSvStartL = frac(ma.refStartL_m, ma.refStartL_t);
    mi.refSvStartR = frac(ma.refStartR_m, ma.refStartR_t);
    mi.refSvRightL = frac(ma.refRightL_m, ma.refRightL_t);
    mi.refSvRightR = frac(ma.refRightR_m, ma.refRightR_t);
  }

  // Accumulate per-read methylation into the MethylAccum for one SV.
  // Collects up to METHYL_WINDOW bp flanking each breakpoint the read overlaps.
  // For ALT insertion reads the internal windows are taken from the inserted
  // sequence itself; for all other SV types internal windows are reference-based
  // and capped at the opposite breakpoint to stay "inside" the SV.
  inline void
  accumulateMethyl(bam1_t* rec,
                   const std::vector<int8_t>& methCall,
                   const StructuralVariantRecord& sv,
                   int32_t refIndex,
                   int32_t chromLen,
                   bool isAlt,
                   const std::vector<int32_t>& candidates,
                   MethylAccum& accum) {
    if (methCall.empty()) return;

    int32_t svStart = sv.svStart;
    int32_t svEnd   = sv.svEnd;
    int32_t insLen  = sv.insLen;
    bool isTrans = _translocation(sv.svt);
    bool isIns   = (sv.svt == 4);

    // Determine which breakpoints this read overlaps on the current chromosome
    bool onStart = false, onEnd = false;
    for (int32_t cand : candidates) {
      if (cand == svStart && sv.chr  == refIndex) onStart = true;
      if (cand == svEnd   && sv.chr2 == refIndex) onEnd   = true;
    }
    if (!onStart && !onEnd) return;

    // Collect reference-coordinate windows; field codes: 0=StartL, 1=StartR, 2=RightL, 3=RightR
    std::vector<std::pair<int32_t,int32_t>> wins;
    std::vector<int32_t> field;
    bool doInsStartR = false, doInsRightL = false;

    if (onStart) {
      // StartL: up to METHYL_WINDOW bp left of svStart
      int32_t beg = std::max(0, svStart - METHYL_WINDOW);
      int32_t end = svStart;
      if (end > beg) { wins.emplace_back(beg, end); field.push_back(0); }

      // StartR: depends on SV type
      if (isAlt && isIns && insLen > 0) {
        doInsStartR = true; // first windowLen bp of inserted sequence
      } else {
        int32_t rbeg = svStart;
        int32_t rend;
        if (!isTrans && !isIns)
          rend = std::min(svStart + METHYL_WINDOW, svEnd); // cap at svEnd (internal only)
        else
          rend = std::min(svStart + METHYL_WINDOW, chromLen);
        if (rend > rbeg) { wins.emplace_back(rbeg, rend); field.push_back(1); }
      }
    }

    if (onEnd) {
      // RightL: depends on SV type
      if (isAlt && isIns && insLen > 0) {
        doInsRightL = true; // last windowLen bp of inserted sequence
      } else {
        int32_t lbeg;
        if (!isTrans && !isIns)
          lbeg = std::max(svStart, svEnd - METHYL_WINDOW); // floor at svStart (internal only)
        else
          lbeg = std::max(0, svEnd - METHYL_WINDOW);
        int32_t lend = svEnd;
        if (lend > lbeg) { wins.emplace_back(lbeg, lend); field.push_back(2); }
      }

      // RightR: up to METHYL_WINDOW bp right of svEnd
      int32_t rbeg = svEnd;
      int32_t rend = std::min(svEnd + METHYL_WINDOW, chromLen);
      if (rend > rbeg) { wins.emplace_back(rbeg, rend); field.push_back(3); }
    }

    // Collect from reference windows in a single CIGAR walk
    if (!wins.empty()) {
      std::vector<uint32_t> methCounts, totCounts;
      collectMethylFromWindows(rec, methCall, wins, methCounts, totCounts);
      for (std::size_t i = 0; i < field.size(); ++i) {
        uint32_t m = methCounts[i], t = totCounts[i];
        if (isAlt) {
          switch (field[i]) {
            case 0: accum.altStartL_m += m; accum.altStartL_t += t; break;
            case 1: accum.altStartR_m += m; accum.altStartR_t += t; break;
            case 2: accum.altRightL_m += m; accum.altRightL_t += t; break;
            case 3: accum.altRightR_m += m; accum.altRightR_t += t; break;
          }
        } else {
          switch (field[i]) {
            case 0: accum.refStartL_m += m; accum.refStartL_t += t; break;
            case 1: accum.refStartR_m += m; accum.refStartR_t += t; break;
            case 2: accum.refRightL_m += m; accum.refRightL_t += t; break;
            case 3: accum.refRightR_m += m; accum.refRightR_t += t; break;
          }
        }
      }
    }

    // Collect from inserted sequence (ALT insertion reads only)
    if (doInsStartR) {
      int32_t wlen = std::min(METHYL_WINDOW, insLen);
      uint32_t m = 0, t = 0;
      collectMethylFromInsertion(rec, methCall, svStart, insLen, true, wlen, m, t);
      accum.altStartR_m += m;
      accum.altStartR_t += t;
    }
    if (doInsRightL) {
      int32_t wlen = std::min(METHYL_WINDOW, insLen);
      uint32_t m = 0, t = 0;
      collectMethylFromInsertion(rec, methCall, svStart, insLen, false, wlen, m, t);
      accum.altRightL_m += m;
      accum.altRightL_t += t;
    }
  }

  template<typename TConfig, typename TPeaks>
  inline void
  peakCalling(TConfig const& c, std::vector<float> const& modfrac, std::vector<uint32_t> const& modpos, TPeaks& peaks) {
    peaks.clear();
    uint32_t vecsize = modfrac.size();
    std::vector<double> avg(vecsize);
    double cumsum = 0;
    double lagging_cumsum = 0;
    for(uint32_t i = 0; i<c.peakwin; ++i) cumsum += modfrac[i];
    uint32_t half_win = (c.peakwin - 1) / 2;
    for(uint32_t k = 0; k <= half_win; ++k) avg[k] = cumsum/(double)(c.peakwin);
    for(uint32_t i = c.peakwin; i < vecsize; ++i) {
      cumsum += modfrac[i];
      lagging_cumsum += modfrac[i - c.peakwin];
      avg[i - half_win] = (cumsum - lagging_cumsum)/(double)(c.peakwin);
    }
    for(uint32_t k = (vecsize - half_win); k < vecsize; ++k) avg[k] = (cumsum - lagging_cumsum)/(double)(c.peakwin);

    // Debug peak signal
    //for(uint32_t i = 0; i < vecsize; ++i) std::cerr << i << '\t' <<  modpos[i] << '\t' << modfrac[i] << '\t' << avg[i] << std::endl;

    // Rolling median
    std::vector<double> medavg(vecsize);
    std::vector<double> rolling_median(c.peakmedwin);
    half_win = (c.peakmedwin - 1) / 2;
    for(uint32_t i = 0; i<c.peakmedwin; ++i) rolling_median[i] = avg[i];
    double medInit = _medianVectorD(rolling_median);
    for(uint32_t k = 0; k <= half_win; ++k) medavg[k] = medInit;
    for(uint32_t i = c.peakmedwin; i < vecsize; ++i) {
      rolling_median[i % c.peakmedwin] = avg[i];
      medavg[i - half_win] = _medianVectorD(rolling_median);
    }
    medInit = _medianVectorD(rolling_median);
    for(uint32_t k = (vecsize - half_win); k < vecsize; ++k) medavg[k] = medInit;

    // Peaks
    std::vector<uint32_t> beginPeak;
    std::vector<uint32_t> endPeak;
    if (medavg[0] > c.peakmethylth) beginPeak.push_back(0);
    for(uint32_t i = 0; i < (vecsize - 1); ++i) {
      if ((medavg[i] <= c.peakmethylth) && (medavg[i+1] > c.peakmethylth)) beginPeak.push_back(i+1);
      if ((medavg[i] > c.peakmethylth) && (medavg[i+1] <= c.peakmethylth)) endPeak.push_back(i+1);
    }
    if (medavg[vecsize - 1] > c.peakmethylth) endPeak.push_back(vecsize);

    // Debug
    /*
    char peaksign = 'x';
    for(uint32_t i = 0; i < vecsize; ++i) {
      if (std::find(beginPeak.begin(), beginPeak.end(), i) != beginPeak.end()) peaksign = 'o';
      if (std::find(endPeak.begin(), endPeak.end(), i) != endPeak.end()) peaksign = 'x';
      std::cerr << i << '\t' <<  modpos[i] << '\t' << modfrac[i] << '\t' << avg[i] << '\t' << medavg[i] << '\t' << peaksign << std::endl;
    }
    */

    if (beginPeak.size()) {
      for(uint32_t i = 0; i < beginPeak.size(); ++i) {
	//std::cerr << beginPeak[i] << ',' << endPeak[i] << std::endl;
	uint32_t peakMods = endPeak[i] - beginPeak[i];
	double avgPeakSig = 0;
	for(uint32_t k = beginPeak[i]; k < endPeak[i]; ++k) avgPeakSig+=modfrac[k];
	avgPeakSig /= (double) peakMods;

	//std::cerr << "Peak\t" << modpos[beginPeak[i]] << '\t' << modpos[endPeak[i]] << '\t' << peakMods << '\t' << avgPeakSig << std::endl;
	std::string idname = "numMods" + boost::lexical_cast<std::string>(peakMods) + "sig" + boost::lexical_cast<std::string>((int)(100 * avgPeakSig));
	peaks.push_back(BedEntry(modpos[beginPeak[i]], modpos[endPeak[i]], idname));
      }
    }
  }
  
  template<typename TConfig>
  inline bool
  screenMods(TConfig const& c) {
    // Open file handles
    samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Output table
    std::streambuf* buf;
    std::ofstream of;
    if(c.outfile.string() != "-") {
      of.open(c.outfile.string().c_str());
      buf = of.rdbuf();
    } else buf = std::cout.rdbuf();
    std::ostream out(buf);
    out << "chr\tstart\tend\tname\tmodbase\tstrand\tmodcount\tunmod\tothermod\tcoverage\tfail\tfracmod" << std::endl;

    // Load FASTA
    char* seq = NULL;
    faidx_t* fai = fai_load(c.genome.string().c_str());
    if (fai == NULL) {
      if (fai_build(c.genome.string().c_str()) == -1) {
	std::cerr << "Fail to open genome fai index for " << c.genome.string() << std::endl;
	return false;
      } else fai = fai_load(c.genome.string().c_str());
    }
    
    // Keep track of avg. coverage
    typedef boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance> > TAccumulator;
    TAccumulator acc;
    
    // Parse genome chr-by-chr
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Methylation scanning" << std::endl;

    // Parse BED
    typedef std::vector<BedEntry> TChrIntervals;
    typedef std::vector<TChrIntervals> TRegionsGenome;
    TRegionsGenome scanRegions;
    if (c.hasBedFile) {
      scanRegions.resize(hdr->n_targets, TChrIntervals());
      if (!c.callPeaks) {
	int32_t nreg = _parseBedIntervals(c, scanRegions);
	if (nreg == 0) {
	  std::cerr << "Error: Couldn't parse BED intervals. Do the chromosome names match?" << std::endl;
	  return 1;
	} else {
	  std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parsed " << nreg << " input regions." << std::endl;
	}
      }
    }
    
    // Calling threshold
    uint8_t probModTh = (uint8_t) ((int) (c.minMod * 256));
    uint8_t probUnmodTh = (uint8_t) ((int) (c.maxUnmod * 256));
    
    // Iterate chromosomes
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      // Load chromosome
      std::string tname = hdr->target_name[refIndex];
      int32_t seqlen = -1;
      seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
      
      // Counters
      typedef uint16_t TCovValue;
      TCovValue maxval = std::numeric_limits<TCovValue>::max();
      std::vector<TCovValue> unmod_plus(hdr->target_len[refIndex], 0);
      std::vector<TCovValue> unmod_minus(hdr->target_len[refIndex], 0);
      std::vector<TCovValue> uncertain_plus(hdr->target_len[refIndex], 0);
      std::vector<TCovValue> uncertain_minus(hdr->target_len[refIndex], 0);
      std::vector<TCovValue> m_plus(hdr->target_len[refIndex], 0);
      std::vector<TCovValue> m_minus(hdr->target_len[refIndex], 0);
      std::vector<TCovValue> h_plus(hdr->target_len[refIndex], 0);
      std::vector<TCovValue> h_minus(hdr->target_len[refIndex], 0);
      
      // Read alignments
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	
	// Keep secondary alignments
	if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	
	  // Parse quality and seq
	bool readRev = (bool) (rec->core.flag & BAM_FREVERSE);
	typedef std::vector<uint8_t> TQuality;
	TQuality quality;
	quality.resize(rec->core.l_qseq);
	std::string sequence;
	sequence.resize(rec->core.l_qseq);
	uint8_t const* seqptr = bam_get_seq(rec);
	uint8_t const* qualptr = bam_get_qual(rec);
	for (int32_t i = 0; i < rec->core.l_qseq; ++i) {
	  quality[i] = qualptr[i];
	  sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	}
	std::string fwdseq = sequence;
	if (readRev) reverseComplement(fwdseq);
	
	// Parse MM and ML tags (if present)
	std::vector<ModHit> modhits;
	uint8_t* mm_aux = bam_aux_get(rec, "MM");
	if (mm_aux && *mm_aux == 'Z') {
	  const char* mmstr = (const char*)(mm_aux + 1);
	  std::string s(mmstr);
	  std::vector<std::string> tokens;
	  boost::split(tokens, s, boost::is_any_of(";"));
	  for(const auto& tok : tokens) {
	    if (tok.empty()) continue;
	    // parse <base><strand><modcode>,<pos>,<pos>,... e.g. C+m,0,5,...
	    std::size_t idx_tok = 0;
	    char base = tok[idx_tok++];
	    if (idx_tok >= tok.size()) continue;
	    char strand = tok[idx_tok++];
	    bool revMod = false;
	    if (strand == '-') revMod = true;
	    std::string mod_codes;
	    while ((idx_tok < tok.size()) && (tok[idx_tok] != ',')) {
	      char ch = tok[idx_tok++];
	      if (std::isalpha(static_cast<unsigned char>(ch))) mod_codes.push_back(ch); // ignore ?
	    }
	    // parse positions
	    if ((idx_tok < tok.size()) && (tok[idx_tok] == ',')) {
	      std::string pos_str = tok.substr(idx_tok+1);
	      if (!pos_str.empty()) {
		std::vector<std::string> pos_tokens;
		boost::split(pos_tokens, pos_str, boost::is_any_of(","));
		int32_t current = -1;
		for(const auto& ptoken : pos_tokens) {
		  if (ptoken.empty()) continue;
		  int32_t delta = std::stoi(ptoken);
		  current += (delta + 1);
		  for(char mc : mod_codes) modhits.push_back(ModHit(current, mc, 255, revMod, base));
		}
	      }
	    }
	  }
	}

	// Parse ML tag
	uint8_t* ml_aux = bam_aux_get(rec, "ML");
	if (ml_aux && *ml_aux == 'B') {
	  char subtype = *(char*)(ml_aux + 1);
	  if (subtype == 'C') {
	    const uint8_t* p = ml_aux + 2;
	    int32_t n = 0;
	    n = p[0] | (p[1]<<8) | (p[2]<<16) | (p[3]<<24);
	    const uint8_t* data = p + 4;
	    int32_t assign = std::min<int32_t>(n, (int32_t) modhits.size());
	    for(int32_t i = 0; i < assign; ++i) modhits[i].prob = data[i];
	  }
	}
	
	// Find read position of As, Cs, ...
	std::unordered_map<char, std::vector<int32_t>> basepos;
	for (int32_t i = 0; i < (int32_t) fwdseq.size(); ++i) basepos[std::toupper(static_cast<unsigned char>(fwdseq[i]))].push_back(i);
	
	// Replace base positions with absolute read positions
	std::unordered_map<int32_t, std::vector<ModHit> > modByPos;
	for (const auto& mh : modhits) {
	  char ub = std::toupper(static_cast<unsigned char>(mh.base));
	  char target_base = (mh.rev) ? complement_base(ub) : ub;
	  auto it = basepos.find(target_base);
	  if (it == basepos.end()) continue;
	  const auto& occs = it->second;
	  if ((mh.pos < 0) || ((std::size_t)mh.pos >= occs.size())) continue;
	  // Absolute read pos is occs[mh.pos]
	  if ( (uint16_t) (quality[occs[mh.pos]]) >= c.minBaseQual) {
	    modByPos[occs[mh.pos]].push_back(ModHit(occs[mh.pos], mh.code, mh.prob, mh.rev, mh.base));
	  }
	}
	
	// Parse cigar
	uint32_t rp = rec->core.pos; // reference pointer
	uint32_t sp = 0; // sequence pointer
	uint32_t const* cigar = bam_get_cigar(rec);
	for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	  uint32_t op = bam_cigar_op(cigar[i]);
	  uint32_t oplen = bam_cigar_oplen(cigar[i]);
	  if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF)) {
	    for(uint32_t k = 0; k < oplen; ++k, ++rp, ++sp) {
	      if (rp >= hdr->target_len[refIndex]) break;
	      // Filter for base quality
	      if (quality[sp] >= c.minBaseQual) {
		uint32_t lookUpPos = sp;
		if (readRev) lookUpPos = rec->core.l_qseq - sp - 1;
		auto it = modByPos.find(lookUpPos);
		if (it != modByPos.end()) {
		  uint32_t unmod = 0;
		  uint32_t uncertain = 0;
		  std::set<char> modmatch;
		  for(const auto& mh : it->second) {
		    if (((mh.rev == readRev) && (sequence[sp] == mh.base)) || ((mh.rev != readRev) && (sequence[sp] == complement_base(mh.base)))) {
		      if (mh.prob >= probModTh) modmatch.insert(mh.code);
		      if (mh.prob > probUnmodTh) ++uncertain;
		      else ++unmod;
		    }
		  }
		  if (modmatch.size() == 1) {
		    char code = *(modmatch.begin());
		    if ((code == 'm') || (code == 'M')) {
		      if (readRev) {
			if (m_minus[rp] < maxval) ++m_minus[rp];
		      } else {
			if (m_plus[rp] < maxval) ++m_plus[rp];
		      }
		    } else if ((code == 'h') || (code == 'H')) {
		      if (readRev) {
			if (h_minus[rp] < maxval) ++h_minus[rp];
			} else {
			if (h_plus[rp] < maxval) ++h_plus[rp];
		      }
		    } else {
		      std::cerr << "Warning: Unknown modification code! " << code << std::endl;
		    }
		  } else if ((modmatch.size() > 1) || (uncertain)) {
		      // Fail
		    if (readRev) {
		      if (uncertain_minus[rp] < maxval) ++uncertain_minus[rp];
		    } else {
		      if (uncertain_plus[rp] < maxval) ++uncertain_plus[rp];
		    }
		  } else {
		    if (readRev) {
		      if (unmod_minus[rp] < maxval) ++unmod_minus[rp];
		    } else {
		      if (unmod_plus[rp] < maxval) ++unmod_plus[rp];
		    }
		  }
		}
	      }
	    }
	  } else if (op == BAM_CDEL) {
	    rp += oplen;
	  } else if (op == BAM_CINS) {
	    sp += oplen;
	  } else if ((op == BAM_CSOFT_CLIP) || (op == BAM_CHARD_CLIP)) {
	    sp += oplen;
	  } else if (op == BAM_CREF_SKIP) {
	    rp += oplen;
	  } else {
	    std::cerr << "Warning: Unknown Cigar operation!" << std::endl;
	  }
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
      
      // Summarize coverage for this chromosome
      for(uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
	uint32_t total_cov = m_plus[pos] + h_plus[pos] + unmod_plus[pos] + m_minus[pos] + h_minus[pos] + unmod_minus[pos] + uncertain_plus[pos] + uncertain_minus[pos];
	if (total_cov > 0) {
	  double failfrac = (double) (uncertain_plus[pos] + uncertain_minus[pos]) / (double) (total_cov);
	  acc(failfrac);
	}
      }
      double avgfail = boost::accumulators::mean(acc);
      std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << c.sampleName << ", " << hdr->target_name[refIndex] << ", running failed fraction of mod calls=" << avgfail << std::endl;
      
      // Identify peaks
      if (c.callPeaks) {
	std::vector<float> modfrac;
	std::vector<uint32_t> modpos;
	for (uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
	  // CpG sites
	  bool fwdCpG = true;
	  bool revCpG = true;
	  if (c.onlyCpG) {
	    fwdCpG = false;
	    revCpG = false;
	    if (pos + 1 < hdr->target_len[refIndex]) {
	      char b1 = std::toupper(seq[pos]);
	      char b2 = std::toupper(seq[pos+1]);
	      if ((b1 == 'C') && (b2 == 'G')) fwdCpG = true;
	    }
	    if (pos > 0) {
	      char b1 = std::toupper(seq[pos-1]);
	      char b2 = std::toupper(seq[pos]);
	      if ((b1 == 'C') && (b2 == 'G')) revCpG = true;
	    }
	  }
	  // Unstranded
	  if (c.peakStrand == '*') {
	    if ((fwdCpG) && (pos + 1 < hdr->target_len[refIndex])) {
	      uint32_t ucov = m_plus[pos] + h_plus[pos] + unmod_plus[pos] + m_minus[pos+1] + h_minus[pos+1] + unmod_minus[pos+1];
	      if (ucov >= c.minCov) {
		modpos.push_back(pos);
		if (c.modifiedRegions) {
		  if (c.modCode == 'h') modfrac.push_back(float(h_plus[pos] + h_minus[pos+1]) / float(ucov));
		  else modfrac.push_back(float(m_plus[pos] + m_minus[pos+1]) / float(ucov));
		} else {
		  if (c.modCode == 'h') modfrac.push_back(1.0 - float(h_plus[pos] + h_minus[pos+1]) / float(ucov));
		  else modfrac.push_back(1.0 - float(m_plus[pos] + m_minus[pos+1]) / float(ucov));
		}
	      }
	    }
	  } else if (c.peakStrand == '+') {
	    if (fwdCpG) {
	      uint32_t cov_plus = m_plus[pos] + h_plus[pos] + unmod_plus[pos];
	      if (cov_plus >= c.minCov) {
		modpos.push_back(pos);
		if (c.modifiedRegions) {
		  if (c.modCode == 'h') modfrac.push_back(float(h_plus[pos]) / float(cov_plus));
		  else modfrac.push_back(float(m_plus[pos]) / float(cov_plus));
		} else {
		  if (c.modCode == 'h') modfrac.push_back(1.0 - float(h_plus[pos]) / float(cov_plus));
		    else modfrac.push_back(1.0 - float(m_plus[pos]) / float(cov_plus));
		}
	      }
	    }
	  } else {
	    if (revCpG) {
	      uint32_t cov_minus = m_minus[pos] + h_minus[pos] + unmod_minus[pos];
	      if (cov_minus >= c.minCov) {
		modpos.push_back(pos);
		if (c.modifiedRegions) {
		  if (c.modCode == 'h') modfrac.push_back(float(h_minus[pos]) / float(cov_minus));
		  else modfrac.push_back(float(m_minus[pos]) / float(cov_minus));
		} else {
		  if (c.modCode == 'h') modfrac.push_back(1.0 - float(h_minus[pos]) / float(cov_minus));
		  else modfrac.push_back(1.0 - float(m_minus[pos]) / float(cov_minus));
		}
	      }
	    }
	  }
	}
	
	// call peaks
	if (modfrac.size() > 2 * std::max(c.peakwin, c.peakmedwin)) peakCalling(c, modfrac, modpos, scanRegions[refIndex]);
      }
      
      // Output percent modified per site
      if (c.hasBedFile) {
	// BED file
	for(uint32_t i = 0; i < scanRegions[refIndex].size(); ++i) {
	  if (scanRegions[refIndex][i].end <= hdr->target_len[refIndex]) {
	    // Aggregate mod counts
	    uint32_t reg_unmod_plus = 0;
	    uint32_t reg_unmod_minus = 0;
	    uint32_t reg_uncertain_plus = 0;
	    uint32_t reg_uncertain_minus = 0;
	    uint32_t reg_m_plus = 0;
	    uint32_t reg_m_minus = 0;
	    uint32_t reg_h_plus = 0;
	    uint32_t reg_h_minus = 0;
	    for (uint32_t pos = scanRegions[refIndex][i].start; pos < scanRegions[refIndex][i].end; ++pos) {
	      // CpG sites
	      bool fwdCpG = true;
	      bool revCpG = true;
	      if (c.onlyCpG) {
		fwdCpG = false;
		revCpG = false;
		if (pos + 1 < hdr->target_len[refIndex]) {
		  char b1 = std::toupper(seq[pos]);
		  char b2 = std::toupper(seq[pos+1]);
		  if ((b1 == 'C') && (b2 == 'G')) fwdCpG = true;
		}
		if (pos > 0) {
		  char b1 = std::toupper(seq[pos-1]);
		  char b2 = std::toupper(seq[pos]);
		  if ((b1 == 'C') && (b2 == 'G')) revCpG = true;
		}
	      }
	      if (fwdCpG) {
		reg_unmod_plus += unmod_plus[pos];
		reg_uncertain_plus += uncertain_plus[pos];
		reg_m_plus += m_plus[pos];
		reg_h_plus += h_plus[pos];
	      }
	      if (revCpG) {
		reg_unmod_minus += unmod_minus[pos];
		reg_uncertain_minus += uncertain_minus[pos];				  
		reg_m_minus += m_minus[pos];
		reg_h_minus += h_minus[pos];
	      }
	    }
	    // Unstranded
	    if (c.combineStrands) {
	      uint32_t unstranded_cov = reg_m_plus + reg_h_plus + reg_unmod_plus + reg_m_minus + reg_h_minus + reg_unmod_minus;
	      double pct_h = 0;
	      double pct_m = 0;
	      if (unstranded_cov > 0) {
		pct_h = double(reg_h_plus + reg_h_minus) / double(unstranded_cov);
		pct_m = double(reg_m_plus + reg_m_minus) / double(unstranded_cov);
	      }
	      out << hdr->target_name[refIndex] << "\t" << scanRegions[refIndex][i].start << "\t" << scanRegions[refIndex][i].end << "\t" << scanRegions[refIndex][i].id << "\th\t*\t" << (reg_h_plus + reg_h_minus) << "\t" << (reg_unmod_plus + reg_unmod_minus) << "\t" << (reg_m_plus + reg_m_minus) << "\t" << unstranded_cov << "\t" << (reg_uncertain_plus + reg_uncertain_minus) << "\t" << (boost::format("%1$.4f") % pct_h) << std::endl;
	      out << hdr->target_name[refIndex] << "\t" << scanRegions[refIndex][i].start << "\t" << scanRegions[refIndex][i].end << "\t" << scanRegions[refIndex][i].id << "\tm\t*\t" << (reg_m_plus + reg_m_minus) << "\t" << (reg_unmod_plus + reg_unmod_minus) << "\t" << (reg_h_plus + reg_h_minus) << "\t" << unstranded_cov << "\t" << (reg_uncertain_plus + reg_uncertain_minus) << "\t" << (boost::format("%1$.4f") % pct_m) << std::endl;
	    } else {
	      // Plus strand
	      uint32_t reg_cov_plus = reg_m_plus + reg_h_plus + reg_unmod_plus;
	      double pct_h_plus = 0;
	      double pct_m_plus = 0;
	      if (reg_cov_plus) {
		pct_h_plus = double(reg_h_plus) / double(reg_cov_plus);
		pct_m_plus = double(reg_m_plus) / double(reg_cov_plus);
	      }
	      out << hdr->target_name[refIndex] << "\t" << scanRegions[refIndex][i].start << "\t" << scanRegions[refIndex][i].end << "\t" << scanRegions[refIndex][i].id << "\th\t+\t" << reg_h_plus << "\t" << reg_unmod_plus << "\t" << reg_m_plus << "\t" << reg_cov_plus << "\t" << reg_uncertain_plus << "\t" << (boost::format("%1$.4f") % pct_h_plus) << std::endl;
	      out << hdr->target_name[refIndex] << "\t" << scanRegions[refIndex][i].start << "\t" << scanRegions[refIndex][i].end << "\t" << scanRegions[refIndex][i].id << "\tm\t+\t" << reg_m_plus << "\t" << reg_unmod_plus << "\t" << reg_h_plus << "\t" << reg_cov_plus << "\t" << reg_uncertain_plus << "\t" << (boost::format("%1$.4f") % pct_m_plus) << std::endl;
	      uint32_t reg_cov_minus = reg_m_minus + reg_h_minus + reg_unmod_minus;
	      double pct_h_minus = 0;
	      double pct_m_minus = 0;
	      if (reg_cov_minus) {
		pct_h_minus = double(reg_h_minus) / double(reg_cov_minus);
		pct_m_minus = double(reg_m_minus) / double(reg_cov_minus);
	      }
	      out << hdr->target_name[refIndex] << "\t" << scanRegions[refIndex][i].start << "\t" << scanRegions[refIndex][i].end << "\t" << scanRegions[refIndex][i].id << "\th\t-\t" << reg_h_minus << "\t" << reg_unmod_minus << "\t" << reg_m_minus << "\t" << reg_cov_minus << "\t" << reg_uncertain_minus << "\t" << (boost::format("%1$.4f") % pct_h_minus) << std::endl;
	      out << hdr->target_name[refIndex] << "\t" << scanRegions[refIndex][i].start << "\t" << scanRegions[refIndex][i].end << "\t" << scanRegions[refIndex][i].id << "\tm\t-\t" << reg_m_minus << "\t" << reg_unmod_minus << "\t" << reg_h_minus << "\t" << reg_cov_minus << "\t" << reg_uncertain_minus << "\t" << (boost::format("%1$.4f") % pct_m_minus) << std::endl;
	    }
	  }
	}
      } else {
	// No BED input
	for (uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
	  // CpG sites
	  bool fwdCpG = true;
	  bool revCpG = true;
	  if (c.onlyCpG) {
	    fwdCpG = false;
	    revCpG = false;
	    if (pos + 1 < hdr->target_len[refIndex]) {
	      char b1 = std::toupper(seq[pos]);
	      char b2 = std::toupper(seq[pos+1]);
	      if ((b1 == 'C') && (b2 == 'G')) fwdCpG = true;
	    }
	    if (pos > 0) {
	      char b1 = std::toupper(seq[pos-1]);
	      char b2 = std::toupper(seq[pos]);
	      if ((b1 == 'C') && (b2 == 'G')) revCpG = true;
	    }
	  }
	  
	  // Unstranded
	  if (c.combineStrands) {
	    if ((c.onlyCpG) && (fwdCpG) && (pos + 1 < hdr->target_len[refIndex])) {
	      uint32_t unstranded_cov = m_plus[pos] + h_plus[pos] + unmod_plus[pos] + m_minus[pos+1] + h_minus[pos+1] + unmod_minus[pos+1];
	      if (unstranded_cov > 0) {
		double pct_h = double(h_plus[pos] + h_minus[pos+1]) / double(unstranded_cov);
		double pct_m = double(m_plus[pos] + m_minus[pos+1]) / double(unstranded_cov);
		out << hdr->target_name[refIndex] << "\t" << pos << "\t" << (pos + 1) << "\t" << c.sampleName << "\th\t*\t" << (h_plus[pos] + h_minus[pos+1]) << "\t" << (unmod_plus[pos] + unmod_minus[pos+1]) << "\t" << (m_plus[pos] + m_minus[pos+1]) << "\t" << unstranded_cov << "\t" << (uncertain_plus[pos] + uncertain_minus[pos+1]) << "\t" << (boost::format("%1$.4f") % pct_h) << std::endl;
		out << hdr->target_name[refIndex] << "\t" << pos << "\t" << (pos + 1) << "\t" << c.sampleName << "\tm\t*\t" << (m_plus[pos] + m_minus[pos+1]) << "\t" << (unmod_plus[pos] + unmod_minus[pos+1]) << "\t" << (h_plus[pos] + h_minus[pos+1]) << "\t" << unstranded_cov << "\t" << (uncertain_plus[pos] + uncertain_minus[pos+1]) << "\t" << (boost::format("%1$.4f") % pct_m) << std::endl;
	      }
	    }
	  } else {
	    // Plus strand
	    if (fwdCpG) {
	      uint32_t cov_plus = m_plus[pos] + h_plus[pos] + unmod_plus[pos];
	      if (cov_plus > 0) {
		double pct_h_plus = double(h_plus[pos]) / double(cov_plus);
		double pct_m_plus = double(m_plus[pos]) / double(cov_plus);
		out << hdr->target_name[refIndex] << "\t" << pos << "\t" << (pos + 1) << "\t" << c.sampleName << "\th\t+\t" << h_plus[pos] << "\t" << unmod_plus[pos] << "\t" << m_plus[pos] << "\t" << cov_plus << "\t" << uncertain_plus[pos] << "\t" << (boost::format("%1$.4f") % pct_h_plus) << std::endl;
		out << hdr->target_name[refIndex] << "\t" << pos << "\t" << (pos + 1) << "\t" << c.sampleName << "\tm\t+\t" << m_plus[pos] << "\t" << unmod_plus[pos] << "\t" << h_plus[pos] << "\t" << cov_plus << "\t" << uncertain_plus[pos] << "\t" << (boost::format("%1$.4f") % pct_m_plus) << std::endl;
	      }
	    }
	      // Minus strand
	    if (revCpG) {
	      uint32_t cov_minus = m_minus[pos] + h_minus[pos] + unmod_minus[pos];
	      if (cov_minus > 0) {
		double pct_h_minus = double(h_minus[pos]) / double(cov_minus);
		double pct_m_minus = double(m_minus[pos]) / double(cov_minus);
		out << hdr->target_name[refIndex] << "\t" << pos << "\t" << (pos + 1) << "\t" << c.sampleName << "\th\t-\t" << h_minus[pos] << "\t" << unmod_minus[pos] << "\t" << m_minus[pos] << "\t" << cov_minus << "\t" << uncertain_minus[pos] << "\t" << (boost::format("%1$.4f") % pct_h_minus) << std::endl;
		out << hdr->target_name[refIndex] << "\t" << pos << "\t" << (pos + 1) << "\t" << c.sampleName << "\tm\t-\t" << m_minus[pos] << "\t" << unmod_minus[pos] << "\t" << h_minus[pos] << "\t" << cov_minus << "\t" << uncertain_minus[pos] << "\t" << (boost::format("%1$.4f") % pct_m_minus) << std::endl;
	      }
	    }
	  }
	}
      }
    
      // Free space
      if (seq != NULL) {
	free(seq);
	seq = NULL;
      }
    }

    // Close file
    if(c.outfile.string() != "-") of.close();

    // Clean-up
    if (fai) fai_destroy(fai);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);

    return true;
  }

  
}

#endif
