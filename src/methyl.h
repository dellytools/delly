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

  static const int32_t METHYL_WINDOW = 2000;
  static const uint8_t METHYL_PROB_THRESHOLD = 128;   // Min. ML prob (128/256, 50%)

  struct MethylInfo {
    float altSvStartL;
    float altSvStartR;
    float altSvRightL;
    float altSvRightR;
    float refSvStartL;
    float refSvStartR;
    float refSvRightL;
    float refSvRightR;

    MethylInfo() : altSvStartL(-1.0f), altSvStartR(-1.0f), altSvRightL(-1.0f), altSvRightR(-1.0f), refSvStartL(-1.0f), refSvStartR(-1.0f), refSvRightL(-1.0f), refSvRightR(-1.0f) {}
  };

  // Methyl struct
  struct MethylAccum {
    uint32_t altStartL_m;
    uint32_t altStartL_t;
    uint32_t altStartR_m;
    uint32_t altStartR_t;
    uint32_t altRightL_m;
    uint32_t altRightL_t;
    uint32_t altRightR_m;
    uint32_t altRightR_t;
    uint32_t refStartL_m;
    uint32_t refStartL_t;
    uint32_t refStartR_m;
    uint32_t refStartR_t;
    uint32_t refRightL_m;
    uint32_t refRightL_t;
    uint32_t refRightR_m;
    uint32_t refRightR_t;

    MethylAccum() : altStartL_m(0), altStartL_t(0), altStartR_m(0), altStartR_t(0), altRightL_m(0), altRightL_t(0), altRightR_m(0), altRightR_t(0), refStartL_m(0), refStartL_t(0), refStartR_m(0), refStartR_t(0), refRightL_m(0), refRightL_t(0), refRightR_m(0), refRightR_t(0) {}
  };


  struct ModHit {
    int32_t pos;
    char code;
    uint8_t prob;
    bool rev;
    char base;

    ModHit(int32_t const p, char const c, uint8_t const r, bool const s, char const b) : pos(p), code(c), prob(r), rev(s), base(b) {}
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

  inline double
  _medianVectorD(std::vector<double> myVector) {
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
  
  //   1  = methylated 5mC (prob >= probModTh)
  //   0  = unmodified C (explicit in MM or implicit)
  //  -1  = position is not a C in fwdseq (no call)
  inline bool
  buildMethylCalls(bam1_t* rec, uint8_t probModTh, std::vector<int8_t>& methCall) {
    int32_t l = rec->core.l_qseq;
    methCall.assign(l, -1);

    uint8_t* mm_aux = bam_aux_get(rec, "MM");
    if ((!mm_aux) || (*mm_aux != 'Z')) return false;
    bool readRev = (rec->core.flag & BAM_FREVERSE);

    // Sequence
    const uint8_t* seqptr = bam_get_seq(rec);
    std::string sequence(l, 'N');
    for (int32_t i = 0; i < l; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

    // Fwd seq
    std::string fwdseq = sequence;
    if (readRev) reverseComplement(fwdseq);
    std::unordered_map<char, std::vector<int32_t>> basepos;
    for (int32_t i = 0; i < l; ++i) basepos[std::toupper(static_cast<unsigned char>(fwdseq[i]))].push_back(i);

    // Parse MM tag
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
      while ((idx < tok.size()) && (tok[idx] != ',')) {
        char ch = tok[idx++];
        if ((ch == '?') || (ch == '.')) hasSkip = true;
        else if (std::isalpha(static_cast<unsigned char>(ch))) mod_codes.push_back(ch);
      }
      if (hasSkip)
        for (char c : mod_codes)
          if ((c == 'm') || (c == 'M')) mCHasSkipStatus = true;
      if ((idx < tok.size()) && (tok[idx] == ',')) {
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

    // Pre-mark all C positions as unmodified
    if (!mCHasSkipStatus) {
      for (int32_t i = 0; i < l; ++i)
        if (std::toupper(static_cast<unsigned char>(fwdseq[i])) == 'C') methCall[i] = 0;
    }

    // Parse ML tag for mod probabilities
    uint8_t* ml_aux = bam_aux_get(rec, "ML");
    if ((ml_aux) && (*ml_aux == 'B')) {
      char subtype = *reinterpret_cast<char*>(ml_aux + 1);
      if (subtype == 'C') {
        const uint8_t* p = ml_aux + 2;
        int32_t n = (int32_t)(p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24));
        const uint8_t* data = p + 4;
        int32_t assign = std::min<int32_t>(n, (int32_t)modhits.size());
        for (int32_t i = 0; i < assign; ++i) modhits[i].prob = data[i];
      }
    }

    // Map 5mC modhits to fwdseq positions
    for (const auto& mh : modhits) {
      if ((mh.code != 'm') && (mh.code != 'M')) continue; // only 5mC
      char ub = std::toupper(static_cast<unsigned char>(mh.base));
      char target_base = mh.rev ? complement_base(ub) : ub;
      auto it = basepos.find(target_base);
      if (it == basepos.end()) continue;
      const auto& occs = it->second;
      if ((mh.pos < 0) || (static_cast<std::size_t>(mh.pos) >= occs.size())) continue;
      int32_t fwdPos = occs[mh.pos];
      methCall[fwdPos] = (mh.prob >= probModTh) ? 1 : 0;
    }

    return true;
  }

  // Collect methylation counts
  inline void
  collectMethylFromWindows(bam1_t* rec, std::vector<int8_t> const& methCall, std::vector<std::pair<int32_t,int32_t>> const& wins, std::vector<uint32_t>& methCounts, std::vector<uint32_t>& totCounts) {
    methCounts.assign(wins.size(), 0);
    totCounts.assign(wins.size(), 0);
    if ((wins.empty()) || (methCall.empty())) return;

    bool readRev = (rec->core.flag & BAM_FREVERSE);
    int32_t l = rec->core.l_qseq;
    // Maximum reference end
    int32_t maxEnd = 0;
    for (const auto& w : wins) maxEnd = std::max(maxEnd, w.second);
    
    int32_t rp = rec->core.pos;
    int32_t sp = 0;
    const uint32_t* cigar = bam_get_cigar(rec);
    for (uint32_t ci = 0; ci < rec->core.n_cigar; ++ci) {
      uint32_t op = bam_cigar_op(cigar[ci]);
      int32_t oplen = (int32_t)bam_cigar_oplen(cigar[ci]);
      if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF)) {
        for (int32_t k = 0; k < oplen; ++k, ++rp, ++sp) {
          if (rp >= maxEnd) return;
          int32_t fwdPos = readRev ? (l - sp - 1) : sp;
          if ((fwdPos < 0) || (fwdPos >= l)) continue;
          int8_t call = methCall[fwdPos];
          if (call < 0) continue; // not a C in fwdseq
          for (std::size_t wi = 0; wi < wins.size(); ++wi) {
            if ((rp >= wins[wi].first) && (rp < wins[wi].second)) {
              ++totCounts[wi];
              if (call == 1) ++methCounts[wi];
            }
          }
        }
      } else if ((op == BAM_CDEL) || (op == BAM_CREF_SKIP)) {
        rp += oplen;
        if (rp >= maxEnd) return;
      } else if ((op == BAM_CINS) || (op == BAM_CSOFT_CLIP)) {
        sp += oplen;
      }
    }
  }

  // Collect methylation counts from the insertion sequence
  inline void
  collectMethylFromInsertion(bam1_t* rec, std::vector<int8_t> const& methCall, int32_t insRefPos, int32_t insLen, bool useFirst, int32_t windowLen, uint32_t& methCount, uint32_t& totCount) {
    methCount = 0;
    totCount = 0;
    if ((methCall.empty()) || (insLen <= 0) || (windowLen <= 0)) return;

    bool readRev = (rec->core.flag & BAM_FREVERSE);
    int32_t l = rec->core.l_qseq;
    int32_t rp = rec->core.pos;
    int32_t sp = 0;
    const uint32_t* cigar = bam_get_cigar(rec);
    for (uint32_t ci = 0; ci < rec->core.n_cigar; ++ci) {
      uint32_t op = bam_cigar_op(cigar[ci]);
      int32_t oplen = (int32_t)bam_cigar_oplen(cigar[ci]);
      if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF)) {
        for (int32_t k = 0; k < oplen; ++k, ++rp, ++sp) {
          if (rp > insRefPos) return; // passed insertion site
        }
      } else if ((op == BAM_CDEL) || (op == BAM_CREF_SKIP)) {
        rp += oplen;
        if (rp > insRefPos) return;
      } else if (op == BAM_CINS) {
        if (rp == insRefPos) {
          // Found the insertion at the right reference position
          int32_t wlen = std::min(windowLen, oplen);
          int32_t spBeg = useFirst ? sp : (sp + oplen - wlen);
          int32_t spEnd = useFirst ? (sp + wlen) : (sp + oplen);
          for (int32_t si = spBeg; ((si < spEnd) && (si < l)); ++si) {
            int32_t fwdPos = readRev ? (l - si - 1) : si;
            if ((fwdPos < 0) || (fwdPos >= l)) continue;
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
    }
  }

  // Compute final methylation fractions
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

  // Accumulate per-read methylation into MethylAccum for one SV
  inline void
  accumulateMethyl(bam1_t* rec, std::vector<int8_t> const& methCall, StructuralVariantRecord const& sv, int32_t refIndex, int32_t chromLen, bool isAlt, std::vector<int32_t> const& candidates, MethylAccum& accum) {
    if (methCall.empty()) return;

    int32_t svStart = sv.svStart;
    int32_t svEnd   = sv.svEnd;
    int32_t insLen  = sv.insLen;
    bool isTrans = _translocation(sv.svt);
    bool isIns   = (sv.svt == 4);

    // Determine which breakpoints this read overlaps on the current chromosome
    bool onStart = false;
    bool onEnd = false;
    for (int32_t cand : candidates) {
      if ((cand == svStart) && (sv.chr  == refIndex)) onStart = true;
      if ((cand == svEnd) && (sv.chr2 == refIndex)) onEnd = true;
    }
    if ((!onStart) && (!onEnd)) return;

    // Collect reference-coordinate windows
    // 0=StartL, 1=StartR, 2=RightL, 3=RightR
    std::vector<std::pair<int32_t,int32_t>> wins;
    std::vector<int32_t> field;
    bool doInsStartR = false;
    bool doInsRightL = false;

    if (onStart) {
      // StartL: up to METHYL_WINDOW bp left of svStart
      int32_t beg = std::max(0, svStart - METHYL_WINDOW);
      int32_t end = svStart;
      if (end > beg) {
	wins.emplace_back(beg, end);
	field.push_back(0);
      }

      // StartR: depends on SV type
      if ((isAlt) && (isIns) && (insLen > 0)) {
        doInsStartR = true; // first windowLen bp of inserted sequence
      } else {
        int32_t rbeg = svStart;
        int32_t rend;
        if ((!isTrans) && (!isIns)) rend = std::min(svStart + METHYL_WINDOW, svEnd); // cap at svEnd (internal sequence only)
        else rend = std::min(svStart + METHYL_WINDOW, chromLen);
        if (rend > rbeg) {
	  wins.emplace_back(rbeg, rend);
	  field.push_back(1);
	}
      }
    }

    if (onEnd) {
      // RightL: depends on SV type
      if ((isAlt) && (isIns) && (insLen > 0)) {
        doInsRightL = true; // last windowLen bp of inserted sequence
      } else {
        int32_t lbeg;
        if ((!isTrans) && (!isIns)) lbeg = std::max(svStart, svEnd - METHYL_WINDOW); // floor at svStart (internal sequence only)
        else lbeg = std::max(0, svEnd - METHYL_WINDOW);
        int32_t lend = svEnd;
        if (lend > lbeg) {
	  wins.emplace_back(lbeg, lend);
	  field.push_back(2);
	}
      }

      // RightR: up to METHYL_WINDOW bp right of svEnd
      int32_t rbeg = svEnd;
      int32_t rend = std::min(svEnd + METHYL_WINDOW, chromLen);
      if (rend > rbeg) {
	wins.emplace_back(rbeg, rend);
	field.push_back(3);
      }
    }

    // Collect from reference windows in a single CIGAR walk
    if (!wins.empty()) {
      std::vector<uint32_t> methCounts, totCounts;
      collectMethylFromWindows(rec, methCall, wins, methCounts, totCounts);
      for (std::size_t i = 0; i < field.size(); ++i) {
        uint32_t m = methCounts[i];
	uint32_t t = totCounts[i];
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

    // Collect from inserted sequence
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
  
}

#endif
