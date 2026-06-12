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
#include <unordered_set>
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

  // Methylation summary per SV (-1 = no call)
  struct MethylInfo {
    int32_t altSvStartL;   // ALT methylation % [0-100]
    int32_t altSvStartR;
    int32_t altSvRightL;
    int32_t altSvRightR;
    int32_t refSvStartL;   // REF methylation % [0-100]
    int32_t refSvStartR;
    int32_t refSvRightL;
    int32_t refSvRightR;
    int32_t mncStartL;     // unique CpG sites observed (combined alt+ref reads)
    int32_t mncStartR;
    int32_t mncRightL;
    int32_t mncRightR;
    int32_t mdpStartL;     // avg. read depth per CpG site
    int32_t mdpStartR;
    int32_t mdpRightL;
    int32_t mdpRightR;

    MethylInfo() : altSvStartL(-1), altSvStartR(-1), altSvRightL(-1), altSvRightR(-1), refSvStartL(-1), refSvStartR(-1), refSvRightL(-1), refSvRightR(-1), mncStartL(-1), mncStartR(-1), mncRightL(-1), mncRightR(-1), mdpStartL(-1), mdpStartR(-1), mdpRightL(-1), mdpRightR(-1) {}
  };

  // Accumulated raw counts per SV from reads
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
    // Unique CpG positions observed across all reads (alt+ref combined)
    std::unordered_set<int32_t> cpgStartL;
    std::unordered_set<int32_t> cpgStartR;
    std::unordered_set<int32_t> cpgRightL;
    std::unordered_set<int32_t> cpgRightR;

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
  collectMethylFromWindows(bam1_t* rec, std::vector<int8_t> const& methCall, std::vector<std::pair<int32_t,int32_t>> const& wins, std::vector<uint32_t>& methCounts, std::vector<uint32_t>& totCounts, std::vector<std::unordered_set<int32_t>>& cpgPos) {
    methCounts.assign(wins.size(), 0);
    totCounts.assign(wins.size(), 0);
    cpgPos.assign(wins.size(), std::unordered_set<int32_t>());
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
              cpgPos[wi].insert(rp);
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

  // Collect methylation from the insertion sequence using edlib alignment.
  template<typename TConfig>
  inline void
  collectMethylFromInsertionEdlib(TConfig const& c, bam1_t* rec, std::string const& readSeq, std::vector<int8_t> const& methCall, int32_t svStart, std::string const& consensus, int32_t consBp, int32_t insLen, uint32_t& methStartR, uint32_t& totStartR, std::unordered_set<int32_t>& cpgStartR, uint32_t& methRightL, uint32_t& totRightL, std::unordered_set<int32_t>& cpgRightL) {
    int32_t l = (int32_t)readSeq.size();
    if ((l == 0) || (insLen <= 0)) return;
    if ((consBp < 0) || ((consBp + insLen) > (int32_t)consensus.size())) return;

    bool readRev = (rec->core.flag & BAM_FREVERSE);

    // Locate the insertion
    int32_t insReadStart = -1;
    int32_t insReadLen   = 0;
    int32_t minOpLen = std::min(insLen / 2, c.methylWindow / 2);
    if (minOpLen < 1) minOpLen = 1;
    bool mapsBeforeBp = (rec->core.pos < svStart);
    {
      int32_t rp = rec->core.pos;
      int32_t sp = 0;
      const uint32_t* cigar = bam_get_cigar(rec);
      int32_t nCigar = (int32_t)rec->core.n_cigar;
      for (int32_t ci = 0; ci < nCigar; ++ci) {
        uint32_t op    = bam_cigar_op(cigar[ci]);
        int32_t  oplen = (int32_t)bam_cigar_oplen(cigar[ci]);
        if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF)) {
          rp += oplen; sp += oplen;
        } else if ((op == BAM_CDEL) || (op == BAM_CREF_SKIP)) {
          rp += oplen;
        } else if (op == BAM_CINS) {
          if ((rp == svStart) && (oplen >= minOpLen)) {
            insReadStart = sp;   
            insReadLen   = oplen;
            break;
          } else {
            sp += oplen;
          }
        } else if (op == BAM_CSOFT_CLIP) {
          if (oplen >= minOpLen) {
            bool isFirst = (ci == 0);
            bool isLast  = (ci == nCigar - 1);
            bool wantFirst = !mapsBeforeBp;
            if ((wantFirst && isFirst) || (!wantFirst && isLast)) {
              if (oplen > insReadLen) {  
                insReadStart = sp;
                insReadLen   = oplen;
              }
            }
          }
          sp += oplen;
        }
      }
    }
    if (insReadStart < 0) return;
    insReadLen = std::min(insReadLen, l - insReadStart);
    if (insReadLen <= 0) return;

    // extract insertion sequence
    std::string readIns = readSeq.substr(insReadStart, insReadLen);

    // consensus insertion sequence
    std::string consIns = consensus.substr(consBp, insLen);

    // align
    EdlibAlignResult aln = edlibAlign(readIns.c_str(), (int)readIns.size(), consIns.c_str(), (int)consIns.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    if ((aln.editDistance < 0) || (!aln.alignment) || (!aln.startLocations)) {
      edlibFreeAlignResult(aln);
      return;
    }

    // map coordinates
    std::vector<int32_t> consToRead(insLen, -1);
    {
      int32_t qi = 0;
      int32_t ti = aln.startLocations[0];
      for (int32_t ai = 0; (ai < aln.alignmentLength) && (qi < (int32_t)readIns.size()) && (ti < insLen); ++ai) {
        unsigned char op = aln.alignment[ai];
        if ((op == EDLIB_EDOP_MATCH) || (op == EDLIB_EDOP_MISMATCH)) {
          consToRead[ti] = qi; ++qi; ++ti;
        } else if (op == EDLIB_EDOP_INSERT) {
          consToRead[ti] = -1; ++ti;
        } else {
          ++qi;
        }
      }
    }
    edlibFreeAlignResult(aln);

    // methylation
    int32_t wlen = std::min(c.methylWindow, insLen);
    for (int32_t window = 0; window < 2; ++window) {
      int32_t winStart = (window == 0) ? 0 : (insLen - wlen);
      int32_t winEnd = (window == 0) ? wlen : insLen;
      uint32_t& meth_r = (window == 0) ? methStartR : methRightL;
      uint32_t& tot_r = (window == 0) ? totStartR : totRightL;
      std::unordered_set<int32_t>& cpg_r = (window == 0) ? cpgStartR : cpgRightL;

      for (int32_t k = winStart; k < (winEnd - 1); ++k) {
        // Identify CpG in consensus
        char ck  = std::toupper(static_cast<unsigned char>(consIns[k]));
        char ck1 = std::toupper(static_cast<unsigned char>(consIns[k + 1]));
        if ((ck != 'C') || (ck1 != 'G')) continue;

        // Fuzzy CpG search
        static const int32_t CPGTOL = 8;
        int8_t call = -1;
        if (!readRev) {
          // Forward read
          int32_t ri_center = consToRead[k];
          if (ri_center < 0) continue;
          for (int32_t delta = 0; delta <= CPGTOL && call < 0; ++delta) {
            for (int32_t sign : {0, -1, 1}) {
              if (delta == 0 && sign != 0) continue;
              if (delta != 0 && sign == 0) continue;
              int32_t ri = ri_center + delta * sign;
              if (ri < 0 || ri + 1 >= (int32_t)readIns.size()) continue;
              if (std::toupper(static_cast<unsigned char>(readIns[ri]))   != 'C') continue;
              if (std::toupper(static_cast<unsigned char>(readIns[ri+1])) != 'G') continue;
              int32_t idx = insReadStart + ri;
              if (idx < 0 || idx >= l) continue;
              if (methCall[idx] < 0) continue;
              call = methCall[idx];
              break;
            }
          }
        } else {
          // Reverse read
          if ((k + 1) >= insLen) continue;
          int32_t ri_g_center = consToRead[k + 1];
          if (ri_g_center < 0) continue;
          for (int32_t delta = 0; delta <= CPGTOL && call < 0; ++delta) {
            for (int32_t sign : {0, -1, 1}) {
              if (delta == 0 && sign != 0) continue;
              if (delta != 0 && sign == 0) continue;
              int32_t ri_g = ri_g_center + delta * sign;
              if (ri_g < 1 || ri_g >= (int32_t)readIns.size()) continue;
              if (std::toupper(static_cast<unsigned char>(readIns[ri_g-1])) != 'C') continue;
              if (std::toupper(static_cast<unsigned char>(readIns[ri_g]))   != 'G') continue;
              int32_t idx = l - 1 - (insReadStart + ri_g);
              if (idx < 0 || idx >= l) continue;
              if (methCall[idx] < 0) continue;
              call = methCall[idx];
              break;
            }
          }
        }

        if (call < 0) continue;

        // CpG key k, stable in consensus
        cpg_r.insert(k);
        ++tot_r;
        if (call == 1) ++meth_r;
      }
    }
  }

  // Compute methylation
  inline void
  finalizeMethylInfo(const MethylAccum& ma, MethylInfo& mi) {
    auto pct = [](uint32_t m, uint32_t t) -> int32_t {
      return (t > 0) ? (int32_t)std::round(100.0f * (float)m / (float)t) : -1;
    };
    auto mdp = [](uint32_t ta, uint32_t tr, int32_t cpg) -> int32_t {
      if (cpg <= 0) return -1;
      return (int32_t)std::round((float)(ta + tr) / (float)cpg);
    };
    mi.altSvStartL = pct(ma.altStartL_m, ma.altStartL_t);
    mi.altSvStartR = pct(ma.altStartR_m, ma.altStartR_t);
    mi.altSvRightL = pct(ma.altRightL_m, ma.altRightL_t);
    mi.altSvRightR = pct(ma.altRightR_m, ma.altRightR_t);
    mi.refSvStartL = pct(ma.refStartL_m, ma.refStartL_t);
    mi.refSvStartR = pct(ma.refStartR_m, ma.refStartR_t);
    mi.refSvRightL = pct(ma.refRightL_m, ma.refRightL_t);
    mi.refSvRightR = pct(ma.refRightR_m, ma.refRightR_t);
    mi.mncStartL = (int32_t)ma.cpgStartL.size();
    mi.mncStartR = (int32_t)ma.cpgStartR.size();
    mi.mncRightL = (int32_t)ma.cpgRightL.size();
    mi.mncRightR = (int32_t)ma.cpgRightR.size();
    mi.mdpStartL = mdp(ma.altStartL_t, ma.refStartL_t, mi.mncStartL);
    mi.mdpStartR = mdp(ma.altStartR_t, ma.refStartR_t, mi.mncStartR);
    mi.mdpRightL = mdp(ma.altRightL_t, ma.refRightL_t, mi.mncRightL);
    mi.mdpRightR = mdp(ma.altRightR_t, ma.refRightR_t, mi.mncRightR);
  }

  // Accumulate per-read methylation into MethylAccum for one SV
  template<typename TConfig>
  inline void
  accumulateMethyl(TConfig const& c, bam1_t* rec, std::vector<int8_t> const& methCall, StructuralVariantRecord const& sv, int32_t refIndex, int32_t chromLen, bool isAlt, std::vector<int32_t> const& candidates, MethylAccum& accum, std::string const& readSeq) {
    if (methCall.empty()) return;

    int32_t svStart = sv.svStart;
    int32_t svEnd = sv.svEnd;
    int32_t insLen = sv.insLen;
    bool isTrans = _translocation(sv.svt);
    bool isIns = (sv.svt == 4);
    bool isDel = (sv.svt == 2);

    // Determine which breakpoints this read overlaps on the current chromosome
    bool onStart = false;
    bool onEnd = false;
    for (int32_t cand : candidates) {
      if ((cand == svStart) && (sv.chr  == refIndex)) onStart = true;
      if ((cand == svEnd) && (sv.chr2 == refIndex)) onEnd = true;
    }
    if ((!onStart) && (!onEnd)) return;

    // Collect reference-coordinate windows
    // field: 0=StartL, 1=StartR, 2=RightL, 3=RightR
    std::vector<std::pair<int32_t,int32_t>> wins;
    std::vector<int32_t> field;

    if (onStart) {
      // StartL
      int32_t beg = std::max(0, svStart - c.methylWindow);
      int32_t end = svStart;
      if (end > beg) {
	wins.emplace_back(beg, end);
	field.push_back(0);
      }

      // StartR
      if (!isIns && !(isAlt && isDel)) {
        int32_t rbeg = svStart;
        int32_t rend = (!isTrans) ? std::min(svStart + c.methylWindow, svEnd) : std::min(svStart + c.methylWindow, chromLen);
        if (rend > rbeg) {
	  wins.emplace_back(rbeg, rend);
	  field.push_back(1);
	}
      }
    }

    if (onEnd) {
      // RightL
      if (!isIns && !(isAlt && isDel)) {
        int32_t lbeg = (!isTrans) ? std::max(svStart, svEnd - c.methylWindow) : std::max(0, svEnd - c.methylWindow);
        int32_t lend = svEnd;
        if (lend > lbeg) {
	  wins.emplace_back(lbeg, lend);
	  field.push_back(2);
	}
      }

      // RightR
      int32_t rbeg = svEnd;
      int32_t rend = std::min(svEnd + c.methylWindow, chromLen);
      if (rend > rbeg) {
	wins.emplace_back(rbeg, rend);
	field.push_back(3);
      }
    }

    // Collect
    if (!wins.empty()) {
      std::vector<uint32_t> methCounts, totCounts;
      std::vector<std::unordered_set<int32_t>> cpgPos;
      collectMethylFromWindows(rec, methCall, wins, methCounts, totCounts, cpgPos);
      for (std::size_t i = 0; i < field.size(); ++i) {
        uint32_t m = methCounts[i];
	uint32_t t = totCounts[i];
        if (isAlt) {
          switch (field[i]) {
            case 0: accum.altStartL_m += m; accum.altStartL_t += t; accum.cpgStartL.insert(cpgPos[i].begin(), cpgPos[i].end()); break;
            case 1: accum.altStartR_m += m; accum.altStartR_t += t; accum.cpgStartR.insert(cpgPos[i].begin(), cpgPos[i].end()); break;
            case 2: accum.altRightL_m += m; accum.altRightL_t += t; accum.cpgRightL.insert(cpgPos[i].begin(), cpgPos[i].end()); break;
            case 3: accum.altRightR_m += m; accum.altRightR_t += t; accum.cpgRightR.insert(cpgPos[i].begin(), cpgPos[i].end()); break;
          }
        } else {
          switch (field[i]) {
            case 0: accum.refStartL_m += m; accum.refStartL_t += t; accum.cpgStartL.insert(cpgPos[i].begin(), cpgPos[i].end()); break;
            case 1: accum.refStartR_m += m; accum.refStartR_t += t; accum.cpgStartR.insert(cpgPos[i].begin(), cpgPos[i].end()); break;
            case 2: accum.refRightL_m += m; accum.refRightL_t += t; accum.cpgRightL.insert(cpgPos[i].begin(), cpgPos[i].end()); break;
            case 3: accum.refRightR_m += m; accum.refRightR_t += t; accum.cpgRightR.insert(cpgPos[i].begin(), cpgPos[i].end()); break;
          }
        }
      }
    }

    // For INS, alignment
    if (isIns && isAlt && (insLen > 0) && (onStart || onEnd) && !readSeq.empty() && !sv.consensus.empty()) {
      collectMethylFromInsertionEdlib(c, rec, readSeq, methCall, svStart, sv.consensus, sv.consBp, insLen, accum.altStartR_m, accum.altStartR_t, accum.cpgStartR, accum.altRightL_m, accum.altRightL_t, accum.cpgRightL);
    }
  }
  
}

#endif
