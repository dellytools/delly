/*
============================================================================
DELLY: Structural variant discovery by integrated PE mapping and SR analysis
============================================================================
Copyright (C) 2012 Tobias Rausch

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Tobias Rausch (rausch@embl.de)
============================================================================
*/

#ifndef SPLIT_H
#define SPLIT_H

#include <iostream>
#include "gotoh.h"
#include "needle.h"

namespace torali
{



  // Deletions
  template<typename TString, typename TSvRecord, typename TAIndex, typename TPosition>
  inline bool
  _coordTransform(TString const&, TSvRecord const& sv, TAIndex, TAIndex, TAIndex rStart, TAIndex rEnd, TPosition& finalGapStart, TPosition& finalGapEnd, SVType<DeletionTag>) {
    TAIndex annealed = sv.svStartEnd - sv.svStartBeg;
    if ((rStart >= annealed) || (rEnd < annealed)) return false;
    finalGapStart = sv.svStartBeg + rStart;
    finalGapEnd = sv.svEndBeg + (rEnd - annealed);
    return true;
  }


  template<typename TString, typename TSvRecord, typename TAIndex, typename TPosition>
  inline bool
  _coordTransform(TString const&, TSvRecord const& sv, TAIndex, TAIndex, TAIndex rStart, TAIndex rEnd, TPosition& finalGapStart, TPosition& finalGapEnd, SVType<DuplicationTag>) {
    TAIndex annealed = sv.svEndEnd - sv.svEndBeg;
    if ((rStart >= annealed) || (rEnd < annealed)) return false;
    finalGapStart = sv.svStartBeg + (rEnd - annealed);
    finalGapEnd = sv.svEndBeg + rStart;
    return true;
  }

  template<typename TString, typename TSvRecord, typename TAIndex, typename TPosition>
  inline bool
  _coordTransform(TString const& ref, TSvRecord const& sv, TAIndex cStart, TAIndex cEnd, TAIndex rStart, TAIndex rEnd, TPosition& finalGapStart, TPosition& finalGapEnd, SVType<InversionTag>) {
    return false;
  }

  template<typename TString, typename TSvRecord, typename TAIndex, typename TPosition>
  inline bool
  _coordTransform(TString const& ref, TSvRecord const& sv, TAIndex cStart, TAIndex cEnd, TAIndex rStart, TAIndex rEnd, TPosition& finalGapStart, TPosition& finalGapEnd, SVType<TranslocationTag>) {
    return false;
  }

  template<typename TString, typename TSvRecord, typename TAIndex, typename TPosition>
  inline bool
  _coordTransform(TString const& ref, TSvRecord const& sv, TAIndex cStart, TAIndex cEnd, TAIndex rStart, TAIndex rEnd, TPosition& finalGapStart, TPosition& finalGapEnd, SVType<InsertionTag>) {
    return false;
  }


  template<typename TAlign, typename TAIndex>
  inline bool
  _findSplit(TAlign const& align, TAIndex& cStart, TAIndex& cEnd, TAIndex& rStart, TAIndex& rEnd) {
    // Initializiation
    cStart=0;
    cEnd=0;
    rStart=0;
    rEnd=0;

    // Find longest internal gap
    TAIndex refIndex=0;
    TAIndex varIndex=0;
    TAIndex gapStartRefIndex=0;
    TAIndex gapEndRefIndex=0;
    TAIndex gapStartVarIndex=0;
    TAIndex gapEndVarIndex=0;
    bool inGap=false;
    for(TAIndex j = 0; j<align.shape()[1]; ++j) {
      if (align[0][j] != '-') ++varIndex;
      if (align[1][j] != '-') ++refIndex;
      // Internal gap?
      if (((align[0][j] == '-') || (align[1][j] == '-')) && (refIndex>0) && (varIndex>0)) {
	if (!inGap) {
	  gapStartRefIndex=refIndex;
	  gapEndRefIndex=refIndex + 1;
	  gapStartVarIndex=varIndex;
	  gapEndVarIndex=varIndex + 1;
	  inGap = true;
	} else {
	  gapEndRefIndex=refIndex + 1;
	  gapEndVarIndex=varIndex + 1;
	}
      } else {
	if ((inGap) && ((gapEndRefIndex - gapStartRefIndex) > (rEnd - rStart))) {
	  rStart=gapStartRefIndex;
	  rEnd=gapEndRefIndex;
	  cStart=gapStartVarIndex;
	  cEnd=gapEndVarIndex;
	}
	inGap=false;
      }
    }
    return (rEnd > rStart);
  }

  template<typename TConfig, typename TStructuralVariantRecord, typename TTag>
  inline bool
  alignConsensus(TConfig const& c, TStructuralVariantRecord& sv, std::string const& svRefStr, SVType<TTag> svType) {

    // Create input alignments
    typedef boost::multi_array<char, 2> TAlign;
    typedef typename TAlign::index TAIndex;
    TAlign cons;
    cons.resize(boost::extents[1][sv.consensus.size()]);
    TAIndex ind = 0;
    for(typename std::string::const_iterator str = sv.consensus.begin(); str != sv.consensus.end(); ++str) cons[0][ind++] = *str;
    TAlign ref;
    ref.resize(boost::extents[1][svRefStr.size()]);
    ind = 0;
    for(typename std::string::const_iterator str = svRefStr.begin(); str != svRefStr.end(); ++str) ref[0][ind++] = *str;

    // Consensus to reference alignment
    TAlign alignFwd;
    int score = gotoh(cons, ref, alignFwd);

    // Scoring - deduce leading, trailing and one internal gap (assuming go=10, ge=1, match=5)
    score += (alignFwd.shape()[1] - sv.consensus.size()) + 3 * 10;
    int maxScore = sv.consensus.size() * 5;
    if (score < 0) score = 0;
    double quality = 1;
    if (score < maxScore) quality = (double) score / (double) maxScore;

    // Debug consensus to reference alignment
    for(TAIndex i = 0; i<alignFwd.shape()[0]; ++i) {
      for(TAIndex j = 0; j<alignFwd.shape()[1]; ++j) {
	std::cerr << alignFwd[i][j];
      }
      std::cerr << std::endl;
    }
    std::cerr << "Alignment score: " << score << " (Quality: " << quality << ")" << std::endl;
    std::cerr << std::endl;

    // Check quality
    if (quality < ((double) c.flankQuality / 100.0)) return false;

    // Check breakpoint
    TAIndex cStart, cEnd, rStart, rEnd;
    if (!_findSplit(alignFwd, cStart, cEnd, rStart, rEnd)) return false;

    // Check flanking alignment length
    if ((cStart < c.minimumFlankSize) || ((sv.consensus.size() - cEnd) < c.minimumFlankSize)) return false;

    // Get the start and end of the structural variant
    unsigned int finalGapStart = 0;
    unsigned int finalGapEnd = 0;
    if (!_coordTransform(svRefStr, sv, cStart, cEnd, rStart, rEnd, finalGapStart, finalGapEnd, svType)) return false;

    std::cerr << "Coordinate transform" << std::endl;
    std::cerr << sv.svStart << ',' << sv.svEnd << std::endl;
    std::cerr << finalGapStart << ',' << finalGapEnd << std::endl;

    // Set breakpoint & quality
    sv.precise=true;
    sv.svStart=finalGapStart;
    sv.svEnd=finalGapEnd;
    sv.srAlignQuality=quality;

    return true;
  }



}

#endif
