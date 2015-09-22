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

  template<typename TBPoint, typename TCT>
  inline void
  _adjustOrientation(std::string&, TBPoint, TCT, SVType<DeletionTag>) 
  {
    //Nop
  }

  template<typename TBPoint, typename TCT>
  inline void
  _adjustOrientation(std::string&, TBPoint, TCT, SVType<InsertionTag>) 
  {
    //Nop
  }

  template<typename TBPoint, typename TCT>
  inline void
  _adjustOrientation(std::string&, TBPoint, TCT, SVType<DuplicationTag>) 
  {
    //Nop
  }

  template<typename TBPoint, typename TCT>
  inline void
  _adjustOrientation(std::string& sequence, TBPoint bpPoint, TCT ct, SVType<InversionTag>) 
  {
    if (((!ct) && (bpPoint)) || ((ct) && (!bpPoint))) reverseComplement(sequence);
  }

  template<typename TBPoint, typename TCT>
  inline void
  _adjustOrientation(std::string& sequence, TBPoint bpPoint, TCT ct, SVType<TranslocationTag>) 
  {
    if (((ct==0) && (bpPoint)) || ((ct==1) && (!bpPoint))) reverseComplement(sequence);
  }


  inline bool 
  _validSoftClip(bam1_t* rec, int& clipSize, int& splitPoint, bool& leadingSC) {
    // Check read-length
    if (rec->core.l_qseq < 35) return false;

    // Check for soft-clips
    bool hasSoftClip = false;
    uint32_t* cigar = bam_get_cigar(rec);
    for (unsigned int i = 0; i < rec->core.n_cigar; ++i) 
      if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) hasSoftClip = true;
    if (!hasSoftClip) return false;
    
    // Get quality vector
    typedef std::vector<uint8_t> TQuality;
    TQuality quality;
    quality.resize(rec->core.l_qseq);
    uint8_t* qualptr = bam_get_qual(rec);
    for (int i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
    
    // Get soft-clips
    unsigned int alen = 0;
    unsigned int numSoftClip = 0;
    unsigned int meanQuality = 0;
    for (unsigned int i = 0; i < rec->core.n_cigar; ++i) {
      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CINS)) alen += bam_cigar_oplen(cigar[i]);
      else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	if (!alen) leadingSC = true;
	else leadingSC = false;
	++numSoftClip;
	clipSize = bam_cigar_oplen(cigar[i]);
	splitPoint = rec->core.pos + alen;
	unsigned int qualSum = 0;
	for(unsigned int i = alen; i < (alen+clipSize); ++i) qualSum += quality[i];
	meanQuality = qualSum / clipSize;
      }
    }
    //std::cerr << numSoftClip << ',' << clipSize << ',' << meanQuality << ',' << splitPoint << std::endl;
    return ((numSoftClip==1) && (meanQuality>=20));
  }

  template<typename TBPoint, typename TCT>
  inline bool
  _validSCOrientation(TBPoint bpPoint, bool leadingSC, TCT, SVType<DeletionTag>) 
  {
    if (((!bpPoint) && (!leadingSC)) || ((bpPoint) && (leadingSC))) return true;
    else return false;
  }

  template<typename TBPoint, typename TCT>
  inline bool
  _validSCOrientation(TBPoint bpPoint, bool leadingSC, TCT, SVType<InsertionTag>) 
  {
    if (((!bpPoint) && (!leadingSC)) || ((bpPoint) && (leadingSC))) return true;
    else return false;
  }

  template<typename TBPoint, typename TCT>
  inline bool
  _validSCOrientation(TBPoint bpPoint, bool leadingSC, TCT, SVType<DuplicationTag>) 
  {
    if (((!bpPoint) && (leadingSC)) || ((bpPoint) && (!leadingSC))) return true;
    else return false;
  }

  template<typename TBPoint, typename TCT>
  inline bool
  _validSCOrientation(TBPoint, bool leadingSC, TCT ct, SVType<InversionTag>) 
  {
    return (ct ? leadingSC : (!leadingSC));
  }

  template<typename TBPoint, typename TCT>
  inline bool
  _validSCOrientation(TBPoint bpPoint, bool leadingSC, TCT ct, SVType<TranslocationTag>) 
  {
    if (ct == 0) return (!leadingSC);
    else if (ct == 1) return leadingSC;
    else if (ct == 2) {
      if (((!bpPoint) && (!leadingSC)) || ((bpPoint) && (leadingSC))) return true;
      else return false;
    } 
    else if (ct == 3) {
      if (((!bpPoint) && (leadingSC)) || ((bpPoint) && (!leadingSC))) return true;
      else return false;
    } else return false;
  }


  // Deletions
  template<typename TSeq, typename TSVRecord, typename TRef>
  inline std::string
  _getSVRef(TSeq const* const ref, TSVRecord const& svRec, TRef const, SVType<DeletionTag>) {
    return boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svStartEnd)) + boost::to_upper_copy(std::string(ref + svRec.svEndBeg, ref + svRec.svEndEnd));
  }

  // Insertions
  template<typename TSeq, typename TSVRecord, typename TRef>
  inline std::string
  _getSVRef(TSeq const* const ref, TSVRecord const& svRec, TRef const, SVType<InsertionTag>) {
    return boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svEndEnd));
  }

  // Duplications
  template<typename TSeq, typename TSVRecord, typename TRef>
  inline std::string
  _getSVRef(TSeq const* const ref, TSVRecord const& svRec, TRef const, SVType<DuplicationTag>) {
    return boost::to_upper_copy(std::string(ref + svRec.svEndBeg, ref + svRec.svEndEnd)) + boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svStartEnd));
  }

  // Inversions
  template<typename TSeq, typename TSVRecord, typename TRef>
  inline std::string
  _getSVRef(TSeq const* const ref, TSVRecord const& svRec, TRef const, SVType<InversionTag>) {
    if (!svRec.ct) {
      std::string strEnd=boost::to_upper_copy(std::string(ref + svRec.svEndBeg, ref + svRec.svEndEnd));
      std::string strRevComp=strEnd;
      std::string::reverse_iterator itR = strEnd.rbegin();
      std::string::reverse_iterator itREnd = strEnd.rend();
      for(unsigned int i = 0; itR!=itREnd; ++itR, ++i) {
	switch (*itR) {
	case 'A': strRevComp[i]='T'; break;
	case 'C': strRevComp[i]='G'; break;
	case 'G': strRevComp[i]='C'; break;
	case 'T': strRevComp[i]='A'; break;
	case 'N': strRevComp[i]='N'; break;
	default: break;
	}
      }
      return boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svStartEnd)) + strRevComp;
    } else {
      std::string strStart=boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svStartEnd));
      std::string strRevComp=strStart;
      std::string::reverse_iterator itR = strStart.rbegin();
      std::string::reverse_iterator itREnd = strStart.rend();
      for(unsigned int i = 0; itR!=itREnd; ++itR, ++i) {
	switch (*itR) {
	case 'A': strRevComp[i]='T'; break;
	case 'C': strRevComp[i]='G'; break;
	case 'G': strRevComp[i]='C'; break;
	case 'T': strRevComp[i]='A'; break;
	case 'N': strRevComp[i]='N'; break;
	default: break;
	}
      }
      return strRevComp + boost::to_upper_copy(std::string(ref + svRec.svEndBeg, ref + svRec.svEndEnd));
    }
  }

  // Translocations
  template<typename TSeq, typename TSVRecord, typename TRef>
  inline std::string
  _getSVRef(TSeq const* const ref, TSVRecord const& svRec, TRef const refIndex, SVType<TranslocationTag>) {
    if (svRec.chr==refIndex) {
      if ((svRec.ct==0) || (svRec.ct == 2)) return boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svStartEnd)) + svRec.consensus;
      else if (svRec.ct == 1) {
	std::string strEnd=boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svStartEnd));
	std::string refPart=strEnd;
	std::string::reverse_iterator itR = strEnd.rbegin();
	std::string::reverse_iterator itREnd = strEnd.rend();
	for(unsigned int i = 0; itR!=itREnd; ++itR, ++i) {
	  switch (*itR) {
	  case 'A': refPart[i]='T'; break;
	  case 'C': refPart[i]='G'; break;
	  case 'G': refPart[i]='C'; break;
	  case 'T': refPart[i]='A'; break;
	  case 'N': refPart[i]='N'; break;
	  default: break;
	  }
	}
	return refPart + svRec.consensus;
      } else return svRec.consensus + boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svStartEnd));
    } else {
      // chr2
      if (svRec.ct==0) {
	std::string strEnd=boost::to_upper_copy(std::string(ref + svRec.svEndBeg, ref + svRec.svEndEnd));
	std::string refPart=strEnd;
	std::string::reverse_iterator itR = strEnd.rbegin();
	std::string::reverse_iterator itREnd = strEnd.rend();
	for(unsigned int i = 0; itR!=itREnd; ++itR, ++i) {
	  switch (*itR) {
	  case 'A': refPart[i]='T'; break;
	  case 'C': refPart[i]='G'; break;
	  case 'G': refPart[i]='C'; break;
	  case 'T': refPart[i]='A'; break;
	  case 'N': refPart[i]='N'; break;
	  default: break;
	  }
	}
	return refPart;
      } else return boost::to_upper_copy(std::string(ref + svRec.svEndBeg, ref + svRec.svEndEnd));
    }
  }


  // Deletions
  template<typename TString, typename TSvRecord, typename TAIndex, typename TPosition>
  inline bool
  _coordTransform(TString const&, TSvRecord const& sv, TAIndex rStart, TAIndex rEnd, TPosition& finalGapStart, TPosition& finalGapEnd, SVType<DeletionTag>) {
    TAIndex annealed = sv.svStartEnd - sv.svStartBeg;
    if ((rStart >= annealed) || (rEnd < annealed)) return false;
    finalGapStart = sv.svStartBeg + rStart - 1;
    finalGapEnd = sv.svEndBeg + (rEnd - annealed);
    return true;
  }

  // Duplications
  template<typename TString, typename TSvRecord, typename TAIndex, typename TPosition>
  inline bool
  _coordTransform(TString const&, TSvRecord const& sv, TAIndex rStart, TAIndex rEnd, TPosition& finalGapStart, TPosition& finalGapEnd, SVType<DuplicationTag>) {
    TAIndex annealed = sv.svEndEnd - sv.svEndBeg;
    if ((rStart >= annealed) || (rEnd < annealed)) return false;
    finalGapStart = sv.svStartBeg + (rEnd - annealed);
    finalGapEnd = sv.svEndBeg + rStart - 1;
    return true;
  }

  // Inversion
  template<typename TString, typename TSvRecord, typename TAIndex, typename TPosition>
  inline bool
  _coordTransform(TString const& ref, TSvRecord const& sv, TAIndex rStart, TAIndex rEnd, TPosition& finalGapStart, TPosition& finalGapEnd, SVType<InversionTag>) {
    TAIndex annealed = sv.svStartEnd - sv.svStartBeg;
    if ((rStart >= annealed) || (rEnd < annealed)) return false;
    if (!sv.ct) {
      finalGapStart = sv.svStartBeg + rStart - 1;
      finalGapEnd = sv.svEndBeg + (ref.size() - rEnd) + 1;
    } else {
      finalGapStart = sv.svStartBeg + (annealed - rStart) + 2;
      finalGapEnd = sv.svEndBeg + (rEnd - annealed);
    } 
    return true;
  }

  // Translocation
  template<typename TString, typename TSvRecord, typename TAIndex, typename TPosition>
  inline bool
  _coordTransform(TString const& ref, TSvRecord const& sv, TAIndex rStart, TAIndex rEnd, TPosition& finalGapStart, TPosition& finalGapEnd, SVType<TranslocationTag>) {
    if (sv.ct == 0) {
      TAIndex annealed = sv.svStartEnd - sv.svStartBeg;
      if ((rStart >= annealed) || (rEnd < annealed)) return false;
      finalGapStart = sv.svStartBeg + rStart -1;
      finalGapEnd = sv.svEndBeg + (ref.size() - rEnd) + 1;
    }
    else if (sv.ct == 1) {
      TAIndex annealed = sv.svStartEnd - sv.svStartBeg;
      finalGapStart = sv.svStartBeg + (annealed - rStart) + 2;
      finalGapEnd = sv.svEndBeg + (rEnd - annealed);
    }
    else if (sv.ct == 2) {
      TAIndex annealed = sv.svStartEnd - sv.svStartBeg;
      if ((rStart >= annealed) || (rEnd < annealed)) return false;
      finalGapStart = sv.svStartBeg + rStart - 1;
      finalGapEnd = sv.svEndBeg + (rEnd - annealed);
    } 
    else if (sv.ct == 3) {
      TAIndex annealed = sv.svEndEnd - sv.svEndBeg;
      if ((rStart >= annealed) || (rEnd < annealed)) return false;
      finalGapStart = sv.svStartBeg + (rEnd - annealed);
      finalGapEnd = sv.svEndBeg + rStart - 1;
    }
    else return false;
    return true;
  }

  template<typename TString, typename TSvRecord, typename TAIndex, typename TPosition>
  inline bool
  _coordTransform(TString const&, TSvRecord const& sv, TAIndex rStart, TAIndex rEnd, TPosition& finalGapStart, TPosition& finalGapEnd, SVType<InsertionTag>) {
    finalGapStart = sv.svStartBeg + rStart;
    finalGapEnd = sv.svStartBeg + rEnd;
    return true;
  }


  template<typename TPos, typename TTag>
  inline bool
  _validSRAlignment(TPos const cStart, TPos const cEnd, TPos const rStart, TPos const rEnd, SVType<TTag>)
  {
    return (((cEnd - cStart) < 5) && ((rEnd - rStart) > 15));
  }

  template<typename TPos>
  inline bool
  _validSRAlignment(TPos const cStart, TPos const cEnd, TPos const rStart, TPos const rEnd, SVType<InsertionTag>)
  {
    return (((rEnd - rStart) < 5) && ((cEnd - cStart) > 15));
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
    for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
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
    if (sv.consensus.size() < 35) return false;

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
    AlignConfig<true, false> semiglobal;
    DnaScore<int> sc(5, -4, -5 * c.minimumFlankSize, 0); // Don't penalize the split, ge=0 and make sure we have aligned segments > c.minimumFlankSize
    int score = gotoh(cons, ref, alignFwd, semiglobal, sc);
    score += 5 * c.minimumFlankSize; // Increase the score by allowing one internal gap

    // Check breakpoint
    TAIndex cStart, cEnd, rStart, rEnd;
    if (!_findSplit(alignFwd, cStart, cEnd, rStart, rEnd)) return false;
    if (!_validSRAlignment(cStart, cEnd, rStart, rEnd, svType)) return false;

    // Check quality
    double quality = (double) ((score < 0) ? 0 : score ) / (double) ( (sv.consensus.size() - (cEnd - cStart - 1)) * 5);
    if (quality > 1) quality = 1;
    if (quality < ((double) c.flankQuality / 100.0)) return false;

    // Debug consensus to reference alignment
    //for(TAIndex i = 0; i<alignFwd.shape()[0]; ++i) {
    //for(TAIndex j = 0; j<alignFwd.shape()[1]; ++j) {
    //std::cerr << alignFwd[i][j];
    //}
    //std::cerr << std::endl;
    //}
    //std::cerr << "Alignment score: " << score << " (Quality: " << quality << ")" << std::endl;
    //std::cerr << std::endl;

    // Check flanking alignment length
    if ((cStart < c.minimumFlankSize) || ((sv.consensus.size() - cEnd) < c.minimumFlankSize)) return false;

    // Get the start and end of the structural variant
    unsigned int finalGapStart = 0;
    unsigned int finalGapEnd = 0;
    if (!_coordTransform(svRefStr, sv, rStart, rEnd, finalGapStart, finalGapEnd, svType)) return false;

    // Set breakpoint & quality
    sv.precise=true;
    sv.svStart=finalGapStart;
    sv.svEnd=finalGapEnd;
    sv.srAlignQuality=quality;
    sv.insLen=cEnd - cStart - 1;
    return true;
  }



}

#endif
