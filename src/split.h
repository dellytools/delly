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
  _validSoftClip(bam1_t* rec, int& clipSize, int& splitPoint, bool& leadingSC, unsigned short qualCut) {
    // Check read-length
    if (rec->core.l_qseq < 35) return false;

    // Check for single soft-clip
    unsigned int numSoftClip = 0;
    uint32_t* cigar = bam_get_cigar(rec);
    for (unsigned int i = 0; i < rec->core.n_cigar; ++i) {
      if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	++numSoftClip;
	clipSize = bam_cigar_oplen(cigar[i]);
      }
    }
    if (numSoftClip != 1) return false;

    // Check clip size
    if (clipSize <= (int32_t) (log10(rec->core.l_qseq) * 10)) return false;
    
    // Get quality vector
    typedef std::vector<uint8_t> TQuality;
    TQuality quality;
    quality.resize(rec->core.l_qseq);
    uint8_t* qualptr = bam_get_qual(rec);
    for (int i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
    
    // Get soft-clips
    unsigned int alen = 0;
    unsigned int lastIns = 0;
    unsigned int meanQuality = 0;
    for (unsigned int i = 0; i < rec->core.n_cigar; ++i) {
      if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
	alen += bam_cigar_oplen(cigar[i]) + lastIns;
	lastIns = 0;
      } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	lastIns = bam_cigar_oplen(cigar[i]);   // Only add if followed by 'M'
      } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	if (!alen) leadingSC = true;
	else leadingSC = false;
	splitPoint = rec->core.pos + alen;
	unsigned int qualSum = 0;
	for(unsigned int i = alen; i < (alen+clipSize); ++i) qualSum += quality[i];
	meanQuality = qualSum / clipSize;
      }
    }
    //std::cerr << clipSize << ',' << meanQuality << ',' << splitPoint << std::endl;
    return (meanQuality >= (unsigned int) qualCut);
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
      if ((rStart >= annealed) || (rEnd < annealed)) return false;
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

  template<typename TGap, typename TTag>
  inline bool
  _checkSVGap(TGap const refGap, TGap const oldRefGap, TGap const, TGap const, SVType<TTag>) {
    return (refGap > oldRefGap);
  }

  template<typename TGap>
  inline bool
  _checkSVGap(TGap const, TGap const, TGap const varGap, TGap const oldVarGap, SVType<InsertionTag>) {
    return (varGap > oldVarGap);
  }


  template<typename TAlign, typename TAIndex, typename TFloat, typename TSVType>
  inline bool
  _findSplit(TAlign const& align, TAIndex& cStart, TAIndex& cEnd, TAIndex& rStart, TAIndex& rEnd, TAIndex& gS, TAIndex& gE, TFloat& percId, TSVType svt) {
    // Initializiation
    gS=0;
    gE=0;
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
    TAIndex a1 = 0;
    TAIndex a2 = 0;
    bool inGap=false;
    for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
      if (align[0][j] != '-') ++varIndex;
      if (align[1][j] != '-') ++refIndex;
      // Internal gap?
      if (((align[0][j] == '-') || (align[1][j] == '-')) && (refIndex>0) && (varIndex>0)) {
	if (!inGap) {
	  gapStartRefIndex=refIndex;
	  gapStartVarIndex=varIndex;
	  a1 = j;
	  inGap = true;
	}
	gapEndRefIndex=refIndex + 1;
	gapEndVarIndex=varIndex + 1;
	a2 = j;
      } else {
	if ((inGap) && (_checkSVGap((gapEndRefIndex - gapStartRefIndex), (rEnd - rStart), (gapEndVarIndex - gapStartVarIndex), (cEnd - cStart), svt))) {
	  rStart=gapStartRefIndex;
	  rEnd=gapEndRefIndex;
	  cStart=gapStartVarIndex;
	  cEnd=gapEndVarIndex;
	  gS = a1;
	  gE = a2;
	}
	inGap=false;
      }
    }

    // Find percent identity
    if (rEnd > rStart) {
      TAIndex rI=0;
      TAIndex vI=0;
      uint32_t gapMM = 0;
      uint32_t mm = 0;
      uint32_t ma = 0;
      bool inGap=false;
      for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
	//std::cerr << j << ',' << ma << ',' << mm << std::endl;
	if (align[0][j] != '-') ++vI;
	if (align[1][j] != '-') ++rI;
	// Internal gap?
	if ((align[0][j] == '-') || (align[1][j] == '-')) {
	  if ((rI>0) && (vI>0)) {
	    if (!inGap) {
	      inGap = true;
	      gapMM = 0;
	    }
	    if ((rI < rStart) || (rI >= rEnd)) gapMM += 1;
	  }
	} else {
	  if (inGap) {
	    mm += gapMM;
	    inGap=false;
	  }
	  if (align[0][j] == align[1][j]) ma += 1;
	  else mm += 1;
	}
      }
      percId = (TFloat) ma / (TFloat) (ma + mm);
    }
    return (rEnd > rStart);
  }

  template<typename TAlign, typename TAIndex, typename TSVType>
  inline bool
  _findSplit(TAlign const& align, TAIndex& cStart, TAIndex& cEnd, TAIndex& rStart, TAIndex& rEnd, TSVType svt) {
    TAIndex alignJ1 = 0;
    TAIndex alignJ2 = 0;
    double percId = 0;
    return _findSplit(align, cStart, cEnd, rStart, rEnd, alignJ1, alignJ2, percId, svt);
  }

  template<typename TAlign, typename TAIndex, typename TLength>
  inline void
  _findHomology(TAlign const& align, TAIndex const gS, TAIndex const gE, TLength& homLeft, TLength& homRight) {
    int32_t mmThres = 1;
    if (align[1][gS] == '-') {
      // Insertion
      int32_t mismatch = 0;
      int32_t offset = 0;
      for(TAIndex i = 0; i < gS; ++i, ++homLeft) {
	if (align[1][gS-i-1] != align[0][gE-i-offset]) ++mismatch;
	if (mismatch > mmThres) {
	  // Try 1bp insertion
	  if (!offset) {
	    if (align[1][gS-i-1] == align[0][gE-i-(++offset)]) {
	      --mismatch;
	      continue;
	    }
	  }
	  break;
	}
      }
      mismatch = 0;
      offset = 0;
      for(TAIndex i = 0; i < (TAIndex) (align.shape()[1] - gE - 1); ++i, ++homRight) {
	if (align[0][gS+i] != align[1][gE+i+1]) ++mismatch;
	if (mismatch > mmThres) {
	  // Try 1bp insertion
	  if (!offset) {
	    if (align[0][gS+i+(++offset)] == align[1][gE+i+1]) {
	      --mismatch;
	      continue;
	    }
	  }
	  break;
	}
      }
    } else if (align[0][gS] == '-') {
      // Deletion
      int32_t mismatch = 0;
      int32_t offset = 0;
      for(TAIndex i = 0; i < gS; ++i, ++homLeft) {
	if (align[0][gS-i-1] != align[1][gE-i-offset]) ++mismatch;
	if (mismatch > mmThres) {
	  // Try 1bp deletion
	  if (!offset) {
	    if (align[0][gS-i-1] == align[1][gE-i-(++offset)]) {
	      --mismatch;
	      continue;
	    }
	  }
	  break;
	}
      }
      mismatch = 0;
      offset = 0;
      for(TAIndex i = 0; i < (TAIndex) (align.shape()[1] - gE - 1); ++i, ++homRight) {
	if (align[1][gS+i] != align[0][gE+i+1]) ++mismatch;
	if (mismatch > mmThres) {
	  // Try 1bp deletion
	  if (!offset) {
	    if (align[1][gS+i+(++offset)] == align[0][gE+i+1]) {
	      --mismatch;
	      continue;
	    }
	  }
	  break;
	}
      }
    }
  }

  template<typename TAlign, typename TTag>
  inline bool
    _consRefAlignment(std::string const& cons, std::string const& svRefStr, TAlign& aln, SVType<TTag>)
  {
    AlignConfig<true, false> semiglobal;
    DnaScore<int> lnsc(5, -4, -4, -4);
    bool reNeedle = longNeedle(cons, svRefStr, aln, semiglobal, lnsc);
    return reNeedle;
  }

  template<typename TAlign>
  inline bool
    _consRefAlignment(std::string const& cons, std::string const& svRefStr, TAlign& aln, SVType<InsertionTag>)
  {
    typedef typename TAlign::index TAIndex;
    AlignConfig<false, true> semiglobal;
    DnaScore<int> lnsc(5, -4, -4, -4);
    bool reNeedle = longNeedle(svRefStr, cons, aln, semiglobal, lnsc);
    for(TAIndex j = 0; j < (TAIndex) aln.shape()[1]; ++j) {
      char tmp = aln[0][j];
      aln[0][j] = aln[1][j];
      aln[1][j] = tmp;
    }	
    return reNeedle;
  }

  template<typename TConfig, typename TStructuralVariantRecord, typename TTag>
  inline bool
  alignConsensus(TConfig const& c, TStructuralVariantRecord& sv, std::string const& svRefStr, SVType<TTag> svType) {
    if (sv.consensus.size() < 35) return false;
    typedef boost::multi_array<char, 2> TAlign;
    typedef typename TAlign::index TAIndex;

    // Consensus to reference alignment
    TAlign alignFwd;
    if (!_consRefAlignment(sv.consensus, svRefStr, alignFwd, svType)) return false;

    // Check breakpoint
    TAIndex cStart, cEnd, rStart, rEnd, gS, gE;
    double quality = 0;
    if (!_findSplit(alignFwd, cStart, cEnd, rStart, rEnd, gS, gE, quality, svType)) return false;

    // Debug consensus to reference alignment
    //for(TAIndex i = 0; i<alignFwd.shape()[0]; ++i) {
    //for(TAIndex j = 0; j<alignFwd.shape()[1]; ++j) {
    //std::cerr << alignFwd[i][j];
    //}
    //std::cerr << std::endl;
    //}
    //std::cerr << "Flanking alignment percent identity: " << quality << std::endl;
    //std::cerr << cStart << ',' << cEnd << ',' << rStart << ',' << rEnd << std::endl;
    //std::cerr << std::endl;

    // Check quality
    if (quality < c.flankQuality) return false;

    // Find homology
    if (!_validSRAlignment(cStart, cEnd, rStart, rEnd, svType)) return false;
    int32_t homLeft = 0;
    int32_t homRight = 0;
    _findHomology(alignFwd, gS, gE, homLeft, homRight);

    // Check flanking alignment length
    if ((homLeft + c.minimumFlankSize > (int32_t) cStart) || ( (int32_t) (sv.consensus.size() - cEnd) < homRight + c.minimumFlankSize)) return false;
    if ((homLeft + c.minimumFlankSize > (int32_t) rStart) || ( (int32_t) (svRefStr.size() - rEnd) < homRight + c.minimumFlankSize)) return false;

    // Get the start and end of the structural variant
    unsigned int finalGapStart = 0;
    unsigned int finalGapEnd = 0;
    if (c.technology == "illumina") {
      if (!_coordTransform(svRefStr, sv, rStart, rEnd, finalGapStart, finalGapEnd, svType)) return false;
    } else if (c.technology == "pacbio") {
      int32_t rs = std::max(0, sv.svStart - (int32_t) (sv.consensus.size()));
      finalGapStart = rs + rStart - 1;
      finalGapEnd = rs + rEnd - 1;
    }

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
