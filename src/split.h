#ifndef SPLIT_H
#define SPLIT_H

#include <iostream>
#include "gotoh.h"
#include "needle.h"

namespace torali
{

  struct AlignDescriptor {
    int32_t cStart;
    int32_t cEnd;
    int32_t rStart;
    int32_t rEnd;
    int32_t homLeft;
    int32_t homRight;
    float percId;
    
    AlignDescriptor() : cStart(0), cEnd(0), rStart(0), rEnd(0), homLeft(0), homRight(0), percId(0) {}
  };

  template<typename TBPoint>
  inline void
  _adjustOrientation(std::string& sequence, TBPoint bpPoint, int32_t const svt) {
    if (_translocation(svt)) {
      uint8_t ct = _getSpanOrientation(svt);
      if (((ct==0) && (bpPoint)) || ((ct==1) && (!bpPoint))) reverseComplement(sequence);
    } else {
      if (svt == 0) {
	if (bpPoint) reverseComplement(sequence);
      } else if (svt == 1) {
	if (!bpPoint) reverseComplement(sequence);
      }
    }
  }

  inline bool
  _largeClipFraction(int32_t const clipSize, int32_t const qlen, int32_t const svt) {
    if (svt == 2) return (((double) clipSize / (double) qlen) > 0.5);
    else return false;
  }

  inline bool 
  _validSoftClip(bam1_t* rec, int32_t& clipSize, int32_t& splitPoint, bool& leadingSC, unsigned short qualCut, int32_t const svt) {
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

    // Check clip fraction
    if (_largeClipFraction(clipSize, rec->core.l_qseq, svt)) return false;
    
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

  template<typename TBPoint>
  inline bool
  _validSCOrientation(TBPoint bpPoint, bool leadingSC, int32_t const svt) {
    if (_translocation(svt)) {
      uint8_t ct = _getSpanOrientation(svt);
      if (ct == 0) return (!leadingSC);
      else if (ct == 1) return leadingSC;
      else if (ct == 2) {
	if (((!bpPoint) && (!leadingSC)) || ((bpPoint) && (leadingSC))) return true;
	else return false;
      } else if (ct == 3) {
	if (((!bpPoint) && (leadingSC)) || ((bpPoint) && (!leadingSC))) return true;
	else return false;
      } else return false;
    } else {
      if (svt == 0) return (!leadingSC);
      else if (svt == 1) return leadingSC;
      else if (svt == 2) {
	if (((!bpPoint) && (!leadingSC)) || ((bpPoint) && (leadingSC))) return true;
	else return false;
      } else if (svt == 3) {
	if (((!bpPoint) && (leadingSC)) || ((bpPoint) && (!leadingSC))) return true;
	else return false;
      } else if (svt == 4) {
	if (((!bpPoint) && (!leadingSC)) || ((bpPoint) && (leadingSC))) return true;
	else return false;
      }
    }
    return false;
  }

  // Deletions
  template<typename TSeq, typename TSVRecord, typename TRef>
  inline std::string
  _getSVRef(TSeq const* const ref, TSVRecord const& svRec, TRef const refIndex, int32_t const svt) {
    if (_translocation(svt)) {
      uint8_t ct = _getSpanOrientation(svt);
      if (svRec.chr==refIndex) {
	if ((ct==0) || (ct == 2)) return boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svStartEnd)) + svRec.part1;
	else if (ct == 1) {
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
	  return refPart + svRec.part1;
	} else return svRec.part1 + boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svStartEnd));
      } else {
	// chr2
	if (ct==0) {
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
    } else {
      if (svt == 2) {
	return boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svStartEnd)) + boost::to_upper_copy(std::string(ref + svRec.svEndBeg, ref + svRec.svEndEnd));
      } else if (svt == 4) {
	return boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svEndEnd));
      } else if (svt == 3) {
	return boost::to_upper_copy(std::string(ref + svRec.svEndBeg, ref + svRec.svEndEnd)) + boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svStartEnd));
      } else if (svt == 0) {
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
      } else if (svt == 1) {
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
    return "";
  }


  template<typename TString, typename TSvRecord, typename TAlignDescriptor, typename TPosition>
  inline bool
  _coordTransform(TString const& ref, TSvRecord const& sv, TAlignDescriptor const& ad, TPosition& finalGapStart, TPosition& finalGapEnd, int32_t svt) {
    if (_translocation(svt)) {
      uint8_t ct = _getSpanOrientation(svt);
      if (ct == 0) {
	int32_t annealed = sv.svStartEnd - sv.svStartBeg;
	if ((ad.rStart >= annealed) || (ad.rEnd < annealed)) return false;
	finalGapStart = sv.svStartBeg + ad.rStart;
	finalGapEnd = sv.svEndBeg + (ref.size() - ad.rEnd) + 1;
      } else if (ct == 1) {
	int32_t annealed = sv.svStartEnd - sv.svStartBeg;
	if ((ad.rStart >= annealed) || (ad.rEnd < annealed)) return false;
	finalGapStart = sv.svStartBeg + (annealed - ad.rStart) + 1;
	finalGapEnd = sv.svEndBeg + (ad.rEnd - annealed);
      } else if (ct == 2) {
	int32_t annealed = sv.svStartEnd - sv.svStartBeg;
	if ((ad.rStart >= annealed) || (ad.rEnd < annealed)) return false;
	finalGapStart = sv.svStartBeg + ad.rStart;
	finalGapEnd = sv.svEndBeg + (ad.rEnd - annealed);
      } else if (ct == 3) {
	int32_t annealed = sv.svEndEnd - sv.svEndBeg;
	if ((ad.rStart >= annealed) || (ad.rEnd < annealed)) return false;
	finalGapStart = sv.svStartBeg + (ad.rEnd - annealed);
	finalGapEnd = sv.svEndBeg + ad.rStart;
      } else return false;
      return true;
    } else {
      if (svt == 2) {
	int32_t annealed = sv.svStartEnd - sv.svStartBeg;
	if ((ad.rStart >= annealed) || (ad.rEnd < annealed)) return false;
	finalGapStart = sv.svStartBeg + ad.rStart;
	finalGapEnd = sv.svEndBeg + (ad.rEnd - annealed);
	return true;
      } else if (svt == 3) {
	int32_t annealed = sv.svEndEnd - sv.svEndBeg;
	if ((ad.rStart >= annealed) || (ad.rEnd < annealed)) return false;
	finalGapStart = sv.svStartBeg + (ad.rEnd - annealed);
	finalGapEnd = sv.svEndBeg + ad.rStart;
	return true;
      } else if (svt == 0) {
	int32_t annealed = sv.svStartEnd - sv.svStartBeg;
	if ((ad.rStart >= annealed) || (ad.rEnd < annealed)) return false;
	finalGapStart = sv.svStartBeg + ad.rStart;
	finalGapEnd = sv.svEndBeg + (ref.size() - ad.rEnd) + 1;
	return true;
      } else if (svt == 1) {
	int32_t annealed = sv.svStartEnd - sv.svStartBeg;
	if ((ad.rStart >= annealed) || (ad.rEnd < annealed)) return false;
	finalGapStart = sv.svStartBeg + (annealed - ad.rStart) + 1;
	finalGapEnd = sv.svEndBeg + (ad.rEnd - annealed);
	return true;
      } else if (svt == 4) {
	finalGapStart = sv.svStartBeg + ad.rStart;
	finalGapEnd = sv.svStartBeg + ad.rEnd;
	return true;
      }
    }
    return true;
  }


  template<typename TPos>
  inline bool
  _validSRAlignment(TPos const cStart, TPos const cEnd, TPos const rStart, TPos const rEnd, int32_t const svt)
  {
    if (svt == 4) return (((rEnd - rStart) < 5) && ((cEnd - cStart) > 15));
    else return (((cEnd - cStart) < 5) && ((rEnd - rStart) > 15));
  }

  template<typename TGap>
  inline bool
  _checkSVGap(TGap const refGap, TGap const oldRefGap, TGap const varGap, TGap const oldVarGap, int32_t const svt) {
    if (svt == 4) return (varGap > oldVarGap);
    else return (refGap > oldRefGap);
  }

  template<typename TAlignDescriptor>
  inline void
  _findHomology(std::string const& consensus, std::string const& svRefStr, TAlignDescriptor& ad, int32_t const svt) {
    if (svt == 4) {
      ad.homRight = longestHomology(consensus.substr(ad.cStart), svRefStr.substr(ad.rEnd -1), -1);
      std::string preC = consensus.substr(0, ad.cEnd - 1);
      std::string preR = svRefStr.substr(0, ad.rStart);
      std::reverse(preC.begin(), preC.end());
      std::reverse(preR.begin(), preR.end());
      ad.homLeft = longestHomology(preC, preR, -1);
    } else {
      ad.homRight = longestHomology(consensus.substr(ad.cEnd - 1), svRefStr.substr(ad.rStart), -1);
      std::string preC = consensus.substr(0, ad.cStart);
      std::string preR = svRefStr.substr(0, ad.rEnd - 1);
      std::reverse(preC.begin(), preC.end());
      std::reverse(preR.begin(), preR.end());
      ad.homLeft = longestHomology(preC, preR, -1);
    }
  }

  template<typename TAlign, typename TAIndex, typename TFloat>
  inline void
  _percentIdentity(TAlign const& align, TAIndex const gS, TAIndex const gE, TFloat& percId) {
    // Find percent identity
    bool varSeen = false;
    bool refSeen = false;
    uint32_t gapMM = 0;
    uint32_t mm = 0;
    uint32_t ma = 0;
    bool inGap=false;
    for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
      if ((j < gS) || (j > gE)) {
	if (align[0][j] != '-') varSeen = true;
	if (align[1][j] != '-') refSeen = true;
	// Internal gap?
	if ((align[0][j] == '-') || (align[1][j] == '-')) {
	  if ((refSeen) && (varSeen)) {
	    if (!inGap) {
	      inGap = true;
	      gapMM = 0;
	    }
	    gapMM += 1;
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
    }
    percId = (TFloat) ma / (TFloat) (ma + mm);
  }


  template<typename TConfig, typename TAlign, typename TAlignDescriptor>
  inline bool
  _findSplit(TConfig const& c, std::string const& consensus, std::string const& svRefStr, TAlign const& align, TAlignDescriptor& ad, int32_t const svt) {
    // Initializiation
    int32_t gS=0;
    int32_t gE=0;

    // Find longest internal gap
    int32_t refIndex=0;
    int32_t varIndex=0;
    int32_t gapStartRefIndex=0;
    int32_t gapStartVarIndex=0;
    int32_t a1 = 0;
    bool inGap=false;
    for(int32_t j = 0; j < (int32_t) align.shape()[1]; ++j) {
      if (align[0][j] != '-') ++varIndex;
      if (align[1][j] != '-') ++refIndex;
      // Internal gap?
      if (((align[0][j] == '-') || (align[1][j] == '-')) && (refIndex>0) && (varIndex>0)) {
	if (!inGap) {
	  gapStartVarIndex = (align[0][j] != '-') ? (varIndex - 1) : varIndex;
	  gapStartRefIndex = (align[1][j] != '-') ? (refIndex - 1) : refIndex;
	  a1 = j;
	  inGap = true;
	}
      } else {
	if ((inGap) && (_checkSVGap((refIndex - gapStartRefIndex), (ad.rEnd - ad.rStart), (varIndex - gapStartVarIndex), (ad.cEnd - ad.cStart), svt))) {
	  ad.rStart=gapStartRefIndex;
	  ad.rEnd=refIndex;
	  ad.cStart=gapStartVarIndex;
	  ad.cEnd=varIndex;
	  gS = a1;
	  gE = j - 1;
	}
	inGap=false;
      }
    }
    if (ad.rEnd <= ad.rStart) return false;

    // Is this a valid split-read alignment?
    if (!_validSRAlignment(ad.cStart, ad.cEnd, ad.rStart, ad.rEnd, svt)) return false;

    // Check percent identity
    _percentIdentity(align, gS, gE, ad.percId);
    if (ad.percId < c.flankQuality) return false;

    // Find homology
    _findHomology(consensus, svRefStr, ad, svt);
    
    // Check flanking alignment length
    if ((ad.homLeft + c.minimumFlankSize > ad.cStart) || ( varIndex < ad.cEnd + ad.homRight + c.minimumFlankSize)) return false;
    if ((ad.homLeft + c.minimumFlankSize > ad.rStart) || ( refIndex < ad.rEnd + ad.homRight + c.minimumFlankSize)) return false;

    // Valid split-read alignment
    return true;
  }

  template<typename TAlign>
  inline bool
  _consRefAlignment(std::string const& cons, std::string const& svRefStr, TAlign& aln, int32_t const svt) {
    AlignConfig<true, false> semiglobal;
    DnaScore<int> lnsc(5, -4, -4, -4);
    bool reNeedle = false;
    if (svt == 4) {
      reNeedle = longNeedle(svRefStr, cons, aln, semiglobal, lnsc);
      for(uint32_t j = 0; j < aln.shape()[1]; ++j) {
	char tmp = aln[0][j];
	aln[0][j] = aln[1][j];
	aln[1][j] = tmp;
      }	
    } else {
      reNeedle = longNeedle(cons, svRefStr, aln, semiglobal, lnsc);
    }
    return reNeedle;
  }

  template<typename TConfig>
  inline bool
  alignConsensus(TConfig const& c, bam_hdr_t* hdr, char const* seq, char const* sndSeq, StructuralVariantRecord& sv) {
    if ( (int32_t) sv.consensus.size() < (2 * c.minimumFlankSize + sv.insLen)) return false;
    
    // Get reference slice
    Breakpoint bp(sv);
    if (sv.svt ==4) {
      int32_t bufferSpace = std::max((int32_t) ((sv.consensus.size() - sv.insLen) / 3), c.minimumFlankSize);
      _initBreakpoint(hdr, bp, bufferSpace, sv.svt);
    } else _initBreakpoint(hdr, bp, sv.consensus.size(), sv.svt);
    if (bp.chr != bp.chr2) bp.part1 = _getSVRef(sndSeq, bp, bp.chr2, sv.svt);
    std::string svRefStr = _getSVRef(seq, bp, bp.chr, sv.svt);

    // Consensus to reference alignment
    typedef boost::multi_array<char, 2> TAlign;
    TAlign align;
    if (!_consRefAlignment(sv.consensus, svRefStr, align, sv.svt)) return false;

    // Debug consensus to reference alignment
    //std::cerr << "Consensus-to-Reference alignment" << std::endl;
    //typedef typename TAlign::index TAIndex;
    //for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
    //for(TAIndex j = 0; j< (TAIndex) align.shape()[1]; ++j) {
    //std::cerr << align[i][j];
    //}
    //std::cerr << std::endl;
    //}
    
    // Check breakpoint
    AlignDescriptor ad;
    if (!_findSplit(c, sv.consensus, svRefStr, align, ad, sv.svt)) return false;

    // Get the start and end of the structural variant
    unsigned int finalGapStart = 0;
    unsigned int finalGapEnd = 0;
    if (!_coordTransform(svRefStr, bp, ad, finalGapStart, finalGapEnd, sv.svt)) return false;
	
    sv.precise=true;
    sv.svStart=finalGapStart;
    sv.svEnd=finalGapEnd;
    sv.srAlignQuality = ad.percId;
    sv.insLen=ad.cEnd - ad.cStart - 1;
    sv.homLen=std::max(0, ad.homLeft + ad.homRight - 2);
    int32_t ci_wiggle = std::max(ad.homLeft, ad.homRight);
    sv.ciposlow = -ci_wiggle;
    sv.ciposhigh = ci_wiggle;
    sv.ciendlow = -ci_wiggle;
    sv.ciendhigh = ci_wiggle;
    return true;
  }



}

#endif
