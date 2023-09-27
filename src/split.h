#ifndef SPLIT_H
#define SPLIT_H

#include <iostream>

//#include "bindings/cpp/WFAligner.hpp"

#include "edlib.h"
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


  /*
  template<typename TWFAligner>
  inline void
  printWfaAlignment(std::string const& pattern, std::string const& text, TWFAligner& aligner) {
    std::string cigar = aligner.getAlignmentCigar();

    std::cerr << "Alignment score: " << aligner.getAlignmentScore() << std::endl;

    // pattern
    uint32_t pidx = 0;
    for (uint32_t j = 0; j < cigar.size(); ++j) {
      if (cigar[j] == 'I') std::cerr << '-';
      else std::cerr << pattern[pidx++];
    }
    std::cerr << std::endl;

    // text
    uint32_t tidx = 0;
    for (uint32_t j = 0; j < cigar.size(); ++j) {
      if (cigar[j] == 'D') std::cerr << '-';
      else std::cerr <<	text[tidx++];
    }
    std::cerr << std::endl;
  }
  */
  
  
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

  template<typename TConfig, typename TSeq, typename TSVRecord, typename TRef>
  inline std::string
  _getSVRef(TConfig const& c, TSeq const* const ref, TSVRecord const& svRec, TRef const refIndex, int32_t const svt) {
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
	if (svRec.svEnd - svRec.svStart <= c.indelsize) return boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svEndEnd));
	else return boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svStartEnd)) + boost::to_upper_copy(std::string(ref + svRec.svEndBeg, ref + svRec.svEndEnd));
      } else if (svt == 4) {
	return boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svEndEnd));
      } else if (svt == 3) {
	return boost::to_upper_copy(std::string(ref + svRec.svEndBeg, ref + svRec.svEndEnd)) + boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svStartEnd));
      } else if (svt == 0) {
	std::string strEnd;
	if ((svRec.svEnd - svRec.svStart) > c.minConsWindow) strEnd = boost::to_upper_copy(std::string(ref + svRec.svEndBeg, ref + svRec.svEndEnd));
	else strEnd=boost::to_upper_copy(std::string(ref + svRec.svStart, ref + svRec.svEndEnd));
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
	if ((svRec.svEnd - svRec.svStart) > c.minConsWindow) return boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svStartEnd)) + strRevComp;
	else return boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svStartEnd)) + strRevComp + boost::to_upper_copy(std::string(ref + svRec.svEnd, ref + svRec.svEndEnd));
      } else if (svt == 1) {
	std::string strStart;
	if ((svRec.svEnd - svRec.svStart) > c.minConsWindow) strStart = boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svStartEnd));
	else strStart = boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svEnd));
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
	if ((svRec.svEnd - svRec.svStart) > c.minConsWindow) return strRevComp + boost::to_upper_copy(std::string(ref + svRec.svEndBeg, ref + svRec.svEndEnd));
	else return boost::to_upper_copy(std::string(ref + svRec.svStartBeg, ref + svRec.svStart)) + strRevComp + boost::to_upper_copy(std::string(ref + svRec.svEndBeg, ref + svRec.svEndEnd));
      }
    }
    return "";
  }


  template<typename TConfig, typename TString, typename TSvRecord, typename TAlignDescriptor, typename TPosition>
  inline bool
  _coordTransform(TConfig const& c, TString const& ref, TSvRecord const& sv, TAlignDescriptor const& ad, TPosition& finalGapStart, TPosition& finalGapEnd, int32_t svt) {
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
	if (sv.svEnd - sv.svStart > c.indelsize) {
	  int32_t annealed = sv.svStartEnd - sv.svStartBeg;
	  if ((ad.rStart >= annealed) || (ad.rEnd < annealed)) return false;
	  finalGapStart = sv.svStartBeg + ad.rStart;
	  finalGapEnd = sv.svEndBeg + (ad.rEnd - annealed);
	} else {
	  finalGapStart = sv.svStartBeg + ad.rStart;
	  finalGapEnd = sv.svStartBeg + ad.rEnd;
	}
	return true;
      } else if (svt == 3) {
	int32_t annealed = sv.svEndEnd - sv.svEndBeg;
	if ((ad.rStart >= annealed) || (ad.rEnd < annealed)) return false;
	finalGapStart = sv.svStartBeg + (ad.rEnd - annealed);
	finalGapEnd = sv.svEndBeg + ad.rStart;
	return true;
      } else if (svt == 0) {
	if ((sv.svEnd - sv.svStart) > c.minConsWindow) {
	  int32_t annealed = sv.svStartEnd - sv.svStartBeg;
	  if ((ad.rStart >= annealed) || (ad.rEnd < annealed)) return false;
	  finalGapStart = sv.svStartBeg + ad.rStart;
	  finalGapEnd = sv.svEndBeg + (ref.size() - ad.rEnd) + 1;
	} else {
	  int32_t annealed = sv.svStartEnd - sv.svStartBeg;
	  if ((ad.rStart >= annealed) || (ad.rEnd < annealed)) return false;
	  finalGapStart = sv.svStartBeg + ad.rStart;
	  finalGapEnd = sv.svEndEnd - (ad.rEnd - annealed);
	}
	return true;
      } else if (svt == 1) {
	if ((sv.svEnd - sv.svStart) > c.minConsWindow) {
	  int32_t annealed = sv.svStartEnd - sv.svStartBeg;
	  if ((ad.rStart >= annealed) || (ad.rEnd < annealed)) return false;
	  finalGapStart = sv.svStartBeg + (annealed - ad.rStart) + 1;
	  finalGapEnd = sv.svEndBeg + (ad.rEnd - annealed);
	} else {
	  int32_t annealed = (sv.svStart - sv.svStartBeg) + (sv.svEnd - sv.svStartBeg);
	  if ((ad.rStart >= annealed) || (ad.rEnd < annealed)) return false;
	  finalGapStart = sv.svStartBeg + (annealed - ad.rStart) + 1;
	  finalGapEnd = sv.svEndBeg + (ad.rEnd - annealed);
	}
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
    //std::cerr << ad.percId << ',' << c.flankQuality << std::endl;
    if (ad.percId < c.flankQuality) return false;

    // Find homology
    _findHomology(consensus, svRefStr, ad, svt);
    
    // Check flanking alignment length
    if ((ad.homLeft + c.minimumFlankSize > ad.cStart) || ( varIndex < ad.cEnd + ad.homRight + c.minimumFlankSize)) return false;
    if ((ad.homLeft + c.minimumFlankSize > ad.rStart) || ( refIndex < ad.rEnd + ad.homRight + c.minimumFlankSize)) return false;

    // Valid split-read alignment
    return true;
  }


  inline void
  editDistanceVec(std::string const& seqI, std::string const& seqJ, EdlibAlignResult& cigar, std::vector<uint32_t>& dist) {
    // Fill edit distance vector
    dist.resize(seqI.size());
    std::fill(dist.begin(), dist.end(), 0);

    // Compute edit distance
    int32_t tIdx = -1;
    int32_t qIdx = -1;
    uint32_t editDist = 0;
    // seqI
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_DELETE) {
	++tIdx;
	++editDist;
      }
      else if (cigar.alignment[j] == EDLIB_EDOP_INSERT) {
	++qIdx;
	++editDist;
	dist[qIdx] = editDist;
      }
      else {
	++tIdx;
	++qIdx;
	if (seqI[qIdx] != seqJ[tIdx]) ++editDist;
	dist[qIdx] = editDist;
      }
    }
  }

  template<typename TAlign>
  inline void
  glueAlignment(std::string const& query, std::string const& target, uint32_t const gaplen, EdlibAlignMode const modeCode, EdlibAlignResult& cigarLeft, EdlibAlignResult& cigarRight, TAlign& align) {
    // Create new alignment
    int32_t tIdx = -1;
    int32_t qIdx = -1;
    uint32_t missingStart = 0;
    if (modeCode == EDLIB_MODE_HW) {
      tIdx = cigarLeft.endLocations[0];
      for (int32_t i = 0; i < cigarLeft.alignmentLength; i++) {
	if (cigarLeft.alignment[i] != EDLIB_EDOP_INSERT) --tIdx;
      }
      if (tIdx >= 0) missingStart = tIdx + 1;
    }
    uint32_t missingEnd = cigarRight.endLocations[0];
    if (missingEnd < target.size()) missingEnd = target.size() - missingEnd - 1;
    align.resize(boost::extents[2][missingStart + cigarLeft.alignmentLength + gaplen + cigarRight.alignmentLength + missingEnd]);

    // fix start
    if (modeCode == EDLIB_MODE_HW) {
      if (missingStart) {
	for (uint32_t j = 0; j < missingStart; ++j) {
	  align[1][j] = target[j];
	  align[0][j] = '-';
	}
      }
    }
    // target
    for (int32_t j = 0; j < cigarLeft.alignmentLength; ++j) {
      if (cigarLeft.alignment[j] == EDLIB_EDOP_INSERT) align[1][j + missingStart] = '-';
      else {
	++tIdx;
	align[1][j + missingStart] = target[tIdx];
      }
    }
    // query
    for (int32_t j = 0; j < cigarLeft.alignmentLength; ++j) {
      if (cigarLeft.alignment[j] == EDLIB_EDOP_DELETE) align[0][j + missingStart] = '-';
      else align[0][j + missingStart] = query[++qIdx];
    }
    // gap
    for (uint32_t j = 0; j < gaplen; ++j) {
      align[0][j + missingStart + cigarLeft.alignmentLength] = '-';
      ++tIdx;
      align[1][j + missingStart + cigarLeft.alignmentLength] = target[tIdx];
    }
    // target
    for (int32_t j = 0; j < cigarRight.alignmentLength; ++j) {
      if (cigarRight.alignment[j] == EDLIB_EDOP_INSERT) align[1][j + missingStart + cigarLeft.alignmentLength + gaplen] = '-';
      else {
	++tIdx;
	align[1][j + missingStart + cigarLeft.alignmentLength + gaplen] = target[tIdx];
      }
    }
    // query
    for (int32_t j = 0; j < cigarRight.alignmentLength; ++j) {
      if (cigarRight.alignment[j] == EDLIB_EDOP_DELETE) align[0][j + missingStart + cigarLeft.alignmentLength + gaplen] = '-';
      else align[0][j + missingStart + cigarLeft.alignmentLength + gaplen] = query[++qIdx];
    }
    // fix end
    if (modeCode == EDLIB_MODE_HW) {
      if (missingEnd) {
	for (uint32_t j = 0; j < missingEnd; ++j) {
	  ++tIdx;
	  align[1][j + missingStart + cigarLeft.alignmentLength + gaplen + cigarRight.alignmentLength] = target[tIdx];
	  align[0][j + missingStart + cigarLeft.alignmentLength + gaplen + cigarRight.alignmentLength] = '-';
	}
      }
    }
  }


  template<typename TAlign>
  inline bool
  splitAlign(std::string const& cons, std::string const& svRefStr, TAlign& align) {
    // Prepare SHW alignment
    std::string prefix = svRefStr.substr(0, svRefStr.size()/3);
    EdlibAlignResult cigarPrefix = edlibAlign(prefix.c_str(), prefix.size(), cons.c_str(), cons.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    uint32_t csStart = infixStart(cigarPrefix);
    edlibFreeAlignResult(cigarPrefix);
    std::string suffix = svRefStr.substr(2 * svRefStr.size()/3);
    EdlibAlignResult cigarSuffix = edlibAlign(suffix.c_str(), suffix.size(), cons.c_str(), cons.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    uint32_t csEnd = infixEnd(cigarSuffix);
    edlibFreeAlignResult(cigarSuffix);
    
    if (csStart >= csEnd) return false;
    std::string cs = cons.substr(csStart, (csEnd - csStart));
    
    // Prefix alignment
    EdlibAlignResult cigar = edlibAlign(svRefStr.c_str(), svRefStr.size(), cs.c_str(), cs.size(), edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
    std::vector<uint32_t> distFwd;
    //printAlignment(svRefStr, cs, EDLIB_MODE_HW, cigar);
    editDistanceVec(svRefStr, cs, cigar, distFwd);
    edlibFreeAlignResult(cigar);
    // Suffix alignment
    std::string svRefStrRev = svRefStr;
    std::string csRev = cs;
    reverseComplement(svRefStrRev);
    reverseComplement(csRev);
    EdlibAlignResult cigarRev = edlibAlign(svRefStrRev.c_str(), svRefStrRev.size(), csRev.c_str(), csRev.size(), edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
    std::vector<uint32_t> distRev;
    //printAlignment(svRefStrRev, csRev, EDLIB_MODE_HW, cigarRev);
    editDistanceVec(svRefStrRev, csRev, cigarRev, distRev);
    edlibFreeAlignResult(cigarRev);
    
    // Best join
    uint32_t bestJoin = 0;
    for(uint32_t i = 1; i < distFwd.size() - 1; ++i) {
      if (distFwd[i] + distRev[distFwd.size() - i - 2] < distFwd[bestJoin] + distRev[distFwd.size() - bestJoin - 2]) bestJoin = i;
    }
    
    // Split reference and do infix alignment
    std::string svRefLeft = svRefStr.substr(0, bestJoin + 1);
    std::string svRefRight = svRefStr.substr(bestJoin + 1);
    
    // Infix alignment
    EdlibAlignResult cigarLeft = edlibAlign(svRefLeft.c_str(), svRefLeft.size(), cons.c_str(), cons.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    uint32_t leftEnd = infixEnd(cigarLeft);
    //printAlignment(svRefLeft, cons, EDLIB_MODE_HW, cigarLeft);
    EdlibAlignResult cigarRight = edlibAlign(svRefRight.c_str(), svRefRight.size(), cons.c_str(), cons.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    //printAlignment(svRefRight, cons, EDLIB_MODE_HW, cigarRight);
    uint32_t rightStart = infixStart(cigarRight);
    
    // Glue alignment
    if (leftEnd + 15 >= rightStart) return false;
    glueAlignment(svRefStr, cons, rightStart - leftEnd - 1, EDLIB_MODE_HW, cigarLeft, cigarRight, align);
    edlibFreeAlignResult(cigarLeft);
    edlibFreeAlignResult(cigarRight);
    return true;
  }
  

  template<typename TAlign>
  inline bool
  _consRefAlignment(std::string const& cons, std::string const& svRefStr, TAlign& aln, int32_t const svt) {
    AlignConfig<true, false> semiglobal;
    DnaScore<int> lnsc(1, -1, -1, -1);
    bool reNeedle = false;
    if (svt == 4) {
      reNeedle = splitAlign(cons, svRefStr, aln);
      //reNeedle = longNeedle(svRefStr, cons, aln, semiglobal, lnsc);
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
  _alignConsensus(TConfig const& c, std::string& consensus, std::string& svRefStr, StructuralVariantRecord& sv, Breakpoint& bp, bool const realign) {
    // Realign?
    if (realign) {
      std::string revc = consensus;
      reverseComplement(revc);
      EdlibAlignResult alignFwd = edlibAlign(svRefStr.c_str(), svRefStr.size(), consensus.c_str(), consensus.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
      EdlibAlignResult alignRev = edlibAlign(svRefStr.c_str(), svRefStr.size(), revc.c_str(), revc.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
      if (alignRev.editDistance < alignFwd.editDistance) consensus = revc;
      edlibFreeAlignResult(alignFwd);
      edlibFreeAlignResult(alignRev);
    }

    //wfa::WFAlignerGapAffine2Pieces aligner(8, 12, 4, 60, 1, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
    //aligner.alignEnd2End(consensus, svRefStr);
    //printWfaAlignment(consensus, svRefStr, aligner);

    // Consensus to reference alignment
    typedef boost::multi_array<char, 2> TAlign;
    TAlign align;
    //std::cerr << "Consensus-to-Reference alignment" << std::endl;
    if (!_consRefAlignment(consensus, svRefStr, align, sv.svt)) return false;
    
    // Debug consensus to reference alignment
    //for(uint32_t i = 0; i < align.shape()[0]; ++i) {
    //for(uint32_t j = 0; j< align.shape()[1]; ++j) {
    //std::cerr << align[i][j];
    //}
    //std::cerr << std::endl;
    //}
    //std::cerr << std::endl;
    //std::cerr << std::endl;

    // Check breakpoint
    AlignDescriptor ad;
    if (!_findSplit(c, consensus, svRefStr, align, ad, sv.svt)) return false;

    // Transform coordinates
    unsigned int finalGapStart = 0;
    unsigned int finalGapEnd = 0;
    if (!_coordTransform(c, svRefStr, bp, ad, finalGapStart, finalGapEnd, sv.svt)) return false;
    
    if ((_translocation(sv.svt)) || (finalGapStart < finalGapEnd)) {

      // Get exact alleles for INS and DEL
      if (sv.svEnd - sv.svStart <= c.indelsize) {
	if ((sv.svt == 2) || (sv.svt == 4)) {
	  std::string refVCF;
	  std::string altVCF;
	  int32_t cpos = 0;
	  bool inSV = false;
	  for(uint32_t j = 0; j<align.shape()[1]; ++j) {
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
	  sv.alleles = _addAlleles(refVCF, altVCF);
	}
      }
      sv.precise=true;
      sv.svStart=finalGapStart;
      sv.svEnd=finalGapEnd;
      sv.srAlignQuality = ad.percId;
      sv.insLen=ad.cEnd - ad.cStart - 1;
      sv.consBp = ad.cStart;
      sv.homLen=std::max(0, ad.homLeft + ad.homRight - 2);
      int32_t ci_wiggle = std::max(ad.homLeft, ad.homRight);
      sv.ciposlow = -ci_wiggle;
      sv.ciposhigh = ci_wiggle;
      sv.ciendlow = -ci_wiggle;
      sv.ciendhigh = ci_wiggle;
      return true;
    } else {
      return false;
    }
  }
  
  template<typename TConfig>
  inline bool
  alignConsensus(TConfig const& c, bam_hdr_t* hdr, char const* seq, char const* sndSeq, StructuralVariantRecord& sv, bool const realign) {
    if ( (int32_t) sv.consensus.size() < (2 * c.minimumFlankSize + sv.insLen)) return false;
    
    // Get reference slice
    Breakpoint bp(sv);
    if (sv.svt ==4) {
      int32_t bufferSpace = std::max((int32_t) ((sv.consensus.size() - sv.insLen) / 3), c.minimumFlankSize);
      _initBreakpoint(hdr, bp, bufferSpace, sv.svt);
    } else _initBreakpoint(hdr, bp, sv.consensus.size(), sv.svt);
    if (bp.chr != bp.chr2) bp.part1 = _getSVRef(c, sndSeq, bp, bp.chr2, sv.svt);
    std::string svRefStr = _getSVRef(c, seq, bp, bp.chr, sv.svt);

    // Debug
    //std::cerr << ">Cons:" << sv.svStart << std::endl;
    //std::cerr << sv.consensus << std::endl;
    //std::cerr << ">Ref:" << sv.svStart << std::endl;
    //std::cerr << svRefStr << std::endl;

    // Generate consensus alignment
    return _alignConsensus(c, sv.consensus, svRefStr, sv, bp, realign);
  }

  template<typename TConfig>
  inline bool
  alignConsensus(TConfig const& c, bam_hdr_t* hdr, char const* seq, char const* sndSeq, StructuralVariantRecord& sv) {
    return alignConsensus(c, hdr, seq, sndSeq, sv, false);
  }


}

#endif
