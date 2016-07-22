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

#ifndef TAGS_H
#define TAGS_H

namespace torali {

  // Tags
  struct DeletionTag;
  struct DuplicationTag;
  struct InversionTag;
  struct TranslocationTag;
  struct InsertionTag;

  template<typename SvTag>
    struct SVType {
    };


  template<typename CoverageTag>
    struct CoverageType {
    };

  struct BpLevelCount;
  struct NoBpLevelCount;
  
  template<typename BpLevelTag>
    struct BpLevelType {
    };

  // F+ 0
  // F- 1
  // R+ 2
  // R- 3

  template<typename TBamRecord>
  inline uint8_t
    getStrandIndependentOrientation(TBamRecord const& al) {
    if (al.flag & BAM_FREAD1) {
      if (!(al.flag & BAM_FREVERSE)) {
	if (!(al.flag & BAM_FMREVERSE)) return (al.pos < al.mpos) ? 0 : 1;
	else return (al.pos < al.mpos) ? 2 : 3;
      } else {
	if (!(al.flag & BAM_FMREVERSE)) return (al.pos > al.mpos) ? 2 : 3;
	else return (al.pos > al.mpos) ? 0 : 1;
      }
    } else {
      if (!(al.flag & BAM_FREVERSE)) {
	if (!(al.flag & BAM_FMREVERSE)) return (al.pos < al.mpos) ? 1 : 0;
	else return (al.pos < al.mpos) ? 2 : 3;
      } else {
	if (!(al.flag & BAM_FMREVERSE)) return (al.pos > al.mpos) ? 2 : 3;
	else return (al.pos > al.mpos) ? 1 : 0;
      }
    }
  }

  //FF+ 0
  //FF- 1
  //FR+ 2
  //FR- 3
  //RF+ 4
  //RF- 5
  //RR+ 6
  //RR- 7

  template<typename TBamRecord>
  inline uint8_t
    getStrandSpecificOrientation(TBamRecord const& al) {
    if (!(al.flag  & BAM_FREVERSE)) {
      if (!(al.flag & BAM_FMREVERSE)) {
        return (al.pos < al.mpos) ? 0 : 1;
      } else {
        return (al.pos < al.mpos) ? 2 : 3;
      }
    } else {
      if (!(al.flag & BAM_FMREVERSE)) {
        return (al.pos > al.mpos) ? 4 : 5;
      } else {
        return (al.pos > al.mpos) ? 6 : 7;
      }
    }
  }

  // Deletions
  inline uint8_t
  _getCT(SVType<DeletionTag>) {
    return 2;
  }

  // Duplications
  inline uint8_t
  _getCT(SVType<DuplicationTag>) {
    return 3;
  }

  // Inversions
  inline uint8_t
  _getCT(SVType<InversionTag>) {
    return 0;
  }
  
  // Translocations
  inline uint8_t
  _getCT(SVType<TranslocationTag>) {
    return 0;
  }
  
  // Insertion
  inline uint8_t
  _getCT(SVType<InsertionTag>) {
    return 4;
  }
  

  template<typename TTag>
  inline bool
  _translocationMode(SVType<TTag>) {
    return false;
  }

  inline bool
  _translocationMode(SVType<TranslocationTag>) {
    return true;
  }

  
  // Deletions
  template<typename TBamRecord>
    inline uint8_t
    _getSpanOrientation(TBamRecord const&, uint8_t const, SVType<DeletionTag>) {
    return 2;
  }

  // Duplications
  template<typename TBamRecord>
    inline uint8_t 
    _getSpanOrientation(TBamRecord const&, uint8_t const, SVType<DuplicationTag>) {
    return 3;
  }

  // Left- or right-spanning
  template<typename TBamRecord>
    inline uint8_t 
    _getSpanOrientation(TBamRecord const& al, uint8_t const defaultOrient, SVType<InversionTag>) {
    uint8_t orient = getStrandIndependentOrientation(al);
    if (al.flag & BAM_FREAD1) {
      if (defaultOrient == 0) {
	if (((orient==2) && (al.pos < al.mpos)) || ((orient == 3) && (al.pos > al.mpos))) return 0;
      } else if (defaultOrient == 1) {
	if (((orient==2) && (al.pos > al.mpos)) || ((orient == 3) && (al.pos < al.mpos))) return 0;
      } else if (defaultOrient == 2) {
	if (((orient==0) && (al.pos < al.mpos)) || ((orient == 1) && (al.pos > al.mpos))) return 0;
      } else if (defaultOrient == 3) {
	if (((orient==0) && (al.pos > al.mpos)) || ((orient == 1) && (al.pos < al.mpos))) return 0;
      }
      return 1;
    } else {
      if (defaultOrient == 0) {
	if (((orient==2) && (al.pos > al.mpos)) || ((orient == 3) && (al.pos < al.mpos))) return 0;
      } else if (defaultOrient == 1) {
	if (((orient==2) && (al.pos < al.mpos)) || ((orient == 3) && (al.pos > al.mpos))) return 0;
      } else if (defaultOrient == 2) {
	if (((orient==0) && (al.pos > al.mpos)) || ((orient == 1) && (al.pos < al.mpos))) return 0;
      } else if (defaultOrient == 3) {
	if (((orient==0) && (al.pos < al.mpos)) || ((orient == 1) && (al.pos > al.mpos))) return 0;
      }
      return 1;
    }
  }




  template<typename TBamRecord>
    inline uint8_t 
    _inOrderAssign(TBamRecord const& al, bool flipped) {
    if (!flipped) {
      if (!(al.flag & BAM_FREVERSE)) {
	if (!(al.flag & BAM_FMREVERSE)) {
	  return (al.flag & BAM_FREAD1) ? 0 : 1;
	} else {
	  return 2;
	}
      } else {
	if (!(al.flag & BAM_FMREVERSE)) {
	  return 3;
	} else {
	  return (al.flag & BAM_FREAD1) ? 1 : 0;
	}
      }
    } else {
      if (!(al.flag & BAM_FREVERSE)) {
	if (!(al.flag & BAM_FMREVERSE)) {
	  return 2;
	} else {
	  return (al.flag & BAM_FREAD1) ? 0 : 1;
	}
      } else {
	if (!(al.flag & BAM_FMREVERSE)) {
	  return (al.flag & BAM_FREAD1) ? 1 : 0;
	} else {
	  return 3;
	}
      }
    }
  }


  template<typename TBamRecord>
    inline uint8_t 
    _getSpanOrientation(TBamRecord const& al, uint8_t const defaultOrient, SVType<TranslocationTag>) {
    uint8_t orient = getStrandIndependentOrientation(al);
    bool flipped = ( ((defaultOrient<2) && (orient>=2)) || ((defaultOrient>=2) && (orient<2)) );
    bool inOrder = (_inOrderAssign(al, flipped) == defaultOrient);
    if (flipped) {
      if (inOrder) return 0;
      else return 1;
    } else {
      if (inOrder) return 2;
      else return 3;
    }
  }

  // Insertion
  template<typename TBamRecord>
    inline uint8_t 
    _getSpanOrientation(TBamRecord const&, uint8_t const, SVType<InsertionTag>) {
    return 4;
  }

  // Reduced structural variant record for cov
  struct CovRecord {
    int32_t chr;
    int32_t svStart;
    int32_t svEnd;
    int32_t peSupport;
    uint32_t id;

  CovRecord() : chr(0), svStart(0), svEnd(0), peSupport(0), id(0) {}
  CovRecord(int32_t const c, int32_t const s, int32_t const e) : chr(c), svStart(s), svEnd(e), peSupport(0), id(0) {}
  };

  // Structural variant record
  struct StructuralVariantRecord {
    int32_t svStart;
    int32_t svEnd;
    int32_t peSupport;
    int32_t srSupport;
    int32_t wiggle;
    int32_t insLen;
    int32_t homLen;
    uint32_t id;
    float srAlignQuality;
    bool precise;
    uint8_t ct;
    uint8_t peMapQuality;
    int32_t chr;
    int32_t chr2;
    std::string alleles;
    std::string consensus;

  StructuralVariantRecord() : svStart(0), svEnd(0), peSupport(0), srSupport(0), wiggle(0), insLen(0), homLen(0), id(0), srAlignQuality(0), precise(false), ct(0), peMapQuality(0), chr(0), chr2(0) {}
  StructuralVariantRecord(int32_t const c, int const s, int const e) : svStart(s), svEnd(e), peSupport(0), srSupport(0), wiggle(0), insLen(0), homLen(0), id(0), srAlignQuality(0), precise(false), ct(0), peMapQuality(0), chr(c), chr2(c) {}
  };

  template<typename TSV>
  struct SortSVs : public std::binary_function<TSV, TSV, bool>
  {
    inline bool operator()(TSV const& sv1, TSV const& sv2) {
      return ((sv1.chr<sv2.chr) || ((sv1.chr==sv2.chr) && (sv1.svStart<sv2.svStart)) || ((sv1.chr==sv2.chr) && (sv1.svStart==sv2.svStart) && (sv1.svEnd<sv2.svEnd)) || ((sv1.chr==sv2.chr) && (sv1.svStart==sv2.svStart) && (sv1.svEnd==sv2.svEnd) && (sv1.peSupport > sv2.peSupport)));
    }
  };


  struct Breakpoint {
    int32_t svStartBeg;
    int32_t svStartEnd;
    int32_t svEndBeg;
    int32_t svEndEnd;
    int32_t svStart;
    int32_t svEnd;
    uint8_t ct;
    int32_t chr;
    int32_t chr2;
    std::string part1;

    Breakpoint() : svStartBeg(0), svStartEnd(0), svEndBeg(0), svEndEnd(0), svStart(0), svEnd(0), ct(0), chr(0), chr2(0) {}
    Breakpoint(StructuralVariantRecord const& sv) : svStartBeg(sv.svStart), svStartEnd(sv.svStart), svEndBeg(sv.svEnd), svEndEnd(sv.svEnd), svStart(sv.svStart), svEnd(sv.svEnd), ct(sv.ct), chr(sv.chr), chr2(sv.chr2) {}
    
  };

  // Initialize breakpoint

  // Deletions, insertions, duplications and inversions
  template<typename TBreakpoint, typename TTag>
  inline void
  _initBreakpoint(bam_hdr_t* hdr, TBreakpoint& bp, int32_t const boundary, SVType<TTag>) {
    bp.svStartBeg = std::max(0, bp.svStart - boundary);
    bp.svStartEnd = std::min(bp.svStart + boundary, (bp.svStart + bp.svEnd)/2);
    bp.svEndBeg = std::max((bp.svStart + bp.svEnd)/2 + 1, bp.svEnd - boundary);
    bp.svEndEnd = std::min((int32_t) (hdr->target_len[bp.chr2]), bp.svEnd + boundary);
  }

  template<typename TBreakpoint>
  inline void
  _initBreakpoint(bam_hdr_t* hdr, TBreakpoint& bp, int32_t const boundary, SVType<TranslocationTag>) {
    bp.svStartBeg = std::max(0, bp.svStart - boundary);
    bp.svStartEnd = std::min((int32_t) (hdr->target_len[bp.chr]), bp.svStart + boundary);
    bp.svEndBeg = std::max(0, bp.svEnd - boundary);
    bp.svEndEnd = std::min((int32_t) (hdr->target_len[bp.chr2]), bp.svEnd + boundary);
  }

  template<typename TBreakpoint>
  inline void
  _initBreakpoint(bam_hdr_t* hdr, TBreakpoint& bp, int32_t const boundary, SVType<InsertionTag>) {
    bp.svStartBeg = std::max(0, bp.svStart - boundary);
    bp.svStartEnd = std::min((int32_t) (hdr->target_len[bp.chr]), bp.svStart + boundary);
    bp.svEndBeg = std::max(0, bp.svEnd - boundary);
    bp.svEndEnd = std::min((int32_t) (hdr->target_len[bp.chr2]), bp.svEnd + boundary);
  }

  
  
  // SV Paired-end checks

  // Deletions, insertions, duplications and inversions
  template<typename TPos, typename TTag>
    inline TPos
    _minCoord(TPos const position, TPos const matePosition, SVType<TTag>) {
    return std::min(position, matePosition);
  }

  // Translocations
  template<typename TPos>
    inline TPos
    _minCoord(TPos const position, TPos const, SVType<TranslocationTag>) {
    return position;
  }

  // Deletions, insertions, duplications and inversions
  template<typename TPos, typename TTag>
    inline TPos
    _maxCoord(TPos const position, TPos const matePosition, SVType<TTag>) {
    return std::max(position, matePosition);
  }

  // Translocations
  template<typename TPos>
    inline TPos
    _maxCoord(TPos const, TPos const matePosition, SVType<TranslocationTag>) {
    return matePosition;
  }

  // Deletions, duplications and inversions
  template<typename TQualities, typename TAlen, typename TTag>
    inline void
    _resetQualities(TQualities& qualities, TAlen& alen, SVType<TTag>) {
    qualities.clear();
    alen.clear();
  }

  // Translocations
  template<typename TQualities, typename TAlen>
    inline void
    _resetQualities(TQualities&, TAlen&, SVType<TranslocationTag>) {
    // Nop
  }

  // Deletions
  template<typename TRef, typename TPos>
    inline bool
    _mappingPos(TRef const refID, TRef const mateRefID, TPos const position, TPos const matePosition, SVType<DeletionTag>) {
    return ((refID!=mateRefID) || (position==matePosition));
  }

  // Duplications
  template<typename TRef, typename TPos>
    inline bool
    _mappingPos(TRef const refID, TRef const mateRefID, TPos const position, TPos const matePosition, SVType<DuplicationTag>) {
    return ((refID!=mateRefID) || (std::abs(position - matePosition) < 100 ));
  }

  // Inversions
  template<typename TRef, typename TPos>
    inline bool
    _mappingPos(TRef const refID, TRef const mateRefID, TPos const position, TPos const matePosition, SVType<InversionTag>) {
    return ((refID!=mateRefID) || (position==matePosition));
  }

  // Translocations
  template<typename TRef, typename TPos>
    inline bool
    _mappingPos(TRef const refID, TRef const mateRefID, TPos const, TPos const, SVType<TranslocationTag>) {
    return (refID==mateRefID);
  }

  // Insertion
  template<typename TRef, typename TPos>
    inline bool
    _mappingPos(TRef const refID, TRef const mateRefID, TPos const position, TPos const matePosition, SVType<InsertionTag>) {
    return ((refID!=mateRefID) || (position==matePosition));
  }

  // Deletions, duplications and inversions
  template<typename TRef, typename TPos, typename TTag>
    inline bool
    _mappingPosGeno(TRef const refID, TRef const mateRefID, TPos const position, TPos const matePosition, SVType<TTag> svType) {
    return _mappingPos(refID, mateRefID, position, matePosition, svType);
  }

  // Translocations
  template<typename TRef, typename TPos>
    inline bool
    _mappingPosGeno(TRef const refID, TRef const mateRefID, TPos const position, TPos const matePosition, SVType<TranslocationTag>) {
    return ((refID==mateRefID) && (position==matePosition));
  }

  // Deletions, duplications and inversions
  template<typename TRef, typename TPos, typename TTag>
    inline bool
    _firstPairObs(TRef const, TRef const, TPos const position, TPos const matePosition, SVType<TTag>) {
    return (position<matePosition);
  }

  // Translocations
  template<typename TRef, typename TPos>
    inline bool
    _firstPairObs(TRef const refID, TRef const mateRefID, TPos const position, TPos const matePosition, SVType<TranslocationTag>) {
    return ((refID<mateRefID) || ((refID==mateRefID) && (position<matePosition)));
  }

  // Deletions
  template<typename TISize, typename TLibInfo>
  inline bool
  _acceptedInsertSize(TLibInfo& libInfo, TISize const iSize, SVType<DeletionTag>) {
    return ((libInfo.maxISizeCutoff > iSize) || (libInfo.median == 0));
  }

  // Insertions
  template<typename TISize, typename TLibInfo>
  inline bool
  _acceptedInsertSize(TLibInfo&, TISize const, SVType<InsertionTag>) {
    return true;
    //return ((libInfo.minISizeCutoff <= iSize) || (libInfo.median == 0));
  }

  // Duplications
  template<typename TISize, typename TLibInfo>
  inline bool
  _acceptedInsertSize(TLibInfo& libInfo, TISize const iSize, SVType<DuplicationTag>) {
    // Exclude the chimeras in mate-pair libraries
    return !(((libInfo.median>0) && (libInfo.median<1000)) || ((libInfo.median>=1000) && (iSize >=1000)));
  }

  // Inversions
  template<typename TISize, typename TLibInfo>
  inline bool
  _acceptedInsertSize(TLibInfo& libInfo, TISize const, SVType<InversionTag>) {
    return (libInfo.median == 0);
  }

  // Translocations
  template<typename TISize, typename TLibInfo>
  inline bool
  _acceptedInsertSize(TLibInfo& libInfo, TISize const, SVType<TranslocationTag>) {
    return (libInfo.median == 0);
  }

  // Deletions
  template<typename TOrientation>
  inline bool
    _acceptedOrientation(TOrientation const def, TOrientation const lib, SVType<DeletionTag>) {
    return (def != lib);
  }

  // Insertions
  template<typename TOrientation>
  inline bool
    _acceptedOrientation(TOrientation const def, TOrientation const lib, SVType<InsertionTag>) {
    return (def != lib);
  }

  // Duplications
  template<typename TOrientation>
    inline bool
    _acceptedOrientation(TOrientation const def, TOrientation const lib, SVType<DuplicationTag>) {
    if (def==0) return (lib != 1);
    else if (def==1) return (lib!=0);
    else if (def==2) return (lib!=3);
    else if (def==3) return (lib!=2);
    else return true;
  }

  // Inversion
  template<typename TOrientation>
    inline bool
    _acceptedOrientation(TOrientation const def, TOrientation const lib, SVType<InversionTag>) {
    return (((def<=1) && (lib!=2) && (lib!=3)) || ((def>=2) && (lib!=0) && (lib!=1)));
  }

  // Translocation
  template<typename TOrientation>
  inline bool
  _acceptedOrientation(TOrientation const, TOrientation const, SVType<TranslocationTag>) {
    return false;
  }


  // Deletions
  template<typename TSize, typename TISize>
    inline bool
    _pairsDisagree(TSize const pair1Min, TSize const pair1Max, TSize const pair1ReadLength, TISize const pair1maxNormalISize, TSize const pair2Min, TSize const pair2Max, TSize const pair2ReadLength, TISize const pair2maxNormalISize, uint8_t const, uint8_t const, SVType<DeletionTag>) {
    //std::cout << pair1Min << ',' << pair1Max << ',' << pair1ReadLength << ',' << pair1maxNormalISize << ',' << pair2Min << ',' << pair2Max << ',' << pair2ReadLength << ',' << pair2maxNormalISize << std::endl;
    if ((pair2Min + pair2ReadLength - pair1Min) > pair1maxNormalISize) return true;
    if ((pair2Max < pair1Max) && ((pair1Max + pair1ReadLength - pair2Max) > pair1maxNormalISize)) return true;
    if ((pair2Max >= pair1Max) && ((pair2Max + pair2ReadLength - pair1Max) > pair2maxNormalISize)) return true;
    if ((pair1Max < pair2Min) || (pair2Max < pair1Min)) return true;
    return false;
  }

  // Insertions
  template<typename TSize, typename TISize>
    inline bool
    _pairsDisagree(TSize const pair1Min, TSize const pair1Max, TSize const pair1ReadLength, TISize const pair1maxNormalISize, TSize const pair2Min, TSize const pair2Max, TSize const pair2ReadLength, TISize const pair2maxNormalISize, uint8_t const, uint8_t const, SVType<InsertionTag>) {
    //std::cout << pair1Min << ',' << pair1Max << ',' << pair1ReadLength << ',' << pair1maxNormalISize << ',' << pair2Min << ',' << pair2Max << ',' << pair2ReadLength << ',' << pair2maxNormalISize << std::endl;
    if ((pair2Min + pair2ReadLength - pair1Min) > pair1maxNormalISize) return true;
    if ((pair2Max < pair1Max) && ((pair1Max + pair1ReadLength - pair2Max) > pair1maxNormalISize)) return true;
    if ((pair2Max >= pair1Max) && ((pair2Max + pair2ReadLength - pair1Max) > pair2maxNormalISize)) return true;
    if ((pair1Max < pair2Min) || (pair2Max < pair1Min)) return true;
    return false;
  }


  // Duplications
  template<typename TSize, typename TISize>
    inline bool
    _pairsDisagree(TSize const pair1Min, TSize const pair1Max, TSize const pair1ReadLength, TISize const pair1maxNormalISize, TSize const pair2Min, TSize const pair2Max, TSize const pair2ReadLength, TISize const pair2maxNormalISize, uint8_t const, uint8_t const, SVType<DuplicationTag>) {
    if ((pair2Min + pair2ReadLength - pair1Min) > pair2maxNormalISize) return true;
    if ((pair2Max < pair1Max) && ((pair1Max + pair1ReadLength - pair2Max) > pair2maxNormalISize)) return true;
    if ((pair2Max >= pair1Max) && ((pair2Max + pair2ReadLength - pair1Max) > pair1maxNormalISize)) return true;
    return false;
  }


  // Inversions
  template<typename TSize, typename TISize>
    inline bool
    _pairsDisagree(TSize const pair1Min, TSize const pair1Max, TSize const pair1ReadLength, TISize const pair1maxNormalISize, TSize const pair2Min, TSize const pair2Max, TSize const pair2ReadLength, TISize const pair2maxNormalISize, uint8_t const ct1, uint8_t const ct2, SVType<InversionTag>) {
    // Do both pairs support the same inversion type (left- or right-spanning)
    if (ct1 != ct2) return true;
    if (!ct1) {
      // Left-spanning inversions
      if ((pair2Min + pair2ReadLength - pair1Min) > pair1maxNormalISize) return true;
      if ((pair2Max < pair1Max) && ((pair1Max + pair1ReadLength - pair2Max) > pair2maxNormalISize)) return true;
      if ((pair2Max >= pair1Max) && ((pair2Max + pair2ReadLength - pair1Max) > pair1maxNormalISize)) return true;
    } else {
      // Right-spanning inversions
      if ((pair2Min + pair2ReadLength - pair1Min) > pair2maxNormalISize) return true;
      if ((pair2Max < pair1Max) && ((pair1Max + pair1ReadLength - pair2Max) > pair1maxNormalISize)) return true;
      if ((pair2Max >= pair1Max) && ((pair2Max + pair2ReadLength - pair1Max) > pair2maxNormalISize)) return true;
    }
    return false;
  }


  // Translocations
  template<typename TSize, typename TISize>
    inline bool
    _pairsDisagree(TSize const pair1Min, TSize const pair1Max, TSize const pair1ReadLength, TISize const pair1maxNormalISize, TSize const pair2Min, TSize const pair2Max, TSize const pair2ReadLength, TISize const pair2maxNormalISize, uint8_t const ct1, uint8_t const ct2, SVType<TranslocationTag>) {
    // Do both pairs support the same translocation type
    if (ct1 != ct2) return true;

    // Check read offsets
    if (ct1%2==0) {
      if ((pair2Min + pair2ReadLength - pair1Min) > pair1maxNormalISize) return true;
      if (ct1>=2) {
	if (pair2Max < pair1Max) {
	  if ((pair1Max + pair1ReadLength - pair2Max) > pair1maxNormalISize) return true;
	} else {
	  if ((pair2Max + pair2ReadLength - pair1Max) > pair2maxNormalISize) return true;
	}
      } else {
	if (pair2Max < pair1Max) {
	  if ((pair1Max + pair1ReadLength - pair2Max) > pair2maxNormalISize) return true;
	} else {
	  if ((pair2Max + pair2ReadLength - pair1Max) > pair1maxNormalISize) return true;
	}
      }
    } else {
      if ((pair2Min + pair2ReadLength - pair1Min) > pair2maxNormalISize) return true;
      if (ct1>=2) {
	if (pair2Max < pair1Max) {
	  if ((pair1Max + pair1ReadLength - pair2Max) > pair2maxNormalISize) return true;
	} else {
	  if ((pair2Max + pair2ReadLength - pair1Max) > pair1maxNormalISize) return true;
	}
      } else {
	if (pair2Max < pair1Max) {
	  if ((pair1Max + pair1ReadLength - pair2Max) > pair1maxNormalISize) return true;
	} else {
	  if ((pair2Max + pair2ReadLength - pair1Max) > pair2maxNormalISize) return true;
	}
      }
    }
    return false;
  }






}

#endif
