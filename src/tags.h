/*
============================================================================
DELLY: Structural variant discovery by integrated PE mapping and SR analysis
============================================================================
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

  #ifndef DELLY_SVT_TRANS
  #define DELLY_SVT_TRANS 5
  #endif

  inline bool
  _translocation(int32_t const svt) {
    return (DELLY_SVT_TRANS <= svt);
  }

  inline bool
  _translocation(bam1_t* rec) {
    return (rec->core.tid != rec->core.mtid);
  }
  
  // Deletions
  inline uint8_t
  _getSpanOrientation(int32_t const svt) {
    if (_translocation(svt)) {
      return svt - DELLY_SVT_TRANS;
    } else {
      return svt;
    }
  }

  // Structural variant record
  struct StructuralVariantRecord {
    int32_t chr;
    int32_t svStart;
    int32_t chr2;
    int32_t svEnd;
    int32_t ciposlow;
    int32_t ciposhigh;
    int32_t ciendlow;
    int32_t ciendhigh;
    int32_t peSupport;
    int32_t srSupport;
    int32_t wiggle;
    int32_t insLen;
    int32_t homLen;
    int32_t svt;
    uint32_t id;
    float srAlignQuality;
    uint8_t peMapQuality;
    bool precise;
    std::string alleles;
    std::string consensus;

    
    StructuralVariantRecord() : svt(-1) {}

  StructuralVariantRecord(int32_t const c, int32_t const s, int32_t const e) : chr(c), svStart(s), chr2(c), svEnd(e), peSupport(0), svt(-1) {}


    StructuralVariantRecord(int32_t const c1, int32_t const s, int32_t const c2, int32_t const e, int32_t const cipl, int32_t const ciph, int32_t const ciel, int32_t const cieh, int32_t const sup, uint8_t mq, int32_t const ilen): chr(c1), svStart(s), chr2(c2), svEnd(e), ciposlow(cipl), ciposhigh(ciph), ciendlow(ciel), ciendhigh(cieh), peSupport(0), srSupport(sup), insLen(ilen), svt(-1), peMapQuality(mq), precise(true) {}
    
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
    int32_t peSupport;
    int32_t svt;
    int32_t chr;
    int32_t chr2;
    std::string part1;

    Breakpoint() : svStartBeg(0), svStartEnd(0), svEndBeg(0), svEndEnd(0), svStart(0), svEnd(0), peSupport(0), svt(-1), chr(0), chr2(0) {}
    explicit Breakpoint(StructuralVariantRecord const& sv) : svStartBeg(sv.svStart), svStartEnd(sv.svStart), svEndBeg(sv.svEnd), svEndEnd(sv.svEnd), svStart(sv.svStart), svEnd(sv.svEnd), peSupport(sv.peSupport), svt(sv.svt), chr(sv.chr), chr2(sv.chr2) {}
    
  };

  // Initialize breakpoint
  template<typename TBreakpoint>
  inline void
  _initBreakpoint(bam_hdr_t* hdr, TBreakpoint& bp, int32_t const boundary, int32_t const svt) {
    if (_translocation(svt)) {
        bp.svStartBeg = std::max(0, bp.svStart - boundary);
	bp.svStartEnd = std::min((int32_t) (hdr->target_len[bp.chr]), bp.svStart + boundary);
	bp.svEndBeg = std::max(0, bp.svEnd - boundary);
	bp.svEndEnd = std::min((int32_t) (hdr->target_len[bp.chr2]), bp.svEnd + boundary);
    } else {
      if (svt == 4) {
	bp.svStartBeg = std::max(0, bp.svStart - boundary);
	bp.svStartEnd = std::min((int32_t) (hdr->target_len[bp.chr]), bp.svStart + boundary);
	bp.svEndBeg = std::max(0, bp.svEnd - boundary);
	bp.svEndEnd = std::min((int32_t) (hdr->target_len[bp.chr2]), bp.svEnd + boundary);
      } else {
	bp.svStartBeg = std::max(0, bp.svStart - boundary);
	bp.svStartEnd = std::min(bp.svStart + boundary, (bp.svStart + bp.svEnd)/2);
	bp.svEndBeg = std::max((bp.svStart + bp.svEnd)/2 + 1, bp.svEnd - boundary);
	bp.svEndEnd = std::min((int32_t) (hdr->target_len[bp.chr2]), bp.svEnd + boundary);
      }
    }
  }
  
  template<typename TPos>
  inline TPos
  _minCoord(TPos const position, TPos const matePosition, int32_t const svt) {
    if (_translocation(svt)) return position;
    else return std::min(position, matePosition);
  }

  template<typename TPos>
  inline TPos
  _maxCoord(TPos const position, TPos const matePosition, int32_t const svt) {
    if (_translocation(svt)) return matePosition;
    else return std::max(position, matePosition);
  }


  // Deletions, duplications and inversions
  template<typename TRef, typename TPos>
  inline bool
  _mappingPosGeno(TRef const refID, TRef const mateRefID, TPos const position, TPos const matePosition, int32_t const svt) {
    if (_translocation(svt)) return ((refID==mateRefID) && (position==matePosition));
    else {
      if (svt == 3) return ((refID!=mateRefID) || (std::abs(position - matePosition) < 100 ));
      else return ((refID!=mateRefID) || (position==matePosition));
    }
  }
  

  template<typename TSize>
  inline bool
  _svSizeCheck(TSize const s, TSize const e, int32_t const svt) {
    if (svt == 0) return (( e - s ) >= 100);
    else if (svt == 1) return (( e - s ) >= 100);
    else if (svt == 2) return (( e - s ) >= 300);
    else if (svt == 3) return (( e - s ) >= 100);
    else return true;
  }


  // 0: Left-spanning inversion
  // 1: Right-spanning inversion
  // 2: Deletion-type
  // 3: Duplication-type

  template<typename TBamRecord>
  inline uint8_t
  getSVType(TBamRecord const& al) {
    if (al.flag & BAM_FREAD1) {
      if (!(al.flag & BAM_FREVERSE)) {
	if (!(al.flag & BAM_FMREVERSE)) return 0;
	else return (al.pos < al.mpos) ? 2 : 3;
      } else {
	if (!(al.flag & BAM_FMREVERSE)) return (al.pos > al.mpos) ? 2 : 3;
	else return 1;
      }
    } else {
      if (!(al.flag & BAM_FREVERSE)) {
	if (!(al.flag & BAM_FMREVERSE)) return 0;
	else return (al.pos < al.mpos) ? 2 : 3;
      } else {
	if (!(al.flag & BAM_FMREVERSE)) return (al.pos > al.mpos) ? 2 : 3;
	else return 1;
      }
    }
  }

  inline int32_t
  _isizeMappingPos(bam1_t* rec, int32_t isize) {
    if (_translocation(rec)) return DELLY_SVT_TRANS + getSVType(rec->core);
    else {
      if (rec->core.pos == rec->core.mpos) return -1; // No SV
      uint8_t orient = getSVType(rec->core);
      if (orient == 0) return 0;
      else if (orient == 1) return 1;
      else if (orient == 2) {
	if (isize > std::abs(rec->core.isize)) return -1;
	else return 2;
      } else {
	if (std::abs(rec->core.pos - rec->core.mpos) < 100) return -1; // Too small
	return 3;
      }
    } 
  }

  // Decode Orientation
  inline int32_t
  _decodeOrientation(std::string const& value, std::string const& svt) {
    if (svt == "BND") {
      if (value=="3to3") return DELLY_SVT_TRANS + 0;
      else if (value=="5to5") return DELLY_SVT_TRANS + 1;
      else if (value=="3to5") return DELLY_SVT_TRANS + 2;
      else if (value=="5to3") return DELLY_SVT_TRANS + 3;
      else return -1;
    } else {
      if (value=="3to3") return 0;
      else if (value=="5to5") return 1;
      else if (value=="3to5") return 2;
      else if (value=="5to3") return 3;
      else return 4;
    }
  }

  inline unsigned hash_string(const char *s) {
    unsigned h = 37;
    while (*s) {
      h = (h * 54059) ^ (s[0] * 76963);
      s++;
    }
    return h;
  }
  
  template<typename TAlignedReads>
  inline bool
  _firstPairObs(bam1_t* rec, TAlignedReads const& lastAlignedPosReads) {
    if (rec->core.tid == rec->core.mtid) return ((rec->core.pos < rec->core.mpos) || ((rec->core.pos == rec->core.mpos) && (lastAlignedPosReads.find(hash_string(bam_get_qname(rec))) == lastAlignedPosReads.end())));
    else return (rec->core.tid < rec->core.mtid);
  }

  // Deletions
  template<typename TSize, typename TISize>
  inline bool
  _pairsDisagree(TSize const pair1Min, TSize const pair1Max, TSize const pair1ReadLength, TISize const pair1maxNormalISize, TSize const pair2Min, TSize const pair2Max, TSize const pair2ReadLength, TISize const pair2maxNormalISize, int32_t const svt) {
    if (_translocation(svt)) {
      uint8_t ct = _getSpanOrientation(svt);

      // Check read offsets
      if (ct%2==0) {
	if ((pair2Min + pair2ReadLength - pair1Min) > pair1maxNormalISize) return true;
	if (ct>=2) {
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
	if (ct>=2) {
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
    } else {
      if (svt < 2) {
	// Inversion
	if (!svt) {
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
      } else if (svt == 2) {
	// Deletion
	if ((pair2Min + pair2ReadLength - pair1Min) > pair1maxNormalISize) return true;
	if ((pair2Max < pair1Max) && ((pair1Max + pair1ReadLength - pair2Max) > pair1maxNormalISize)) return true;
	if ((pair2Max >= pair1Max) && ((pair2Max + pair2ReadLength - pair1Max) > pair2maxNormalISize)) return true;
	if ((pair1Max < pair2Min) || (pair2Max < pair1Min)) return true;
	return false;
      } else if (svt == 3) {
	if ((pair2Min + pair2ReadLength - pair1Min) > pair2maxNormalISize) return true;
	if ((pair2Max < pair1Max) && ((pair1Max + pair1ReadLength - pair2Max) > pair2maxNormalISize)) return true;
	if ((pair2Max >= pair1Max) && ((pair2Max + pair2ReadLength - pair1Max) > pair1maxNormalISize)) return true;
	return false;
      }
    }
    return false;
  }


}

#endif
