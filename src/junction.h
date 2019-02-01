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

#ifndef JUNCTION_H
#define JUNCTION_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/math/distributions/binomial.hpp>

#include <htslib/sam.h>

#include "util.h"

namespace torali
{

  struct SRBamRecord {
    int32_t chr;
    int32_t pos;
    int32_t chr2;
    int32_t pos2;
    int32_t qual;
    int32_t inslen;
    int32_t svid;
    std::size_t id;
        
  SRBamRecord(int32_t const c, int32_t const p, int32_t const c2, int32_t const p2, int32_t const qval, int32_t const il, std::size_t const idval) : chr(c), pos(p), chr2(c2), pos2(p2), qual(qval), inslen(il), svid(-1), id(idval) {}
  };

  template<typename TSRBamRecord>
  struct SortSRBamRecord : public std::binary_function<TSRBamRecord, TSRBamRecord, bool>
  {
    inline bool operator()(TSRBamRecord const& sv1, TSRBamRecord const& sv2) {
      return ((sv1.chr<sv2.chr) || ((sv1.chr==sv2.chr) && (sv1.pos<sv2.pos)) || ((sv1.chr==sv2.chr) && (sv1.pos==sv2.pos) && (sv1.chr2<sv2.chr2)) || ((sv1.chr==sv2.chr) && (sv1.pos==sv2.pos) && (sv1.chr2==sv2.chr2) && (sv1.pos2 < sv2.pos2)));
    }
  };
  
  
  struct Junction {
    bool forward;
    bool scleft;
    int32_t refidx;
    int32_t refpos;
    int32_t seqpos;
    uint16_t qual;
    
    Junction(bool const fw, bool const cl, int32_t const idx, int32_t const r, int32_t const s, uint16_t const qval) : forward(fw), scleft(cl), refidx(idx), refpos(r), seqpos(s), qual(qval) {}
  };


  template<typename TReadBp>
  inline void
  _insertJunction(TReadBp& readBp, unsigned const seed, bam1_t* rec, int32_t const rp, int32_t const sp, bool const scleft) {
    bool fw = true;
    if (rec->core.flag & BAM_FREVERSE) fw = false;
    typedef typename TReadBp::mapped_type TJunctionVector;
    typename TReadBp::iterator it = readBp.find(seed);
    int32_t seqlen = sequenceLength(rec);
    if (sp <= seqlen) {
      if (rec->core.flag & BAM_FREVERSE) {
	if (it != readBp.end()) it->second.push_back(Junction(fw, scleft, rec->core.tid, rp, seqlen - sp, rec->core.qual));
	else readBp.insert(std::make_pair(seed, TJunctionVector(1, Junction(fw, scleft, rec->core.tid, rp, seqlen - sp, rec->core.qual))));
      } else {
	if (it != readBp.end()) it->second.push_back(Junction(fw, scleft, rec->core.tid, rp, sp, rec->core.qual));
	else readBp.insert(std::make_pair(seed, TJunctionVector(1, Junction(fw, scleft, rec->core.tid, rp, sp, rec->core.qual))));
      }
    }
  }

  template<typename TJunction>
  struct SortJunction : public std::binary_function<TJunction, TJunction, bool>
  {
    inline bool operator()(TJunction const& j1, TJunction const& j2) {
      return ((j1.seqpos<j2.seqpos) || ((j1.seqpos==j2.seqpos) && (j1.refidx<j2.refidx)) || ((j1.seqpos==j2.seqpos) && (j1.refidx==j2.refidx) && (j1.refpos<j2.refpos)) || ((j1.seqpos==j2.seqpos) && (j1.refidx==j2.refidx) && (j1.refpos==j2.refpos) && (j1.scleft < j2.scleft)));
    }
  };

  // Deletion junctions
  template<typename TConfig, typename TReadBp, typename TSRRecords>
  inline void
  selectDeletions(TConfig const& c, TReadBp const& readBp, TSRRecords& br) {
    for(typename TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
      if (it->second.size() > 1) {
	for(uint32_t i = 0; i < it->second.size(); ++i) {
	  for(uint32_t j = i+1; j < it->second.size(); ++j) {
	    if ((uint32_t) (it->second[j].seqpos - it->second[i].seqpos) > c.maxReadSep) break;
	    // Same chr, same direction, opposing soft-clips
	    if ((it->second[j].refidx == it->second[i].refidx) && (it->second[j].forward == it->second[i].forward) && (it->second[i].scleft != it->second[j].scleft)) {
	      // Min. deletion size
	      if ( (uint32_t) std::abs(it->second[j].refpos - it->second[i].refpos) > c.minRefSep) {
		// Avg. qval
		int32_t qval = (int32_t) (((int32_t) it->second[i].qual + (int32_t) it->second[j].qual) / 2);
		// Correct clipping architecture, note: soft-clipping of error-prone reads can lead to switching left/right breakpoints
		if (it->second[i].refpos <= it->second[j].refpos) {
		  if ((!it->second[i].scleft) && (it->second[j].scleft)) {
		    br[2].push_back(SRBamRecord(it->second[i].refidx, it->second[i].refpos, it->second[j].refidx, it->second[j].refpos, qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  }
		} else {
		  if ((it->second[i].scleft) && (!it->second[j].scleft)) {
		    br[2].push_back(SRBamRecord(it->second[j].refidx, it->second[j].refpos, it->second[i].refidx, it->second[i].refpos, qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }


  // Duplication junctions
  template<typename TConfig, typename TReadBp, typename TSRRecords>
  inline void
  selectDuplications(TConfig const& c, TReadBp const& readBp, TSRRecords& br) {
    for(typename TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
      if (it->second.size() > 1) {
	for(uint32_t i = 0; i < it->second.size(); ++i) {
	  for(uint32_t j = i+1; j < it->second.size(); ++j) {
	    if ((uint32_t) (it->second[j].seqpos - it->second[i].seqpos) > c.maxReadSep) break;
	    // Same chr, same direction, opposing soft-clips
	    if ((it->second[j].refidx == it->second[i].refidx) && (it->second[j].forward == it->second[i].forward) && (it->second[i].scleft != it->second[j].scleft)) {
	      // Min. duplication size
	      if ( (uint32_t) std::abs(it->second[j].refpos - it->second[i].refpos) > c.minRefSep) {
		// Avg. qval
		int32_t qval = (int32_t) (((int32_t) it->second[i].qual + (int32_t) it->second[j].qual) / 2);
		// Correct clipping architecture, note: soft-clipping of error-prone reads can lead to switching left/right breakpoints
		if (it->second[i].refpos <= it->second[j].refpos) {
		  if ((it->second[i].scleft) && (!it->second[j].scleft)) {
		    br[3].push_back(SRBamRecord(it->second[i].refidx, it->second[i].refpos, it->second[j].refidx, it->second[j].refpos, qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  }
		} else {
		  if ((!it->second[i].scleft) && (it->second[j].scleft)) {
		    br[3].push_back(SRBamRecord(it->second[j].refidx, it->second[j].refpos, it->second[i].refidx, it->second[i].refpos, qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  // Inversion junctions
  template<typename TConfig, typename TReadBp, typename TSRRecords>
  inline void
  selectInversions(TConfig const& c, TReadBp const& readBp, TSRRecords& br) {
    for(typename TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
      if (it->second.size() > 1) {
	for(uint32_t i = 0; i < it->second.size(); ++i) {
	  for(uint32_t j = i+1; j < it->second.size(); ++j) {
	    if ((uint32_t) (it->second[j].seqpos - it->second[i].seqpos) > c.maxReadSep) break;
	    // Same chr, different direction, agreeing soft-clips
	    if ((it->second[j].refidx == it->second[i].refidx) && (it->second[j].forward != it->second[i].forward) && (it->second[i].scleft == it->second[j].scleft)) {
	      // Min. inversion size
	      if ( (uint32_t) std::abs(it->second[j].refpos - it->second[i].refpos) > c.minRefSep) {
		// Avg. qval
		int32_t qval = (int32_t) (((int32_t) it->second[i].qual + (int32_t) it->second[j].qual) / 2);
		//std::cout << it->second[i].refidx << ':' <<  it->second[i].refpos << "(S:" << it->second[i].seqpos << ',' << it->second[i].forward << ',' << it->second[i].scleft << ')' << it->second[j].refidx << ':' << it->second[j].refpos << "(S:" << it->second[j].seqpos << ',' << it->second[j].forward << ',' << it->second[j].scleft << ')' << std::endl;
		if (it->second[i].refpos <= it->second[j].refpos) {
		  // Need to differentiate 3to3 and 5to5
		  if (it->second[i].scleft) br[1].push_back(SRBamRecord(it->second[i].refidx, it->second[i].refpos, it->second[j].refidx, it->second[j].refpos, qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  else br[0].push_back(SRBamRecord(it->second[i].refidx, it->second[i].refpos, it->second[j].refidx, it->second[j].refpos, qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		} else {
		  // Need to differentiate 3to3 and 5to5
		  if (it->second[i].scleft) br[1].push_back(SRBamRecord(it->second[j].refidx, it->second[j].refpos, it->second[i].refidx, it->second[i].refpos, qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  else br[0].push_back(SRBamRecord(it->second[j].refidx, it->second[j].refpos, it->second[i].refidx, it->second[i].refpos, qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		}
	      }
	    }
	  }
	}
      }
    }
  }

  // Insertion junctions
  template<typename TConfig, typename TReadBp, typename TSRRecords>
  inline void
  selectInsertions(TConfig const& c, TReadBp const& readBp, TSRRecords& br) {
    for(typename TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
      if (it->second.size() > 1) {
	for(uint32_t i = 0; i < it->second.size(); ++i) {
	  for(uint32_t j = i+1; j < it->second.size(); ++j) {
	    // Same chr, same direction, opposing soft-clips
	    if ((it->second[j].refidx == it->second[i].refidx) && (it->second[j].forward == it->second[i].forward) && (it->second[i].scleft != it->second[j].scleft)) {
	      // Reference insertion footprint should be small
	      if ( (uint32_t) std::abs(it->second[j].refpos - it->second[i].refpos) < c.maxReadSep) {
		// Large separation in sequence space
		if ((uint32_t) (it->second[j].seqpos - it->second[i].seqpos) > c.minRefSep) {
		  // Avg. qval
		  int32_t qval = (int32_t) (((int32_t) it->second[i].qual + (int32_t) it->second[j].qual) / 2);
		  if (it->second[i].refpos <= it->second[j].refpos) {
		    br[4].push_back(SRBamRecord(it->second[i].refidx, it->second[i].refpos, it->second[j].refidx, it->second[j].refpos, qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  } else {
		    br[4].push_back(SRBamRecord(it->second[j].refidx, it->second[j].refpos, it->second[i].refidx, it->second[i].refpos, qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }


  // Translocation junctions
  template<typename TConfig, typename TReadBp, typename TSRRecords>
  inline void
  selectTranslocations(TConfig const& c, TReadBp const& readBp, TSRRecords& br) {
    for(typename TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
      if (it->second.size() > 1) {
	for(uint32_t i = 0; i < it->second.size(); ++i) {
	  for(uint32_t j = i+1; j < it->second.size(); ++j) {
	    if ((uint32_t) (it->second[j].seqpos - it->second[i].seqpos) > c.maxReadSep) break;
	    // Different chr
	    if (it->second[j].refidx != it->second[i].refidx) {
	      int32_t chr1ev = j;
	      int32_t chr2ev = i;
	      if (it->second[i].refidx < it->second[j].refidx) {
		chr1ev = i;
		chr2ev = j;
	      }
	      // Avg. qval
	      int32_t qval = (int32_t) (((int32_t) it->second[i].qual + (int32_t) it->second[j].qual) / 2);
	      if (it->second[chr1ev].forward == it->second[chr2ev].forward) {
		// Same direction, opposing soft-clips
		if (it->second[chr1ev].scleft != it->second[chr2ev].scleft) {
		  if (it->second[chr1ev].scleft) {
		    // 5to3
		    br[DELLY_SVT_TRANS + 2].push_back(SRBamRecord(it->second[chr2ev].refidx, it->second[chr2ev].refpos, it->second[chr1ev].refidx, it->second[chr1ev].refpos, qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  } else {
		    // 3to5
		    br[DELLY_SVT_TRANS + 3].push_back(SRBamRecord(it->second[chr2ev].refidx, it->second[chr2ev].refpos, it->second[chr1ev].refidx, it->second[chr1ev].refpos, qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  }
		}
	      } else {
		// Opposing direction, same soft-clips
		if (it->second[chr1ev].scleft == it->second[chr2ev].scleft) {
		  if (it->second[chr1ev].scleft) {
		    // 3to3
		    br[DELLY_SVT_TRANS + 1].push_back(SRBamRecord(it->second[chr2ev].refidx, it->second[chr2ev].refpos, it->second[chr1ev].refidx, it->second[chr1ev].refpos, qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  } else {
		    // 5to5
		    br[DELLY_SVT_TRANS + 0].push_back(SRBamRecord(it->second[chr2ev].refidx, it->second[chr2ev].refpos, it->second[chr1ev].refidx, it->second[chr1ev].refpos, qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }


  template<typename TConfig, typename TSRBamRecords>
  inline void
  outputSRBamRecords(TConfig const& c, TSRBamRecords const& br) {
    samFile* samfile = sam_open(c.files[0].string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Header
    std::cerr << "chr1\tpos1\tchr2\tpos2\tsvtype\tct" << std::endl;

    // SVs
    for(uint32_t svt = 0; svt < br.size(); ++svt) {
      for(uint32_t i = 0; i < br[svt].size(); ++i) {
	std::cerr << hdr->target_name[br[svt][i].chr] << '\t' << br[svt][i].pos << '\t' << hdr->target_name[br[svt][i].chr2] << '\t' << br[svt][i].pos2 << '\t' << _addID(svt) << '\t' << _addOrientation(svt) << std::endl;
      }
    }

    // Clean-up
    bam_hdr_destroy(hdr);
    sam_close(samfile);
  }

  template<typename TConfig>
  inline void
  outputStructuralVariants(TConfig const& c, std::vector<StructuralVariantRecord> const& svs) {
    samFile* samfile = sam_open(c.files[0].string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Header
    std::cerr << "chr1\tpos1\tchr2\tpos2\tsvtype\tct\tpeSupport\tsrSupport" << std::endl;

    // SVs
    for(uint32_t i = 0; i < svs.size(); ++i) {
      std::cerr << hdr->target_name[svs[i].chr] << '\t' << svs[i].svStart << '\t' << hdr->target_name[svs[i].chr2] << '\t' << svs[i].svEnd << '\t' << _addID(svs[i].svt) << '\t' << _addOrientation(svs[i].svt) << '\t' << svs[i].peSupport << '\t' << svs[i].srSupport << std::endl;
    }

    // Clean-up
    bam_hdr_destroy(hdr);
    sam_close(samfile);
  }
  

}

#endif
