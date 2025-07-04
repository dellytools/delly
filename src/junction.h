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

#include <htslib/sam.h>

#include "util.h"
#include "assemble.h"
#include "pangenome.h"

namespace torali
{

  template<typename TReadBp>
  inline void
    _insertJunction(TReadBp& readBp, std::size_t const seed, bam1_t* rec, int32_t const rp, int32_t const sp, bool const scleft) {
    bool fw = true;
    if (rec->core.flag & BAM_FREVERSE) fw = false;
    int32_t readStart = rec->core.pos;
    if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) readStart = -1;
    typedef typename TReadBp::mapped_type TJunctionVector;
    typename TReadBp::iterator it = readBp.find(seed);
    int32_t seqlen = readLength(rec);
    if (sp <= seqlen) {
      if (rec->core.flag & BAM_FREVERSE) {
	if (it != readBp.end()) it->second.push_back(Junction(fw, scleft, rec->core.tid, readStart, rp, seqlen - sp, rec->core.qual));
	else readBp.insert(std::make_pair(seed, TJunctionVector(1, Junction(fw, scleft, rec->core.tid, readStart, rp, seqlen - sp, rec->core.qual))));
      } else {
	if (it != readBp.end()) it->second.push_back(Junction(fw, scleft, rec->core.tid, readStart, rp, sp, rec->core.qual));
	else readBp.insert(std::make_pair(seed, TJunctionVector(1, Junction(fw, scleft, rec->core.tid, readStart, rp, sp, rec->core.qual))));
      }
    }
  }

  inline int32_t
  _selectReadStart(std::vector<Junction> const& jcvec) {
    for(uint32_t i = 0; i < jcvec.size(); ++i) {
      if (jcvec[i].rstart != -1) return jcvec[i].rstart;
    }
    return -1;
  }  
  
  // Deletion junctions
  template<typename TConfig, typename TReadBp>
  inline void
  selectDeletions(TConfig const& c, TReadBp const& readBp, std::vector<std::vector<SRBamRecord> >& br) {
    for(typename TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
      if (it->second.size() > 1) {
	int32_t rst = _selectReadStart(it->second);
	for(uint32_t i = 0; i < it->second.size(); ++i) {
	  for(uint32_t j = i+1; j < it->second.size(); ++j) {
	    if ((uint32_t) (it->second[j].seqpos - it->second[i].seqpos) > c.maxReadSep) break;
	    // Same chr, same direction, opposing soft-clips
	    if ((it->second[j].refidx == it->second[i].refidx) && (it->second[j].forward == it->second[i].forward) && (it->second[i].scleft != it->second[j].scleft)) {
	      // Min. deletion size
	      int32_t dellen = 0;
	      if (it->second[i].forward) {
		if (!it->second[i].scleft) {
		  if (it->second[i].refpos <= it->second[j].refpos) dellen = (it->second[j].refpos - it->second[i].refpos) - (it->second[j].seqpos - it->second[i].seqpos);
		  else dellen = 0;
		} else {
		  if (it->second[i].refpos <= it->second[j].refpos) dellen = 0;
		  else dellen = (it->second[i].refpos - it->second[j].refpos) + (it->second[j].seqpos - it->second[i].seqpos);
		}
	      } else {
		if (it->second[i].scleft) {
		  if (it->second[i].refpos <= it->second[j].refpos) dellen = 0;
		  else dellen = (it->second[i].refpos - it->second[j].refpos) - (it->second[j].seqpos - it->second[i].seqpos);
		} else {
		  if (it->second[i].refpos <= it->second[j].refpos) dellen = (it->second[j].refpos - it->second[i].refpos) + (it->second[j].seqpos - it->second[i].seqpos); 
		  else dellen = 0;
		}
	      }
	      if (dellen > (int32_t) c.minRefSep) {
		// Avg. qval
		int32_t qval = (int32_t) (((int32_t) it->second[i].qual + (int32_t) it->second[j].qual) / 2);
		// Correct clipping architecture, note: soft-clipping of error-prone reads can lead to switching left/right breakpoints
		if (it->second[i].refpos <= it->second[j].refpos) {
		  if ((!it->second[i].scleft) && (it->second[j].scleft)) {
		    br[2].push_back(SRBamRecord(it->second[i].refidx, it->second[i].refpos, it->second[j].refidx, it->second[j].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  }
		} else {
		  if ((it->second[i].scleft) && (!it->second[j].scleft)) {
		    br[2].push_back(SRBamRecord(it->second[j].refidx, it->second[j].refpos, it->second[i].refidx, it->second[i].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  }
		}
		break; // Only report first SV split to avoid ambiguities
	      }
	    }
	  }
	}
      }
    }
  }


  // Duplication junctions
  template<typename TConfig, typename TReadBp>
  inline void
  selectDuplications(TConfig const& c, TReadBp const& readBp, std::vector<std::vector<SRBamRecord> >& br) {
    for(typename TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
      if (it->second.size() > 1) {
	int32_t rst = _selectReadStart(it->second);
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
		    br[3].push_back(SRBamRecord(it->second[i].refidx, it->second[i].refpos, it->second[j].refidx, it->second[j].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  }
		} else {
		  if ((!it->second[i].scleft) && (it->second[j].scleft)) {
		    br[3].push_back(SRBamRecord(it->second[j].refidx, it->second[j].refpos, it->second[i].refidx, it->second[i].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
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
  template<typename TConfig, typename TReadBp>
  inline void
  selectInversions(TConfig const& c, TReadBp const& readBp, std::vector<std::vector<SRBamRecord> >& br) {
    for(typename TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
      if (it->second.size() > 1) {
	int32_t rst = _selectReadStart(it->second);
	for(uint32_t i = 0; i < it->second.size(); ++i) {
	  for(uint32_t j = i+1; j < it->second.size(); ++j) {
	    if ((uint32_t) (it->second[j].seqpos - it->second[i].seqpos) > c.maxReadSep) break;
	    // Same chr, different direction, agreeing soft-clips
	    if ((it->second[j].refidx == it->second[i].refidx) && (it->second[j].forward != it->second[i].forward) && (it->second[i].scleft == it->second[j].scleft)) {
	      // Min. inversion size
	      if ( (uint32_t) std::abs(it->second[j].refpos - it->second[i].refpos) > c.minRefSep) {
		// Avg. qval
		int32_t qval = (int32_t) (((int32_t) it->second[i].qual + (int32_t) it->second[j].qual) / 2);
		if (it->second[i].refpos <= it->second[j].refpos) {
		  // Need to differentiate 3to3 and 5to5
		  if (it->second[i].scleft) br[1].push_back(SRBamRecord(it->second[i].refidx, it->second[i].refpos, it->second[j].refidx, it->second[j].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  else br[0].push_back(SRBamRecord(it->second[i].refidx, it->second[i].refpos, it->second[j].refidx, it->second[j].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		} else {
		  // Need to differentiate 3to3 and 5to5
		  if (it->second[i].scleft) br[1].push_back(SRBamRecord(it->second[j].refidx, it->second[j].refpos, it->second[i].refidx, it->second[i].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  else br[0].push_back(SRBamRecord(it->second[j].refidx, it->second[j].refpos, it->second[i].refidx, it->second[i].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		}
	      }
	    }
	  }
	}
      }
    }
  }

  // Insertion junctions
  template<typename TConfig, typename TReadBp>
  inline void
  selectInsertions(TConfig const& c, TReadBp const& readBp, std::vector<std::vector<SRBamRecord> >& br) {
    for(typename TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
      if (it->second.size() > 1) {
	int32_t rst = _selectReadStart(it->second);
	for(uint32_t i = 0; i < it->second.size(); ++i) {
	  for(uint32_t j = i+1; j < it->second.size(); ++j) {
	    // Same chr, same direction, opposing soft-clips
	    if ((it->second[j].refidx == it->second[i].refidx) && (it->second[j].forward == it->second[i].forward) && (it->second[i].scleft != it->second[j].scleft)) {
	      // Reference insertion footprint should be small
	      if ( (uint32_t) std::abs(it->second[j].refpos - it->second[i].refpos) < c.maxReadSep) {
		// Large separation in sequence space
		int32_t isizelen = 0;
		if (it->second[i].forward) {
		  if (!it->second[i].scleft) {
		    if (it->second[i].refpos <= it->second[j].refpos) isizelen = (it->second[j].seqpos - it->second[i].seqpos) - (it->second[j].refpos - it->second[i].refpos);
		    else isizelen = (it->second[j].seqpos - it->second[i].seqpos) + (it->second[i].refpos - it->second[j].refpos);
		  } else isizelen = 0;
		} else {
		  if (it->second[i].scleft) {
		    if (it->second[i].refpos <= it->second[j].refpos) isizelen = (it->second[j].seqpos - it->second[i].seqpos) + (it->second[j].refpos - it->second[i].refpos);
		    else isizelen = (it->second[j].seqpos - it->second[i].seqpos) - (it->second[i].refpos - it->second[j].refpos);
		  } else isizelen = 0;
		}
		if ((isizelen > (int32_t) c.minRefSep) && (isizelen < std::max(it->second[i].seqpos, it->second[j].seqpos))) {
		  // Avg. qval
		  int32_t qval = (int32_t) (((int32_t) it->second[i].qual + (int32_t) it->second[j].qual) / 2);
		  if (it->second[i].refpos <= it->second[j].refpos) {
		    br[4].push_back(SRBamRecord(it->second[i].refidx, it->second[i].refpos, it->second[j].refidx, it->second[j].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, isizelen, it->first));
		  } else {
		    br[4].push_back(SRBamRecord(it->second[j].refidx, it->second[j].refpos, it->second[i].refidx, it->second[i].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, isizelen, it->first));
		  }
		  break; // Only report first SV split to avoid ambiguities
		}
	      }
	    }
	  }
	}
      }
    }
  }


  // Translocation junctions
  template<typename TConfig, typename TReadBp>
  inline void
  selectTranslocations(TConfig const& c, TReadBp const& readBp, std::vector<std::vector<SRBamRecord> >& br) {
    for(typename TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
      if (it->second.size() > 1) {
	int32_t rst = _selectReadStart(it->second);
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
		    // 3to5
		    br[DELLY_SVT_TRANS + 2].push_back(SRBamRecord(it->second[chr2ev].refidx, it->second[chr2ev].refpos, it->second[chr1ev].refidx, it->second[chr1ev].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  } else {
		    // 5to3
		    br[DELLY_SVT_TRANS + 3].push_back(SRBamRecord(it->second[chr2ev].refidx, it->second[chr2ev].refpos, it->second[chr1ev].refidx, it->second[chr1ev].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  }
		}
	      } else {
		// Opposing direction, same soft-clips
		if (it->second[chr1ev].scleft == it->second[chr2ev].scleft) {
		  if (it->second[chr1ev].scleft) {
		    // 5to5
		    br[DELLY_SVT_TRANS + 1].push_back(SRBamRecord(it->second[chr2ev].refidx, it->second[chr2ev].refpos, it->second[chr1ev].refidx, it->second[chr1ev].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  } else {
		    // 3to3
		    br[DELLY_SVT_TRANS + 0].push_back(SRBamRecord(it->second[chr2ev].refidx, it->second[chr2ev].refpos, it->second[chr1ev].refidx, it->second[chr1ev].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }


  // Insertion junctions
  template<typename TReadBp>
  inline void
  bridgeInsertions(TReadBp const& readBp, std::vector<std::vector<SRBamRecord> >& br) {
    // Insertion map
    std::set<std::size_t> readIds;
    typedef std::pair<uint32_t, uint32_t> TRefCoord;
    typedef std::map<TRefCoord, int32_t> TPosInsMap;
    TPosInsMap pins;

    // Iterate all split-reads with insertions
    for(uint32_t i = 0; i < br[4].size(); ++i) {
      readIds.insert(br[4][i].id);
      for(int32_t k = br[4][i].pos; k <= br[4][i].pos2; ++k) {
	TRefCoord rc = std::make_pair(br[4][i].chr, k);
	typename TPosInsMap::iterator it = pins.find(rc);
	if (it == pins.end()) pins.insert(std::make_pair(rc, br[4][i].inslen));
	else it->second = (it->second + br[4][i].inslen) / 2;
      }
    }
    if (!pins.empty()) {
      for(typename TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
	if (readIds.find(it->first) == readIds.end()) {
	  for(uint32_t i = 0; i < it->second.size(); ++i) {
	    // Any insertion?
	    TRefCoord rc = std::make_pair(it->second[i].refidx, it->second[i].refpos);
	    if (pins.find(rc) != pins.end()) {
	      int32_t qval = (int32_t) it->second[i].qual;
	      int32_t rst = _selectReadStart(it->second);	      
	      br[4].push_back(SRBamRecord(it->second[i].refidx, it->second[i].refpos, it->second[i].refidx, it->second[i].refpos + 1, rst, it->second[i].seqpos, qval, pins[rc], it->first));
	    }
	  }
	}
      }
    }
  }
  

  template<typename TConfig, typename TValidRegion, typename TReadBp>
  inline void
  findJunctions(TConfig const& c, TValidRegion const& validRegions, TReadBp& readBp, std::set<size_t>& validSR) {
    typedef typename TValidRegion::value_type TChrIntervals;

    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile[0]);
    
    // Parse genome chr-by-chr
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Split-read scanning" << std::endl;

    // Iterate chromosomes
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      if (validRegions[refIndex].empty()) continue;

      // Collect reads from all samples
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	// Read alignments
	for(typename TChrIntervals::const_iterator vRIt = validRegions[refIndex].begin(); vRIt != validRegions[refIndex].end(); ++vRIt) {
	  hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, vRIt->lower(), vRIt->upper());
	  bam1_t* rec = bam_init1();
	  while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {

	    // Keep secondary alignments
	    if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
	    if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;

	    std::size_t seed = hash_lr(rec);
	    if ((!validSR.empty()) && (validSR.find(seed) == validSR.end())) continue;
	    //std::cerr << bam_get_qname(rec) << '\t' << seed << std::endl;
	    uint32_t rp = rec->core.pos; // reference pointer
	    uint32_t sp = 0; // sequence pointer
	    
	    // Parse the CIGAR
	    const uint32_t* cigar = bam_get_cigar(rec);
	    for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
		sp += bam_cigar_oplen(cigar[i]);
		rp += bam_cigar_oplen(cigar[i]);
	      } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
		if (bam_cigar_oplen(cigar[i]) > c.minRefSep) _insertJunction(readBp, seed, rec, rp, sp, false);
		rp += bam_cigar_oplen(cigar[i]);
		if (bam_cigar_oplen(cigar[i]) > c.minRefSep) { // Try look-ahead
		  uint32_t spOrig = sp;
		  uint32_t rpTmp = rp;
		  uint32_t spTmp = sp;
		  uint32_t dlen = bam_cigar_oplen(cigar[i]);
		  for (std::size_t j = i + 1; j < rec->core.n_cigar; ++j) {
		    if ((bam_cigar_op(cigar[j]) == BAM_CMATCH) || (bam_cigar_op(cigar[j]) == BAM_CEQUAL) || (bam_cigar_op(cigar[j]) == BAM_CDIFF)) {
		      spTmp += bam_cigar_oplen(cigar[j]);
		      rpTmp += bam_cigar_oplen(cigar[j]);
		      if ((double) (spTmp - sp) / (double) (dlen + (rpTmp - rp)) > c.indelExtension) break;
		    } else if (bam_cigar_op(cigar[j]) == BAM_CDEL) {
		      rpTmp += bam_cigar_oplen(cigar[j]);
		      if (bam_cigar_oplen(cigar[j]) > c.minRefSep) {
			// Extend deletion
			dlen += (rpTmp - rp);
			rp = rpTmp;
			sp = spTmp;
			i = j;
		      }
		    } else if (bam_cigar_op(cigar[j]) == BAM_CINS) {
		      if (bam_cigar_oplen(cigar[j]) > c.minRefSep) break; // No extension
		      spTmp += bam_cigar_oplen(cigar[j]);
		    } else break; // No extension
		  }
		  _insertJunction(readBp, seed, rec, rp, spOrig, true);
		}
	      } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
		if (bam_cigar_oplen(cigar[i]) > c.minRefSep) _insertJunction(readBp, seed, rec, rp, sp, false);
		sp += bam_cigar_oplen(cigar[i]);
		if (bam_cigar_oplen(cigar[i]) > c.minRefSep) { // Try look-ahead
		  uint32_t rpOrig = rp;
		  uint32_t rpTmp = rp;
		  uint32_t spTmp = sp;
		  uint32_t ilen = bam_cigar_oplen(cigar[i]);
		  for (std::size_t j = i + 1; j < rec->core.n_cigar; ++j) {
		    if ((bam_cigar_op(cigar[j]) == BAM_CMATCH) || (bam_cigar_op(cigar[j]) == BAM_CEQUAL) || (bam_cigar_op(cigar[j]) == BAM_CDIFF)) {
		      spTmp += bam_cigar_oplen(cigar[j]);
		      rpTmp += bam_cigar_oplen(cigar[j]);
		      if ((double) (rpTmp - rp) / (double) (ilen + (spTmp - sp)) > c.indelExtension) break;
		    } else if (bam_cigar_op(cigar[j]) == BAM_CDEL) {
		      if (bam_cigar_oplen(cigar[j]) > c.minRefSep) break; // No extension
		      rpTmp += bam_cigar_oplen(cigar[j]);
		    } else if (bam_cigar_op(cigar[j]) == BAM_CINS) {
		      spTmp += bam_cigar_oplen(cigar[j]);
		      if (bam_cigar_oplen(cigar[j]) > c.minRefSep) {
			// Extend insertion
			ilen += (spTmp - sp);
			rp = rpTmp;
			sp = spTmp;
			i = j;
		      }
		    } else {
		      break; // No extension
		    }
		  }
		  _insertJunction(readBp, seed, rec, rpOrig, sp, true);
		}
	      } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
		rp += bam_cigar_oplen(cigar[i]);
	      } else if ((bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) {
		int32_t finalsp = sp;
		bool scleft = false;
		if (sp == 0) {
		  finalsp += bam_cigar_oplen(cigar[i]); // Leading soft-clip / hard-clip
		  scleft = true;
		}
		sp += bam_cigar_oplen(cigar[i]);
		//std::cerr << bam_get_qname(rec) << ',' << rp << ',' << finalsp << ',' << scleft << std::endl;
		if (bam_cigar_oplen(cigar[i]) > c.minClip) _insertJunction(readBp, seed, rec, rp, finalsp, scleft);
	      } else {
		std::cerr << "Unknown Cigar options" << std::endl;
	      }
	    }
	  }
	  bam_destroy1(rec);
	  hts_itr_destroy(iter);
	}
      }
    }

    // Sort junctions
    for(typename TReadBp::iterator it = readBp.begin(); it != readBp.end(); ++it) std::sort(it->second.begin(), it->second.end());

    // Clean-up
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }


  template<typename TConfig, typename TReadBp>
  inline void
  findL1(TConfig const& c, TReadBp& readBp) {
    std::string line1 = boost::to_upper_copy(std::string("gggggaggagccaagatggccgaataggaacagctccggtctacagctcccagcgtgagcgacgcagaagacggtgatttctgcatttccatctgaggtaccgggttcatctcactagggagtgccagacagtgggcgcaggccagtgtgtgtgcgcaccgtgcgcgagccgaagcagggcgaggcattgcctcacctgggaagcgcaaggggtcagggagttccctttctgagtcaaagaaaggggtgacggtcgcacctggaaaatcgggtcactcccacccgaatattgcgcttttcagaccggcttaagaaacggcgcaccacgagactatatcccacacctggctcggagggtcctacgcccacggaatctcgctgattgctagcacagcagtctgagatcaaactgcaaggcggcaacgaggctgggggaggggcgcccgccattgcccaggcttgcttaggtaaacaaagcagccgggaagctcgaactgggtggagcccaccacagctcaaggaggcctgcctgcctctgtaggctccacctctgggggcagggcacagacaaacaaaaagacagcagtaacctctgcagacttaagtgtccctgtctgacagctttgaagagagcagtggttctcccagcacgcagctggagatctgagaacgggcagacagactgcctcctcaagtgggtccctgactcctgacccccgagcagcctaactgggaggcaccccccagcaggggcacactgacacctcacacggcagggtattccaacagacctgcagctgagggtcctgtctgttagaaggaaaactaacaaccagaaaggacatctacaccgaaaacccatctgtacatcaccatcatcaaagaccaaaagtagataaaaccacaaagatggggaaaaaacagaacagaaaaactggaaactctaaaacgcagagcgcctctcctcctccaaaggaacgcagttcctcaccagcaacggaacaaagctggatggagaatgattttgacgagctgagagaagaaggcttcagacgatcaaattactctgagctacgggaggacattcaaaccaaaggcaaagaagttgaaaactttgaaaaaaatttagaagaatgtataactagaataaccaatacagagaagtgcttaaaggagctgatggagctgaaaaccaaggctcgagaactacgtgaagaatgcagaagcctcaggagccgatgcgatcaactggaagaaagggtatcagcaatggaagatgaaatgaatgaaatgaagcgagaagggaagtttagagaaaaaagaataaaaagaaatgagcaaagcctccaagaaatatgggactatgtgaaaagaccaaatctacgtctgattggtgtacctgaaagtgatgtggagaatggaaccaagttggaaaacactctgcaggatattatccaggagaacttccccaatctagcaaggcaggccaacgttcagattcaggaaatacagagaacgccacaaagatactcctcgagaagagcaactccaagacacataattgtcagattcaccaaagttgaaatgaaggaaaaaatgttaagggcagccagagagaaaggtcgggttaccctcaaaggaaagcccatcagactaacagtggatctctcggcagaaaccctacaagccagaagagagtgggggccaatattcaacattcttaaagaaaagaattttcaacccagaatttcatatccagccaaactaagcttcataagtgaaggagaaataaaatactttatagacaagcaaatgttgagagattttgtcaccaccaggcctgccctaaaagagctcctgaaggaagcgctaaacatggaaaggaacaaccggtaccagccgctgcaaaatcatgccaaaatgtaaagaccatcgagactaggaagaaactgcatcaactaatgagcaaaatcaccagctaacatcataatgacaggatcaaattcacacataacaatattaactttaaatataaatggactaaattctgcaattaaaagacacagactggcaagttggataaagagtcaagacccatcagtgtgctgtattcaggaaacccatctcacgtgcagagacacacataggctcaaaataaaaggatggaggaagatctaccaagccaatggaaaacaaaaaaaggcaggggttgcaatcctagtctctgataaaacagactttaaaccaacaaagatcaaaagagacaaagaaggccattacataatggtaaagggatcaattcaacaagaggagctaactatcctaaatatttatgcacccaatacaggagcacccagattcataaagcaagtcctcagtgacctacaaagagacttagactcccacacattaataatgggagactttaacaccccactgtcaacattagacagatcaacgagacagaaagtcaacaaggatacccaggaattgaactcagctctgcaccaagcagacctaatagacatctacagaactctccaccccaaatcaacagaatatacctttttttcagcaccacaccacacctattccaaaattgaccacatagttggaagtaaagctctcctcagcaaatgtaaaagaacagaaattataacaaactatctctcagaccacagtgcaatcaaactagaactcaggattaagaatctcactcaaagccgctcaactacatggaaactgaacaacctgctcctgaatgactactgggtacataacgaaatgaaggcagaaataaagatgttctttgaaaccaacgagaacaaagacaccacataccagaatctctgggacgcattcaaagcagtgtgtagagggaaatttatagcactaaatgcctacaagagaaagcaggaaagatccaaaattgacaccctaacatcacaattaaaagaactagaaaagcaagagcaaacacattcaaaagctagcagaaggcaagaaataactaaaatcagagcagaactgaaggaaatagagacacaaaaaacccttcaaaaaatcaatgaatccaggagctggttttttgaaaggatcaacaaaattgatagaccgctagcaagactaataaagaaaaaaagagagaagaatcaaatagacacaataaaaaatgataaaggggatatcaccaccgatcccacagaaatacaaactaccatcagagaatactacaaacacctctacgcaaataaactagaaaatctagaagaaatggatacattcctcgacacatacactctcccaagactaaaccaggaagaagttgaatctctgaatagaccaataacaggctctgaaattgtggcaataatcaatagtttaccaaccaaaaagagtccaggaccagatggattcacagccgaattctaccagaggtacatggaggaactggtaccattccttctgaaactattccaatcaatagaaaaagagggaatcctccctaactcattttatgaggccagcatcattctgataccaaagccgggcagagacacaaccaaaaaagagaattttagaccaatatccttgatgaacattgatgcaaaaatcctcaataaaatactggcaaaccgaatccagcagcacatcaaaaagcttatccaccatgatcaagtgggcttcatccctgggatgcaaggctggttcaatatacgcaaatcaataaatgtaatccagcatataaacagagccaaagacaaaaaccacatgattatctcaatagatgcagaaaaagcctttgacaaaattcaacaacccttcatgctaaaaactctcaataaattaggtattgatgggacgtatttcaaaataataagagctatctatgacaaacccacagccaatatcatactgaatgggcaaaaactggaagcattccctttgaaaaccggcacaagacagggatgccctctctcaccgctcctattcaacatagtgttggaagttctggccagggcaatcaggcaggagaaggaaataaagggtattcaattaggaaaagaggaagtcaaattgtccctgtttgcagacgacatgattgtatatctagaaaaccccatcgtctcagcccaaaatctccttaagctgataagcaacttcagcaaagtctcaggatacaaaatcaatgtacaaaaatcacaagcattcttatacaccaacaacagacaaacagagagccaaatcatgggtgaactcccattcgtaattgcttcaaagagaataaaatacctaggaatccaacttacaagggatgtgaaggacctcttcaaggagaactacaaaccactgctcaaggaaataaaagaggacacaaacaaatggaagaacattccatgctcatgggtaggaagaatcaatatcgtgaaaatggccatactgcccaaggtaatttacagattcaatgccatccccatcaagctaccaatgactttcttcacagaattggaaaaaactactttaaagttcatatggaaccaaaaaagagcccgcattgccaagtcaatcctaagccaaaagaacaaagctggaggcatcacactacctgacttcaaactatactacaaggctacagtaaccaaaacagcatggtactggtaccaaaacagagatatagatcaatggaacagaacagagccctcagaaataatgccgcatatctacaactatctgatctttgacaaacctgagaaaaacaagcaatggggaaaggattccctatttaataaatggtgctgggaaaactggctagccatatgtagaaagctgaaactggatcccttccttacaccttatacaaaaatcaattcaagatggattaaagatttaaacgttaaacctaaaaccataaaaaccctagaagaaaacctaggcattaccattcaggacataggcgtgggcaaggacttcatgtccaaaacaccaaaagcaatggcaacaaaagacaaaattgacaaatgggatctaattaaactaaagagcttctgcacagcaaaagaaactaccatcagagtgaacaggcaacctacaacatgggagaaaattttcgcaacctactcatctgacaaagggctaatatccagaatctacaatgaacttaaacaaatttacaagaaaaaaacaaacaaccccatcaaaaagtgggcgaaggacatgaacagacacttctcaaaagaagacatttatgcagccaaaaaacacatgaagaaatgctcatcatcactggccatcagagaaatgcaaatcaaaaccactatgagatatcatctcacaccagttagaatggcaatcattaaaaagtcaggaaacaacaggtgctggagaggatgcggagaaataggaacacttttacactgttggtgggactgtaaactagttcaaccattgtggaagtcagtgtggcgattcctcagggatctagaactagaaataccatttgacccagccatcccattactgggtatatacccaaatgagtataaatcatgctgctataaagacacatgcacacgtatgtttattgcggcactattcacaatagcaaagacttggaaccaacccaaatgtccaacaatgatagactggattaagaaaatgtggcacatatacaccatggaatactatgcagccataaaaaatgatgagttcatatcctttgtagggacatggatgaaattggaaaccatcattctcagtaaactatcgcaagaacaaaaaaccaaacaccgcatattctcactcataggtgggaattgaacaatgagatcacatggacacaggaaggggaatatcacactctggggactgtggtggggtcgggggaggggggagggatagcattgggagatatacctaatgctagatgacacattagtgggtgcagcgcaccagcatggcacatgtatacatatgtaactaacctgcacaatgtgcacatgtaccctaaaacttagagtat"));
    std::string line1Rev = line1;
    reverseComplement(line1Rev);

    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile[0]);
    
    // Parse genome chr-by-chr
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "L1 serach" << std::endl;

    // L1 candidate read search
    std::set<std::size_t> clusteredReads;
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet bp(hdr->target_len[refIndex], false);
      TBitSet firstHit(hdr->target_len[refIndex], false);
      TBitSet secondHit(hdr->target_len[refIndex], false);

      // Cluster breakpoints +/-5bp
      for(typename TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
	if (it->second.size() > 1) {
	  for(uint32_t i = 0; i < it->second.size(); ++i) {
	    if (it->second[i].refidx == refIndex) {
	      int32_t refStart = 0;
	      if (it->second[i].refpos > 5) refStart = it->second[i].refpos - 5;
	      int32_t refEnd = hdr->target_len[refIndex];
	      if (it->second[i].refpos + 5 < refEnd) refEnd = it->second[i].refpos + 5;
	      for(int32_t k = refStart; k < refEnd; ++k) {
		if (secondHit[k]) bp[k] = true;
		else {
		  if (firstHit[k]) secondHit[k] = true;
		  else firstHit[k] = true;
		}
	      }
	    }
	  }
	}
      }
      
      // 2nd pass to fetch reads
      for(typename TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
	if (it->second.size() > 1) {
	  for(uint32_t i = 0; i < it->second.size(); ++i) {
	    if (it->second[i].refidx == refIndex) {
	      int32_t refStart = 0;
	      if (it->second[i].refpos > 5) refStart = it->second[i].refpos - 5;
	      int32_t refEnd = hdr->target_len[refIndex];
	      if (it->second[i].refpos + 5 < refEnd) refEnd = it->second[i].refpos + 5;
	      for(int32_t k = refStart; k < refEnd; ++k) {
		if (bp[k]) clusteredReads.insert(it->first);
	      }
	    }
	  }
	}
      }
      //std::cerr << refIndex << ',' << clusteredReads.size() << std::endl;
    }

    // Find L1
    if (!clusteredReads.empty()) {
      for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
	// Collect reads from all samples
	for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	  // Read alignments
	  hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr->target_len[refIndex]);
	  bam1_t* rec = bam_init1();
	  while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	    // Keep only primary alignments
	    if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
	    std::size_t seed = hash_lr(rec);

	    if (clusteredReads.find(seed) != clusteredReads.end()) {
	      if (readBp[seed].size()>1) {
		std::string sequence;
		
		// Get all sequences embedded in breakpoints
		for(uint32_t k = 1; k < readBp[seed].size(); ++k) {
		  int32_t fragsize = readBp[seed][k].seqpos - readBp[seed][k-1].seqpos;
		  if ((fragsize > 500) && (fragsize < 7000)) {
		    // L1 fragment?
		    if (sequence.empty()) {
		      sequence.resize(rec->core.l_qseq);
		      const uint8_t* seqptr = bam_get_seq(rec);
		      for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
		      // Reverse complement if necessary
		      if (rec->core.flag & BAM_FREVERSE) reverseComplement(sequence);
		    }

		    // Get subsequence
		    std::string subseq = sequence.substr(readBp[seed][k-1].seqpos, fragsize);
		    //std::cerr << ">" << bam_get_qname(rec) << "_" << readBp[seed][k-1].seqpos << "_" << subseq.size() << std::endl;
		    //std::cerr << subseq << std::endl;

		    // Align to L1
		    EdlibAlignResult cigarFwd = edlibAlign(subseq.c_str(), subseq.size(), line1.c_str(), line1.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
		    int32_t alignStartFwd = infixStart(cigarFwd);
		    int32_t alignEndFwd = infixEnd(cigarFwd);
		    double percentIdentityFwd = 1.0 - (double) (cigarFwd.editDistance) / (double) (alignEndFwd - alignStartFwd);
		  
		    // Align to L1 (reverse)
		    EdlibAlignResult cigarRev = edlibAlign(subseq.c_str(), subseq.size(), line1Rev.c_str(), line1Rev.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
		    int32_t alignStartRev = infixStart(cigarRev);
		    int32_t alignEndRev = infixEnd(cigarRev);
		    double percentIdentityRev = 1.0 - (double) (cigarRev.editDistance) / (double) (alignEndRev - alignStartRev);
		  
		    if ( ((percentIdentityFwd > 0.9) && ((alignEndFwd - alignStartFwd) > 500)) || ((percentIdentityRev > 0.9) && ((alignEndRev - alignStartRev) > 500)) ) {
		      std::string l1Type = "L1 insertion";
		      if (readBp[seed][k-1].refidx != readBp[seed][k].refidx) {
			l1Type = "Inter-chromosomal SV with L1 fragment";
		      } else {
			int32_t offset = std::abs(readBp[seed][k].refpos - readBp[seed][k-1].refpos);
			if (offset > 1000) {
			  l1Type = "Intra-chromosomal SV with L1 fragment";
			}
		      }
		      std::cout << l1Type << '\t';
		      std::cout << bam_get_qname(rec) << '\t';
		      std::cout << hdr->target_name[readBp[seed][k-1].refidx] << ':' << readBp[seed][k-1].refpos << " (ReadPos: " << readBp[seed][k-1].seqpos << ')' << '\t';
		      std::cout << hdr->target_name[readBp[seed][k].refidx] << ':' << readBp[seed][k].refpos << " (ReadPos: " << readBp[seed][k].seqpos << ')' << '\t';
		      if ((percentIdentityFwd > 0.9) && ((alignEndFwd - alignStartFwd) > 500)) {
			std::cout << alignStartFwd << '-' << alignEndFwd << '\t';
			std::cout << (alignEndFwd - alignStartFwd) << '\t';
			std::cout << "fwd" << '\t';
			std::cout << percentIdentityFwd << std::endl;
			//printAlignment(subseq, line1, EDLIB_MODE_HW, cigarFwd);
		      } else {
			std::cout << alignStartRev << '-' << alignEndRev << '\t';
			std::cout << (alignEndRev - alignStartRev) << '\t';
			std::cout << "rev" << '\t';
			std::cout << percentIdentityRev << std::endl;
			//printAlignment(subseq, line1Rev, EDLIB_MODE_HW, cigarRev);
		      }
		    }
		    edlibFreeAlignResult(cigarFwd);
		    edlibFreeAlignResult(cigarRev);
		  }
		}
	      }
	    }
	  }
	  bam_destroy1(rec);
	  hts_itr_destroy(iter);
	}
      }
    }

    // Clean-up
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }   
  
  template<typename TConfig, typename TReadBp>
  inline void
  fetchSVs(TConfig const& c, TReadBp& readBp, std::vector<std::vector<SRBamRecord> >& br) {
    //std::cout << "L1InsType\tReadName\tRefCoordBefore\tRefCoordAfter\tL1Coord\tL1FragmentSize\tL1Direction\tL1Similarity" << std::endl;
    //findL1(c, readBp);
    // Extract BAM records
    if ((c.svtset.empty()) || (c.svtset.find(2) != c.svtset.end())) selectDeletions(c, readBp, br);
    if ((c.svtset.empty()) || (c.svtset.find(3) != c.svtset.end())) selectDuplications(c, readBp, br);
    if ((c.svtset.empty()) || (c.svtset.find(0) != c.svtset.end()) || (c.svtset.find(1) != c.svtset.end())) selectInversions(c, readBp, br);
    if ((c.svtset.empty()) || (c.svtset.find(4) != c.svtset.end())) {
      selectInsertions(c, readBp, br);
      bridgeInsertions(readBp, br);
    }
    if ((c.svtset.empty()) || (c.svtset.find(5) != c.svtset.end()) || (c.svtset.find(6) != c.svtset.end()) || (c.svtset.find(7) != c.svtset.end()) || (c.svtset.find(8) != c.svtset.end())) selectTranslocations(c, readBp, br);
  }

  template<typename TConfig, typename TValidRegions, typename TSvtSRBamRecord>
  inline void
  _findSRBreakpoints(TConfig const& c, TValidRegions const& validRegions, TSvtSRBamRecord& srBR, std::set<size_t>& validSR) {
    // Breakpoints
    typedef std::vector<Junction> TJunctionVector;
    typedef std::map<std::size_t, TJunctionVector> TReadBp;
    TReadBp readBp;
    findJunctions(c, validRegions, readBp, validSR);
    fetchSVs(c, readBp, srBR);
  }

  template<typename TConfig, typename TValidRegions, typename TSvtSRBamRecord>
  inline void
  _findSRBreakpoints(TConfig const& c, TValidRegions const& validRegions, TSvtSRBamRecord& srBR) {
    std::set<size_t> validSR;
    _findSRBreakpoints(c, validRegions, srBR, validSR);
  }

  template<typename TConfig, typename TValidRegions, typename TSVs, typename TSRStore>
  inline void
  _clusterSRReads(TConfig const& c, TValidRegions const& validRegions, TSVs& svc, TSRStore& srStore) {
    // Split-reads
    typedef std::vector<SRBamRecord> TSRBamRecord;
    typedef std::vector<TSRBamRecord> TSvtSRBamRecord;
    TSvtSRBamRecord srBR(2 * DELLY_SVT_TRANS, TSRBamRecord());
    _findSRBreakpoints(c, validRegions, srBR);

    // Alternate alignments
    if (c.hasAltFile) {

      // Insert all split-reads
      std::set<std::size_t> validSR;
      for(uint32_t svt = 0; svt < srBR.size(); ++svt) {
	for(uint32_t i = 0; i < srBR[svt].size(); ++i) validSR.insert(srBR[svt][i].id);
      }
      
      // Parse alternate alignments
      std::vector<boost::filesystem::path> align;
      std::vector<boost::filesystem::path> genome;
      _alternateAlignments(c, align, genome);
      for(uint32_t i = 0; i < align.size(); ++i) {
	std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Parsing alternate alignment " << align[i] << std::endl;
	
	TSvtSRBamRecord altSR(2 * DELLY_SVT_TRANS, TSRBamRecord());
	TConfig altConfig(c);
	altConfig.hasExcludeFile = false;
	altConfig.exclude.clear();
	altConfig.genome = genome[i];
	altConfig.files.clear();
	altConfig.files.push_back(align[i]);
	if (isBamCram(align[i].string())) {
	  // Alternate alignment in BAM/CRAM format
	  samFile* samfile = sam_open(altConfig.files[0].string().c_str(), "r");
	  bam_hdr_t* hdr = sam_hdr_read(samfile);
	  TValidRegions vR;
	  _parseExcludeIntervals(altConfig, hdr, vR);
	  bam_hdr_destroy(hdr);
	  sam_close(samfile);
	  _findSRBreakpoints(altConfig, vR, altSR, validSR);

	  // Debug
	  //std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "AltLinearAlign: " << align[i] << std::endl;
	  //outputSRBamRecords(altConfig, altSR, true);
	} else {
	  Graph g;
	  parseGfa(altConfig, g);
	  altConfig.nchr = g.smap.size();

	  // Minigraph lacks alignment extension so larger BP separation expected?
	  altConfig.maxReadSep = 2000;
	  
	  // Alternate alignment in GAF format
	  _findGraphSRBreakpoints(altConfig, g, altSR, validSR);

	  // Debug
	  //std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "AltGraphAlign: " << align[i] << std::endl;
	  //outputGraphSRBamRecords(altConfig, g, altSR);
	}

	// Update valid SR set
	std::set<std::size_t> newValidSR;
	for(uint32_t svt = 0; svt < altSR.size(); ++svt) {
	  //std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "SVT Size: " << svt << ',' << altSR[svt].size() << std::endl;
	  if (altSR[svt].empty()) continue;
	  for(uint32_t ik = 0; ik < altSR[svt].size(); ++ik) {
	    if (validSR.find(altSR[svt][ik].id) != validSR.end()) newValidSR.insert(altSR[svt][ik].id);
	  }
	}
	validSR = newValidSR;

	// Debug
	//std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "SR size: " << validSR.size() << std::endl;
      }

      // Filter primary SR records
      uint32_t origSRsize = 0;
      uint32_t newSRsize = 0;
      TSvtSRBamRecord copySR(2 * DELLY_SVT_TRANS, TSRBamRecord());
      for(uint32_t svt = 0; svt < srBR.size(); ++svt) {
	origSRsize += srBR[svt].size();
	copySR[svt] = srBR[svt];
	srBR[svt].clear();
	for(uint32_t i = 0; i < copySR[svt].size(); ++i) {
	  if (validSR.find(copySR[svt][i].id) != validSR.end()) srBR[svt].push_back(copySR[svt][i]);
	}
	newSRsize += srBR[svt].size();
      }
      uint32_t filtSR = origSRsize - newSRsize;
      double filtRatio = 0;
      if (origSRsize) filtRatio = (double) filtSR / (double) origSRsize;
      std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Filtered " << filtSR << " out of " << origSRsize << " split-read records (" << (filtRatio * 100) << "%)" << std::endl;      
    }


    // Debug
    //std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "MainAlign: " << std::endl;
    //outputSRBamRecords(c, srBR, true);

    // Cluster BAM records
    for(uint32_t svt = 0; svt < srBR.size(); ++svt) {
      if (srBR[svt].empty()) continue;
      
      // Sort
      std::sort(srBR[svt].begin(), srBR[svt].end());
      
      // Cluster
      cluster(c, srBR[svt], svc, svt);

      // Debug
      //outputStructuralVariants(c, svc, srBR, svt, true);
      // Track split-reads
      samFile* samfile = sam_open(c.files[0].string().c_str(), "r");
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      for(uint32_t i = 0; i < srBR[svt].size(); ++i) {
	// Read assigned?
	if ((srBR[svt][i].svid != -1) && (srBR[svt][i].rstart != -1)) {
	  if (srBR[svt][i].rstart < (int32_t) hdr->target_len[srBR[svt][i].chr]) {
	    if (srStore[srBR[svt][i].chr].find(std::make_pair(srBR[svt][i].rstart, srBR[svt][i].id)) == srStore[srBR[svt][i].chr].end()) srStore[srBR[svt][i].chr].insert(std::make_pair(std::make_pair(srBR[svt][i].rstart, srBR[svt][i].id), std::vector<SeqSlice>()));
	    srStore[srBR[svt][i].chr][std::make_pair(srBR[svt][i].rstart, srBR[svt][i].id)].push_back(SeqSlice(srBR[svt][i].svid, srBR[svt][i].sstart, srBR[svt][i].inslen, srBR[svt][i].qual));
	  }
	  if (srBR[svt][i].chr != srBR[svt][i].chr2) {
	    // Unclear which chr was primary alignment so insert both if and only if rstart < reference length
	    if (srBR[svt][i].rstart < (int32_t) hdr->target_len[srBR[svt][i].chr2]) {
	      if (srStore[srBR[svt][i].chr2].find(std::make_pair(srBR[svt][i].rstart, srBR[svt][i].id)) == srStore[srBR[svt][i].chr2].end()) srStore[srBR[svt][i].chr2].insert(std::make_pair(std::make_pair(srBR[svt][i].rstart, srBR[svt][i].id), std::vector<SeqSlice>()));
	      srStore[srBR[svt][i].chr2][std::make_pair(srBR[svt][i].rstart, srBR[svt][i].id)].push_back(SeqSlice(srBR[svt][i].svid, srBR[svt][i].sstart, srBR[svt][i].inslen, srBR[svt][i].qual));
	    }
	  }
	}
      }
      bam_hdr_destroy(hdr);
      sam_close(samfile);
    }
  }


  template<typename TConfig>
  inline void
  outputSRBamRecords(TConfig const& c, std::vector<std::vector<SRBamRecord> > const& br, bool const longread) {
    // Header
    std::cerr << "qname\tid\tchr1\tpos1\tchr2\tpos2\tsvtype\tct\tqual\tseqpos\tinslen" << std::endl;

    // Hash reads
    typedef std::map<std::size_t, std::string> THashMap;
    THashMap hm;
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile[0]);
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
	  std::size_t seed;
	  if (longread) seed = hash_lr(rec);
	  else seed = hash_sr(rec);
	  std::string qname = bam_get_qname(rec);
	  if (hm.find(seed) == hm.end()) hm.insert(std::make_pair(seed, qname));
	  else {
	    if (hm[seed] != qname) {
	      std::cerr << "Warning: Hash collision! " << seed << ',' << hm[seed] << ',' << qname << std::endl;
	    }
	  }
	}
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }
    }
    
    // SVs
    for(uint32_t svt = 0; svt < br.size(); ++svt) {
      for(uint32_t i = 0; i < br[svt].size(); ++i) {
	std::cerr << hm[br[svt][i].id] << '\t' << br[svt][i].id << '\t' << hdr->target_name[br[svt][i].chr] << '\t' << br[svt][i].pos << '\t' << hdr->target_name[br[svt][i].chr2] << '\t' << br[svt][i].pos2 << '\t' << _addID(svt) << '\t' << _addOrientation(svt) << '\t' << br[svt][i].qual << '\t' << br[svt][i].sstart << '\t' << br[svt][i].inslen << std::endl;
      }
    }

    // Clean-up
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }

  template<typename TConfig, typename TSvtSRBamRecord>
  inline void
  outputStructuralVariants(TConfig const& c, std::vector<StructuralVariantRecord> const& svs, TSvtSRBamRecord const& srBR, int32_t const svt, bool const longread) {
    // Header
    std::cerr << "chr1\tpos1\tchr2\tpos2\tsvtype\tct\tinslen\tpeSupport\tsrSupport" << std::endl;
    
    // Hash reads
    typedef std::map<std::size_t, std::string> THashMap;
    THashMap hm;
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile[0]);
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
	  std::size_t seed;
	  if (longread) seed = hash_lr(rec);
	  else seed = hash_sr(rec);
	  std::string qname = bam_get_qname(rec);
	  if (hm.find(seed) == hm.end()) hm.insert(std::make_pair(seed, qname));
	  else {
	    if (hm[seed] != qname) {
	      std::cerr << "Warning: Hash collision! " << seed << ',' << hm[seed] << ',' << qname << std::endl;
	    }
	  }
	}
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }
    }
    
    // Track split-reads
    typedef std::vector<std::string> TReadNameVector;
    typedef std::vector<TReadNameVector> TSVReadNames;
    TSVReadNames svReadNames(svs.size(), TReadNameVector());
    for(uint32_t i = 0; i < srBR[svt].size(); ++i) {
      if (srBR[svt][i].svid != -1) {
	svReadNames[srBR[svt][i].svid].push_back(hm[srBR[svt][i].id]);
      }
    }
  
    // SVs
    for(uint32_t i = 0; i < svs.size(); ++i) {
      if (svs[i].svt != svt) continue;
      std::cerr << hdr->target_name[svs[i].chr] << '\t' << svs[i].svStart << '\t' << hdr->target_name[svs[i].chr2] << '\t' << svs[i].svEnd << '\t' << _addID(svs[i].svt) << '\t' << _addOrientation(svs[i].svt) << '\t' << svs[i].insLen << '\t' << svs[i].peSupport << '\t' << svs[i].srSupport << '\t';
      for(uint32_t k = 0; k < svReadNames[svs[i].id].size(); ++k) std::cerr << svReadNames[svs[i].id][k] << ',';
      std::cerr << std::endl;
    }

    // Clean-up
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }


}

#endif
