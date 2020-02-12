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
#include <boost/progress.hpp>

#include <htslib/sam.h>

#include "util.h"

namespace torali
{

  struct SRBamRecord {
    int32_t chr;
    int32_t pos;
    int32_t chr2;
    int32_t pos2;
    int32_t rstart;
    int32_t sstart;
    int32_t qual;
    int32_t inslen;
    int32_t svid;
    std::size_t id;
        
    SRBamRecord(int32_t const c, int32_t const p, int32_t const c2, int32_t const p2, int32_t const rst, int32_t const sst, int32_t const qval, int32_t const il, std::size_t const idval) : chr(c), pos(p), chr2(c2), pos2(p2), rstart(rst), sstart(sst), qual(qval), inslen(il), svid(-1), id(idval) {}
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
    int32_t rstart;
    int32_t refpos;
    int32_t seqpos;
    uint16_t qual;
    
    Junction(bool const fw, bool const cl, int32_t const idx, int32_t const rst, int32_t const r, int32_t const s, uint16_t const qval) : forward(fw), scleft(cl), refidx(idx), rstart(rst), refpos(r), seqpos(s), qual(qval) {}
  };


  template<typename TReadBp>
  inline void
    _insertJunction(TReadBp& readBp, std::size_t const seed, bam1_t* rec, int32_t const rp, int32_t const sp, bool const scleft) {
    bool fw = true;
    if (rec->core.flag & BAM_FREVERSE) fw = false;
    int32_t readStart = rec->core.pos;
    if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) readStart = -1;
    typedef typename TReadBp::mapped_type TJunctionVector;
    typename TReadBp::iterator it = readBp.find(seed);
    int32_t seqlen = sequenceLength(rec);
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

  template<typename TJunction>
  struct SortJunction : public std::binary_function<TJunction, TJunction, bool>
  {
    inline bool operator()(TJunction const& j1, TJunction const& j2) {
      return ((j1.seqpos<j2.seqpos) || ((j1.seqpos==j2.seqpos) && (j1.refidx<j2.refidx)) || ((j1.seqpos==j2.seqpos) && (j1.refidx==j2.refidx) && (j1.refpos<j2.refpos)) || ((j1.seqpos==j2.seqpos) && (j1.refidx==j2.refidx) && (j1.refpos==j2.refpos) && (j1.scleft < j2.scleft)));
    }
  };

  // Deletion junctions
  template<typename TConfig, typename TReadBp>
  inline void
  selectDeletions(TConfig const& c, TReadBp const& readBp, std::vector<std::vector<SRBamRecord> >& br) {
    for(typename TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
      if (it->second.size() > 1) {
	for(uint32_t i = 0; i < it->second.size(); ++i) {
	  for(uint32_t j = i+1; j < it->second.size(); ++j) {
	    if ((uint32_t) (it->second[j].seqpos - it->second[i].seqpos) > c.maxReadSep) break;
	    // Same chr, same direction, opposing soft-clips
	    if ((it->second[j].refidx == it->second[i].refidx) && (it->second[j].forward == it->second[i].forward) && (it->second[i].scleft != it->second[j].scleft)) {
	      // Min. deletion size
	      if ( (uint32_t) std::abs(it->second[j].refpos - it->second[i].refpos) > c.minRefSep) {
		int32_t rst = it->second[i].rstart;
		if (rst == -1) rst = it->second[j].rstart;
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
	for(uint32_t i = 0; i < it->second.size(); ++i) {
	  for(uint32_t j = i+1; j < it->second.size(); ++j) {
	    if ((uint32_t) (it->second[j].seqpos - it->second[i].seqpos) > c.maxReadSep) break;
	    // Same chr, same direction, opposing soft-clips
	    if ((it->second[j].refidx == it->second[i].refidx) && (it->second[j].forward == it->second[i].forward) && (it->second[i].scleft != it->second[j].scleft)) {
	      // Min. duplication size
	      if ( (uint32_t) std::abs(it->second[j].refpos - it->second[i].refpos) > c.minRefSep) {
		int32_t rst = it->second[i].rstart;
		if (rst == -1) rst = it->second[j].rstart;
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
	for(uint32_t i = 0; i < it->second.size(); ++i) {
	  for(uint32_t j = i+1; j < it->second.size(); ++j) {
	    if ((uint32_t) (it->second[j].seqpos - it->second[i].seqpos) > c.maxReadSep) break;
	    // Same chr, different direction, agreeing soft-clips
	    if ((it->second[j].refidx == it->second[i].refidx) && (it->second[j].forward != it->second[i].forward) && (it->second[i].scleft == it->second[j].scleft)) {
	      // Min. inversion size
	      if ( (uint32_t) std::abs(it->second[j].refpos - it->second[i].refpos) > c.minRefSep) {
		int32_t rst = it->second[i].rstart;
		if (rst == -1) rst = it->second[j].rstart;
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
	for(uint32_t i = 0; i < it->second.size(); ++i) {
	  for(uint32_t j = i+1; j < it->second.size(); ++j) {
	    // Same chr, same direction, opposing soft-clips
	    if ((it->second[j].refidx == it->second[i].refidx) && (it->second[j].forward == it->second[i].forward) && (it->second[i].scleft != it->second[j].scleft)) {
	      // Reference insertion footprint should be small
	      if ( (uint32_t) std::abs(it->second[j].refpos - it->second[i].refpos) < c.maxReadSep) {
		// Large separation in sequence space
		if ((uint32_t) (it->second[j].seqpos - it->second[i].seqpos) > c.minRefSep) {
		  int32_t rst = it->second[i].rstart;
		  if (rst == -1) rst = it->second[j].rstart;
		  // Avg. qval
		  int32_t qval = (int32_t) (((int32_t) it->second[i].qual + (int32_t) it->second[j].qual) / 2);
		  if (it->second[i].refpos <= it->second[j].refpos) {
		    br[4].push_back(SRBamRecord(it->second[i].refidx, it->second[i].refpos, it->second[j].refidx, it->second[j].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  } else {
		    br[4].push_back(SRBamRecord(it->second[j].refidx, it->second[j].refpos, it->second[i].refidx, it->second[i].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
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
  template<typename TConfig, typename TReadBp>
  inline void
  selectTranslocations(TConfig const& c, TReadBp const& readBp, std::vector<std::vector<SRBamRecord> >& br) {
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
	      int32_t rst = it->second[i].rstart;
	      if (rst == -1) rst = it->second[j].rstart;
	      // Avg. qval
	      int32_t qval = (int32_t) (((int32_t) it->second[i].qual + (int32_t) it->second[j].qual) / 2);
	      if (it->second[chr1ev].forward == it->second[chr2ev].forward) {
		// Same direction, opposing soft-clips
		if (it->second[chr1ev].scleft != it->second[chr2ev].scleft) {
		  if (it->second[chr1ev].scleft) {
		    // 5to3
		    br[DELLY_SVT_TRANS + 2].push_back(SRBamRecord(it->second[chr2ev].refidx, it->second[chr2ev].refpos, it->second[chr1ev].refidx, it->second[chr1ev].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  } else {
		    // 3to5
		    br[DELLY_SVT_TRANS + 3].push_back(SRBamRecord(it->second[chr2ev].refidx, it->second[chr2ev].refpos, it->second[chr1ev].refidx, it->second[chr1ev].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  }
		}
	      } else {
		// Opposing direction, same soft-clips
		if (it->second[chr1ev].scleft == it->second[chr2ev].scleft) {
		  if (it->second[chr1ev].scleft) {
		    // 3to3
		    br[DELLY_SVT_TRANS + 1].push_back(SRBamRecord(it->second[chr2ev].refidx, it->second[chr2ev].refpos, it->second[chr1ev].refidx, it->second[chr1ev].refpos, rst, std::min(it->second[j].seqpos, it->second[i].seqpos), qval, std::abs(it->second[j].seqpos - it->second[i].seqpos), it->first));
		  } else {
		    // 5to5
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


  template<typename TConfig>
  inline int32_t
  scoreThreshold(std::string const& s1, std::string const& s2, TConfig const& c) {
    int32_t l1 = 0;
    int32_t l2 = 0;
    for(uint32_t i = 0; i < s1.size(); ++i) {
      if (s1[i] != '-') ++l1;
    }
    for(uint32_t i = 0; i < s2.size(); ++i) {
      if (s2[i] != '-') ++l2;
    }
    int32_t minSize = std::min(l1, l2);
    return (int32_t) ((c.flankQuality * minSize * c.aliscore.match) + ((1.0 - c.flankQuality) * minSize * c.aliscore.mismatch));
  }
  
  template<typename TScore>
  inline int32_t
  scoreAlign(std::string const& s1, std::string const& s2, TScore const& sc) {
    int32_t score = 0;
    bool inGap = false;
    for(uint32_t i = 0; i < s1.size(); ++i) {
      if ((s1[i] == '-') || (s2[i] == '-')) {
	if (inGap) score += sc.ge;
	else {
	  inGap = true;
	  score += (sc.go + sc.ge);
	}
      } else {
	inGap = false;
	if (s1[i] == s2[i]) score += sc.match;
	else score += sc.mismatch;
      }
    }
    return score;
  }
  

  template<typename TConfig, typename TValidRegion, typename TJunctionMap, typename TReadCountMap>
  inline void
  trackRef(TConfig const& c, TValidRegion const& validRegions, std::vector<StructuralVariantRecord>& svs, TJunctionMap& jctMap, TReadCountMap& covMap) {
    typedef typename TValidRegion::value_type TChrIntervals;
    if (svs.empty()) return;

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
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "SV annotation" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );

    // Collect ref. supporting reads
    typedef std::map<uint32_t, uint32_t> TQualMap;
    typedef std::vector<TQualMap> TSvSupport;
    TSvSupport svSup(svs.size(), TQualMap());

    // Coverage distribution
    typedef uint16_t TMaxCoverage;
    uint32_t maxCoverage = std::numeric_limits<TMaxCoverage>::max();
    typedef std::vector<uint32_t> TCovDist;
    TCovDist covDist(maxCoverage, 0);

    // Error rates
    uint64_t matchCount = 0;
    uint64_t mismatchCount = 0;
    uint64_t delCount = 0;
    uint64_t insCount = 0;

    // Read length distribution
    typedef uint16_t TMaxReadLength;
    uint32_t maxReadLength = std::numeric_limits<TMaxReadLength>::max();
    uint32_t rlBinSize = 100;
    typedef std::vector<uint32_t> TReadLengthDist;
    TReadLengthDist rlDist(maxReadLength * rlBinSize, 0);
    
    // Track reference reads
    char* seq = NULL;
    faidx_t* fai = fai_load(c.genome.string().c_str());
    
    // Iterate chromosomes
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      ++show_progress;
      if (validRegions[refIndex].empty()) continue;

      // Coverage track
      typedef std::vector<TMaxCoverage> TBpCoverage;
      TBpCoverage covBases(hdr->target_len[refIndex], 0);

      // Load sequence
      std::string tname(hdr->target_name[refIndex]);
      int32_t seqlen = -1;
      seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);

      // N-content
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet nrun(hdr->target_len[refIndex], false);
      for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	if ((seq[i] == 'n') || (seq[i] == 'N')) nrun[i] = 1;
      }

      // Flag breakpoints
      typedef std::set<int32_t> TIdSet;
      typedef std::map<uint32_t, TIdSet> TBpToIdMap;
      TBpToIdMap bpid;
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet bpOccupied(hdr->target_len[refIndex], false);
      for(uint32_t i = 0; i < svs.size(); ++i) {
	if ((svs[i].chr == refIndex) || (svs[i].chr2 == refIndex)) {
	  if (svs[i].chr == refIndex) {
	    svs[i].alleles = _addAlleles(boost::to_upper_copy(std::string(seq + svs[i].svStart - 1, seq + svs[i].svStart)), std::string(hdr->target_name[svs[i].chr2]), svs[i], svs[i].svt);
	    bpOccupied[svs[i].svStart] = 1;
	    if (bpid.find(svs[i].svStart) == bpid.end()) bpid.insert(std::make_pair(svs[i].svStart, TIdSet()));
	    bpid[svs[i].svStart].insert(svs[i].id);
	  }
	  if (svs[i].chr2 == refIndex) {
	    bpOccupied[svs[i].svEnd] = 1;
	    if (bpid.find(svs[i].svEnd) == bpid.end()) bpid.insert(std::make_pair(svs[i].svEnd, TIdSet()));
	    bpid[svs[i].svEnd].insert(svs[i].id);
	  }
	}
      }

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

	    // Get sequence
	    std::string sequence;
	    sequence.resize(rec->core.l_qseq);
	    uint8_t* seqptr = bam_get_seq(rec);
	    for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	    // Read length
	    if (!(rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP))) {
	      uint32_t slen = rec->core.l_qseq;
	      if (!slen) slen = sequenceLength(rec);
	      if (slen < maxReadLength * rlBinSize) ++rlDist[(int32_t) (slen / rlBinSize)];
	    }

	    // Hashing
	    std::size_t seed = hash_lr(rec);
	    uint32_t rp = rec->core.pos; // reference pointer
	    uint32_t sp = 0; // sequence pointer
      
	    // Parse the CIGAR
	    uint32_t* cigar = bam_get_cigar(rec);
	    for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
		// Fetch reference alignments
		int32_t bpHit = 0;
		for(uint32_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k) {
		  if ((rp < hdr->target_len[refIndex]) && (covBases[rp] < maxCoverage - 1)) ++covBases[rp];
		  if (bpOccupied[rp]) bpHit = rp;
		  if (sequence[sp] == seq[rp]) ++matchCount;
		  else ++mismatchCount;
		  ++sp;
		  ++rp;
		}
		if (bpHit) {
		  // Only spanning alignments
		  int32_t bg = std::max((int32_t) 0, bpHit - c.minimumFlankSize);
		  int32_t ed = std::min((int32_t) hdr->target_len[refIndex], bpHit + c.minimumFlankSize);
		  if ((bg < rec->core.pos) || ((uint32_t) ed > rec->core.pos + alignmentLength(rec))) continue;
		  
		  // Re-iterate cigar to track alignment		
		  int32_t rit = rec->core.pos; // reference pointer
		  int32_t sit = 0; // sequence pointer
		  std::string refAlign;
		  std::string altAlign;
		  for (std::size_t ij = 0; ij < rec->core.n_cigar; ++ij) {
		    if (rit >= ed) break;
		    if ((bam_cigar_op(cigar[ij]) == BAM_CMATCH) || (bam_cigar_op(cigar[ij]) == BAM_CEQUAL) || (bam_cigar_op(cigar[ij]) == BAM_CDIFF)) {
		      for(uint32_t k = 0; k < bam_cigar_oplen(cigar[ij]); ++k) {
			if ((rit >= bg) && (rit < ed)) {
			  refAlign.push_back(seq[rit]);
			  altAlign.push_back(sequence[sit]);
			}
			++rit;
			++sit;
		      }
		    } else if ((bam_cigar_op(cigar[ij]) == BAM_CDEL) || (bam_cigar_op(cigar[ij]) == BAM_CREF_SKIP)) {
		      for(uint32_t k = 0; k < bam_cigar_oplen(cigar[ij]); ++k) {
			if ((rit >= bg) && (rit < ed)) {
			  refAlign.push_back(seq[rit]);
			  altAlign.push_back('-');
			}
			++rit;
		      }
		    } else if (bam_cigar_op(cigar[ij]) == BAM_CINS) {
		      for(uint32_t k = 0; k < bam_cigar_oplen(cigar[ij]); ++k) {
			if ((rit >= bg) && (rit < ed)) {
			  refAlign.push_back('-');
			  altAlign.push_back(sequence[sit]);
			}
			++sit;
		      }
		    } else if (bam_cigar_op(cigar[ij]) == BAM_CSOFT_CLIP) {
		      sit += bam_cigar_oplen(cigar[ij]);
		    } else if (bam_cigar_op(cigar[ij]) == BAM_CHARD_CLIP) {
		      // Do nothing
		    }
		  }
		  int32_t score = scoreAlign(boost::to_upper_copy<std::string>(refAlign), boost::to_upper_copy<std::string>(altAlign), c.aliscore);
		  int32_t scth = scoreThreshold(refAlign, altAlign, c);
		  if ((score > 0) && (scth > 0)) {
		    int32_t qual = (int32_t) ((((float) score / (float) (scth)) - 1.0) * 100.0);
		    if ((qual > 0) && (qual > c.minMapQual)) {
		      // Debug
		      //std::cerr << "Qual: " << qual << std::endl;
		      //std::cerr << refAlign << std::endl;
		      //std::cerr << altAlign << std::endl;
		      //std::cerr << std::endl;
		  
		      // Assign read to all SVs at this position
		      for(typename TIdSet::const_iterator it = bpid[bpHit].begin(); it != bpid[bpHit].end(); ++it) {
			if (svSup[*it].find(seed) == svSup[*it].end()) svSup[*it].insert(std::make_pair(seed, qual));
			else {
			  if (svSup[*it][seed] < (uint32_t) qual) svSup[*it][seed] = (uint32_t) qual;
			}
		      }
		    }
		  }
		}
	      } else if ((bam_cigar_op(cigar[i]) == BAM_CDEL) || (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP)) {
		++delCount;
		rp += bam_cigar_oplen(cigar[i]);
	      } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
		++insCount;
		sp += bam_cigar_oplen(cigar[i]);
	      } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
		sp += bam_cigar_oplen(cigar[i]);
	      } else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
		// Do nothing
	      } else {
		std::cerr << "Unknown Cigar options" << std::endl;
	      }
	    }
	  }
	  bam_destroy1(rec);
	  hts_itr_destroy(iter);
	}
      }
      if (seq != NULL) free(seq);
      
      // Summarize coverage for this chromosome
      for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	if (!nrun[i]) ++covDist[covBases[i]];
      }
      
      // Assign SV support
      for(uint32_t i = 0; i < svs.size(); ++i) {
	if (svs[i].chr == refIndex) {
	  int32_t halfSize = (svs[i].svEnd - svs[i].svStart)/2;
	  if ((_translocation(svs[i].svt)) || (svs[i].svt == 4)) halfSize = 500;

	  // Left region
	  int32_t lstart = std::max(svs[i].svStart - halfSize, 0);
	  int32_t lend = svs[i].svStart;
	  int32_t covbase = 0;
	  for(uint32_t k = lstart; ((k < (uint32_t) lend) && (k < hdr->target_len[refIndex])); ++k) covbase += covBases[k];
	  covMap[0][svs[i].id].leftRC = covbase;

	  // Actual SV
	  covbase = 0;
	  int32_t mstart = svs[i].svStart;
	  int32_t mend = svs[i].svEnd;
	  if ((_translocation(svs[i].svt)) || (svs[i].svt == 4)) {
	    mstart = std::max(svs[i].svStart - halfSize, 0);
	    mend = std::min(svs[i].svStart + halfSize, (int32_t) hdr->target_len[refIndex]);
	  }
	  for(uint32_t k = mstart; ((k < (uint32_t) mend) && (k < hdr->target_len[refIndex])); ++k) covbase += covBases[k];
	  covMap[0][svs[i].id].rc = covbase;

	  // Right region
	  covbase = 0;
	  int32_t rstart = svs[i].svEnd;
	  int32_t rend = std::min(svs[i].svEnd + halfSize, (int32_t) hdr->target_len[refIndex]);
	  if ((_translocation(svs[i].svt)) || (svs[i].svt == 4)) {
	    rstart = svs[i].svStart;
	    rend = std::min(svs[i].svStart + halfSize, (int32_t) hdr->target_len[refIndex]);
	  }
	  for(uint32_t k = rstart; ((k < (uint32_t) rend) && (k < hdr->target_len[refIndex])); ++k) covbase += covBases[k];
	  covMap[0][svs[i].id].rightRC = covbase;
	}
      }
    }
    fai_destroy(fai);

    // Assign SV support
    for(uint32_t i = 0; i < svs.size(); ++i) {
      for (typename TQualMap::const_iterator it = svSup[svs[i].id].begin(); it != svSup[svs[i].id].end(); ++it) {
	jctMap[0][svs[i].id].ref.push_back(it->second);
      }
    }

    // Output coverage info
    uint64_t totalCovCount = 0;
    for (uint32_t i = 0; i < covDist.size(); ++i) totalCovCount += covDist[i];
    std::vector<uint32_t> covPercentiles(5, 0);  // 5%, 25%, 50%, 75%, 95%
    uint64_t cumCovCount = 0;
    for (uint32_t i = 0; i < covDist.size(); ++i) {
      cumCovCount += covDist[i];
      double frac = (double) cumCovCount / (double) totalCovCount;
      if (frac < 0.05) covPercentiles[0] = i + 1;
      if (frac < 0.25) covPercentiles[1] = i + 1;
      if (frac < 0.5) covPercentiles[2] = i + 1;
      if (frac < 0.75) covPercentiles[3] = i + 1;
      if (frac < 0.95) covPercentiles[4] = i + 1;
    }
    std::cout << "Coverage distribution (^COV)" << std::endl;
    std::cout << "COV\t95% of bases are >=" << covPercentiles[0] << "x" << std::endl;
    std::cout << "COV\t75% of bases are >=" << covPercentiles[1] << "x" << std::endl;
    std::cout << "COV\t50% of bases are >=" << covPercentiles[2] << "x" << std::endl;
    std::cout << "COV\t25% of bases are >=" << covPercentiles[3] << "x" << std::endl;
    std::cout << "COV\t5% of bases are >=" << covPercentiles[4] << "x" << std::endl;

    
    // Output read length info
    uint64_t totalRlCount = 0;
    for (uint32_t i = 0; i < rlDist.size(); ++i) totalRlCount += rlDist[i];
    std::vector<uint32_t> rlPercentiles(5, 0);  // 5%, 25%, 50%, 75%, 95%
    uint64_t cumRlCount = 0;
    for (uint32_t i = 0; i < rlDist.size(); ++i) {
      cumRlCount += rlDist[i];
      double frac = (double) cumRlCount / (double) totalRlCount;
      if (frac < 0.05) rlPercentiles[0] = (i + 1) * rlBinSize;
      if (frac < 0.25) rlPercentiles[1] = (i + 1) * rlBinSize;
      if (frac < 0.5) rlPercentiles[2] = (i + 1) * rlBinSize;
      if (frac < 0.75) rlPercentiles[3] = (i + 1) * rlBinSize;
      if (frac < 0.95) rlPercentiles[4] = (i + 1) * rlBinSize;
    }
    std::cout << "Read-length distribution (^RL)" << std::endl;
    std::cout << "RL\t95% of reads are >=" << rlPercentiles[0] << "bp" << std::endl;
    std::cout << "RL\t75% of reads are >=" << rlPercentiles[1] << "bp" << std::endl;
    std::cout << "RL\t50% of reads are >=" << rlPercentiles[2] << "bp" << std::endl;
    std::cout << "RL\t25% of reads are >=" << rlPercentiles[3] << "bp" << std::endl;
    std::cout << "RL\t5% of reads are >=" << rlPercentiles[4] << "bp" << std::endl;

    // Output sequencing error rates
    uint64_t alignedbases = matchCount + mismatchCount + delCount + insCount;
    std::cout << "Sequencing error rates (^ERR)" << std::endl;
    std::cout << "ERR\tMatchRate\t" << (double) matchCount / (double) alignedbases << std::endl;
    std::cout << "ERR\tMismatchRate\t" << (double) mismatchCount / (double) alignedbases << std::endl;
    std::cout << "ERR\tDeletionRate\t" << (double) delCount / (double) alignedbases << std::endl;
    std::cout << "ERR\tInsertionRate\t" << (double) insCount / (double) alignedbases << std::endl;

    // Clean-up
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }
     
  
  template<typename TConfig, typename TValidRegion, typename TReadBp>
  inline void
  findJunctions(TConfig const& c, TValidRegion const& validRegions, TReadBp& readBp) {
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
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Split-read scanning" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );

    // Iterate chromosomes
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      ++show_progress;
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
	    uint32_t rp = rec->core.pos; // reference pointer
	    uint32_t sp = 0; // sequence pointer
	    
	    // Parse the CIGAR
	    uint32_t* cigar = bam_get_cigar(rec);
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
	      } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
		int32_t finalsp = sp;
		bool scleft = false;
		if (sp == 0) {
		  finalsp += bam_cigar_oplen(cigar[i]); // Leading soft-clip / hard-clip
		  scleft = true;
		}
		sp += bam_cigar_oplen(cigar[i]);
		if (bam_cigar_oplen(cigar[i]) > c.minClip) _insertJunction(readBp, seed, rec, rp, finalsp, scleft);
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
    for(typename TReadBp::iterator it = readBp.begin(); it != readBp.end(); ++it) {
      std::sort(it->second.begin(), it->second.end(), SortJunction<Junction>());
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
    // Extract BAM records
    if ((!c.svtcmd) || (c.svtset.find(2) != c.svtset.end())) selectDeletions(c, readBp, br);
    if ((!c.svtcmd) || (c.svtset.find(3) != c.svtset.end())) selectDuplications(c, readBp, br);
    if ((!c.svtcmd) || (c.svtset.find(0) != c.svtset.end()) || (c.svtset.find(1) != c.svtset.end())) selectInversions(c, readBp, br);
    if ((!c.svtcmd) || (c.svtset.find(4) != c.svtset.end())) selectInsertions(c, readBp, br);
    if ((!c.svtcmd) || (c.svtset.find(5) != c.svtset.end()) || (c.svtset.find(6) != c.svtset.end()) || (c.svtset.find(7) != c.svtset.end()) || (c.svtset.find(8) != c.svtset.end())) selectTranslocations(c, readBp, br);
  }


  template<typename TConfig>
  inline void
  outputSRBamRecords(TConfig const& c, std::vector<std::vector<SRBamRecord> > const& br) {
    samFile* samfile = sam_open(c.files[0].string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Header
    std::cerr << "chr1\tpos1\tchr2\tpos2\tsvtype\tct\tinslen" << std::endl;

    // SVs
    for(uint32_t svt = 0; svt < br.size(); ++svt) {
      for(uint32_t i = 0; i < br[svt].size(); ++i) {
	std::cerr << hdr->target_name[br[svt][i].chr] << '\t' << br[svt][i].pos << '\t' << hdr->target_name[br[svt][i].chr2] << '\t' << br[svt][i].pos2 << '\t' << _addID(svt) << '\t' << _addOrientation(svt) << '\t' << br[svt][i].inslen << std::endl;
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
