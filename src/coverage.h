/*
============================================================================
DELLY: Structural variant discovery by integrated PE mapping and SR analysis
============================================================================
Copyright (C) 2012-2018 Tobias Rausch

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

#ifndef COVERAGE_H
#define COVERAGE_H

#include <boost/container/flat_set.hpp>
#include <htslib/sam.h>
#include "tags.h"
#include "util.h"

namespace torali {

  
  // Reduced structural variant record for cov
  struct CovRecord {
    int32_t svStart;
    int32_t svEnd;
    uint32_t id;

    CovRecord() : svStart(0), svEnd(0), id(0) {}
  };

  template<typename TRecord>
  struct SortCovRecord : public std::binary_function<TRecord, TRecord, bool> {
    inline bool operator()(TRecord const& s1, TRecord const& s2) const {
      return ((s1.svStart<s2.svStart) || ((s1.svStart==s2.svStart) && (s1.svEnd<s2.svEnd)) || ((s1.svStart==s2.svStart) && (s1.svEnd==s2.svEnd) && (s1.id < s2.id)));
    }
  };
  

template<typename TFiles, typename TSVs, typename TCountMap>
inline void
annotateCoverage(TFiles const& files, uint16_t minMapQual, TSVs& svs, TCountMap& countMap)
{
  // Open file handles
  typedef std::vector<samFile*> TSamFile;
  typedef std::vector<hts_idx_t*> TIndex;
  typedef std::vector<bam_hdr_t*> THeader;
  TSamFile samfile(files.size());
  TIndex idx(files.size());
  THeader hdr(files.size());
  int32_t totalTarget = 0;
  for(unsigned int file_c = 0; file_c < files.size(); ++file_c) {
    samfile[file_c] = sam_open(files[file_c].string().c_str(), "r");
    idx[file_c] = sam_index_load(samfile[file_c], files[file_c].string().c_str());
    hdr[file_c] = sam_hdr_read(samfile[file_c]);
    totalTarget += hdr[file_c]->n_targets;
  }

  // Sort Structural Variants
  uint32_t totalSVs = 0;
  for(int32_t refIndex=0; refIndex < hdr[0]->n_targets; ++refIndex) {
    sort(svs[refIndex].begin(), svs[refIndex].end(), SortCovRecord<CovRecord>());
    totalSVs += svs[refIndex].size();
  }

  // Initialize count maps
  countMap.resize(files.size());
  for(uint32_t file_c = 0; file_c < files.size(); ++file_c) countMap[file_c].resize(totalSVs);

  // Iterate all samples
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Read-depth annotation" << std::endl;
  boost::progress_display show_progress( totalTarget );

#pragma omp parallel for default(shared)
  for(unsigned int file_c = 0; file_c < files.size(); ++file_c) {
    // Pair qualities and features
    typedef boost::unordered_map<std::size_t, uint8_t> TQualities;
    TQualities qualities;
  
    // Iterate chromosomes
    for(int32_t refIndex=0; refIndex < (int32_t) hdr[file_c]->n_targets; ++refIndex) {
      ++show_progress;

      // Any SVs on this chromosome?
      if (svs[refIndex].empty()) continue;

      // Check we have mapped reads on this chromosome
      bool nodata = true;
      std::string suffix("cram");
      std::string str(files[file_c].string());
      if ((str.size() >= suffix.size()) && (str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0)) nodata = false;
      uint64_t mapped = 0;
      uint64_t unmapped = 0;
      hts_idx_get_stat(idx[file_c], refIndex, &mapped, &unmapped);
      if (mapped) nodata = false;
      if (nodata) continue;
      
      // Coverage track
      typedef uint16_t TCount;
      uint32_t maxCoverage = std::numeric_limits<TCount>::max();
      typedef std::vector<TCount> TCoverage;
      TCoverage covFragment(hdr[file_c]->target_len[refIndex], 0);
      TCoverage covBases(hdr[file_c]->target_len[refIndex], 0);

      // Count reads
      hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr[file_c]->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      int32_t lastAlignedPos = 0;
      std::set<std::size_t> lastAlignedPosReads;
      while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	if ((rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP | BAM_FMUNMAP)) || (rec->core.tid != rec->core.mtid) || (!(rec->core.flag & BAM_FPAIRED))) continue;
	if (rec->core.qual < minMapQual) continue;

	// Clean-up the read store for identical alignment positions
	if (rec->core.pos > lastAlignedPos) {
	  lastAlignedPosReads.clear();
	  lastAlignedPos = rec->core.pos;
	}
	
	if ((rec->core.pos < rec->core.mpos) || ((rec->core.pos == rec->core.mpos) && (lastAlignedPosReads.find(hash_string(bam_get_qname(rec))) == lastAlignedPosReads.end()))) {
	  // First read
	  lastAlignedPosReads.insert(hash_string(bam_get_qname(rec)));
	  std::size_t hv = hash_pair(rec);
	  qualities[hv] = rec->core.qual;
	} else {
	  // Second read
	  std::size_t hv = hash_pair_mate(rec);
	  if (qualities.find(hv) == qualities.end()) continue; // Mate discarded
	  uint8_t pairQuality = std::min((uint8_t) qualities[hv], (uint8_t) rec->core.qual);
	  qualities[hv] = 0;

	  // Pair quality
	  if (pairQuality < minMapQual) continue; // Low quality pair

	  // Count mid point (fragment counting)
	  int32_t midPoint = rec->core.pos + halfAlignmentLength(rec);
	  if ((midPoint < (int32_t) hdr[file_c]->target_len[refIndex]) && (covFragment[midPoint] < maxCoverage - 1)) ++covFragment[midPoint];

	  // Count basepair (small InDels)
	  uint32_t rp = 0; // reference pointer
	  uint32_t* cigar = bam_get_cigar(rec);
	  for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	    if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
	      for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k) {
		if ((rec->core.pos + rp < hdr[file_c]->target_len[refIndex]) && (covBases[rec->core.pos + rp] < maxCoverage - 1)) ++covBases[rec->core.pos + rp];
		++rp;
	      }
	    } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	      rp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	      rp += bam_cigar_oplen(cigar[i]);
	    }
	  }
	}
      }
      // Clean-up
      bam_destroy1(rec);
      hts_itr_destroy(iter);
      qualities.clear();

      // Assign fragment and base counts to SVs
      for(uint32_t i = 0; i < svs[refIndex].size(); ++i) {
	int32_t covfrag = 0;
	int32_t covbase = 0;
	for(int32_t k = svs[refIndex][i].svStart; ((k < svs[refIndex][i].svEnd) && (k < (int32_t) hdr[file_c]->target_len[refIndex])); ++k) {
	  covfrag += covFragment[k]; 
	  covbase += covBases[k];
	}
	// Store counts
	countMap[file_c][svs[refIndex][i].id].first = covbase;
	countMap[file_c][svs[refIndex][i].id].second = covfrag;
      }
    }
  }
  // Clean-up
  for(unsigned int file_c = 0; file_c < files.size(); ++file_c) {
    bam_hdr_destroy(hdr[file_c]);
    hts_idx_destroy(idx[file_c]);
    sam_close(samfile[file_c]);    
  }
}

}

#endif
