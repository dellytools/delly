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

#ifndef COVERAGE_H
#define COVERAGE_H

#include <boost/container/flat_set.hpp>
#include <htslib/sam.h>
#include "tags.h"
#include "util.h"

namespace torali {


template<typename TWindow, typename TCount>
inline void
_addBpCounts(bam1_t* rec, TWindow posBeg, TWindow posEnd, TCount& bp_sum, BpLevelType<BpLevelCount>)
{  
  if (rec->core.pos >= posEnd) return;
  int32_t bpPos = rec->core.pos;
  uint32_t* cigar = bam_get_cigar(rec);
  for (unsigned int i = 0; i < rec->core.n_cigar; ++i) {
    int op = bam_cigar_op(cigar[i]);
    int ol = bam_cigar_oplen(cigar[i]);
    if (op == BAM_CMATCH) 
      for(int k = 0; k<ol; ++k, ++bpPos) {
	if ((bpPos>=posBeg) && (bpPos<posEnd)) ++bp_sum;
      }
    else if ((op == BAM_CREF_SKIP) || (op == BAM_CDEL)) bpPos += ol;
  }
}

template<typename TWindow, typename TCount>
inline void
_addBpCounts(bam1_t*, TWindow, TWindow, TCount, BpLevelType<NoBpLevelCount>)
{  
  //Nop
}

template<typename TFiles, typename TSVs, typename TCountMap, typename TBpLevelType>
inline void
annotateCoverage(TFiles const& files, uint16_t minMapQual, TSVs& svs, TCountMap& countMap, TBpLevelType bpLevel)
{
  typedef typename TSVs::value_type TSV;

  // Open file handles
  typedef std::vector<samFile*> TSamFile;
  typedef std::vector<hts_idx_t*> TIndex;
  TSamFile samfile(files.size());
  TIndex idx(files.size());
  for(unsigned int file_c = 0; file_c < files.size(); ++file_c) {
    samfile[file_c] = sam_open(files[file_c].string().c_str(), "r");
    idx[file_c] = sam_index_load(samfile[file_c], files[file_c].string().c_str());
  }

  // Sort Structural Variants
  sort(svs.begin(), svs.end(), SortSVs<TSV>());

  // Initialize count maps
  countMap.resize(files.size());
  for(uint32_t file_c = 0; file_c < files.size(); ++file_c) countMap[file_c].resize(svs.size());

  // Iterate all samples
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Read-depth annotation" << std::endl;
  boost::progress_display show_progress( (svs.end() - svs.begin()) );
#pragma omp parallel for default(shared)
  for(unsigned int file_c = 0; file_c < files.size(); ++file_c) {
    // Read alignments
    int32_t oldChr = -1;
    typedef std::vector< std::pair<int32_t, int32_t> > TInterval;
    TInterval intervals;
    typedef std::vector<unsigned int> TIntervalSum;
    TIntervalSum readSum;
    TIntervalSum bpSum;

    typename TSVs::const_iterator itSV = svs.begin();
    typename TSVs::const_iterator itSVEnd = svs.end();
    for(;itSV!=itSVEnd;++itSV) {
      if (file_c==(files.size()-1)) ++show_progress;
      
      // Find all non-overlapping sub-intervals
      if (itSV->chr!=oldChr) {
	oldChr=itSV->chr;
	intervals.clear();
	bpSum.clear();
	readSum.clear();
	typedef std::vector<int32_t> TBreaks;
	TBreaks iBounds;
	iBounds.push_back(itSV->svStart);
	iBounds.push_back(itSV->svEnd);
	int32_t maxPos = itSV->svEnd;
	typename TSVs::const_iterator itSVFor = itSV;
	for(++itSVFor;itSVFor!=itSVEnd; ++itSVFor) {
	  if (oldChr!=itSVFor->chr) break;
	  if (itSVFor->svStart>maxPos) {
	    std::sort(iBounds.begin(), iBounds.end());
	    typename TBreaks::const_iterator itBreak = iBounds.begin();
	    typename TBreaks::const_iterator itBreakNext = iBounds.begin();
	    for(++itBreakNext; itBreakNext!=iBounds.end(); ++itBreakNext, ++itBreak)
	      if (*itBreak < *itBreakNext) intervals.push_back(std::make_pair(*itBreak, *itBreakNext));
	    maxPos=itSVFor->svEnd;
	    iBounds.clear();
	  }
	  iBounds.push_back(itSVFor->svStart);
	  iBounds.push_back(itSVFor->svEnd);
	  if (itSVFor->svEnd > maxPos) maxPos=itSVFor->svEnd;
	}
	std::sort(iBounds.begin(), iBounds.end());
	typename TBreaks::const_iterator itBreak = iBounds.begin();
	typename TBreaks::const_iterator itBreakNext = iBounds.begin();
	for(++itBreakNext; itBreakNext!=iBounds.end(); ++itBreakNext, ++itBreak) 
	  if (*itBreak < *itBreakNext) intervals.push_back(std::make_pair(*itBreak, *itBreakNext));
	iBounds.clear();

	// Process sub-intervals
	for(typename TInterval::const_iterator itInt = intervals.begin(); itInt!=intervals.end(); ++itInt) {
	  // Count reads / aligned base-pairs
	  unsigned int bp_sum = 0;
	  unsigned int read_sum = 0;

	  hts_itr_t* iter = sam_itr_queryi(idx[file_c], oldChr, itInt->first, itInt->second);
	  bam1_t* rec = bam_init1();
	  while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	    if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	    if (rec->core.qual < minMapQual) continue;
	  
	    int32_t midPoint = rec->core.pos + halfAlignmentLength(rec);
	    if ((midPoint >= itInt->first) && (midPoint < itInt->second)) ++read_sum;
	    _addBpCounts(rec, itInt->first, itInt->second, bp_sum, bpLevel);
	  }
	  bam_destroy1(rec);
	  hts_itr_destroy(iter);
	  bpSum.push_back(bp_sum);
	  readSum.push_back(read_sum);
	}
      }

      // Sum up sub-intervals
      unsigned int cumBpSum = 0;
      unsigned int cumReadSum = 0;
      typename TInterval::const_iterator itInt = std::lower_bound(intervals.begin(), intervals.end(), std::make_pair(itSV->svStart, 0));
      typename TIntervalSum::const_iterator itRead = readSum.begin();
      itRead += (itInt - intervals.begin());
      typename TIntervalSum::const_iterator itBp = bpSum.begin();
      itBp += (itInt -intervals.begin());
      for(;itInt!=intervals.end(); ++itInt, ++itRead, ++itBp) {
	if (itInt->first >= itSV->svEnd) break;
	cumBpSum+=(*itBp);
	cumReadSum+=(*itRead);
      }

      // Store counts
#pragma omp critical
      {
	countMap[file_c][itSV->id] = std::make_pair(cumBpSum, cumReadSum);
      }
    }
  }
  // Clean-up
  for(unsigned int file_c = 0; file_c < files.size(); ++file_c) {
    hts_idx_destroy(idx[file_c]);
    sam_close(samfile[file_c]);
  }
}

}

#endif
