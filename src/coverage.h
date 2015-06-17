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

namespace torali {


template<typename TCigarVec>
inline std::string cigarString(TCigarVec const& cigarOperations) {
  std::stringstream cigar;
  typename TCigarVec::const_iterator coIter = cigarOperations.begin();
  if (coIter == cigarOperations.end()) cigar << "*";
  else
    for(; coIter != cigarOperations.end(); ++coIter) cigar << coIter->Length << coIter->Type;
  return cigar.str();
}

inline unsigned int halfAlignmentLength(bam1_t* rec) {
  uint32_t* cigar = bam_get_cigar(rec);
  unsigned int alen = 0;
  for (unsigned int i = 0; i < rec->core.n_cigar; ++i) 
    if (bam_cigar_op(cigar[i]) == BAM_CMATCH) alen+=bam_cigar_oplen(cigar[i]);
  return (alen/2);
}


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

template<typename TPos, typename TUniquePairs>
inline bool 
_redundancyFilter(TPos matePos, TUniquePairs& uRead, CoverageType<RedundancyFilterTag>) {
  return uRead.insert(matePos).second;
}

template<typename TPos, typename TUniquePairs>
inline bool
_redundancyFilter(TPos, TUniquePairs&, CoverageType<NoRedundancyFilterTag>) {
  return true;
}



template<typename TFiles, typename TSampleLibrary, typename TSVs, typename TCountMap, typename TBpLevelType, typename TCoverageType>
inline void
annotateCoverage(TFiles const& files, uint16_t minMapQual, TSampleLibrary& sampleLib, TSVs& svs, TCountMap& countMap, TBpLevelType bpLevel, TCoverageType covType)
{
  typedef typename TSVs::value_type TSV;

  // Open file handles
  typedef std::vector<samFile*> TSamFile;
  typedef std::vector<hts_idx_t*> TIndex;
  TSamFile samfile;
  TIndex idx;
  samfile.resize(files.size());
  idx.resize(files.size());
  for(unsigned int file_c = 0; file_c < files.size(); ++file_c) {
    samfile[file_c] = sam_open(files[file_c].string().c_str(), "r");
    idx[file_c] = sam_index_load(samfile[file_c], files[file_c].string().c_str());
  }

  // For alignment midpoint, maximum read-length
  int32_t maxReadLen = 1000;

  // Sort Structural Variants
  sort(svs.begin(), svs.end(), SortSVs<TSV>());

  // Initialize count maps
  for(typename TSampleLibrary::iterator sIt = sampleLib.begin(); sIt!=sampleLib.end(); ++sIt) {
    for(typename TSVs::const_iterator itSV = svs.begin(); itSV!=svs.end(); ++itSV) {
      countMap.insert(std::make_pair(std::make_pair(sIt->first, itSV->id), std::make_pair(0,0)));
    }
  }

  // Iterate all samples
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Read-depth annotation" << std::endl;
  boost::progress_display show_progress( (svs.end() - svs.begin()) );
#pragma omp parallel for default(shared)
  for(unsigned int file_c = 0; file_c < files.size(); ++file_c) {
    // Get a sample name
    std::string sampleName(files[file_c].stem().string());

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
	typedef std::set<int32_t> TBreaks;
	TBreaks iBounds;
	iBounds.insert(itSV->svStart);
	iBounds.insert(itSV->svEnd);
	int32_t maxPos = itSV->svEnd;
	typename TSVs::const_iterator itSVFor = itSV;
	for(++itSVFor;itSVFor!=itSVEnd; ++itSVFor) {
	  if (oldChr!=itSVFor->chr) break;
	  if (itSVFor->svStart>maxPos) {
	    typename TBreaks::const_iterator itBreak = iBounds.begin();
	    typename TBreaks::const_iterator itBreakNext = iBounds.begin();
	    for(++itBreakNext; itBreakNext!=iBounds.end(); ++itBreakNext, ++itBreak) intervals.push_back(std::make_pair(*itBreak, *itBreakNext));
	    maxPos=itSVFor->svEnd;
	    iBounds.clear();
	  }
	  iBounds.insert(itSVFor->svStart);
	  iBounds.insert(itSVFor->svEnd);
	  if (itSVFor->svEnd > maxPos) maxPos=itSVFor->svEnd;
	}
	typename TBreaks::const_iterator itBreak = iBounds.begin();
	typename TBreaks::const_iterator itBreakNext = iBounds.begin();
	for(++itBreakNext; itBreakNext!=iBounds.end(); ++itBreakNext, ++itBreak) intervals.push_back(std::make_pair(*itBreak, *itBreakNext));

	// Process sub-intervals
	typename TInterval::const_iterator itInt = intervals.begin();
	for(;itInt!=intervals.end(); ++itInt) {
	  // Unique pairs for the given interval
	  typedef boost::container::flat_set<int32_t> TUniquePairs;
	  TUniquePairs unique_pairs_read;

	  // Count reads / aligned base-pairs
	  unsigned int bp_sum = 0;
	  unsigned int read_sum = 0;
	  int32_t oldPos=-1;

	  // Read alignments
	  hts_itr_t* iter = sam_itr_queryi(idx[file_c], oldChr, std::max(0, (int32_t) itInt->first - maxReadLen), itInt->second);
	  bam1_t* rec = bam_init1();
	  while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	    if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	    if (rec->core.qual < minMapQual) continue;
	  
	    // Is it a unique pair
	    if (rec->core.pos!=oldPos) {
	      oldPos=rec->core.pos;
	      unique_pairs_read.clear();
	    }
	    if (_redundancyFilter(rec->core.mpos, unique_pairs_read, covType)) {
	      int32_t midPoint = rec->core.pos + halfAlignmentLength(rec);
	      if ((midPoint >= itInt->first) && (midPoint < itInt->second)) ++read_sum;
	      _addBpCounts(rec, itInt->first, itInt->second, bp_sum, bpLevel);
	    }
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
      typename TInterval::const_iterator itInt = intervals.begin();
      typename TIntervalSum::const_iterator itRead = readSum.begin();
      typename TIntervalSum::const_iterator itBp = bpSum.begin();
      for(;itInt!=intervals.end(); ++itInt, ++itRead, ++itBp) {
	//std::cerr << itInt->first << ',' << itInt->second << ':' << *itRead << ',' << *itBp << '(' << itSV->svStart << ',' << itSV->svEnd << ')' << cumReadSum << std::endl;
	if (itInt->first > itSV->svEnd) break;
	if ((itInt->first >= itSV->svStart) && (itInt->second <= itSV->svEnd)) {
	  cumBpSum+=(*itBp);
	  cumReadSum+=(*itRead);
	}
      }

      // Store counts
#pragma omp critical
      {
	typedef typename TCountMap::key_type TSampleSVPair;
	TSampleSVPair svSample = std::make_pair(sampleName, itSV->id);
	typename TCountMap::iterator countMapIt=countMap.find(svSample);
	//std::cerr << itSV->id << ':' << cumBpSum << ',' << cumReadSum << std::endl;
	countMapIt->second.first=cumBpSum;
	countMapIt->second.second=cumReadSum;
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
