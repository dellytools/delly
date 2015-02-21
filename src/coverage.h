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
#include "tags.h"
#include "api/BamReader.h"
#include "api/BamIndex.h"
#include "api/BamAlignment.h"

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

template<typename TCigarVec>
inline unsigned int alignmentLength(TCigarVec const& cigarOperations) {
  typename TCigarVec::const_iterator coIter = cigarOperations.begin();
  if (coIter == cigarOperations.end()) return 0;
  unsigned int alen = 0;
  unsigned int offset = 0;
  if (coIter->Type == 'S') offset+=coIter->Length;
  for(; coIter != cigarOperations.end(); ++coIter) {
    if (coIter->Type == 'M') alen+=coIter->Length;
  }
  return (offset + alen/2);
}



template<typename TPos, typename TCigar, typename TWindow, typename TCount>
inline void
_addBpCounts(TPos pos, TCigar cigar, TWindow posBeg, TWindow posEnd, TCount& bp_sum, BpLevelType<BpLevelCount>)
{  
  if (pos >= posEnd) return;
  int32_t bpPos = pos;
  int num_start = 0;
  int num_end = 0;
  std::string::const_iterator cig = cigar.begin();
  std::string::const_iterator cigEnd = cigar.end();
  for(;cig!=cigEnd; ++cig, ++num_end) {
    if (((int) *cig >=48) && ((int) *cig <= 57)) continue;
    unsigned int len = atoi(cigar.substr(num_start, (num_end-num_start)).c_str());
    if (*cig == 'M') 
      for(unsigned int i = 0; i<len; ++i, ++bpPos) {
	if ((bpPos>=posBeg) && (bpPos<posEnd)) ++bp_sum;
      }
    else if ((*cig == 'N') || (*cig=='D')) bpPos += len;
    num_start = num_end + 1;
  }
}

template<typename TPos, typename TCigar, typename TWindow, typename TCount>
inline void
_addBpCounts(TPos, TCigar, TWindow, TWindow, TCount, BpLevelType<NoBpLevelCount>)
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
  typedef typename TCountMap::key_type TSampleSVPair;

   // Get the references
  BamTools::BamReader readerRef;
  if ( ! readerRef.Open(files[0].string()) ) return;
  BamTools::RefVector references = readerRef.GetReferenceData();

  // Sort Structural Variants
  sort(svs.begin(), svs.end(), SortSVs<TSV>());

  // Initialize count maps
  int32_t maxReadLen = 1000;
  for(typename TSampleLibrary::iterator sIt = sampleLib.begin(); sIt!=sampleLib.end();++sIt) {
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

    // Initialize bam file
    BamTools::BamReader reader;
    reader.Open(files[file_c].string());
    reader.LocateIndex();

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
	  BamTools::BamAlignment al;
	  if (reader.SetRegion(oldChr, std::max(0, (int32_t) itInt->first - maxReadLen), oldChr, itInt->second)) {
	    while( reader.GetNextAlignmentCore(al) ) {
	      if ((al.AlignmentFlag & 0x0004) || (al.AlignmentFlag & 0x0200) || (al.AlignmentFlag & 0x0400) || (al.MapQuality < minMapQual) || (al.AlignmentFlag & 0x0100) || (al.AlignmentFlag & 0x0800)) continue;
	  
	      // Is it a unique pair
	      if (al.Position!=oldPos) {
		oldPos=al.Position;
		unique_pairs_read.clear();
	      }
	      if (_redundancyFilter(al.MatePosition, unique_pairs_read, covType)) {
		int32_t midPoint = al.Position + alignmentLength(al.CigarData);
		if ((midPoint >= itInt->first) && (midPoint < itInt->second)) ++read_sum;
		_addBpCounts(al.Position, cigarString(al.CigarData), itInt->first, itInt->second, bp_sum, bpLevel);
	      }
	    }
	  }
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
	TSampleSVPair svSample = std::make_pair(sampleName, itSV->id);
	typename TCountMap::iterator countMapIt=countMap.find(svSample);
	//std::cerr << itSV->id << ':' << cumBpSum << ',' << cumReadSum << std::endl;
	countMapIt->second.first=cumBpSum;
	countMapIt->second.second=cumReadSum;
      }
    }
  }
}

}

#endif
