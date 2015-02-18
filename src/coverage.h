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

template<typename TBamRecord, typename TUniquePairs>
inline bool 
_redundancyFilter(TBamRecord const& al, TUniquePairs& uRead1, TUniquePairs& uRead2, CoverageType<RedundancyFilterTag>) {
  Hit hitPos(al);
  bool inserted;
  typename TUniquePairs::const_iterator pos;
  if (al.AlignmentFlag & 0x0040) boost::tie(pos, inserted) = uRead1.insert(hitPos);
  else boost::tie(pos, inserted) = uRead2.insert(hitPos);
  return inserted;
}

template<typename TBamRecord, typename TUniquePairs>
inline bool
_redundancyFilter(TBamRecord const&, TUniquePairs&, TUniquePairs&, CoverageType<NoRedundancyFilterTag>) {
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
    typename TSVs::const_iterator itSV = svs.begin();
    typename TSVs::const_iterator itSVEnd = svs.end();
    for(;itSV!=itSVEnd;++itSV) {
      if (file_c==(files.size()-1)) ++show_progress;

      // Count reads / aligned base-pairs
      unsigned int bp_sum = 0;
      unsigned int read_sum = 0;

      // Unique pairs for the given sample and SV
      typedef std::set<Hit> TUniquePairs;
      TUniquePairs unique_pairs_read1;
      TUniquePairs unique_pairs_read2;

      // Read alignments
      BamTools::BamAlignment al;
      if (reader.SetRegion(itSV->chr, itSV->svStart, itSV->chr, itSV->svEnd)) {
	while( reader.GetNextAlignmentCore(al) ) {
	  if ((al.AlignmentFlag & 0x0004) || (al.AlignmentFlag & 0x0200) || (al.AlignmentFlag & 0x0400) || (al.MapQuality < minMapQual) || (al.AlignmentFlag & 0x0100) || (al.AlignmentFlag & 0x0800)) continue;
	  
	  // Is it a unique pair
	  bool inserted = _redundancyFilter(al, unique_pairs_read1, unique_pairs_read2, covType);
	  if ((inserted) || !(al.AlignmentFlag & 0x0001)) {
	    int32_t midPoint = al.Position + (int32_t) al.Length/2;
	    if ((midPoint >= itSV->svStart) && (midPoint < itSV->svEnd)) ++read_sum;
	    _addBpCounts(al.Position, cigarString(al.CigarData), itSV->svStart, itSV->svEnd, bp_sum, bpLevel);
	  }
	}
      }

      // Store counts
#pragma omp critical
      {
	TSampleSVPair svSample = std::make_pair(sampleName, itSV->id);
	typename TCountMap::iterator countMapIt=countMap.find(svSample);
	countMapIt->second.first=bp_sum;
	countMapIt->second.second=read_sum;
      }
    }
  }
}


}

#endif
