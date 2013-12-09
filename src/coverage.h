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
inline std::string  cigarString(TCigarVec const& cigarOperations) {
  std::stringstream cigar;
  typename TCigarVec::const_iterator coIter = cigarOperations.begin();
  if (coIter == cigarOperations.end()) cigar << "*";
  else
    for(; coIter != cigarOperations.end(); ++coIter) cigar << coIter->Length << coIter->Type;
  return cigar.str();
}



// Single hit data structure for alignment position of a single read
template<typename TPos, typename TCigar>
struct SingleHit {
  TPos pos;
  TCigar cigar;

  SingleHit() {}

  SingleHit(BamTools::BamAlignment const& al) : pos(al.Position+1), cigar(cigarString(al.CigarData)) {}
};

template<typename TPos>
struct SingleHit<TPos, void> {
  TPos pos;

  SingleHit() {}

  SingleHit(BamTools::BamAlignment const& al) : pos(al.Position+1) {}
};


template<typename TPos, typename TCigar, typename TArrayType>
inline int
_addReadAndBpCounts(std::vector<SingleHit<TPos, TCigar> > const& hit_vector, TArrayType* read_count, TArrayType* bp_count)
{
  typedef std::vector<SingleHit<TPos, TCigar> > THits;
  typename THits::const_iterator vecBeg = hit_vector.begin();
  typename THits::const_iterator vecEnd = hit_vector.end();
  
  // Add read and bp counts
  int lastCoveredPos=0;
  for(;vecBeg!=vecEnd; ++vecBeg) {
    ++read_count[vecBeg->pos];
    int num_start = 0;
    int num_end = 0;
    TArrayType* bpPoint = &bp_count[vecBeg->pos];
    std::string::const_iterator cig = vecBeg->cigar.begin();
    std::string::const_iterator cigEnd = vecBeg->cigar.end();
    for(;cig!=cigEnd; ++cig, ++num_end) {
      if (((int) *cig >=48) && ((int) *cig <= 57)) continue;
      unsigned int len = atoi(vecBeg->cigar.substr(num_start, (num_end-num_start)).c_str());
      if (*cig == 'M') 
	for(unsigned int i = 0; i<len; ++i, ++bpPoint) ++(*bpPoint);
      else if ((*cig == 'N') || (*cig=='D')) bpPoint += len;
      num_start = num_end + 1;
    }
    if ((bpPoint - bp_count) > lastCoveredPos) lastCoveredPos=(bpPoint - bp_count);
  }
  return lastCoveredPos;
}

template<typename TPos, typename TArrayType>
inline int
_addReadAndBpCounts(std::vector<SingleHit<TPos, void> > const& hit_vector, TArrayType* read_count, TArrayType*)
{
  typedef std::vector<SingleHit<TPos, void> > THits;
  typename THits::const_iterator vecBeg = hit_vector.begin();
  typename THits::const_iterator vecEnd = hit_vector.end();
  
  // Add read and bp counts
  int lastCoveredPos=0;
  for(;vecBeg!=vecEnd; ++vecBeg) {
    ++read_count[vecBeg->pos];
    if (vecBeg->pos > lastCoveredPos) lastCoveredPos=vecBeg->pos;
  }
  return lastCoveredPos;
}

template<typename THit>
struct SortSingleHits : public std::binary_function<THit, THit, bool>
{
  inline bool operator()(THit const& hit1, THit const& hit2) {
    return (hit1.pos < hit2.pos);
  }
};



template<typename TFiles, typename TSampleLibrary, typename TSVs, typename TCountMap, typename TSingleHit>
inline void
annotateCoverage(TFiles const& files, uint16_t minMapQual, bool inclCigar, TSampleLibrary& sampleLib, TSVs& svs, TCountMap& countMap, TSingleHit)
{
  typedef typename TCountMap::key_type TSampleSVPair;
  typedef typename TSampleLibrary::mapped_type TLibraryMap;

   // Get the references
  BamTools::BamReader readerRef;
  if ( ! readerRef.Open(files[0].string()) ) return;
  BamTools::RefVector references = readerRef.GetReferenceData();

  // Reset duplicate counters
  typename TSampleLibrary::iterator sampleIt=sampleLib.begin();
  for(;sampleIt!=sampleLib.end();++sampleIt) {
    typename TLibraryMap::iterator libIt=sampleIt->second.begin();
    for(;libIt!=sampleIt->second.end();++libIt) {
      libIt->second.non_unique_pairs=0;
      libIt->second.unique_pairs=0;
      libIt->second.mappedReads=0;
    }
  }

  // Sort Structural Variants
  sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());

  // Declare the chromosome array
  typedef unsigned short TArrayType;
  TArrayType* read_count = new TArrayType[MAX_CHROM_SIZE]();
  TArrayType* bp_count = new TArrayType[MAX_CHROM_SIZE]();
  int lastCoveredPos=MAX_CHROM_SIZE;

  // Process chromosome by chromosome
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Read-depth annotation" << std::endl;
  boost::progress_display show_progress( (references.end() - references.begin()) );
  typename BamTools::RefVector::const_iterator  itRef = references.begin();
  for(int refIndex=0;itRef!=references.end();++itRef, ++refIndex) {
    ++show_progress;

    // Iterate all samples
    for(unsigned int file_c = 0; file_c < files.size(); ++file_c) {

      // Iterate all SAM files and collect all hits
      typedef std::vector<TSingleHit> THits;
      THits hit_vector;  

      // Get a sample name
      std::string sampleName(files[file_c].stem().string());
      typename TSampleLibrary::iterator sampleIt=sampleLib.find(sampleName);

      // Initialize bam file
      BamTools::BamReader reader;
      if ( ! reader.Open(files[file_c].string()) ) return;
      reader.LocateIndex();
      if ( !reader.HasIndex() ) return;
	
      // Unique pairs for the given sample
      typedef std::set<Hit> TUniquePairs;
      TUniquePairs unique_pairs_read1;
      TUniquePairs unique_pairs_read2;

      // Read alignments
      BamTools::BamAlignment al;
      if (reader.Jump(refIndex, 0)) {
	while( reader.GetNextAlignmentCore(al) ) {
	  if (al.RefID != refIndex) break;
	  if ((al.AlignmentFlag & 0x0004) || (al.AlignmentFlag & 0x0200) || (al.AlignmentFlag & 0x0400) || (al.MapQuality < minMapQual) || (!(inclCigar) && (al.AlignmentFlag & 0x0100)) ) continue;
	  // Get the library information
	  al.BuildCharData();
	  std::string rG = "DefaultLib";
	  al.GetTag("RG", rG);
	  typename TLibraryMap::iterator libIt=sampleIt->second.find(rG);

	  // Is it a unique pair
	  Hit hitPos(al);
	  typename TUniquePairs::const_iterator pos;
	  bool inserted = true;
	  if (al.AlignmentFlag & 0x0040) boost::tie(pos, inserted) = unique_pairs_read1.insert(hitPos);
	  else boost::tie(pos, inserted) = unique_pairs_read2.insert(hitPos);
	  if ((inserted) || !(al.AlignmentFlag & 0x0001)) {
	    hit_vector.push_back(TSingleHit(al));
	    ++libIt->second.mappedReads;
	    if (al.AlignmentFlag & 0x0040) ++libIt->second.unique_pairs;
	  } else {
	    if (al.AlignmentFlag & 0x0040) ++libIt->second.non_unique_pairs;
	  }
	}
      }
    
      // Sort SAM records by chromosome and position
      sort(hit_vector.begin(), hit_vector.end(), SortSingleHits<TSingleHit>());

      // Reset arrays
      memset(read_count, 0, std::min(lastCoveredPos + 1, MAX_CHROM_SIZE) * sizeof(TArrayType));
      memset(bp_count, 0, std::min(lastCoveredPos + 1, MAX_CHROM_SIZE) * sizeof(TArrayType));

      // Iterate all reads of that chromosome
      lastCoveredPos=_addReadAndBpCounts(hit_vector, read_count, bp_count);

      // Store coverage for all input intervals
      StructuralVariantRecord findSV;
      findSV.svStart=0;
      findSV.chr=references[refIndex].RefName;
      typename TSVs::const_iterator itSV = std::lower_bound(svs.begin(), svs.end(), findSV, SortSVs<StructuralVariantRecord>());
      typename TSVs::const_iterator itSVEnd = svs.end();
      for(;itSV!=itSVEnd;++itSV) {
	if (itSV->chr == references[refIndex].RefName) {
	  TSampleSVPair svSample = std::make_pair(sampleName, itSV->id);
	  typename TCountMap::iterator countMapIt=countMap.find(svSample);
	  if (countMapIt==countMap.end()) {
	    countMapIt=countMap.insert(std::make_pair(svSample, std::make_pair(0,0))).first;
	  }

	  int32_t pos = itSV->svStart;
	  int32_t window_len = itSV->svEnd;
	  unsigned int bp_sum = 0;
	  unsigned int read_sum = 0;
	  TArrayType* bpPoint = &bp_count[pos+1];
	  TArrayType* readPoint = &read_count[pos+1];
	  for(int32_t i=pos; i<window_len; ++i, ++bpPoint, ++readPoint) { 
	    bp_sum += *bpPoint;
	    read_sum += *readPoint;
	  }
	  countMapIt->second.first=bp_sum;
	  countMapIt->second.second=read_sum;
	} else break;
      }
    }
  }
  delete[] read_count;
  delete[] bp_count;

  // Get normalization factors
  typedef std::vector<unsigned int> TReadCounts;
  typedef std::map<std::string, double> TSampleNorm;
  TReadCounts readCounts;
  TSampleNorm sampleNorm;
  sampleIt=sampleLib.begin();
  for(;sampleIt!=sampleLib.end();++sampleIt) {
    unsigned int mappedReads = 0;
    typename TLibraryMap::const_iterator libIt=sampleIt->second.begin();
    for(;libIt!=sampleIt->second.end();++libIt) {
      mappedReads+=libIt->second.mappedReads;
    }
    sampleNorm.insert(std::make_pair(sampleIt->first, mappedReads));
    readCounts.push_back(mappedReads);
  }
  std::sort(readCounts.begin(), readCounts.end());
  unsigned int median=readCounts[readCounts.size()/2];
  sampleIt=sampleLib.begin();
  for(;sampleIt!=sampleLib.end();++sampleIt) {
    typename TSampleNorm::iterator sampleNormIt = sampleNorm.find(sampleIt->first);
    if (sampleNormIt->second!=0) sampleNormIt->second = median/sampleNormIt->second;
    else sampleNormIt->second = 1.0;
  }

  // Normalize read counts
  for(unsigned int file_c = 0; file_c < files.size(); ++file_c) {
    std::string sampleName(files[file_c].stem().string());
    typename TSVs::const_iterator itSV = svs.begin();
    typename TSVs::const_iterator itSVEnd = svs.end();
    for(;itSV!=itSVEnd;++itSV) {
      TSampleSVPair svSample = std::make_pair(sampleName, itSV->id);
      typename TCountMap::iterator countMapIt=countMap.find(svSample);
      typename TSampleNorm::iterator sampleNormIt=sampleNorm.find(sampleName);
      countMapIt->second.first *= sampleNormIt->second;
      countMapIt->second.second *= sampleNormIt->second;
    }
  }

}


}

#endif


