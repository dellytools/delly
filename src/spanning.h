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

#ifndef SPANNING_H
#define SPANNING_H

#include <boost/unordered_map.hpp>
#include "tags.h"

namespace torali {


  template<typename TPos, typename TQual>
    struct HitInterval {
      TPos start;
      TPos end;
      TQual qual;
  
      HitInterval() {}
      
    HitInterval(TPos const s, TPos const e, TQual const q) : start(s+1), end(e+1), qual(q) {}
    };


  template<typename TPos>
    struct HitInterval<TPos, void> {
    TPos start;
    TPos end;
    
    HitInterval() {}
    
  HitInterval(TPos const s, TPos const e) : start(s+1), end(e+1) {}

    template<typename TQual>
      HitInterval(TPos const s, TPos const e, TQual) : start(s+1), end(e+1) {}
  };


  template<typename THitInterval>
    struct SortHitInterval : public std::binary_function<THitInterval, THitInterval, bool>
    {
      inline bool operator()(THitInterval const& hit1, THitInterval const& hit2) {
	return (hit1.start < hit2.start) || ((hit1.start == hit2.start) && (hit1.end < hit2.end));
      }
    };


  template<typename TPos, typename TQual, typename TArrayType>
    inline void
    _addReadAndBpCounts(std::vector<HitInterval<TPos, TQual> > const& hit_vector, TArrayType* bp_count)
  {
    typedef std::vector<HitInterval<TPos, TQual> > THits;
    typename THits::const_iterator vecBeg = hit_vector.begin();
    typename THits::const_iterator vecEnd = hit_vector.end();
  
    // Add bp counts
    for(;vecBeg!=vecEnd; ++vecBeg) {
      TArrayType* bpPoint = &bp_count[vecBeg->start - 1];
      TArrayType* bpPointEnd = &bp_count[vecBeg->end];
      for(;bpPoint!=bpPointEnd; ++bpPoint) ++(*bpPoint);
    }
  }


  template<typename TPos, typename TString>
    inline void
    _buildMAPQString(std::vector< HitInterval<TPos, void> > const&, TPos const, TPos const, std::vector<TString>&)
  {
    // Nothing to do
  }

  template<typename TQual, typename TPos, typename TMapqVector>
    inline void
    _buildMAPQString(std::vector< HitInterval<TPos, TQual> > const& hit_vector, TPos const posStart, TPos const posEnd, std::vector<TMapqVector>& str)
  {
    typedef HitInterval<TPos, TQual> THit;
    typedef std::vector< THit > THits;

    // Initialize result vector
    str.resize(posEnd-posStart);
    std::fill(str.begin(), str.end(), TMapqVector());
  
    // Add mapq counts
    int searchRange = posStart - 10000;
    if (searchRange < 0) searchRange=0;
    typename THits::const_iterator vecIt = std::lower_bound(hit_vector.begin(), hit_vector.end(), THit(searchRange,searchRange,0), SortHitInterval<THit>());
    for(;vecIt!=hit_vector.end(); ++vecIt) {
      if (vecIt->end < posStart) continue;
      if (vecIt->start > posEnd) break;
      for (int i = (vecIt->start - 1); i<vecIt->end; ++i) {
	if (i >= posStart && i<posEnd) str[i-posStart].push_back(vecIt->qual);
      }
    }
  }


  template<typename TDefaultOrientation>
    inline bool
    _mateIsUpstream(TDefaultOrientation defOrient, bool firstRead, bool reverse) {
    if (firstRead) {
      if (reverse) {
	if (defOrient % 2 == 0) return false;
	else return true;
      } else {
	if (defOrient % 2 == 0) return true;
	else return false;
      }
    } else {
      if (reverse) {
	if ((defOrient==1) || (defOrient==2)) return false;
	else return true;
      } else {
	if ((defOrient==1) || (defOrient==2)) return true;
	else return false;
      }
    }
  }   

  template<typename TArrayType, typename THits, typename TCountMapIterator>
    inline void
    _addCounts(TArrayType* normalCount, TArrayType* missingCount, THits const&, THits const&, TCountMapIterator& countMapIt, TCountMapIterator& abCountMapIt, int const posStart, int const posEnd, int) {
    TArrayType* normalCountPoint = &normalCount[posStart];
    TArrayType* missingCountPoint = &missingCount[posStart];
    for(int i=posStart; i<posEnd; ++i, ++normalCountPoint, ++missingCountPoint) {
      countMapIt->second.push_back(*normalCountPoint);
      abCountMapIt->second.push_back(*missingCountPoint);
    }
  }

  template<typename TArrayType, typename THits, typename TCountMapIterator>
    inline void
    _addCounts(TArrayType*, TArrayType*, THits const& normalSpan, THits const& missingSpan, TCountMapIterator& countMapIt, TCountMapIterator& abCountMapIt, int const posStart, int const posEnd, std::vector<uint16_t>) {
    std::vector<std::vector<uint16_t> > normalStr;
    std::vector<std::vector<uint16_t> > missingStr;
    _buildMAPQString(normalSpan, posStart, posEnd, normalStr);
    _buildMAPQString(missingSpan, posStart, posEnd, missingStr);
    for(int i=posStart; i<posEnd; ++i) {
      countMapIt->second.push_back(normalStr[i - posStart]);
      abCountMapIt->second.push_back(missingStr[i-posStart]);
    }
  }

  template<typename TRef, typename TQuality>
    inline TQuality
    _pairQuality(TRef const RefID, TRef const MateRefID, TQuality const q1, TQuality const q2, SVType<TranslocationTag>) {
    if (RefID==MateRefID) return q2;
    else return std::min(q1, q2);
  }

  template<typename TRef, typename TQuality, typename TTag>
    inline TQuality
    _pairQuality(TRef const, TRef const, TQuality const q1, TQuality const q2, SVType<TTag>) {
    return std::min(q1, q2);
  }

  template<typename TFiles, typename TSampleLibrary, typename TSVs, typename TCountMap, typename THitInterval, typename TSVType>
    inline void
    annotateSpanningCoverage(TFiles const& files, int const bpWindowOffset, uint16_t const minMapQual, TSampleLibrary& sampleLib, TSVs& svs, TCountMap& normalCountMap, TCountMap& abnormalCountMap, THitInterval, TSVType svType)
  {
    typedef typename TCountMap::key_type TSampleSVPair;
    typedef typename TCountMap::mapped_type TCountRange;
    typedef typename TCountRange::value_type TCount;
    typedef typename TSampleLibrary::mapped_type TLibraryMap;

    // References
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
      }
    }

    // Sort Structural Variants
    std::sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());

    // Qualities
    typedef boost::unordered_map<unsigned int, uint16_t> TQualities;
    std::vector<TQualities> qualities;
    qualities.resize(files.size());

    // Process chromosome by chromosome
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Breakpoint spanning coverage annotation" << std::endl;
    boost::progress_display show_progress( (references.end() - references.begin()) );
    typename BamTools::RefVector::const_iterator  itRef = references.begin();
    for(int refIndex=0;itRef!=references.end();++itRef, ++refIndex) {
      ++show_progress;

      // Collect all SV sizes on this chromosome
      typedef std::vector<std::pair<int, int> > TSVSizes;
      TSVSizes svSizes;
      typename TSVs::const_iterator itSV = svs.begin();
      typename TSVs::const_iterator itSVEnd = svs.end();
      for(;itSV!=itSVEnd;++itSV) {
	if (itSV->chr == references[refIndex].RefName) svSizes.push_back(std::make_pair(itSV->svStart, itSV->svEnd));
	else if (itSV->chr2 == references[refIndex].RefName) svSizes.push_back(std::make_pair(itSV->svEnd, itSV->svStart));
      }
      std::sort(svSizes.begin(), svSizes.end());
      if (svSizes.empty()) continue;

      // Iterate all samples
      #pragma omp parallel for default(shared)
      for(unsigned int file_c = 0; file_c < files.size(); ++file_c) {
	// Store all spanning ranges
	typedef std::vector<THitInterval> THits;
	THits normalSpan;  
	THits missingSpan;  

	// Get a sample name
	std::string sampleName(files[file_c].stem().string());
	typename TSampleLibrary::iterator sampleIt=sampleLib.find(sampleName);

	// Initialize bam file
	BamTools::BamReader reader;
	reader.Open(files[file_c].string());
	reader.LocateIndex();
	
	// Unique pairs for the given sample
	typedef std::set<Hit> TUniquePairs;
	TUniquePairs unique_pairs;

	// Read alignments
	BamTools::BamAlignment al;
	if (reader.Jump(refIndex, 0)) {
	  while( reader.GetNextAlignmentCore(al) ) {
	    if (al.RefID != refIndex) break;
	    if (!(al.AlignmentFlag & 0x0001) || (al.AlignmentFlag & 0x0004) || (al.AlignmentFlag & 0x0008) || (al.AlignmentFlag & 0x0100) || (al.AlignmentFlag & 0x0200) || (al.AlignmentFlag & 0x0400) || (al.MapQuality < minMapQual)) continue;

	    // Mapping positions valid?
	    if (_mappingPosGeno(al.RefID, al.MateRefID, al.Position, al.MatePosition, svType)) continue;

	    // Get the library information
	    al.BuildCharData();
	    std::string rG = "DefaultLib";
	    al.GetTag("RG", rG);
	    typename TLibraryMap::iterator libIt=sampleIt->second.find(rG);
	    if (libIt->second.median == 0) continue; // Single-end library
	    int outerISize = std::abs(al.Position - al.MatePosition) + al.Length;

	    if ((getStrandIndependentOrientation(al) != libIt->second.defaultOrient) || (outerISize > libIt->second.maxNormalISize) || (al.RefID!=al.MateRefID)) {
	      if (_acceptedInsertSize(libIt->second.maxNormalISize, libIt->second.median, abs(al.InsertSize), svType)) continue;  // Normal paired-end (for deletions only)
	      if (_acceptedOrientation(libIt->second.defaultOrient, getStrandIndependentOrientation(al), svType)) continue;  // Orientation disagrees with SV type

	      // Does the pair confirm a valid SV size
	      int32_t const minPos = _minCoord(al.Position, al.MatePosition, svType);
	      int32_t const maxPos = _maxCoord(al.Position, al.MatePosition, svType);

	      typename TSVSizes::const_iterator itSize = std::lower_bound(svSizes.begin(), svSizes.end(), std::make_pair(minPos, maxPos));
	      bool validSize = false;
	      if (itSize != svSizes.end()) {
		++itSize;
		if (itSize != svSizes.end()) validSize = (!_pairsDisagree(minPos, maxPos, al.Length, libIt->second.median, itSize->first, itSize->second, al.Length, libIt->second.median, _getSpanOrientation(al, libIt->second.defaultOrient, svType), _getSpanOrientation(al, libIt->second.defaultOrient, svType), svType));
		--itSize;
	      }
	      if ((!validSize) && (itSize != svSizes.end())) {
		validSize = (!_pairsDisagree(minPos, maxPos, al.Length, libIt->second.median, itSize->first, itSize->second, al.Length, libIt->second.median, _getSpanOrientation(al, libIt->second.defaultOrient, svType), _getSpanOrientation(al, libIt->second.defaultOrient, svType), svType));
	      }
	      if ((!validSize) && (!svSizes.empty()) && (itSize != svSizes.begin())) {
		--itSize;
		validSize = (!_pairsDisagree(itSize->first, itSize->second, al.Length, libIt->second.median, minPos, maxPos, al.Length, libIt->second.median, _getSpanOrientation(al, libIt->second.defaultOrient, svType), _getSpanOrientation(al, libIt->second.defaultOrient, svType), svType));
	      }
	      if (!validSize) continue;
	    }

	    // Get or store the mapping quality for the partner
	    if (_firstPairObs(al.RefID, al.MateRefID, al.Position, al.MatePosition, svType)) {
	      // Hash the quality
	      unsigned int index=((al.Position % (int)boost::math::pow<14>(2))<<14) + (al.MatePosition % (int)boost::math::pow<14>(2));
	      qualities[file_c][index]=al.MapQuality;
	    } else {
	      // Get the two mapping qualities
	      unsigned int index=((al.MatePosition % (int)boost::math::pow<14>(2))<<14) + (al.Position % (int)boost::math::pow<14>(2));
	      uint16_t pairQuality = _pairQuality(al.RefID, al.MateRefID, qualities[file_c][index], al.MapQuality, svType);
	      qualities[file_c][index]=0;

	      // Pair quality
	      if (pairQuality < minMapQual) continue;

	      // Is it a unique pair
	      Hit hitPos(al);
	      typename TUniquePairs::const_iterator pos = unique_pairs.begin();
	      bool inserted;
	      boost::tie(pos, inserted) = unique_pairs.insert(hitPos);
	      if (inserted) {
		// Insert the interval
		if ((getStrandIndependentOrientation(al) == libIt->second.defaultOrient) && (outerISize >= libIt->second.minNormalISize) && (outerISize <= libIt->second.maxNormalISize) && (al.RefID==al.MateRefID)) {
		  // Normal spanning coverage
		  normalSpan.push_back(THitInterval(std::min(al.Position, al.MatePosition), std::max(al.Position, al.MatePosition) + al.Length, pairQuality));
		} else if ((getStrandIndependentOrientation(al) != libIt->second.defaultOrient) || (outerISize > libIt->second.maxNormalISize) || (al.RefID!=al.MateRefID)) {
		  // Missing spanning coverage
		  if (_mateIsUpstream(libIt->second.defaultOrient, (al.AlignmentFlag & 0x0040), (al.AlignmentFlag & 0x0010))) 
		    missingSpan.push_back(THitInterval(al.Position, al.Position + libIt->second.median, pairQuality));
		  //missingSpan.push_back(THitInterval(al.Position, al.Position + libIt->second.maxNormalISize, pairQuality));
		  else
		    missingSpan.push_back(THitInterval(std::max(0, al.Position + al.Length - libIt->second.median), al.Position + al.Length, pairQuality));
		  //missingSpan.push_back(THitInterval(std::max(0, al.Position + al.Length - libIt->second.maxNormalISize), al.Position + al.Length, pairQuality));
		  if (al.RefID==al.MateRefID) {
		    if (_mateIsUpstream(libIt->second.defaultOrient, !(al.AlignmentFlag & 0x0040), (al.AlignmentFlag & 0x0020)))
		      missingSpan.push_back(THitInterval(al.MatePosition, al.MatePosition + libIt->second.median, pairQuality));
		    //missingSpan.push_back(THitInterval(al.MatePosition, al.MatePosition + libIt->second.maxNormalISize, pairQuality));
		    else
		      missingSpan.push_back(THitInterval(std::max(0, al.MatePosition + al.Length - libIt->second.median), al.MatePosition + al.Length, pairQuality));
		    //missingSpan.push_back(THitInterval(std::max(0, al.MatePosition + al.Length - libIt->second.maxNormalISize), al.MatePosition + al.Length, pairQuality));
		  }
		}
		++libIt->second.unique_pairs;
	      } else {
		++libIt->second.non_unique_pairs;
	      }
	    }
	  }
	}

	// Sort SAM records by start position
	sort(normalSpan.begin(), normalSpan.end(), SortHitInterval<THitInterval>());
	sort(missingSpan.begin(), missingSpan.end(), SortHitInterval<THitInterval>());

	// Reset the chromosome array
	typedef unsigned short TArrayType;
	unsigned int arrayLen=itRef->RefLength + 1000000;
	TArrayType* normalCount = new TArrayType[arrayLen];
	TArrayType* missingCount = new TArrayType[arrayLen];
	memset(normalCount, 0, arrayLen * sizeof(TArrayType));
	memset(missingCount, 0, arrayLen * sizeof(TArrayType));
	_addReadAndBpCounts(normalSpan, normalCount);
	_addReadAndBpCounts(missingSpan, missingCount);

	// Store spanning coverage for all input intervals
	typename TSVs::const_iterator itSV = svs.begin();
	typename TSVs::const_iterator itSVEnd = svs.end();
	#pragma omp critical
	{
	  for(;itSV!=itSVEnd;++itSV) {
	    if (itSV->chr == references[refIndex].RefName) {
	      // First breakpoint
	      TSampleSVPair svSample = std::make_pair(sampleName, -itSV->id);
	      typename TCountMap::iterator countMapIt=normalCountMap.find(svSample);
	      typename TCountMap::iterator abCountMapIt=abnormalCountMap.find(svSample);
	      if (countMapIt==normalCountMap.end()) {
		normalCountMap.insert(std::make_pair(svSample, TCountRange()));
		abnormalCountMap.insert(std::make_pair(svSample, TCountRange()));
		countMapIt=normalCountMap.find(svSample);
		abCountMapIt=abnormalCountMap.find(svSample);
	      }
	      int posStart = (itSV->svStart - bpWindowOffset < 0) ? 0 : (itSV->svStart - bpWindowOffset);
	      int posEnd = (bpWindowOffset) ? (itSV->svStart + bpWindowOffset) : (itSV->svStart + 1);
	      _addCounts(normalCount, missingCount, normalSpan, missingSpan, countMapIt, abCountMapIt, posStart, posEnd, TCount());
	    }

	    if (itSV->chr2 == references[refIndex].RefName) {
	      // Second breakpoint
	      TSampleSVPair svSample = std::make_pair(sampleName, itSV->id);
	      typename TCountMap::iterator countMapIt=normalCountMap.find(svSample);
	      typename TCountMap::iterator abCountMapIt=abnormalCountMap.find(svSample);
	      if (countMapIt==normalCountMap.end()) {
		normalCountMap.insert(std::make_pair(svSample, TCountRange()));
		abnormalCountMap.insert(std::make_pair(svSample, TCountRange()));
		countMapIt=normalCountMap.find(svSample);
		abCountMapIt=abnormalCountMap.find(svSample);
	      }
	      int posStart = (itSV->svEnd - bpWindowOffset < 0) ? 0 : (itSV->svEnd - bpWindowOffset);
	      int posEnd = (bpWindowOffset) ? (itSV->svEnd + bpWindowOffset) : (itSV->svEnd + 1);
	      _addCounts(normalCount, missingCount, normalSpan, missingSpan, countMapIt, abCountMapIt, posStart, posEnd, TCount());
	    }
	  }
	}

	// Clean-up
	delete[] normalCount;
	delete[] missingCount;

	// Clean-up qualities
	_resetQualities(qualities[file_c], svType);
      }
    }
  }

}

#endif
