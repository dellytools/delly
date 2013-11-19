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
    struct SortHits : public std::binary_function<THitInterval, THitInterval, bool>
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
    _buildMAPQString(std::vector< HitInterval<TPos, void> >&, TPos const, TPos const, std::vector<TString>&)
  {
    // Nothing to do
  }

  template<typename TQual, typename TPos, typename TMapqVector>
    inline void
    _buildMAPQString(std::vector< HitInterval<TPos, TQual> >& hit_vector, TPos const posStart, TPos const posEnd, std::vector<TMapqVector>& str)
  {
    typedef HitInterval<TPos, TQual> THit;
    typedef std::vector< THit > THits;
    typename THits::const_iterator vecBeg = hit_vector.begin();
    typename THits::const_iterator vecEnd = hit_vector.end();

    // Initialize result vector
    str.resize(posEnd-posStart);
    std::fill(str.begin(), str.end(), TMapqVector());
  
    // Add mapq counts
    int searchRange = posStart - 10000;
    if (searchRange < 0) searchRange=0;
    THit hit;
    hit.start=searchRange;
    hit.end=searchRange;
    typename THits::const_iterator vecIt = std::lower_bound(vecBeg, vecEnd, hit, SortHits<THit>());
    for(;vecIt!=vecEnd; ++vecIt) {
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
    _addCounts(TArrayType* normalCount, TArrayType* missingCount, THits&, THits&, TCountMapIterator& countMapIt, TCountMapIterator& abCountMapIt, int posStart, int posEnd, int) {
    TArrayType* normalCountPoint = &normalCount[posStart];
    TArrayType* missingCountPoint = &missingCount[posStart];
    for(int i=posStart; i<posEnd; ++i, ++normalCountPoint, ++missingCountPoint) {
      countMapIt->second.push_back(*normalCountPoint);
      abCountMapIt->second.push_back(*missingCountPoint);
    }
  }

  template<typename TArrayType, typename THits, typename TCountMapIterator>
    inline void
    _addCounts(TArrayType*, TArrayType*, THits& normalSpan, THits& missingSpan, TCountMapIterator& countMapIt, TCountMapIterator& abCountMapIt, int posStart, int posEnd, std::vector<uint16_t>) {
    std::vector<std::vector<uint16_t> > normalStr;
    std::vector<std::vector<uint16_t> > missingStr;
    _buildMAPQString(normalSpan, posStart, posEnd, normalStr);
    _buildMAPQString(missingSpan, posStart, posEnd, missingStr);
    for(int i=posStart; i<posEnd; ++i) {
      countMapIt->second.push_back(normalStr[i - posStart]);
      abCountMapIt->second.push_back(missingStr[i-posStart]);
    }
  }


  template<typename TFiles, typename TSampleLibrary, typename TSVs, typename TCountMap, typename THitInterval>
    inline void
    annotateSpanningCoverage(TFiles const& files, int bpWindowOffset, uint16_t minMapQual, TSampleLibrary& sampleLib, TSVs const& svs, TCountMap& normalCountMap, TCountMap& abnormalCountMap, THitInterval)
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

    // Process chromosome by chromosome
    std::cout << "Breakpoint spanning coverage annotation" << std::endl;
    boost::progress_display show_progress( (references.end() - references.begin()) );
    typename BamTools::RefVector::const_iterator  itRef = references.begin();
    for(int refIndex=0;itRef!=references.end();++itRef, ++refIndex) {
      ++show_progress;

      // Iterate all samples
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
	if ( ! reader.Open(files[file_c].string()) ) return;
	reader.LocateIndex();
	if ( !reader.HasIndex() ) return;
	

	// Unique pairs for the given sample
	typedef std::set<Hit> TUniquePairs;
	TUniquePairs unique_pairs;

	// Read alignments and hash qualities
	uint16_t* qualities = new uint16_t[(int)boost::math::pow<28>(2)];
	uint16_t* qualitiesEnd = qualities + (int) boost::math::pow<28>(2);
	std::fill(qualities, qualitiesEnd, 0);
	BamTools::BamAlignment al;
	if (reader.Jump(refIndex, 0)) {
	  while( reader.GetNextAlignment(al) ) {
	    if (al.RefID != refIndex) break;
	    if (!(al.AlignmentFlag & 0x0001) || (al.AlignmentFlag & 0x0004) || (al.AlignmentFlag & 0x0008) || (al.AlignmentFlag & 0x0100) || (al.AlignmentFlag & 0x0200) || (al.AlignmentFlag & 0x0400) || (al.Position==al.MatePosition) || (al.RefID!=al.MateRefID) || (al.MapQuality < minMapQual)) continue;
	    // Get the library information
	    std::string rG = "DefaultLib";
	    al.GetTag("RG", rG);
	    typename TLibraryMap::iterator libIt=sampleIt->second.find(rG);
	    if (libIt->second.median == 0) continue; // Single-end library

	    // Get or store the mapping quality for the partner
	    if (al.Position<al.MatePosition) {
	      // Hash the quality
	      boost::hash<std::string> hashStr;
	      unsigned int index=((hashStr(al.Name) % (int)boost::math::pow<4>(2))<<24) + ((al.Position % (int)boost::math::pow<12>(2))<<12) + (al.MatePosition % (int)boost::math::pow<12>(2));

	      qualities[index]=al.MapQuality;
	    } else {
	      // Get the two mapping qualities
	      boost::hash<std::string> hashStr;
	      unsigned int index=((hashStr(al.Name) % (int)boost::math::pow<4>(2))<<24) + ((al.MatePosition % (int)boost::math::pow<12>(2))<<12) + (al.Position % (int)boost::math::pow<12>(2));
	      uint16_t pairQuality = std::min(qualities[index], al.MapQuality);
	      if (pairQuality < minMapQual) continue;

	      // Is it a unique pair
	      Hit hitPos(al);
	      typename TUniquePairs::const_iterator pos = unique_pairs.begin();
	      bool inserted;
	      boost::tie(pos, inserted) = unique_pairs.insert(hitPos);
	      if (inserted) {
		// Insert the interval
		unsigned int outerISize = (al.Position + al.Length) - al.MatePosition;
		if ((getStrandIndependentOrientation(al) == libIt->second.defaultOrient) && (outerISize >= libIt->second.minNormalISize) && (outerISize <= libIt->second.maxNormalISize)) {
		  // Normal spanning coverage
		  normalSpan.push_back(THitInterval(al.MatePosition, al.Position+al.Length, pairQuality));
		} else if ((getStrandIndependentOrientation(al) != libIt->second.defaultOrient) || (outerISize > libIt->second.maxNormalISize)) {
		  // Missing spanning coverage
		  if (_mateIsUpstream(libIt->second.defaultOrient, (al.AlignmentFlag & 0x0040), (al.AlignmentFlag & 0x0010))) 
		    missingSpan.push_back(THitInterval(al.Position, al.Position + libIt->second.median, pairQuality));
		  else
		    missingSpan.push_back(THitInterval(std::max(0, al.Position + al.Length - libIt->second.median), al.Position + al.Length, pairQuality));
		  if (_mateIsUpstream(libIt->second.defaultOrient, !(al.AlignmentFlag & 0x0040), (al.AlignmentFlag & 0x0020)))
		    missingSpan.push_back(THitInterval(al.MatePosition, al.MatePosition + libIt->second.median, pairQuality));
		  else
		    missingSpan.push_back(THitInterval(std::max(0, al.MatePosition + al.Length - libIt->second.median), al.MatePosition + al.Length, pairQuality));
		}
		++libIt->second.unique_pairs;
	      } else {
		++libIt->second.non_unique_pairs;
	      }
	    }
	  }
	}
	delete [] qualities;

	// Sort SAM records by start position
	sort(normalSpan.begin(), normalSpan.end(), SortHits<THitInterval>());
	sort(missingSpan.begin(), missingSpan.end(), SortHits<THitInterval>());

	// Declare the chromosome array
	typedef unsigned short TArrayType;
	TArrayType* normalCount = new TArrayType[MAX_CHROM_SIZE];
	TArrayType* missingCount = new TArrayType[MAX_CHROM_SIZE];
	std::fill(normalCount, normalCount + MAX_CHROM_SIZE, 0);
	std::fill(missingCount, missingCount + MAX_CHROM_SIZE, 0);
	_addReadAndBpCounts(normalSpan, normalCount);
	_addReadAndBpCounts(missingSpan, missingCount);

	// Write spanning coverage for all input intervals
	typename TSVs::const_iterator itSV = svs.begin();
	typename TSVs::const_iterator itSVEnd = svs.end();
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
	    
	    // Second breakpoint
	    svSample = std::make_pair(sampleName, itSV->id);
	    countMapIt=normalCountMap.find(svSample);
	    abCountMapIt=abnormalCountMap.find(svSample);
	    if (countMapIt==normalCountMap.end()) {
	      normalCountMap.insert(std::make_pair(svSample, TCountRange()));
	      abnormalCountMap.insert(std::make_pair(svSample, TCountRange()));
	      countMapIt=normalCountMap.find(svSample);
	      abCountMapIt=abnormalCountMap.find(svSample);
	    }
	    posStart = (itSV->svEnd - bpWindowOffset < 0) ? 0 : (itSV->svEnd - bpWindowOffset);
	    posEnd = (bpWindowOffset) ? (itSV->svEnd + bpWindowOffset) : (itSV->svEnd + 1);
	    _addCounts(normalCount, missingCount, normalSpan, missingSpan, countMapIt, abCountMapIt, posStart, posEnd, TCount());
	  }
	}

	// Clean-up
	delete[] normalCount;
	delete[] missingCount;
      }
    }
  }


}

#endif
