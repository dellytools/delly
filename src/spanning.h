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

#include <boost/container/flat_set.hpp>
#include <boost/unordered_map.hpp>
#include "tags.h"

namespace torali {


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

  template<typename TFiles, typename TSampleLibrary, typename TSVs, typename TCountMap, typename TSVType>
    inline void
    annotateSpanningCoverage(TFiles const& files, uint16_t const minMapQual, TSampleLibrary& sampleLib, TSVs& svs, TCountMap& spanCountMap, TSVType svType)
  {
    typedef typename TCountMap::key_type TSampleSVPair;
    typedef typename TCountMap::mapped_type TCountPair;
    typedef typename TSampleLibrary::mapped_type TLibraryMap;

    // References
    BamTools::BamReader readerRef;
    if ( ! readerRef.Open(files[0].string()) ) return;
    BamTools::RefVector references = readerRef.GetReferenceData();

    // Get maximum insert size
    int maxInsertSize = 0;
    typename TSampleLibrary::iterator sampleIt=sampleLib.begin();
    for(;sampleIt!=sampleLib.end();++sampleIt) {
      typename TLibraryMap::iterator libIt=sampleIt->second.begin();
      for(;libIt!=sampleIt->second.end();++libIt) {
	if (libIt->second.median > maxInsertSize) maxInsertSize=libIt->second.median;
      }
    }

    // Sort Structural Variants
    std::sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());

    // Initialize count map
    for(typename TSampleLibrary::iterator sIt = sampleLib.begin(); sIt!=sampleLib.end();++sIt) {
      for(typename TSVs::const_iterator itSV = svs.begin(); itSV!=svs.end(); ++itSV) {
	// Left breakpoint
	spanCountMap.insert(std::make_pair(std::make_pair(sIt->first, -itSV->id), TCountPair()));
	// Right breakpoint
	spanCountMap.insert(std::make_pair(std::make_pair(sIt->first, itSV->id), TCountPair()));
      }
    }

    // Iterate all samples
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Breakpoint spanning coverage annotation" << std::endl;
    boost::progress_display show_progress( (svs.end() - svs.begin()) );
#pragma omp parallel for default(shared)
    for(unsigned int file_c = 0; file_c < files.size(); ++file_c) {
      // Get a sample name
      std::string sampleName(files[file_c].stem().string());
      typename TSampleLibrary::iterator sampleIt=sampleLib.find(sampleName);

      // Initialize bam file
      BamTools::BamReader reader;
      reader.Open(files[file_c].string());
      reader.LocateIndex();
	
      // Read alignments
      typename TSVs::const_iterator itSV = svs.begin();
      typename TSVs::const_iterator itSVEnd = svs.end();
      for(;itSV!=itSVEnd;++itSV) {
	if (file_c==(files.size()-1)) ++show_progress;

	// Set up the mapping quality iterators
	TSampleSVPair svSample = std::make_pair(sampleName, -itSV->id);
	typename TCountMap::iterator leftIt = spanCountMap.find(svSample);
	svSample = std::make_pair(sampleName, itSV->id);
	typename TCountMap::iterator rightIt = spanCountMap.find(svSample);

	// Processed reads
	typedef std::set<std::string> TProcessedReads;
	TProcessedReads procReads;

	// Qualities
	typedef boost::unordered_map<unsigned int, uint16_t> TQualities;
	TQualities qualities;

	// Unique pairs for the given sample
	typedef boost::container::flat_set<int32_t> TUniquePairs;
	TUniquePairs unique_pairs;

	// Scan left and right breakpoint
	BamTools::BamAlignment al;
	for (unsigned int bpPoint = 0; bpPoint<2; ++bpPoint) {
	  int32_t regionChr;
	  int regionStart;
	  int regionEnd;
	  if (bpPoint==(unsigned int)(itSV->chr==itSV->chr2)) {
	    regionChr = itSV->chr2;
	    regionStart = std::max(0, (int) itSV->svEnd - (int) maxInsertSize);
	    regionEnd = itSV->svEnd + maxInsertSize;
	  } else {
	    regionChr = itSV->chr;
	    regionStart = std::max(0, (int) itSV->svStart - (int) maxInsertSize);
	    regionEnd = itSV->svStart + maxInsertSize;
	  }
	  int32_t oldAlignPos=-1;
	  if (reader.SetRegion(regionChr, regionStart, regionChr, regionEnd)) {
	    while( reader.GetNextAlignmentCore(al) ) {
	      if (!(al.AlignmentFlag & 0x0001) || (al.AlignmentFlag & 0x0004) || (al.AlignmentFlag & 0x0008) || (al.AlignmentFlag & 0x0100) || (al.AlignmentFlag & 0x0200) || (al.AlignmentFlag & 0x0400) || (al.AlignmentFlag & 0x0800) || (al.MapQuality < minMapQual)) continue;

	      // Mapping positions valid?
	      if (_mappingPosGeno(al.RefID, al.MateRefID, al.Position, al.MatePosition, svType)) continue;

	      // Get the library information
	      al.BuildCharData();
	      std::string rG = "DefaultLib";
	      al.GetTag("RG", rG);
	      typename TLibraryMap::iterator libIt=sampleIt->second.find(rG);
	      if (libIt->second.median == 0) continue; // Single-end library
	      int outerISize = std::abs(al.Position - al.MatePosition) + al.Length;
	      
	      // Have we processed this read already?
	      std::string readId = al.Name;
	      if (al.AlignmentFlag & 0x0040) readId.append("First");
	      else readId.append("Second");
	      typename TProcessedReads::const_iterator procReadPos = procReads.begin();
	      bool inserted;
	      boost::tie(procReadPos, inserted) = procReads.insert(readId);
	      if (!inserted) continue;
	      
	      // Abnormal paired-end
	      if ((getStrandIndependentOrientation(al) != libIt->second.defaultOrient) || (outerISize > libIt->second.maxNormalISize) || (al.RefID!=al.MateRefID)) {
		if (_acceptedInsertSize(libIt->second, abs(al.InsertSize), svType)) continue;  // Normal paired-end (for deletions, insertions only)
		if (_acceptedOrientation(libIt->second.defaultOrient, getStrandIndependentOrientation(al), svType)) continue;  // Orientation disagrees with SV type
		
		// Does the pair confirm the SV
		int32_t const minPos = _minCoord(al.Position, al.MatePosition, svType);
		int32_t const maxPos = _maxCoord(al.Position, al.MatePosition, svType);
		
		bool validSize = true;
		if (al.RefID==itSV->chr) {
		  if (minPos < itSV->svStart) {
		    validSize = (!_pairsDisagree(minPos, maxPos, al.Length, libIt->second.maxNormalISize, itSV->svStart, itSV->svEnd, al.Length, libIt->second.maxNormalISize, _getSpanOrientation(al, libIt->second.defaultOrient, svType), itSV->ct, svType));
		  } else {
		    validSize = (!_pairsDisagree(itSV->svStart, itSV->svEnd, al.Length, libIt->second.maxNormalISize, minPos, maxPos, al.Length, libIt->second.maxNormalISize, itSV->ct, _getSpanOrientation(al, libIt->second.defaultOrient, svType), svType));
		  }
		}
		if (!validSize) continue;
	      }
		    
	      // Get or store the mapping quality for the partner
	      if (_firstPairObs(al.RefID, al.MateRefID, al.Position, al.MatePosition, svType)) {
		// Hash the quality
		unsigned int index=((al.Position % (int)boost::math::pow<14>(2))<<14) + (al.MatePosition % (int)boost::math::pow<14>(2));
		qualities[index]=al.MapQuality;
	      } else {
		// Get the two mapping qualities
		unsigned int index=((al.MatePosition % (int)boost::math::pow<14>(2))<<14) + (al.Position % (int)boost::math::pow<14>(2));
		uint16_t pairQuality = _pairQuality(al.RefID, al.MateRefID, qualities[index], al.MapQuality, svType);
		qualities[index]=0;
		
		// Pair quality
		if (pairQuality < minMapQual) continue;

		// Is it a unique pair
		if (al.Position!=oldAlignPos) {
		  oldAlignPos=al.Position;
		  unique_pairs.clear();
		}
		if (unique_pairs.insert(al.MatePosition).second) {
		  // Insert the interval
		  if ((getStrandIndependentOrientation(al) == libIt->second.defaultOrient) && (outerISize >= libIt->second.minNormalISize) && (outerISize <= libIt->second.maxNormalISize) && (al.RefID==al.MateRefID)) {
		    // Normal spanning coverage
		    int32_t sPos = std::min(al.Position, al.MatePosition);
		    int32_t ePos = std::max(al.Position, al.MatePosition) + al.Length;
		    int32_t midPoint = sPos+(ePos-sPos)/2;
		    sPos=std::max(sPos, midPoint - al.Length);
		    ePos=std::min(ePos, midPoint + al.Length);
		    int32_t innerSPos = std::min(al.Position, al.MatePosition) + al.Length;
		    int32_t innerEPos = std::max(al.Position, al.MatePosition);
		    if ((innerSPos<innerEPos) && ((innerEPos - innerSPos) > (ePos-sPos))) {
		      if ((itSV->svStart>=innerSPos) && (itSV->svStart<=innerEPos)) leftIt->second.first.push_back(pairQuality);
		      if ((itSV->svEnd>=innerSPos) && (itSV->svEnd<=innerEPos)) rightIt->second.first.push_back(pairQuality);
		    } else {
		      if ((itSV->svStart>=sPos) && (itSV->svStart<=ePos)) leftIt->second.first.push_back(pairQuality);
		      if ((itSV->svEnd>=sPos) && (itSV->svEnd<=ePos)) rightIt->second.first.push_back(pairQuality);
		    }		      
		  } else if ((getStrandIndependentOrientation(al) != libIt->second.defaultOrient) || (outerISize > libIt->second.maxNormalISize) || (al.RefID!=al.MateRefID)) {
		    // Missing spanning coverage
		    if (_mateIsUpstream(libIt->second.defaultOrient, (al.AlignmentFlag & 0x0040), (al.AlignmentFlag & 0x0010))) {
		      if ((itSV->chr==al.RefID) && (itSV->svStart>=al.Position) && (itSV->svStart<=(al.Position + libIt->second.maxNormalISize))) leftIt->second.second.push_back(pairQuality);
		      if ((itSV->chr2==al.RefID) && (itSV->svEnd>=al.Position) && (itSV->svEnd<=(al.Position + libIt->second.maxNormalISize))) rightIt->second.second.push_back(pairQuality);
		    } else {
		      if ((itSV->chr==al.RefID) && (itSV->svStart>=std::max(0, al.Position + al.Length - libIt->second.maxNormalISize)) && (itSV->svStart<=(al.Position + al.Length))) leftIt->second.second.push_back(pairQuality);
		      if ((itSV->chr2==al.RefID) && (itSV->svEnd>=std::max(0, al.Position + al.Length - libIt->second.maxNormalISize)) && (itSV->svEnd<=(al.Position + al.Length))) rightIt->second.second.push_back(pairQuality);
		    }
		    if (_mateIsUpstream(libIt->second.defaultOrient, !(al.AlignmentFlag & 0x0040), (al.AlignmentFlag & 0x0020))) {
		      if ((itSV->chr==al.MateRefID) && (itSV->svStart>=al.MatePosition) && (itSV->svStart<=(al.MatePosition + libIt->second.maxNormalISize))) leftIt->second.second.push_back(pairQuality);
		      if ((itSV->chr2==al.MateRefID) && (itSV->svEnd>=al.MatePosition) && (itSV->svEnd<=(al.MatePosition + libIt->second.maxNormalISize))) rightIt->second.second.push_back(pairQuality);
		    } else {
		      if ((itSV->chr==al.MateRefID) && (itSV->svStart>=std::max(0, al.MatePosition + al.Length - libIt->second.maxNormalISize)) && (itSV->svStart<=(al.MatePosition + al.Length))) leftIt->second.second.push_back(pairQuality);
		      if ((itSV->chr2==al.MateRefID) && (itSV->svEnd>=std::max(0, al.MatePosition + al.Length - libIt->second.maxNormalISize)) && (itSV->svEnd<=(al.MatePosition + al.Length))) rightIt->second.second.push_back(pairQuality);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

}

#endif
