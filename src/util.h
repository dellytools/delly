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

#ifndef UTIL_H
#define UTIL_H

#include <math.h>
#include "tags.h"


namespace torali
{


  struct LibraryInfo {
    int median;
    int mad;
    int percentileCutoff;
    int minNormalISize;
    int maxNormalISize;
    int defaultOrient;
    unsigned int mappedReads;
    unsigned int non_unique_abnormal_pairs;
    unsigned int abnormal_pairs;
    unsigned int non_unique_pairs;
    unsigned int unique_pairs;

    
  LibraryInfo() : median(0), mad(0), percentileCutoff(0), maxNormalISize(0), defaultOrient(0), mappedReads(0), non_unique_pairs(0), unique_pairs(0) {}
  };


  struct _LibraryParams {
    unsigned int processedNumPairs;
    unsigned int orient[4];
    std::vector<unsigned int> vecISize;
  };

  
  template<typename TAlphabet>
    inline void
    reverseComplement(std::vector<TAlphabet>& seq) 
    {
      typedef std::vector<TAlphabet> TSequence;
      typename TSequence::iterator itF = seq.begin();
      typename TSequence::iterator itB = seq.end();
      if (itF != itB) {
	--itB;
	do {
	  TAlphabet tmp = *itF;
	  switch((unsigned int) *itB) {
	  case 0: *itF = 3; break;
	  case 1: *itF = 2; break;
	  case 2: *itF = 1; break;
	  case 3: *itF = 0; break;
	  }
	  switch((unsigned int) tmp) {
	  case 0: *itB = 3; break;
	  case 1: *itB = 2; break;
	  case 2: *itB = 1; break;
	  case 3: *itB = 0; break;
	  }
	} while (++itF <= --itB);
      }
    }

  template<typename TIterator, typename TValue>
  inline
    void getMedian(TIterator begin, TIterator end, TValue& median) 
    {
      std::nth_element(begin, begin + (end - begin) / 2, end);
      median = *(begin + (end - begin) / 2);
    }

  template<typename TIterator, typename TPercent, typename TValue>
  inline
    void getPercentileCutoff(TIterator begin, TIterator end, TPercent percentile, TValue& cutoff) 
    {
      if (percentile < 1) {
	std::nth_element(begin, begin + (int) ((end - begin) * percentile), end);
	cutoff = *(begin + (int) ((end - begin) * percentile));
      } else {
	cutoff = *(std::max_element(begin, end));
      }
    }

  template<typename TIterator, typename TValue>
  inline
    void getMAD(TIterator begin, TIterator end, TValue median, TValue& mad) 
    {
      std::vector<TValue> absDev;
      for(;begin<end;++begin) 
	absDev.push_back(std::abs((TValue)*begin - median));
      getMedian(absDev.begin(), absDev.end(), mad);
    }

  template<typename TIterator, typename TValue>
  inline
    void getMean(TIterator begin, TIterator end, TValue& mean) 
    {
      mean = 0;
      unsigned int count = 0;
      for(; begin<end; ++begin,++count) mean += *begin;
      mean /= count;
    }

  template<typename TIterator, typename TValue>
  inline
    void getStdDev(TIterator begin, TIterator end, TValue mean, TValue& stdDev) 
    {
      stdDev = 0;
      unsigned int count = 0;
      for(;begin<end;++begin,++count) stdDev += ((TValue)*begin - mean) * ((TValue)*begin - mean);
      stdDev = sqrt(stdDev / (TValue) count);
    }

  template<typename TIterator, typename TPercentile, typename TValue1, typename TValue2, typename TValue3>
  inline
    void getLibraryStats(TIterator begin, TIterator end, TPercentile percentile, TValue1& median, TValue2& mad, TValue3& percentileCutoff) 
    {
      getMedian(begin,end,median);
      getMAD(begin,end,median,mad);
      getPercentileCutoff(begin, end, 1.0 - percentile, percentileCutoff);
    }


  template<typename TLibraryMap>
  inline void getLibraryParams(boost::filesystem::path const& path, TLibraryMap& libInfo, double const percentile, unsigned short const madCutoff) {
    // Maximum number of pairs used to estimate library parameters for each RG library
    unsigned int maxNumPairs=5000000;

    // Store the counts in an object for each RG librar
    typedef std::map<std::string, _LibraryParams> TParams;
    TParams params;

    // Create SAM Object
    BamTools::BamReader reader;
    if ( ! reader.Open(path.string()) ) {
      std::cerr << "Could not open input bam file: " << path.c_str() << std::endl;
      return;
    }

    // Get read groups
    BamTools::SamHeader samHeader = reader.GetHeader();
    if (samHeader.HasReadGroups()) {
      BamTools::SamReadGroupConstIterator rgIter = samHeader.ReadGroups.ConstBegin();
      BamTools::SamReadGroupConstIterator rgIterEnd = samHeader.ReadGroups.ConstEnd();
      for(;rgIter!=rgIterEnd;++rgIter) {
	libInfo.insert(std::make_pair(rgIter->ID, LibraryInfo()));
	params.insert(std::make_pair(rgIter->ID, _LibraryParams()));
      }
    } else {
      libInfo.insert(std::make_pair("DefaultLib", LibraryInfo()));
      params.insert(std::make_pair("DefaultLib", _LibraryParams()));
    }

    // Initialize arrays
    for(TParams::iterator paramIt = params.begin(); paramIt!=params.end(); ++paramIt) {
      paramIt->second.processedNumPairs=0;
      for(unsigned int i=0;i<4;++i) paramIt->second.orient[i]=0;
      paramIt->second.vecISize.clear();
    }

    // Collect insert sizes
    bool missingPairs=true;
    BamTools::BamAlignment al;
    while ((reader.GetNextAlignmentCore(al)) && (missingPairs)) {
      if ((al.AlignmentFlag & 0x0001) && !(al.AlignmentFlag & 0x0004) && !(al.AlignmentFlag & 0x0008) && (al.AlignmentFlag & 0x0040) && (al.RefID==al.MateRefID) && !(al.AlignmentFlag & 0x0100) && !(al.AlignmentFlag & 0x0200)  && !(al.AlignmentFlag & 0x0400) && !(al.AlignmentFlag & 0x0800)) {
	al.BuildCharData();
	std::string rG = "DefaultLib";
	al.GetTag("RG", rG);
	TParams::iterator paramIt= params.find(rG);
	++paramIt->second.processedNumPairs;
	if (paramIt->second.processedNumPairs>=maxNumPairs) {
	  // Check every now and then if we have a good estimate of the library complexity
	  if (paramIt->second.processedNumPairs % 10000 == 0) {
	    TParams::const_iterator paramIter = params.begin();
	    missingPairs=false;
	    for(;paramIter!=params.end();++paramIter) {
	      if (paramIter->second.processedNumPairs<maxNumPairs) {
		missingPairs=true;
		break;
	      }
	    }
	  }
	  continue;
	}
	paramIt->second.vecISize.push_back(abs(al.InsertSize));
	++paramIt->second.orient[getStrandIndependentOrientation(al)];
      }
    }

    // Set library parameters
    typename TLibraryMap::iterator libInfoIt = libInfo.begin();
    for(;libInfoIt!=libInfo.end();++libInfoIt) {
      TParams::iterator paramIt= params.find(libInfoIt->first);

      // Get default orientation
      libInfoIt->second.defaultOrient=0;
      unsigned int maxOrient=paramIt->second.orient[0];
      for(unsigned int i=1;i<4;++i) {
	if (paramIt->second.orient[i]>maxOrient) {
	  maxOrient=paramIt->second.orient[i];
	  libInfoIt->second.defaultOrient=i;
	}
      }
      
      // Mate-pair library (If yes, trim off the chimera peak < 1000bp)
      if (libInfoIt->second.defaultOrient==3) {
	typedef std::vector<unsigned int> TVecISize;
	TVecISize vecISizeTmp;
	typename TVecISize::const_iterator iSizeBeg = paramIt->second.vecISize.begin();
	typename TVecISize::const_iterator iSizeEnd = paramIt->second.vecISize.end();
	for(;iSizeBeg<iSizeEnd;++iSizeBeg)
	  if (*iSizeBeg >= 1000) vecISizeTmp.push_back(*iSizeBeg);
	paramIt->second.vecISize = vecISizeTmp;
      }
      
      // Check that this is a proper paired-end library and require that there are at least 1000 pairs to estimate the insert size
      if (paramIt->second.vecISize.size()>=10000) {
	// Get library stats
	double median;
	double mad;
	double percentileCutoff;
	getLibraryStats(paramIt->second.vecISize.begin(), paramIt->second.vecISize.end(), percentile, median, mad, percentileCutoff);
	if ((median >= 50) && (median<=100000)) {
	  libInfoIt->second.median = (int) median;
	  libInfoIt->second.mad = (int) mad;
	  libInfoIt->second.percentileCutoff = (int) percentileCutoff;
	  if (percentile!=0) libInfoIt->second.maxNormalISize = libInfoIt->second.percentileCutoff;
	  else libInfoIt->second.maxNormalISize = libInfoIt->second.median + (madCutoff * libInfoIt->second.mad);
	  libInfoIt->second.minNormalISize = libInfoIt->second.median - (madCutoff * libInfoIt->second.mad);
	  if (libInfoIt->second.minNormalISize < 1) libInfoIt->second.minNormalISize=1;
	}
      }
    }
  }

}

#endif
