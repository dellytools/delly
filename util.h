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

#include <fstream>
#include <math.h>

#include "sam.h"
#include "fasta_reader.h"


namespace torali
{
  
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

  // File size in MB
  template<typename TFilePath>
    inline unsigned int
    fileSize(TFilePath const& fpath) 
    {
      std::ifstream f;
      f.open(fpath.c_str(), std::ios_base::binary | std::ios_base::in);
      if (!f.good() || f.eof() || !f.is_open()) { return 0; }
      f.seekg(0, std::ios_base::beg);
      std::ifstream::pos_type begin_pos = f.tellg();
      f.seekg(0, std::ios_base::end);
      return static_cast<unsigned int>(((f.tellg() - begin_pos) / 1000000));
    }

  template<typename TFilePath>
    inline bool
    isBinary(TFilePath const& fpath) 
    {
      int readChar;
      bool binaryFile = false;
      std::ifstream alignFile(fpath.c_str());
      unsigned int count = 0;
      while (((readChar = alignFile.get()) != EOF) && (count++ < 1000)) {
	if (readChar>127) {
	  binaryFile=true;
	  break;
	}
      }
      alignFile.close();
      return binaryFile;
    }


  template<typename TFilePath, typename TSize>
    inline bool
    isValidFile(TFilePath const& fpath, TSize s) 
  {
    std::ifstream empty_check(fpath.c_str(), std::ios::binary);
    if (!empty_check) { 
      //Off#std::cerr << "File does not exist: " << fpath << std::endl; 
      return false;
    }
    empty_check.seekg(0, std::ios::end);
    if ((int) empty_check.tellg() == 0) {
      empty_check.close();
      //Off#std::cerr << "File is empty: " << fpath << std::endl;
      return false;   
    }
    empty_check.close();
    if ((s != (TSize) 0) && (fileSize(fpath) > (long unsigned int) s)) { 
      //Off#std::cerr << "File too large: " << fpath << std::endl; 
      return false;
    }
    return true;
  }

  template<typename TFilePath>
    inline bool
    isValidFile(TFilePath const& fpath) 
  {
    return isValidFile(fpath, 0);
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

  template<typename TIterator, typename TValue1, typename TValue2, typename TValue3>
  inline
    void getStats(TIterator begin, TIterator end, TValue1& mean, TValue2& median, TValue3& stdDev) 
    {
      getMean(begin,end,mean);
      getMedian(begin,end,median);
      getStdDev(begin,end,mean,stdDev);
    }


  template<typename TIterator, typename TValue>
  inline
    void depecated_getLibraryStdDev(TIterator begin, TIterator end, TValue median, TValue& stdDev) 
    {
      stdDev = 0;
      TValue cutoffMax = median + 7 * 0.1 * median;
      TValue cutoffMin = median - 7 * 0.1 * median; //SD calculation cutoffs are 7 SDs to the left and right assuming 10% library deviation 
      if ((cutoffMin < 0) || (cutoffMax < cutoffMin)) cutoffMin = 0; 
      unsigned int count = 0;
      for(;begin<end;++begin) 
	if ((*begin >= cutoffMin) && (*begin <= cutoffMax)) {
	  stdDev += ((TValue)*begin - median) * ((TValue)*begin - median);
	  ++count;
	}
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


  struct LibraryInfo {
    int median;
    int mad;
    int percentileCutoff;
    int minNormalISize;
    int maxNormalISize;
    int defaultOrient;
    unsigned int non_unique_pairs;
    unsigned int unique_pairs;


  LibraryInfo() : median(0), mad(0), percentileCutoff(0), maxNormalISize(0), defaultOrient(0), non_unique_pairs(0), unique_pairs(0) {}
  };

  struct LibraryParams {
    unsigned int processedNumPairs;
    unsigned int orient[4];
    std::vector<unsigned int> vecISize;
  };

  inline void getLibraryParams(boost::filesystem::path const& path, std::map<std::string, LibraryInfo>& libInfo, double const percentile, unsigned short madCutoff) {
    typedef  std::map<std::string, LibraryInfo> TLibraryMap;

    // Maximum number of pairs used to estimate library parameters for each RG library
    unsigned int maxNumPairs=5000000;

    // Store the counts in an object for each RG librar
    typedef std::map<std::string, LibraryParams> TParams;
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
	params.insert(std::make_pair(rgIter->ID, LibraryParams()));
      }
    } else {
      libInfo.insert(std::make_pair("DefaultLib", LibraryInfo()));
      params.insert(std::make_pair("DefaultLib", LibraryParams()));
    }
    bool missingPairs=true;
    BamTools::BamAlignment al;
    while ((reader.GetNextAlignment(al)) && (missingPairs)) {
      if ((al.AlignmentFlag & 0x0001) && !(al.AlignmentFlag & 0x0004) && !(al.AlignmentFlag & 0x0008) && (al.AlignmentFlag & 0x0040) && (al.RefID==al.MateRefID)) {
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
    TLibraryMap::iterator libInfoIt = libInfo.begin();
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
      if (paramIt->second.vecISize.size()>=1000) {
	// Get library stats
	double median;
	double mad;
	double percentileCutoff;
	getLibraryStats(paramIt->second.vecISize.begin(), paramIt->second.vecISize.end(), percentile, median, mad, percentileCutoff);
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

  class Centromere {
  public:
    int centr_ext;
    typedef std::pair<int, int> TInterval;
    typedef std::map<std::string, TInterval> TCentr;
    TCentr centr;
    
    Centromere() : centr_ext(0) {
      Centromere(centr_ext);
    }
      
    Centromere(unsigned int ext) : centr_ext(ext) 
      {
      centr.insert(std::make_pair("chr1.fa", std::make_pair(std::max(121535434 - centr_ext, 0), 124535434 + centr_ext)));
      centr.insert(std::make_pair("chr2.fa", std::make_pair(std::max(92326171 - centr_ext, 0), 95326171 + centr_ext)));
      centr.insert(std::make_pair("chr3.fa", std::make_pair(std::max(90504854 - centr_ext, 0), 93504854 + centr_ext)));
      centr.insert(std::make_pair("chr4.fa", std::make_pair(std::max(49660117 - centr_ext, 0), 52660117 + centr_ext)));
      centr.insert(std::make_pair("chr5.fa", std::make_pair(std::max(46405641 - centr_ext, 0), 49405641 + centr_ext)));
      centr.insert(std::make_pair("chr6.fa", std::make_pair(std::max(58830166 - centr_ext, 0), 61830166 + centr_ext)));
      centr.insert(std::make_pair("chr7.fa", std::make_pair(std::max(58054331 - centr_ext, 0), 61054331 + centr_ext)));
      centr.insert(std::make_pair("chr8.fa", std::make_pair(std::max(43838887 - centr_ext, 0), 46838887 + centr_ext)));
      centr.insert(std::make_pair("chr9.fa", std::make_pair(std::max(47367679 - centr_ext, 0), 50367679 + centr_ext)));
      centr.insert(std::make_pair("chr10.fa", std::make_pair(std::max(39254935 - centr_ext, 0), 42254935 + centr_ext)));
      centr.insert(std::make_pair("chr11.fa", std::make_pair(std::max(51644205 - centr_ext, 0), 54644205 + centr_ext)));
      centr.insert(std::make_pair("chr12.fa", std::make_pair(std::max(34856694 - centr_ext, 0), 37856694 + centr_ext)));
      centr.insert(std::make_pair("chr13.fa", std::make_pair(std::max(16000000 - centr_ext, 0), 19000000 + centr_ext)));
      centr.insert(std::make_pair("chr14.fa", std::make_pair(std::max(16000000 - centr_ext, 0), 19000000 + centr_ext)));
      centr.insert(std::make_pair("chr15.fa", std::make_pair(std::max(17000000 - centr_ext, 0), 20000000 + centr_ext)));
      centr.insert(std::make_pair("chr16.fa", std::make_pair(std::max(35335801 - centr_ext, 0), 38335801 + centr_ext)));
      centr.insert(std::make_pair("chr17.fa", std::make_pair(std::max(22263006 - centr_ext, 0), 25263006 + centr_ext)));
      centr.insert(std::make_pair("chr18.fa", std::make_pair(std::max(15460898 - centr_ext, 0), 18460898 + centr_ext)));
      centr.insert(std::make_pair("chr19.fa", std::make_pair(std::max(24681782 - centr_ext, 0), 27681782 + centr_ext)));
      centr.insert(std::make_pair("chr20.fa", std::make_pair(std::max(26369569 - centr_ext, 0), 29369569 + centr_ext)));
      centr.insert(std::make_pair("chr21.fa", std::make_pair(std::max(11288129 - centr_ext, 0), 14288129 + centr_ext)));
      centr.insert(std::make_pair("chr22.fa", std::make_pair(std::max(13000000 - centr_ext, 0), 16000000 + centr_ext)));
      centr.insert(std::make_pair("chrX.fa", std::make_pair(std::max(58632012 - centr_ext, 0), 61632012 + centr_ext)));
      centr.insert(std::make_pair("chrY.fa", std::make_pair(std::max(10104553 - centr_ext, 0), 13104553 + centr_ext)));
    }

      inline bool 
      inCentromere(std::string const& chr, int const pos) 
	{
	  TCentr::const_iterator find = centr.find(chr);
	  if (find != centr.end()) {
	    if ((find->second.first < pos) && (find->second.second > pos)) return true;
	  }
	  return false;
	}

  };



  class Telomere {
  public:
    int chrStart_telemore;
    int tele_ext;
    typedef std::pair<int, int> TInterval;
    typedef std::map<std::string, TInterval> TTele;
    TTele tele;
    
    Telomere() : chrStart_telemore(10000), tele_ext(0) {
      Telomere(tele_ext);
    }
      
    Telomere(unsigned int ext) : tele_ext(ext) 
      {
      chrStart_telemore += tele_ext;
      tele.insert(std::make_pair("chr1.fa", std::make_pair(std::max(249240621 - tele_ext, 0), 249250621 + tele_ext)));
      tele.insert(std::make_pair("chr2.fa", std::make_pair(std::max(243189373 - tele_ext, 0), 243199373 + tele_ext)));
      tele.insert(std::make_pair("chr3.fa", std::make_pair(std::max(198012430 - tele_ext, 0), 198022430 + tele_ext)));
      tele.insert(std::make_pair("chr4.fa", std::make_pair(std::max(191144276 - tele_ext, 0), 191154276 + tele_ext)));
      tele.insert(std::make_pair("chr5.fa", std::make_pair(std::max(180905260 - tele_ext, 0), 180915260 + tele_ext)));
      tele.insert(std::make_pair("chr6.fa", std::make_pair(std::max(171105067 - tele_ext, 0), 171115067 + tele_ext)));
      tele.insert(std::make_pair("chr7.fa", std::make_pair(std::max(159128663 - tele_ext, 0), 159138663 + tele_ext)));
      tele.insert(std::make_pair("chr8.fa", std::make_pair(std::max(146354022 - tele_ext, 0), 146364022 + tele_ext)));
      tele.insert(std::make_pair("chr9.fa", std::make_pair(std::max(141203431 - tele_ext, 0), 141213431 + tele_ext)));
      tele.insert(std::make_pair("chr10.fa", std::make_pair(std::max(135524747 - tele_ext, 0), 135534747 + tele_ext)));
      tele.insert(std::make_pair("chr11.fa", std::make_pair(std::max(134996516 - tele_ext, 0), 135006516 + tele_ext)));
      tele.insert(std::make_pair("chr12.fa", std::make_pair(std::max(133841895 - tele_ext, 0), 133851895 + tele_ext)));
      tele.insert(std::make_pair("chr13.fa", std::make_pair(std::max(115159878 - tele_ext, 0), 115169878 + tele_ext)));
      tele.insert(std::make_pair("chr14.fa", std::make_pair(std::max(107339540 - tele_ext, 0), 107349540 + tele_ext)));
      tele.insert(std::make_pair("chr15.fa", std::make_pair(std::max(102521392 - tele_ext, 0), 102531392 + tele_ext)));
      tele.insert(std::make_pair("chr16.fa", std::make_pair(std::max(90344753 - tele_ext, 0), 90354753 + tele_ext)));
      tele.insert(std::make_pair("chr17.fa", std::make_pair(std::max(81185210 - tele_ext, 0), 81195210 + tele_ext)));
      tele.insert(std::make_pair("chr18.fa", std::make_pair(std::max(78067248 - tele_ext, 0), 78077248 + tele_ext)));
      tele.insert(std::make_pair("chr19.fa", std::make_pair(std::max(59118983 - tele_ext, 0), 59128983 + tele_ext)));
      tele.insert(std::make_pair("chr20.fa", std::make_pair(std::max(63015520 - tele_ext, 0), 63025520 + tele_ext)));
      tele.insert(std::make_pair("chr21.fa", std::make_pair(std::max(48119895 - tele_ext, 0), 48129895 + tele_ext)));
      tele.insert(std::make_pair("chr22.fa", std::make_pair(std::max(51294566 - tele_ext, 0), 51304566 + tele_ext)));
      tele.insert(std::make_pair("chrX.fa", std::make_pair(std::max(155260560 - tele_ext, 0), 155270560 + tele_ext)));
      tele.insert(std::make_pair("chrY.fa", std::make_pair(std::max(59363566 - tele_ext, 0), 59373566 + tele_ext)));
    }


      inline bool 
	inTelomere(std::string const& chr, int const pos) 
	{
	  if (pos < chrStart_telemore) return true;
	  TTele::const_iterator find = tele.find(chr);
	  if (find != tele.end()) {
	    if ((find->second.first < pos) && (find->second.second > pos)) return true;
	  }
	  return false;
	}

  };
}

#endif



