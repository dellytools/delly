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

#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <htslib/sam.h>
#include <math.h>
#include "tags.h"


namespace torali
{

  struct LibraryInfo {
    int median;
    int mad;
    int percentileCutoff;
    int minNormalISize;
    int minISizeCutoff;
    int maxNormalISize;
    int maxISizeCutoff;
    int defaultOrient;
    unsigned int non_unique_abnormal_pairs;
    unsigned int abnormal_pairs;
  };


  struct _LibraryParams {
    unsigned int processedNumPairs;
    unsigned int orient[4];
    int defaultOrient;
    double median;
    double mad;
    double percentileCutoff;
    std::vector<int32_t> vecISize;
  };

  
  inline unsigned int alignmentLength(bam1_t* rec) {
    uint32_t* cigar = bam_get_cigar(rec);
    unsigned int alen = 0;
    for (unsigned int i = 0; i < rec->core.n_cigar; ++i)
      if (bam_cigar_op(cigar[i]) == BAM_CMATCH) alen+=bam_cigar_oplen(cigar[i]);
    return alen;
  }

  inline unsigned int halfAlignmentLength(bam1_t* rec) {
    return (alignmentLength(rec) / 2);
  }
  
  inline std::size_t hash_pair(bam1_t* rec) {
    std::size_t seed = 0;
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    boost::hash_combine(seed, rec->core.mtid);
    boost::hash_combine(seed, rec->core.mpos);
    return seed;
  }

  inline std::size_t hash_pair_mate(bam1_t* rec) {
    std::size_t seed = 0;
    boost::hash_combine(seed, rec->core.mtid);
    boost::hash_combine(seed, rec->core.mpos);
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    return seed;
  }

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


  template<typename TSVId>
    inline unsigned int
    parseSVid(TSVId id) 
    {
      for(unsigned int i=3; i<id.size();++i) {
	if (id[i]!='0') {
	  id=id.substr(i);
	  break;
	}
      }
      return boost::lexical_cast<unsigned int>(id);
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


  template<typename TFiles, typename TSampleLibrary>
    inline void getLibraryParams(TFiles const& files, TSampleLibrary& sampleLib, double const percentile, unsigned short const madCutoff) {
    typedef typename TSampleLibrary::mapped_type TLibraryMap;

    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    TSamFile samfile;
    samfile.resize(files.size());
    for(unsigned int file_c = 0; file_c < files.size(); ++file_c) {
      samfile[file_c] = sam_open(files[file_c].string().c_str(), "r");
      std::string sampleName(files[file_c].stem().string());
      sampleLib.insert(std::make_pair(sampleName, TLibraryMap()));
    }

#pragma omp parallel for default(shared)
    for(unsigned int file_c = 0; file_c < files.size(); ++file_c) {
      // Minimum and maximum number of pairs used to estimate library parameters for each RG library
      unsigned int maxNumPairs=5000000;
      unsigned int minNumPairs=1000;

      // Store the counts in an object for each RG librar
      typedef boost::unordered_map<std::string, _LibraryParams> TParams;
      TParams params;

      // Create SAM Object
      bam_hdr_t* hdr = sam_hdr_read(samfile[file_c]);

      // Get read groups
      std::string header(hdr->text);
      std::string delimiters("\n");
      typedef std::vector<std::string> TStrParts;
      TStrParts lines;
      boost::split(lines, header, boost::is_any_of(delimiters));
      TStrParts::const_iterator itH = lines.begin();
      TStrParts::const_iterator itHEnd = lines.end();
      bool rgPresent = false;
      for(;itH!=itHEnd; ++itH) {
	if (itH->find("@RG")==0) {
	  std::string delim("\t ");
	  TStrParts keyval;
	  boost::split(keyval, *itH, boost::is_any_of(delim));
	  TStrParts::const_iterator itKV = keyval.begin();
	  TStrParts::const_iterator itKVEnd = keyval.end();
	  for(;itKV != itKVEnd; ++itKV) {
	    size_t sp = itKV->find(":");
	    if (sp != std::string::npos) {
	      std::string field = itKV->substr(0, sp);
	      if (field == "ID") {
		rgPresent = true;
		std::string rgID = itKV->substr(sp+1);
		params.insert(std::make_pair(rgID, _LibraryParams()));
	      }
	    }
	  }
	}
      }
      if (!rgPresent) params.insert(std::make_pair("DefaultLib", _LibraryParams()));

      // Initialize arrays
      for(TParams::iterator paramIt = params.begin(); paramIt!=params.end(); ++paramIt) {
	paramIt->second.processedNumPairs=0;
	for(unsigned int i=0;i<4;++i) paramIt->second.orient[i]=0;
	paramIt->second.vecISize.clear();
      }
      
      // Collect insert sizes
      bool missingPairs=true;
      bam1_t* rec = bam_init1();
      while ((sam_read1(samfile[file_c], hdr, rec) >=0) && (missingPairs)) {
	if ((rec->core.flag & BAM_FPAIRED) && !(rec->core.flag & BAM_FUNMAP) && !(rec->core.flag & BAM_FMUNMAP) && (rec->core.flag & BAM_FREAD1) && (rec->core.tid==rec->core.mtid) && !(rec->core.flag & BAM_FSECONDARY) && !(rec->core.flag & BAM_FQCFAIL)  && !(rec->core.flag & BAM_FDUP) && !(rec->core.flag & BAM_FSUPPLEMENTARY)) {
	  std::string rG = "DefaultLib";
	  uint8_t *rgptr = bam_aux_get(rec, "RG");
	  if (rgptr) {
	    char* rg = (char*) (rgptr + 1);
	    rG = std::string(rg);
	  }
	  TParams::iterator paramIt= params.find(rG);
	  ++paramIt->second.processedNumPairs;
	  if (paramIt->second.processedNumPairs>=maxNumPairs) {
	    // Check every now and then if we have a good estimate of the library complexity
	    if (paramIt->second.processedNumPairs % minNumPairs == 0) {
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
	  paramIt->second.vecISize.push_back(abs(rec->core.isize));
	  ++paramIt->second.orient[getStrandIndependentOrientation(rec->core)];
	}
      }
      
      // Clean-up
      bam_destroy1(rec);

      // Get default library orientation
      for(TParams::iterator paramIt=params.begin(); paramIt != params.end(); ++paramIt) {
	paramIt->second.defaultOrient=0;
	unsigned int maxOrient=paramIt->second.orient[0];
	for(unsigned int i=1;i<4;++i) {
	  if (paramIt->second.orient[i]>maxOrient) {
	    maxOrient=paramIt->second.orient[i];
	    paramIt->second.defaultOrient=i;
	  }
	}
	
	// Mate-pair library (If yes, trim off the chimera peak < 1000bp)
	if (paramIt->second.defaultOrient==3) {
	  typedef std::vector<int32_t> TVecISize;
	  TVecISize vecISizeTmp;
	  typename TVecISize::const_iterator iSizeBeg = paramIt->second.vecISize.begin();
	  typename TVecISize::const_iterator iSizeEnd = paramIt->second.vecISize.end();
	  for(;iSizeBeg<iSizeEnd;++iSizeBeg)
	    if (*iSizeBeg >= 1000) vecISizeTmp.push_back(*iSizeBeg);
	  paramIt->second.vecISize = vecISizeTmp;
	}
	
	// Check that this is a proper paired-end library
	if (paramIt->second.vecISize.size()>=minNumPairs) {
	  // Get library stats
	  getLibraryStats(paramIt->second.vecISize.begin(), paramIt->second.vecISize.end(), percentile, paramIt->second.median, paramIt->second.mad, paramIt->second.percentileCutoff);
	}
      }
      bam_hdr_destroy(hdr);

#pragma omp critical
      {
	std::string sampleName(files[file_c].stem().string());
	typename TSampleLibrary::iterator sampleIt=sampleLib.find(sampleName);
	for(TParams::iterator paramIt=params.begin(); paramIt != params.end(); ++paramIt) {
	  typename TLibraryMap::iterator libInfoIt = sampleIt->second.insert(std::make_pair(paramIt->first, LibraryInfo())).first;
	  if ((paramIt->second.median >= 50) && (paramIt->second.median<=100000)) {
	    libInfoIt->second.defaultOrient = paramIt->second.defaultOrient;
	    libInfoIt->second.median = (int) paramIt->second.median;
	    libInfoIt->second.mad = (int) paramIt->second.mad;
	    libInfoIt->second.percentileCutoff = (int) paramIt->second.percentileCutoff;
	    libInfoIt->second.maxNormalISize = libInfoIt->second.median + (5 * libInfoIt->second.mad);
	    libInfoIt->second.minNormalISize = libInfoIt->second.median - (5 * libInfoIt->second.mad);
	    if (libInfoIt->second.minNormalISize < 0) libInfoIt->second.minNormalISize=0;
	    if (percentile!=0) libInfoIt->second.maxISizeCutoff = libInfoIt->second.percentileCutoff;
	    else libInfoIt->second.maxISizeCutoff = libInfoIt->second.median + (madCutoff * libInfoIt->second.mad);
	    libInfoIt->second.minISizeCutoff = libInfoIt->second.median - (madCutoff * libInfoIt->second.mad);
	    if (libInfoIt->second.minISizeCutoff < 0) libInfoIt->second.minISizeCutoff=0;
	  }
	}
      }
    }

    // Clean-up
    for(unsigned int file_c = 0; file_c < files.size(); ++file_c) sam_close(samfile[file_c]);
  }



}

#endif
