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

#include <boost/multi_array.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <htslib/sam.h>
#include <sstream>
#include <math.h>
#include "tags.h"


namespace torali
{

  struct LibraryInfo {
    int32_t rs;
    int32_t avgdist;
    int32_t median;
    int32_t mad;
    int32_t minNormalISize;
    int32_t minISizeCutoff;
    int32_t maxNormalISize;
    int32_t maxISizeCutoff;
    uint8_t defaultOrient;
    uint32_t abnormal_pairs;
    float mappedAsPair;

    LibraryInfo() : rs(0), avgdist(0), median(0), mad(0), minNormalISize(0), minISizeCutoff(0), maxNormalISize(0), maxISizeCutoff(0), defaultOrient(0), abnormal_pairs(0), mappedAsPair(0) {}
  };


  struct LibraryParams {
    typedef int32_t TValue;
    typedef std::vector<TValue> TSizeVector;
    
    uint32_t processedNumPairs;
    uint32_t processedNumReads;
    uint32_t orient[4];
    uint8_t defaultOrient;
    int32_t median;
    int32_t mad;
    int32_t rs;
    int32_t avgdist;
    TSizeVector vecISize;
    TSizeVector readSize;
    TSizeVector readDist;

    LibraryParams() : processedNumPairs(0), processedNumReads(0), defaultOrient(0), median(0), mad(0), rs(0), avgdist(0) {
      for(uint32_t i=0;i<4;++i) orient[i]=0;
    }
    LibraryParams(int32_t maxPSize, int32_t maxRSize) : processedNumPairs(0), processedNumReads(0), defaultOrient(0), median(0), mad(0), rs(0), avgdist(0) {
      for(uint32_t i=0;i<4;++i) orient[i]=0;
      vecISize.resize(maxPSize, 0);
      readSize.resize(maxRSize, 0);
      readDist.resize(maxRSize, 0);
    }
    
  };

  // Read count struct
  struct ReadCount {
    int32_t leftRC;
    int32_t rc;
    int32_t rightRC;

    ReadCount() {}
    ReadCount(int32_t l, int32_t m, int32_t r) : leftRC(l), rc(m), rightRC(r) {}
  };


  // Deletions
  inline std::string
    _addID(SVType<DeletionTag>) {
    return "DEL";
  }
  
  // Duplications
  inline std::string
    _addID(SVType<DuplicationTag>) {
    return "DUP";
  }
  
  // Inversions
  inline std::string
    _addID(SVType<InversionTag>) {
    return "INV";
  }
  
  // Translocations
  inline std::string
    _addID(SVType<TranslocationTag>) {
    return "TRA";
  }
  
  // Insertion
  inline std::string
    _addID(SVType<InsertionTag>) {
    return "INS";
  }


  // Decode Orientation
  inline uint8_t
    _decodeOrientation(std::string const& value) {
    if (value=="3to3") return 0;
    else if (value=="5to5") return 1;
    else if (value=="3to5") return 2;
    else if (value=="5to3") return 3;
    else return 4;
  }
  
  // Add Orientation
  inline std::string
    _addOrientation(uint8_t const ct) {
    if (ct==0) return "3to3";
    else if (ct==1) return "5to5";
    else if (ct==2) return "3to5";
    else if (ct==3) return "5to3";
    else return "NtoN";
  }
  
  
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

  template<typename TLibraryInfo>
  inline int32_t
  getHalfFragmentLength(TLibraryInfo const& libInfo) {
    if (libInfo.median == 0) return libInfo.rs / 2; // Single-end library
    if ((libInfo.median > (libInfo.rs / 2)) && (libInfo.median < 1000)) return libInfo.median / 2;  // Paired-end library
    return 150;
  }
  
  template<typename TLibraryInfo>
  inline bool
  fragmentMidPoint(bam1_t* rec, TLibraryInfo const& libInfo, int32_t& midPoint) {
    if (libInfo.hbin != (libInfo.median / 2)) midPoint = rec->core.pos + halfAlignmentLength(rec); // Single-end library or mate-pairs
    else {
      // Paired-end library
      if (!(rec->core.flag & BAM_FPAIRED)) return false;
      if (rec->core.flag & BAM_FMUNMAP) return false;
      if (rec->core.pos < rec->core.mpos) return false;
      int32_t outerISize = rec->core.pos - rec->core.mpos + rec->core.l_qseq;
      if ((rec->core.tid!=rec->core.mtid) || (getStrandIndependentOrientation(rec->core) != libInfo.defaultOrient) || (outerISize < libInfo.minNormalISize) || (outerISize > libInfo.maxNormalISize)) return false;
      midPoint = rec->core.pos + outerISize / 2;
    }
    return true;
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

  inline void
  reverseComplement(std::string& sequence) 
  {
    std::string rev = boost::to_upper_copy(std::string(sequence.rbegin(), sequence.rend()));
    std::size_t i = 0;
    for(std::string::iterator revIt = rev.begin(); revIt != rev.end(); ++revIt, ++i) {
      switch (*revIt) {
      case 'A': sequence[i]='T'; break;
      case 'C': sequence[i]='G'; break;
      case 'G': sequence[i]='C'; break;
      case 'T': sequence[i]='A'; break;
      case 'N': sequence[i]='N'; break;
      default: break;
      }
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
      for(unsigned int i=0; i<id.size();++i) {
	if (((int) id[i] < 48) || ((int) id[i] > 57)) {
	  id=id.substr(0,i);
	  break;
	}
      }
      return boost::lexical_cast<unsigned int>(id);
    }


  inline std::string
  compressStr(std::string const& data) {
    std::stringstream compressed;
    std::stringstream origin(data);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> out;
    out.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_speed)));
    out.push(origin);
    boost::iostreams::copy(out, compressed);
    return compressed.str();
  }

  inline std::string
  decompressStr(std::string const& data) {
    std::stringstream compressed(data);
    std::stringstream decompressed;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> out;
    out.push(boost::iostreams::gzip_decompressor());
    out.push(compressed);
    boost::iostreams::copy(out, decompressed);
    return decompressed.str();
  }

  inline double
  entropy(std::string const& st) {
    typedef double TPrecision;
    std::vector<char> stvec(st.begin(), st.end());
    std::set<char> alphabet(stvec.begin(), stvec.end());
    TPrecision ent = 0;
    TPrecision freq = 0;
    for(std::set<char>::const_iterator c = alphabet.begin(); c != alphabet.end(); ++c) {
      int ctr = 0;
      for (std::vector<char>::const_iterator s = stvec.begin(); s != stvec.end(); ++s)
	if (*s == *c) ++ctr;
      freq = (TPrecision) ctr / (TPrecision) stvec.size();
      ent += (freq) * log(freq)/log(2);
    }
    return -ent;
  }


  inline bool
  getSMTag(std::string const& header, std::string const& fileName, std::string& sampleName) {
    std::set<std::string> smIdentifiers;
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
	    if (field == "SM") {
	      rgPresent = true;
	      std::string rgSM = itKV->substr(sp+1);
	      smIdentifiers.insert(rgSM);
	    }
	  }
	}
      }
    }
    if (!rgPresent) {
      sampleName = fileName;
      return true;
    } else if (smIdentifiers.size() == 1) {
      sampleName = *(smIdentifiers.begin());
      return true;
    } else {
      sampleName = "";
      return false;
    }
  }


  template<typename TIterator, typename TValue>
  inline void
  getMedian(TIterator begin, TIterator end, TValue& median) 
  {
    std::nth_element(begin, begin + (end - begin) / 2, end);
    median = *(begin + (end - begin) / 2);
  }
  
  template<typename TVector, typename TPercentile, typename TValue>
  inline void
  getPercentile(TVector& vec, TPercentile p, TValue& percentile) 
  {
    std::nth_element(vec.begin(), vec.begin() + int((vec.size() * p)), vec.end());
    percentile = *(vec.begin() + int(vec.size() * p));
  }

  template<typename TConfig, typename TLibraryMap>
  inline int32_t
  getVariability(TConfig const& c, std::vector<TLibraryMap> const& sampleLib)
  {
    int32_t overallVariability = 0;
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      for(typename TLibraryMap::const_iterator libIt=sampleLib[file_c].begin();libIt!=sampleLib[file_c].end();++libIt) {
	if (libIt->second.maxNormalISize > overallVariability) overallVariability = libIt->second.maxNormalISize;
	if (libIt->second.rs > overallVariability) overallVariability = libIt->second.rs;
      }
    }
    return overallVariability;
  }

  template<typename TConfig>
  inline int32_t
  getVariability(TConfig const&, boost::unordered_map<std::string, LibraryInfo> const& sampleLib)
  {
    typedef boost::unordered_map<std::string, LibraryInfo> TLibraryMap;
    int32_t overallVariability = 0;
    for(typename TLibraryMap::const_iterator libIt=sampleLib.begin();libIt!=sampleLib.end();++libIt) {
      if (libIt->second.maxNormalISize > overallVariability) overallVariability = libIt->second.maxNormalISize;
      if (libIt->second.rs > overallVariability) overallVariability = libIt->second.rs;
    }
    return overallVariability;
  }
  
  
  template<typename TIterator, typename TValue>
  inline void
  getMAD(TIterator begin, TIterator end, TValue median, TValue& mad) 
  {
    std::vector<TValue> absDev;
    for(;begin<end;++begin) 
      absDev.push_back(std::abs((TValue)*begin - median));
    getMedian(absDev.begin(), absDev.end(), mad);
  }

  template<typename TIterator, typename TValue>
  inline void
  getMean(TIterator begin, TIterator end, TValue& mean) 
  {
    mean = 0;
    unsigned int count = 0;
    for(; begin<end; ++begin,++count) mean += *begin;
    mean /= count;
  }

  template<typename TIterator, typename TValue>
  inline void
  getStdDev(TIterator begin, TIterator end, TValue mean, TValue& stdDev) 
  {
    stdDev = 0;
    unsigned int count = 0;
    for(;begin<end;++begin,++count) stdDev += ((TValue)*begin - mean) * ((TValue)*begin - mean);
    stdDev = sqrt(stdDev / (TValue) count);
  }

  template<typename TConfig, typename TValidRegion, typename TSampleLibrary>
    inline void getLibraryParams(TConfig const& c, TValidRegion const& validRegions, TSampleLibrary& sampleLib) {
    typedef typename TValidRegion::value_type TChrIntervals;
    typedef typename TSampleLibrary::value_type TLibraryMap;

    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    typedef std::vector<bam_hdr_t*> TSamHeader;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    TSamHeader hdr(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
      hdr[file_c] = sam_hdr_read(samfile[file_c]);
      sampleLib[file_c] = TLibraryMap();
    }

#pragma omp parallel for default(shared)
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      // Minimum and maximum number of pairs used to estimate library parameters for each RG library
      unsigned int maxNumAlignments=1000000;
      unsigned int minNumAlignments=10000;

      // Store the counts in an object for each RG library
      typedef boost::unordered_map<std::string, LibraryParams> TParams;
      TParams params;

      // Get read groups
      std::string header(hdr[file_c]->text);
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
		params.insert(std::make_pair(rgID, LibraryParams(maxNumAlignments, minNumAlignments)));
	      }
	    }
	  }
	}
      }
      if (!rgPresent) params.insert(std::make_pair("DefaultLib", LibraryParams(maxNumAlignments, minNumAlignments)));

      // Collect insert sizes
      uint32_t libCount = 0;
      uint32_t allLibs = params.size();
      for(int32_t refIndex=0; ((refIndex < (int32_t) hdr[0]->n_targets) && (libCount < allLibs)); ++refIndex) {
	if (validRegions[refIndex].empty()) continue;
	for(typename TChrIntervals::const_iterator vRIt = validRegions[refIndex].begin(); ((vRIt != validRegions[refIndex].end()) && (libCount < allLibs)); ++vRIt) {
	  int32_t lastReadPos = vRIt->lower();
	  hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, vRIt->lower(), vRIt->upper());
	  bam1_t* rec = bam_init1();
	  while ((sam_itr_next(samfile[file_c], iter, rec) >= 0) && (libCount < allLibs)) {
	    if (!(rec->core.flag & BAM_FREAD2) && (rec->core.l_qseq < 65000)) {
	      if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	      std::string rG = "DefaultLib";
	      uint8_t *rgptr = bam_aux_get(rec, "RG");
	      if (rgptr) {
		char* rg = (char*) (rgptr + 1);
		rG = std::string(rg);
	      }
	      TParams::iterator paramIt = params.find(rG);
	      paramIt->second.readSize[paramIt->second.processedNumReads % minNumAlignments] = rec->core.l_qseq;
	      paramIt->second.readDist[paramIt->second.processedNumReads % minNumAlignments] = rec->core.pos - lastReadPos;
	      lastReadPos = rec->core.pos;
	      ++paramIt->second.processedNumReads;
	      if ((paramIt->second.processedNumReads >= minNumAlignments) && (paramIt->second.processedNumPairs == 0)) {
		// Single-end library
		if (paramIt->second.processedNumReads == minNumAlignments) ++libCount;
		continue;
	      }
	      if ((rec->core.flag & BAM_FPAIRED) && !(rec->core.flag & BAM_FMUNMAP) && (rec->core.tid==rec->core.mtid)) {
		paramIt->second.vecISize[paramIt->second.processedNumPairs % maxNumAlignments] = abs(rec->core.isize);
		++paramIt->second.orient[getStrandIndependentOrientation(rec->core)];
		++paramIt->second.processedNumPairs;
		if (paramIt->second.processedNumPairs == maxNumAlignments) ++libCount;
	      }
	    }
	  }
	  bam_destroy1(rec);
	  hts_itr_destroy(iter);
	}
      }

      // Get library parameters
      for(TParams::iterator paramIt=params.begin(); paramIt != params.end(); ++paramIt) {
	// Check single-end library parameters
	if (paramIt->second.processedNumReads > 0) {
	  if (paramIt->second.processedNumReads < minNumAlignments) {
	    paramIt->second.readSize.resize(paramIt->second.processedNumReads);
	    paramIt->second.readDist.resize(paramIt->second.processedNumReads);
	  }
	  std::sort(paramIt->second.readSize.begin(), paramIt->second.readSize.end());
	  paramIt->second.rs = paramIt->second.readSize[paramIt->second.readSize.size() / 2];
	  std::sort(paramIt->second.readDist.begin(), paramIt->second.readDist.end());
	  paramIt->second.avgdist = paramIt->second.readDist[paramIt->second.readDist.size() / 2];
	}
	
	// Check that this is a proper paired-end library
	if (paramIt->second.processedNumPairs>=minNumAlignments) {
	  if (paramIt->second.processedNumPairs < maxNumAlignments) paramIt->second.vecISize.resize(paramIt->second.processedNumPairs);

	  // Get default library orientation
	  paramIt->second.defaultOrient=0;
	  unsigned int maxOrient=paramIt->second.orient[0];
	  for(unsigned int i=1;i<4;++i) {
	    if (paramIt->second.orient[i]>maxOrient) {
	      maxOrient=paramIt->second.orient[i];
	      paramIt->second.defaultOrient=(uint8_t) i;
	    }
	  }
	  
	  typedef typename LibraryParams::TSizeVector TVecISize;
	  std::sort(paramIt->second.vecISize.begin(), paramIt->second.vecISize.end());

	  // Mate-pair library (If yes, trim off the chimera peak < 1000bp)
	  if ((paramIt->second.defaultOrient==3) && (paramIt->second.vecISize[paramIt->second.vecISize.size() / 2] >= 1000)) {
	    TVecISize vecISizeTmp;
	    for(typename TVecISize::const_iterator iSizeBeg = paramIt->second.vecISize.begin(); iSizeBeg<paramIt->second.vecISize.end();++iSizeBeg)
	      if (*iSizeBeg >= 1000) vecISizeTmp.push_back(*iSizeBeg);
	    paramIt->second.vecISize = vecISizeTmp;
	  }
	  paramIt->second.median = paramIt->second.vecISize[paramIt->second.vecISize.size() / 2];
	  std::vector<int32_t> absDev;
	  for(typename TVecISize::const_iterator iSizeBeg = paramIt->second.vecISize.begin(); iSizeBeg != paramIt->second.vecISize.end(); ++iSizeBeg) absDev.push_back(std::abs(*iSizeBeg - paramIt->second.median));
	  std::sort(absDev.begin(), absDev.end());
	  paramIt->second.mad = absDev[absDev.size() / 2];
	}
      }

#pragma omp critical
      {
	for(TParams::iterator paramIt=params.begin(); paramIt != params.end(); ++paramIt) {
	  typename TLibraryMap::iterator libInfoIt = sampleLib[file_c].insert(std::make_pair(paramIt->first, LibraryInfo())).first;
	  if (paramIt->second.processedNumReads > 0) {
	    libInfoIt->second.rs = paramIt->second.rs;
	    libInfoIt->second.avgdist = paramIt->second.avgdist;
	    libInfoIt->second.mappedAsPair = (float) paramIt->second.processedNumPairs / (float) paramIt->second.processedNumReads;
	  }
	  if ((paramIt->second.median >= 50) && (paramIt->second.median<=100000)) {
	    libInfoIt->second.defaultOrient = paramIt->second.defaultOrient;
	    libInfoIt->second.median = paramIt->second.median;
	    libInfoIt->second.mad = paramIt->second.mad;
	    libInfoIt->second.maxNormalISize = libInfoIt->second.median + (5 * libInfoIt->second.mad);
	    libInfoIt->second.minNormalISize = libInfoIt->second.median - (5 * libInfoIt->second.mad);
	    if (libInfoIt->second.minNormalISize < 0) libInfoIt->second.minNormalISize=0;
	    libInfoIt->second.maxISizeCutoff = libInfoIt->second.median + (c.madCutoff * libInfoIt->second.mad);
	    libInfoIt->second.minISizeCutoff = libInfoIt->second.median - (c.madCutoff * libInfoIt->second.mad);
	    if (libInfoIt->second.minISizeCutoff < 0) libInfoIt->second.minISizeCutoff=0;
	  }
	}
      }
    }

    // Clean-up
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      bam_hdr_destroy(hdr[file_c]);
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }



}

#endif
