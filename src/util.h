/*
============================================================================
DELLY: Structural variant discovery by integrated PE mapping and SR analysis
============================================================================
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
    int32_t median;
    int32_t mad;
    int32_t minNormalISize;
    int32_t minISizeCutoff;
    int32_t maxNormalISize;
    int32_t maxISizeCutoff;
    uint32_t abnormal_pairs;

    LibraryInfo() : rs(0), median(0), mad(0), minNormalISize(0), minISizeCutoff(0), maxNormalISize(0), maxISizeCutoff(0), abnormal_pairs(0) {}
  };


  // Read count struct
  struct ReadCount {
    int32_t leftRC;
    int32_t rc;
    int32_t rightRC;

    ReadCount() {}
    ReadCount(int32_t l, int32_t m, int32_t r) : leftRC(l), rc(m), rightRC(r) {}
  };


  // Decode Orientation
  inline int32_t
    _decodeOrientation(std::string const& value) {
    if (value=="3to3") return 0;
    else if (value=="5to5") return 1;
    else if (value=="3to5") return 2;
    else if (value=="5to3") return 3;
    else return 4;
  }
  
  // Decode Orientation
  inline int32_t
  _decodeOrientation(std::string const& value, std::string const& svt) {
    if (svt == "BND") {
      if (value=="3to3") return DELLY_SVT_TRANS + 0;
      else if (value=="5to5") return DELLY_SVT_TRANS + 1;
      else if (value=="3to5") return DELLY_SVT_TRANS + 2;
      else if (value=="5to3") return DELLY_SVT_TRANS + 3;
      else return -1;
    } else {
      if (value=="3to3") return 0;
      else if (value=="5to5") return 1;
      else if (value=="3to5") return 2;
      else if (value=="5to3") return 3;
      else return 4;
    }
  }
  
  
  // Deletions
  inline std::string
  _addID(int32_t const svt) {
    if (svt == 0) return "INV";
    else if (svt == 1) return "INV";
    else if (svt == 2) return "DEL";
    else if (svt == 3) return "DUP";
    else if (svt == 4) return "INS";
    else return "BND";
  }

  inline std::string
  _addAlleles(std::string const& ref, std::string const& chr2, StructuralVariantRecord const& sv, int32_t const svt) {
    if (_translocation(svt)) {
      uint8_t ct = _getSpanOrientation(svt);
      if (ct == 0) {
	return ref + "," + ref + "]" + chr2 + ":" + boost::lexical_cast<std::string>(sv.svEnd) + "]";
      } else if (ct == 1) {
	return ref + "," + "[" + chr2 + ":" + boost::lexical_cast<std::string>(sv.svEnd) + "[" + ref;
      } else if (ct == 2) {
	return ref + "," + ref + "[" + chr2 + ":" + boost::lexical_cast<std::string>(sv.svEnd) + "[";
      } else if (ct == 3) {
	return ref + "," + "]" + chr2 + ":" + boost::lexical_cast<std::string>(sv.svEnd) + "]" + ref;
      } else {
	return ref + ",<" + _addID(svt) + ">";
      }
    } else return ref + ",<" + _addID(svt) + ">";
  }

  
  // Add Orientation
  inline std::string
  _addOrientation(int32_t const svt) {
    uint8_t ct = _getSpanOrientation(svt);
    if (ct==0) return "3to3";
    else if (ct==1) return "5to5";
    else if (ct==2) return "3to5";
    else if (ct==3) return "5to3";
    else return "NtoN";
  }

  inline uint32_t sequenceLength(bam1_t const* rec) {
    uint32_t* cigar = bam_get_cigar(rec);
    uint32_t slen = 0;
    for (uint32_t i = 0; i < rec->core.n_cigar; ++i)
      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF) || (bam_cigar_op(cigar[i]) == BAM_CINS) || (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) slen += bam_cigar_oplen(cigar[i]);
    return slen;
  }
  
  inline uint32_t alignmentLength(bam1_t* rec) {
    uint32_t* cigar = bam_get_cigar(rec);
    uint32_t alen = 0;
    for (uint32_t i = 0; i < rec->core.n_cigar; ++i)
      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF) || (bam_cigar_op(cigar[i]) == BAM_CDEL) || (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP)) alen += bam_cigar_oplen(cigar[i]);
    return alen;
  }

  inline uint32_t halfAlignmentLength(bam1_t* rec) {
    return (alignmentLength(rec) / 2);
  }

  inline std::size_t hash_lr(bam1_t* rec) {
    boost::hash<std::string> string_hash;
    std::string qname = bam_get_qname(rec);
    std::size_t seed = hash_string(qname.c_str());
    boost::hash_combine(seed, string_hash(qname));
    return seed;
  }
      
  
  inline std::size_t hash_pair(bam1_t* rec) {
    std::size_t seed = hash_string(bam_get_qname(rec));
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    boost::hash_combine(seed, rec->core.mtid);
    boost::hash_combine(seed, rec->core.mpos);
    return seed;
  }

  inline std::size_t hash_pair_mate(bam1_t* rec) {
    std::size_t seed = hash_string(bam_get_qname(rec));
    boost::hash_combine(seed, rec->core.mtid);
    boost::hash_combine(seed, rec->core.mpos);
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    return seed;
  }

  inline void
  reverseComplement(std::string& sequence) {
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
    for(std::set<char>::const_iterator c = alphabet.begin(); c != alphabet.end(); ++c) {
      int ctr = 0;
      for (std::vector<char>::const_iterator s = stvec.begin(); s != stvec.end(); ++s)
	if (*s == *c) ++ctr;
      TPrecision freq = (TPrecision) ctr / (TPrecision) stvec.size();
      ent += (freq) * log(freq)/log(2);
    }
    return -ent;
  }


  inline void
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
	std::string delim("\t");
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
    } else if (smIdentifiers.size() == 1) {
      sampleName = *(smIdentifiers.begin());
    } else if (smIdentifiers.size() > 1) {
      sampleName = *(smIdentifiers.begin());
      std::cerr << "Warning: Multiple sample names (@RG:SM) present in the BAM file!" << std::endl;
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


  template<typename TConfig>
  inline int32_t
  getVariability(TConfig const&, std::vector<LibraryInfo> const& lib) {
    int32_t overallVariability = 0;
    for(uint32_t libIdx = 0; libIdx < lib.size(); ++libIdx) {
      if (lib[libIdx].maxNormalISize > overallVariability) overallVariability = lib[libIdx].maxNormalISize;
      if (lib[libIdx].rs > overallVariability) overallVariability = lib[libIdx].rs;
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
  inline void
  getLibraryParams(TConfig const& c, TValidRegion const& validRegions, TSampleLibrary& sampleLib) {
    typedef typename TValidRegion::value_type TChrIntervals;

    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    typedef std::vector<bam_hdr_t*> TSamHeader;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    TSamHeader hdr(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
      hdr[file_c] = sam_hdr_read(samfile[file_c]);
    }

    // Iterate all samples
    for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
      uint32_t maxAlignmentsScreened=10000000;
      uint32_t maxNumAlignments=1000000;
      uint32_t minNumAlignments=1000;
      uint32_t alignmentCount=0;
      uint32_t processedNumPairs = 0;
      uint32_t processedNumReads = 0;
      uint32_t rplus = 0;
      uint32_t nonrplus = 0;
      typedef std::vector<uint32_t> TSizeVector;
      TSizeVector vecISize;
      TSizeVector readSize;

      // Collect insert sizes
      bool libCharacterized = false;
      for(uint32_t refIndex=0; refIndex < (uint32_t) hdr[0]->n_targets; ++refIndex) {
	if (validRegions[refIndex].empty()) continue;
	for(typename TChrIntervals::const_iterator vRIt = validRegions[refIndex].begin(); ((vRIt != validRegions[refIndex].end()) && (!libCharacterized)); ++vRIt) {
	  hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, vRIt->lower(), vRIt->upper());
	  bam1_t* rec = bam_init1();
	  while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	    if (!(rec->core.flag & BAM_FREAD2) && (rec->core.l_qseq < 65000)) {
	      if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	      if ((alignmentCount > maxAlignmentsScreened) || ((processedNumReads >= maxNumAlignments) && (processedNumPairs == 0)) || (processedNumPairs >= maxNumAlignments)) {
		  // Paired-end library with enough pairs
		  libCharacterized = true;
		  break;
	      }
	      ++alignmentCount;
	      
	      // Single-end
	      if (processedNumReads < maxNumAlignments) {
		readSize.push_back(rec->core.l_qseq);
		++processedNumReads;
	      }
	      
	      // Paired-end
	      if ((rec->core.flag & BAM_FPAIRED) && !(rec->core.flag & BAM_FMUNMAP) && (rec->core.tid==rec->core.mtid)) {
		if (processedNumPairs < maxNumAlignments) {
		  vecISize.push_back(abs(rec->core.isize));
		  if (getSVType(rec->core) == 2) ++rplus;
		  else ++nonrplus;
		  ++processedNumPairs;
		}
	      }
	    }
	  }
	  bam_destroy1(rec);
	  hts_itr_destroy(iter);
	  if (libCharacterized) break;
	}
	if (libCharacterized) break;
      }
    
      // Get library parameters
      if (processedNumReads >= minNumAlignments) {
	std::sort(readSize.begin(), readSize.end());
	sampleLib[file_c].rs = readSize[readSize.size() / 2];
      }
      if (processedNumPairs >= minNumAlignments) {
	std::sort(vecISize.begin(), vecISize.end());
	int32_t median = vecISize[vecISize.size() / 2];
	std::vector<uint32_t> absDev;
	for(uint32_t i = 0; i < vecISize.size(); ++i) absDev.push_back(std::abs((int32_t) vecISize[i] - median));
	std::sort(absDev.begin(), absDev.end());
	int32_t mad = absDev[absDev.size() / 2];

	// Get default library orientation
	if ((median >= 50) && (median<=100000)) {
	  if (rplus < nonrplus) {
	    std::cerr << "Warning: Sample has a non-default paired-end layout! File: " << c.files[file_c].string() << std::endl;
	    std::cerr << "The expected paired-end orientation is   ---Read1--->      <---Read2---  which is the default illumina paired-end layout." << std::endl;
	    
	  } else {
	    sampleLib[file_c].median = median;
	    sampleLib[file_c].mad = mad;
	    sampleLib[file_c].maxNormalISize = median + (5 * mad);
	    sampleLib[file_c].minNormalISize = median - (5 * mad);
	    if (sampleLib[file_c].minNormalISize < 0) sampleLib[file_c].minNormalISize=0;
	    sampleLib[file_c].maxISizeCutoff = median + (c.madCutoff * mad);
	    sampleLib[file_c].minISizeCutoff = median - (c.madCutoff * mad);

	    // Deletion insert-size sanity checks
	    sampleLib[file_c].maxISizeCutoff = std::max(sampleLib[file_c].maxISizeCutoff, 2*sampleLib[file_c].rs);
	    sampleLib[file_c].maxISizeCutoff = std::max(sampleLib[file_c].maxISizeCutoff, 500);

	    if (sampleLib[file_c].minISizeCutoff < 0) sampleLib[file_c].minISizeCutoff=0;
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
