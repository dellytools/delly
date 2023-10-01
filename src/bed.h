#ifndef BED_H
#define BED_H

#include <boost/filesystem.hpp>
#include <boost/multi_array.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/zlib.hpp>

#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/sam.h>

namespace torali
{

  // Flattens overlapping intervals
  template<typename TRegionsGenome>
  inline int32_t
    _parseBedIntervals(std::string const& filename, bool const filePresent, bam_hdr_t* hdr, TRegionsGenome& bedRegions) {
    typedef typename TRegionsGenome::value_type TChrIntervals;
    typedef typename TChrIntervals::interval_type TIVal;

    int32_t intervals = 0;
    if (filePresent) {
      bedRegions.resize(hdr->n_targets, TChrIntervals());
      std::ifstream chrFile;
      boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
      if (is_gz(filename)) {
	chrFile.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
	dataIn.push(boost::iostreams::gzip_decompressor(), 16*1024);
      } else chrFile.open(filename.c_str(), std::ifstream::in);
      dataIn.push(chrFile);
      std::istream instream(&dataIn);
      std::string chrFromFile;
      while(std::getline(instream, chrFromFile)) {
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t,;");
	Tokenizer tokens(chrFromFile, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter!=tokens.end()) {
	  std::string chrName = *tokIter++;
	  int32_t tid = bam_name2id(hdr, chrName.c_str());
	  if (tid >= 0) {
	    if (tokIter!=tokens.end()) {
	      int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
	      int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
	      bedRegions[tid].insert(TIVal::right_open(start, end));
	      ++intervals;
	    }
	  }
	}
      }
      dataIn.pop();
      if (is_gz(filename)) dataIn.pop();
      chrFile.close();
    }
    return intervals;
  }


  // Keeps overlapping intervals
  template<typename TRegionsGenome>
  inline int32_t
  _parsePotOverlappingIntervals(std::string const& filename, bool const filePresent, bam_hdr_t* hdr, TRegionsGenome& bedRegions) {
    typedef typename TRegionsGenome::value_type TChrIntervals;
	
    int32_t intervals = 0;
    if (filePresent) {
      bedRegions.resize(hdr->n_targets, TChrIntervals());
      std::ifstream chrFile;
      boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
      if (is_gz(filename)) {
	chrFile.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
	dataIn.push(boost::iostreams::gzip_decompressor(), 16*1024);
      } else chrFile.open(filename.c_str(), std::ifstream::in);
      dataIn.push(chrFile);
      std::istream instream(&dataIn);
      std::string chrFromFile;
      while(std::getline(instream, chrFromFile)) {      
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t,;");
	Tokenizer tokens(chrFromFile, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter!=tokens.end()) {
	  std::string chrName = *tokIter++;
	  int32_t tid = bam_name2id(hdr, chrName.c_str());
	  if (tid >= 0) {
	    if (tokIter!=tokens.end()) {
	      int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
	      int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
	      bedRegions[tid].insert(std::make_pair(start, end));
	      ++intervals;
	    }
	  }
	}
      }
      dataIn.pop();
      if (is_gz(filename)) dataIn.pop();
      chrFile.close();
    }
    return intervals;
  }


  template<typename TChrIntervals>
  inline void
  _mergeOverlappingBedEntries(TChrIntervals const& bedRegions, TChrIntervals& citv) {
    typedef boost::icl::interval_set<uint32_t> TUniqueIntervals;
    typedef typename TUniqueIntervals::interval_type TIVal;
    TUniqueIntervals uitv;

    // Insert intervals
    for(typename TChrIntervals::const_iterator it = bedRegions.begin(); it != bedRegions.end(); ++it) uitv.insert(TIVal::right_open(it->first, it->second));
    
    // Fetch unique intervals
    for(typename TUniqueIntervals::iterator it = uitv.begin(); it != uitv.end(); ++it) citv.insert(std::make_pair(it->lower(), it->upper()));
  }
}

#endif
