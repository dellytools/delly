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

#include <iostream>
#include <fstream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "api/BamMultiReader.h"
#include "api/BamReader.h"

#include "memory_mapped_file.h"
#include "bam_file_adaptor.h"
#include "version.h"
#include "util.h"
#include "sam.h"
#include "intervaltree.h"
#include "fasta_reader.h"


using namespace torali;

#define MAX_CHROM_SIZE 250000000


struct Config {
  bool mapq;
  unsigned short minMapQual;
  int bpWindowOffset;
  boost::filesystem::path int_file;
  boost::filesystem::path outfile;
  std::vector<boost::filesystem::path> files;
};

template<typename TChr, typename TPos, typename TQual>
struct Hit {
  TChr chr;
  TPos start;
  TPos end;
  TQual qual;
  
  Hit() {}

  Hit(BamTools::BamAlignment const& al) : chr(al.RefID), start(al.Position+1), end(al.MatePosition+1), qual(al.MapQuality) {}
};


template<typename TChr, typename TPos>
struct Hit<TChr, TPos, void> {
  TChr chr;
  TPos start;
  TPos end;
  
  Hit() {}

  Hit(BamTools::BamAlignment const& al) : chr(al.RefID), start(al.Position+1), end(al.MatePosition+1) {}
};


template<typename THit>
struct SortHits : public std::binary_function<THit, THit, bool>
{
  inline bool operator()(THit const& hit1, THit const& hit2) {
    if (hit1.chr != hit2.chr) return (hit1.chr < hit2.chr);
    if (hit1.start == hit2.start) return (hit1.end < hit2.end);
    else return (hit1.start < hit2.start);
  }
};


template<typename TChr, typename TPos, typename TChar, typename TIterator, typename TArrayType>
inline void
_addReadAndBpCounts(std::vector<Hit<TChr, TPos, TChar> > const& hit_vector, TIterator const chrSEit, TArrayType* bp_count)
{
  typedef std::vector<Hit<TChr, TPos, TChar> > THits;
  typename THits::const_iterator vecBeg = hit_vector.begin() + chrSEit->second.first;
  typename THits::const_iterator vecEnd = hit_vector.begin() + chrSEit->second.second;
  
  // Add bp counts
  for(;vecBeg!=vecEnd; ++vecBeg) {
    TArrayType* bpPoint = &bp_count[vecBeg->start - 1];
    TArrayType* bpPointEnd = &bp_count[vecBeg->end];
    for(;bpPoint!=bpPointEnd; ++bpPoint) ++(*bpPoint);
  }
}


template<typename TChr, typename TIterator, typename TPos, typename TString>
inline void
  _buildMAPQString(std::vector< Hit<TChr, TPos, void> >&, TIterator const, TPos const, TPos const, std::vector<TString>&)
{
  // Nothing to do
}

template<typename TChr, typename TQual, typename TIterator, typename TPos, typename TString>
inline void
  _buildMAPQString(std::vector< Hit<TChr, TPos, TQual> >& hit_vector, TIterator const chrSEit, TPos const posStart, TPos const posEnd, std::vector<TString>& str)
{
  typedef Hit<TChr, TPos, TQual> THit;
  typedef std::vector<THit> THits;
  typename THits::const_iterator vecBeg = hit_vector.begin() + chrSEit->second.first;
  typename THits::const_iterator vecEnd = hit_vector.begin() + chrSEit->second.second;

  // Add mapq counts
  int searchRange = posStart - 10000;
  if (searchRange < 0) searchRange=0;
  THit hit;
  hit.chr=vecBeg->chr;
  hit.start=searchRange;
  hit.end=searchRange;
  typename THits::const_iterator vecIt = std::lower_bound(vecBeg, vecEnd, hit, SortHits<THit>());
  for(;vecIt!=vecEnd; ++vecIt) {
    if (vecIt->end < posStart) continue;
    if (vecIt->start > posEnd) break;
    for (int i = (vecIt->start - 1); i<vecIt->end; ++i) {
      if (i >= posStart && i<posEnd) {
	std::stringstream s;
	s << (unsigned short) vecIt->qual << ',';
	str[i-posStart].append(s.str());
      }
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

template<typename TChr, typename TInterval, typename TAlignRecord, typename THit>
inline void
_insertInterval(std::map<TChr, IntervalTree<TInterval>* > const& chrIntervals, TAlignRecord const& al, std::vector<THit>& hit_vector) {
  typedef IntervalTree<TInterval> TIntervalTree;
  typedef std::map<TChr, TIntervalTree*> TChrIntervalTrees;
  typename TChrIntervalTrees::const_iterator findTree = chrIntervals.find( al.RefID );
  if (findTree != chrIntervals.end()) {
    TIntervalTree* iTree = findTree->second;
    typedef std::vector<TInterval> TResultVec;
    TResultVec results;
    TInterval searchInt(al.Position, al.MatePosition);
    iTree->enumOverlapInterval(searchInt, results);
    if (!results.empty()) {
      hit_vector.push_back(THit(al));
    }
  }
}

template<typename TChr, typename TPos, typename TChar, typename TChrSE>
inline void 
_dissectHitVector(std::vector<Hit<TChr, TPos, TChar> > const& hit_vector, TChrSE& chrSE) {
  typedef std::vector<Hit<TChr, TPos, TChar> > THits;
  typename THits::const_iterator vecBeg = hit_vector.begin();
  typename THits::const_iterator vecEnd = hit_vector.end();
  if (vecBeg != vecEnd) {
    TChr curChr = vecBeg->chr;
    TPos posStart = 0;
    TPos posEnd = 0;
    for(;vecBeg!=vecEnd; ++vecBeg, ++posEnd) {
      if (curChr != vecBeg->chr) {
	chrSE.insert(std::make_pair(curChr, std::make_pair(posStart, posEnd)));
	curChr = vecBeg->chr;
	posStart=posEnd;
      }
    }
    if (posStart != posEnd) chrSE.insert(std::make_pair(curChr, std::make_pair(posStart, posEnd)));
  }
}

template<typename TChr, typename TPos, typename TQual>
inline void
run(Config const& c)
{
  // Valid interval file?
  if (!(boost::filesystem::exists(c.int_file) && boost::filesystem::is_regular_file(c.int_file) && boost::filesystem::file_size(c.int_file))) {
    std::cerr << "Error: " << c.int_file.string() << " does not exist or it is empty." << std::endl;
    return;
  }

  // Read bam files
  BamTools::BamReader readerFirst;
  if ( ! readerFirst.Open(c.files[0].string()) ) {
    std::cerr << "Could not open input bam file!" << std::endl;
    return;
  }

  // Built chromosome map
  typedef std::map<std::string, TChr> TChrToChar;
  TChrToChar chrToChar;
  const BamTools::RefVector references = readerFirst.GetReferenceData();
  BamTools::RefVector::const_iterator refIt = references.begin();
  BamTools::RefVector::const_iterator refItEnd = references.end();
  for (TChr fa_count=0;refIt!=refItEnd; ++refIt, ++fa_count) {
    chrToChar.insert(std::make_pair(refIt->RefName, fa_count));
  }

  // Put the SV intervals into an interval tree (one for each chromosome)
  typedef Interval<unsigned int> TInterval;
  typedef IntervalTree<TInterval> TIntervalTree;
  typedef std::map<TChr, TIntervalTree*> TChrIntervalTrees;
  TChrIntervalTrees chrIntervals;
  typedef Record<std::string, unsigned int, unsigned int, void, void, void, void, void, void, void, void, void> TRecord;
  unsigned int line_counter = 0;
  Memory_mapped_file map_file(c.int_file.c_str());
  char buffer[Memory_mapped_file::MAX_LINE_LENGTH];
  while (map_file.left_bytes() > 0) {
    map_file.read_line(buffer);
    Tokenizer token(buffer, Memory_mapped_file::MAX_LINE_LENGTH);
    TRecord line;
    addF0(token, line);
    addF1(token, line);
    addF2(token, line);
    typename TChrToChar::const_iterator chrEx = chrToChar.find( line.f0 );
    if (chrEx == chrToChar.end()) {
      std::cerr << "Warning: " << line.f0 << ":" << line.f1 << "-" << line.f2 << " does not exist in bam file. Interval will be ignored!" << std::endl;
      continue;
    }
    typename TChrIntervalTrees::iterator findTree = chrIntervals.find( chrEx->second );
    if (findTree == chrIntervals.end()) findTree = chrIntervals.insert(std::make_pair(chrEx->second, new TIntervalTree())).first;
    TIntervalTree* iTree = findTree->second;
    TInterval newInt;
    newInt.low = line.f1 - c.bpWindowOffset;
    newInt.high = (c.bpWindowOffset) ? (line.f1 + c.bpWindowOffset) : (line.f1 + 1);
    newInt.cargo = line_counter++;
    iTree->insertInterval(newInt);
    newInt.low = line.f2 - c.bpWindowOffset;
    newInt.high = (c.bpWindowOffset) ? (line.f2 + c.bpWindowOffset) : (line.f2 + 1);
    newInt.cargo = line_counter++;
    iTree->insertInterval(newInt);
  }

  // Create hit vector
  typedef Hit<TChr, TPos, TQual> THit;
  typedef std::vector<THit> THits;
  THits hit_vector;  

  // Create library objects
  typedef std::map<std::string, LibraryInfo> TLibraryMap;
  TLibraryMap libInfo;

  // Store all spanning ranges
  typedef std::vector<THit> THits;
  THits normalSpan;  
  THits missingSpan;  

  // Read all input alignments
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    getLibraryParams(c.files[file_c], libInfo, 0, 3);
    unsigned int maxNormalISize = libInfo[file_c].median + 3 * libInfo[file_c].mad;
    unsigned int minNormalISize = (libInfo[file_c].median - 3 * libInfo[file_c].mad > 0) ? (libInfo[file_c].median - 3 * libInfo[file_c].mad) : 1;
    std::cout << "Library: " << c.files[file_c] << " (Median: " << libInfo[file_c].median << ", MAD: " << libInfo[file_c].mad << ", Orientation: " << (int) libInfo[file_c].defaultOrient << ", Insert size cutoffs: [" << minNormalISize << "," << maxNormalISize << "] )" << std::endl;
    int defOrient = libInfo[file_c].defaultOrient;
    if (libInfo[file_c].median == 0) continue; // Single-end library
    BamTools::BamReader reader;
    if ( ! reader.Open(c.files[file_c].string()) ) {
      std::cerr << "Could not open " << c.files[file_c].string() << " !" << std::endl;
      return;
    }
    BamTools::BamAlignment al;
    while( reader.GetNextAlignmentCore(al) ) {
      if ((al.AlignmentFlag & 0x0004) || (al.AlignmentFlag & 0x0200) || (al.AlignmentFlag & 0x0400) || (al.AlignmentFlag & 0x0100)) continue;
      if (al.MapQuality < c.minMapQual) continue;
      if (al.RefID < (TPos) references.size()) {
	// ToDo: Single-anchored reads still missing!!!
	if ((al.AlignmentFlag & 0x0001) && !(al.AlignmentFlag & 0x0004) && !(al.AlignmentFlag & 0x0008)) {
	  if ((al.RefID == al.MateRefID) && (al.Position>al.MatePosition)) {
	    // Get the proper length of the fragment based on the CIGAR string
	    unsigned int num_start = 0;
	    unsigned int num_end = 0;
	    unsigned int readLen = 0;
	    std::string cigar = cigarString(al.CigarData);
	    std::string::const_iterator cig = cigar.begin();
	    std::string::const_iterator cigEnd = cigar.end();
	    for(;cig!=cigEnd; ++num_end, ++cig) {
	      if (((int) *cig >=48) && ((int) *cig <= 57)) continue;
	      unsigned int len = atoi(cigar.substr(num_start, (num_end-num_start)).c_str());
	      if (*cig == 'M') readLen+=len;
	      else if ((*cig == 'D') && (len <= 5)) readLen+=len;
	      else if ((*cig == 'D') && (len > 5)) break;
	      else if ((*cig == 'I') && (len > 5)) break;
	      else if ((*cig == 'N') || (*cig == 'S') || (*cig == 'H')) break;
	      num_start = num_end + 1;
	    }	    
	    if (readLen > 0) {
	      BamTools::BamAlignment alBound;
	      alBound.RefID=al.RefID;
	      alBound.MapQuality=al.MapQuality;
	      // Normal spanning coverage
	      unsigned int outerISize = (al.Position + readLen) - al.MatePosition;
	      if ((getStrandIndependentOrientation(al) == defOrient) && (outerISize >= minNormalISize) && (outerISize <= maxNormalISize)) {
		// Abuse Position and MatePosition as (low,high) fragment boundaries
		alBound.Position=al.MatePosition;
		alBound.MatePosition=al.Position+readLen;
		_insertInterval(chrIntervals, alBound, normalSpan);
	      }
	      // Missing spanning coverage
	      if ((getStrandIndependentOrientation(al) != defOrient) || (outerISize >= (unsigned int) libInfo[file_c].median + 5 * libInfo[file_c].mad)) {
		if (_mateIsUpstream(defOrient, (al.AlignmentFlag & 0x0040), (al.AlignmentFlag & 0x0010))) {
		  alBound.Position=al.Position;
		  alBound.MatePosition=al.Position+minNormalISize;
		  _insertInterval(chrIntervals, alBound, missingSpan);
		} else {
		  alBound.Position=al.Position + readLen - minNormalISize;
		  alBound.MatePosition=al.Position+readLen;
		  _insertInterval(chrIntervals, alBound, missingSpan);
		}
		if (_mateIsUpstream(defOrient, !(al.AlignmentFlag & 0x0040), (al.AlignmentFlag & 0x0020))) {
		  alBound.Position=al.MatePosition;
		  alBound.MatePosition=al.MatePosition+minNormalISize;
		  _insertInterval(chrIntervals, alBound, missingSpan);
		} else { 
		  alBound.Position=al.MatePosition + readLen - minNormalISize;
		  alBound.MatePosition=al.MatePosition + readLen;
		  _insertInterval(chrIntervals, alBound, missingSpan);
		}
	      }
	    }
	  }
	}
      }
    }
  }

  // Free all allocated interval trees
  typename TChrIntervalTrees::iterator treeBeg = chrIntervals.begin();
  typename TChrIntervalTrees::iterator treeEnd = chrIntervals.end();
  for(;treeBeg != treeEnd; ++treeBeg) delete treeBeg->second;
  
  // Sort SAM records by chromosome and start position
  sort(normalSpan.begin(), normalSpan.end(), SortHits<THit>());
  sort(missingSpan.begin(), missingSpan.end(), SortHits<THit>());

  // Get all chromosome ranges in the hit vectors
  typedef std::pair<unsigned int, unsigned int> TStartEnd;
  typedef std::map<TChr, TStartEnd> TChrSE;
  TChrSE chrNormalSE;
  _dissectHitVector(normalSpan, chrNormalSE);
  TChrSE chrMissingSE;
  _dissectHitVector(missingSpan, chrMissingSE);

  // Declare the chromosome array
  typedef unsigned short TArrayType;
  TArrayType* normalCount = new TArrayType[MAX_CHROM_SIZE];
  TArrayType* missingCount = new TArrayType[MAX_CHROM_SIZE];


  // Output file
  boost::iostreams::filtering_ostream dataOut;
  dataOut.push(boost::iostreams::gzip_compressor());
  dataOut.push(boost::iostreams::file_sink(c.outfile.c_str(), std::ios_base::out | std::ios_base::binary));

  // Process each chromosome
  std::cout << "Processing chromosomes..." << std::endl;
  typename TChrToChar::const_iterator chrToCharIt = chrToChar.begin();
  typename TChrToChar::const_iterator chrToCharItEnd = chrToChar.end();
  for(;chrToCharIt != chrToCharItEnd; ++chrToCharIt) {
    const TPos chrLen = references[(int) chrToCharIt->second].RefLength;
    const std::string chrName = references[(int) chrToCharIt->second].RefName;
    std::cout << chrName << " (Length: " << chrLen << ")" << std::endl;

    // Set all values to zero
    std::fill(normalCount, normalCount + MAX_CHROM_SIZE, 0);
    std::fill(missingCount, missingCount + MAX_CHROM_SIZE, 0);

    // Iterate all reads of that chromosome
    typename TChrSE::const_iterator chrSEit = chrNormalSE.find(chrToCharIt->second);
    if (chrSEit != chrNormalSE.end()) 
      _addReadAndBpCounts(normalSpan, chrSEit, normalCount);
    chrSEit = chrMissingSE.find(chrToCharIt->second);
    if (chrSEit != chrMissingSE.end())
      _addReadAndBpCounts(missingSpan, chrSEit, missingCount);

    // Write spanning coverage for all input intervals
    if (isValidFile(c.int_file.string())) {
        Memory_mapped_file interval_file(c.int_file.string().c_str());
	char interval_buffer[Memory_mapped_file::MAX_LINE_LENGTH];
	while (interval_file.left_bytes() > 0) {
	  interval_file.read_line(interval_buffer);
	  // Read single interval line
	  Tokenizer token(interval_buffer, Memory_mapped_file::MAX_LINE_LENGTH);
	  std::string interval_rname;
	  token.getString(interval_rname);
	  if (interval_rname.compare(chrToCharIt->first)) continue;
	  int intervalStart = token.getUInt();
	  int intervalEnd = token.getUInt();
	  std::string id;
	  token.getString(id);
	  // First breakpoint
	  int posStart = (intervalStart - c.bpWindowOffset < 0) ? 0 : (intervalStart - c.bpWindowOffset);
	  int posEnd = (c.bpWindowOffset) ? (intervalStart + c.bpWindowOffset) : (intervalStart + 1);
	  if (!c.mapq) {
	    TArrayType* normalCountPoint = &normalCount[posStart];
	    TArrayType* missingCountPoint = &missingCount[posStart];
	    for(int i=posStart; i<posEnd; ++i, ++normalCountPoint, ++missingCountPoint) dataOut << id << "\t" << 0 << "\t" << chrToCharIt->first << "\t" << i << "\t" << *normalCountPoint << "\t" << *missingCountPoint << std::endl;
	  } else {
	    std::vector<std::string> normalStr;
	    normalStr.resize(posEnd-posStart);
	    std::fill(normalStr.begin(), normalStr.end(), "");
	    chrSEit = chrNormalSE.find(chrToCharIt->second);
	    if (chrSEit != chrNormalSE.end())
	      _buildMAPQString(normalSpan, chrSEit, posStart, posEnd, normalStr);
	    std::vector<std::string> missingStr;
	    missingStr.resize(posEnd-posStart);
	    std::fill(missingStr.begin(), missingStr.end(), "");
	    chrSEit = chrMissingSE.find(chrToCharIt->second);
	    if (chrSEit != chrMissingSE.end())
	      _buildMAPQString(missingSpan, chrSEit, posStart, posEnd, missingStr);
	    for(int i=posStart; i<posEnd; ++i) {
	      if (normalStr[i - posStart].size()) normalStr[i - posStart].erase(normalStr[i - posStart].size()-1);
	      if (missingStr[i - posStart].size()) missingStr[i - posStart].erase(missingStr[i - posStart].size()-1);
	      dataOut << id << "\t" << 0 << "\t" << chrToCharIt->first << "\t" << i << "\t" << normalStr[i - posStart] << "\t" << missingStr[i-posStart] << std::endl;
	    }
	  }
	  
	  // Second breakpoint
	  posStart = (intervalEnd - c.bpWindowOffset < 0) ? 0 : (intervalEnd - c.bpWindowOffset);
	  posEnd = (c.bpWindowOffset) ? (intervalEnd + c.bpWindowOffset) : (intervalEnd + 1);
	  if (!c.mapq) {
	    TArrayType* normalCountPoint = &normalCount[posStart];
	    TArrayType* missingCountPoint = &missingCount[posStart];
	    for(int i=posStart; i<posEnd; ++i, ++normalCountPoint, ++missingCountPoint) dataOut << id << "\t" << 1 << "\t" << chrToCharIt->first << "\t" << i << "\t" << *normalCountPoint << "\t" << *missingCountPoint << std::endl;
	  } else {
	    std::vector<std::string> normalStr;
	    normalStr.resize(posEnd-posStart);
	    std::fill(normalStr.begin(), normalStr.end(), "");
	    chrSEit = chrNormalSE.find(chrToCharIt->second);
	    if (chrSEit != chrNormalSE.end())
	      _buildMAPQString(normalSpan, chrSEit, posStart, posEnd, normalStr);
	    std::vector<std::string> missingStr;
	    missingStr.resize(posEnd-posStart);
	    std::fill(missingStr.begin(), missingStr.end(), "");
	    chrSEit = chrMissingSE.find(chrToCharIt->second);
	    if (chrSEit != chrMissingSE.end())
	      _buildMAPQString(missingSpan, chrSEit, posStart, posEnd, missingStr);
	    for(int i=posStart; i<posEnd; ++i) {
	      if (normalStr[i - posStart].size()) normalStr[i - posStart].erase(normalStr[i - posStart].size()-1);
	      if (missingStr[i - posStart].size()) missingStr[i - posStart].erase(missingStr[i - posStart].size()-1);
	      dataOut << id << "\t" << 1 << "\t" << chrToCharIt->first << "\t" << i << "\t" << normalStr[i-posStart] << "\t" << missingStr[i-posStart] << std::endl;
	    }
	  }
	}
	interval_file.close();
    }
  }
  
  // Clean-up
  delete[] normalCount;
  delete[] missingCount;

  // End
  std::cout << '[' << boost::posix_time::second_clock::local_time() << "] Done." << std::endl;
}


int main(int argc, char **argv) {
  Config c;
  c.mapq=false;

  // Define required options
  boost::program_options::options_description required("Required options");
  required.add_options()
    ("intervals,i", boost::program_options::value<boost::filesystem::path>(&c.int_file), "SV interval file")
    ;

  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("show-mapq,s", "use list of PE MAPQ instead of PE counts")
    ("quality-cut,q", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(0), "min. paired-end mapping quality")
    ("bp-offset,b", boost::program_options::value<int>(&c.bpWindowOffset)->default_value(1000), "breakpoint offset")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("span.gz"), "spanning coverage output file")
    ;

  // Define hidden options
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
    ;
  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  // Set the visibility
  boost::program_options::options_description cmdline_options;
  cmdline_options.add(required).add(generic).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(required).add(generic);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);


  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("intervals")) || (!vm.count("input-file"))) { 
    printTitle("Spanning coverage calculation");
    std::cout << "Usage: " << argv[0] << " [OPTIONS] <aligned1.sam/bam> <aligned2.sam/bam> ..." << std::endl;
    std::cout << visible_options << "\n"; 
    return 1; 
  }

  // Show cmd
  std::cout << '[' << boost::posix_time::second_clock::local_time() << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Do we need to include the qualities
  if (vm.count("show-mapq")) c.mapq=true;

  // Clear the outfile because we use append later
  unlink(c.outfile.c_str());

  // Run spanning coverage
  if (c.mapq) run<int32_t, int32_t, unsigned short>(c);
  else run<int32_t, int32_t, void>(c);

  return 0;
}
