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
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/tokenizer.hpp>
#include <boost/functional/hash.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include "api/BamIndex.h"
#include "api/BamReader.h"

#ifdef OPENMP
#include <omp.h>
#endif

#include "version.h"
#include "util.h"
#include "tags.h"
#include "spanning.h"
#include "intervaltree.h"


using namespace torali;

struct Config {
  bool mapq;
  uint16_t minMapQual;
  int bpWindowOffset;
  std::string svType;
  boost::filesystem::path int_file;
  boost::filesystem::path outfile;
  std::vector<boost::filesystem::path> files;
};

template<typename TDataOut>
inline void
_outputQualities(TDataOut& dataOut, int normal, int abnormal) {
  dataOut << "\t" << normal << "\t" << abnormal;
}

template<typename TDataOut>
inline void
_outputQualities(TDataOut& dataOut, std::vector<uint16_t> const& normal, std::vector<uint16_t> const& abnormal) {
  typedef std::vector<uint16_t> TQualVector;
  
  dataOut << "\t";
  typename TQualVector::const_iterator itQ = normal.begin();
  typename TQualVector::const_iterator itQEnd = normal.end();
  for(;itQ!=itQEnd;++itQ) {
    if (itQ != normal.begin()) dataOut << ',';
    dataOut << *itQ;
  }
  dataOut << "\t";
  itQ = abnormal.begin();
  itQEnd = abnormal.end();
  for(;itQ!=itQEnd;++itQ) {
    if (itQ != abnormal.begin()) dataOut << ',';
    dataOut << *itQ;
  }
}


template<typename THitInterval, typename TCount, typename TSVType>
inline int
run(Config const& c, THitInterval, TCount, TSVType svType)
{
  // Valid interval file?
  if (!(boost::filesystem::exists(c.int_file) && boost::filesystem::is_regular_file(c.int_file) && boost::filesystem::file_size(c.int_file))) {
    std::cerr << "Error: " << c.int_file.string() << " does not exist or it is empty." << std::endl;
    return -1;
  }

  // Create library objects
  typedef std::map<std::string, LibraryInfo> TLibraryMap;
  typedef std::map<std::string, TLibraryMap> TSampleLibrary;
  TSampleLibrary sampleLib;

  // Scan libraries
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    // Get a sample name
    std::string sampleName(c.files[file_c].stem().string());

    // Check that all input bam files exist
    BamTools::BamReader reader;
    if ( ! reader.Open(c.files[file_c].string()) ) {
      std::cerr << "Could not open input bam file: " << c.files[file_c].string() << std::endl;
      reader.Close();
      return -1;
    }
    
    // Check that all input bam files are indexed
    reader.LocateIndex();
    if ( !reader.HasIndex() ) {
      std::cerr << "Missing bam index file: " << c.files[file_c].string() << std::endl;
      reader.Close();
      return -1;
    }

    // Get library parameters and overall maximum insert size
    TLibraryMap libInfo;
    getLibraryParams(c.files[file_c], libInfo, 0, 5);
    sampleLib.insert(std::make_pair(sampleName, libInfo));
  }

  // Get references
  BamTools::BamReader readerRef;
  if ( ! readerRef.Open(c.files[0].string()) ) return -1;
  BamTools::RefVector references = readerRef.GetReferenceData();

  // Read all SV intervals
  typedef std::vector<StructuralVariantRecord> TSVs;
  TSVs svs;
  std::map<unsigned int, std::string> idToName;
  unsigned int intervalCount=1;
  if (boost::filesystem::exists(c.int_file) && boost::filesystem::is_regular_file(c.int_file) && boost::filesystem::file_size(c.int_file)) {
    typedef boost::unordered_map<std::string, unsigned int> TMapChr;
    TMapChr mapChr;
    typename BamTools::RefVector::const_iterator itRef = references.begin();
    for(unsigned int i = 0;itRef!=references.end();++itRef, ++i) mapChr[ itRef->RefName ] = i;
    std::ifstream interval_file(c.int_file.string().c_str(), std::ifstream::in);
    if (interval_file.is_open()) {
      while (interval_file.good()) {
	std::string intervalLine;
	getline(interval_file, intervalLine);
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t,;");
	Tokenizer tokens(intervalLine, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter!=tokens.end()) {
	  std::string chrName=*tokIter++;
	  TMapChr::const_iterator mapChrIt = mapChr.find(chrName);
	  if (mapChrIt != mapChr.end()) {
	    if (tokIter!=tokens.end()) {
	      StructuralVariantRecord sv;
	      sv.chr = mapChrIt->second;
	      sv.chr2 = mapChrIt->second;
	      sv.svStart = boost::lexical_cast<int32_t>(*tokIter++);
	      sv.svEnd = boost::lexical_cast<int32_t>(*tokIter++);
	      std::string svName = *tokIter;
	      idToName.insert(std::make_pair(intervalCount, svName));
	      sv.id = intervalCount++;
	      svs.push_back(sv);
	    }
	  }
	}
      }
      interval_file.close();
    }
  }

  // Output data types
  typedef std::pair<std::string, int> TSampleSVPair;
  typedef std::vector<TCount> TCountRange;
  typedef std::map<TSampleSVPair, TCountRange> TCountMap;
  TCountMap normalCountMap;
  TCountMap abnormalCountMap;

  // Annotate spanning coverage
  annotateSpanningCoverage(c.files, c.bpWindowOffset, c.minMapQual, sampleLib, svs, normalCountMap, abnormalCountMap, THitInterval(), svType);

  // Output library statistics
  std::cout << "Library statistics" << std::endl;
  TSampleLibrary::const_iterator sampleIt=sampleLib.begin();
  for(;sampleIt!=sampleLib.end();++sampleIt) {
    std::cout << "Sample: " << sampleIt->first << std::endl;
    TLibraryMap::const_iterator libIt=sampleIt->second.begin();
    for(;libIt!=sampleIt->second.end();++libIt) {
      std::cout << "RG: ID=" << libIt->first << ",Median=" << libIt->second.median << ",MAD=" << libIt->second.mad << ",Orientation=" << (int) libIt->second.defaultOrient << ",MinInsertSize=" << libIt->second.minNormalISize << ",MaxInsertSize=" << libIt->second.maxNormalISize << ",DuplicatePairs=" << libIt->second.non_unique_pairs << ",UniquePairs=" << libIt->second.unique_pairs << std::endl;
    }
  }

  // Output file
  boost::iostreams::filtering_ostream dataOut;
  dataOut.push(boost::iostreams::gzip_compressor());
  dataOut.push(boost::iostreams::file_sink(c.outfile.c_str(), std::ios_base::out | std::ios_base::binary));

  // Iterate all SVs
  typename TSVs::const_iterator itSV = svs.begin();
  typename TSVs::const_iterator itSVEnd = svs.end();
  for(;itSV!=itSVEnd;++itSV) {
    // Iterate left and right breakpoint
    for(int bpOrder=-1;bpOrder<=1;bpOrder+=2) {
      // Iterate the breakpoint range
      int posStart = (itSV->svStart - c.bpWindowOffset < 0) ? 0 : (itSV->svStart - c.bpWindowOffset);
      int posEnd = (c.bpWindowOffset) ? (itSV->svStart + c.bpWindowOffset) : (itSV->svStart + 1);
      if (bpOrder==1) {
	posStart = (itSV->svEnd - c.bpWindowOffset < 0) ? 0 : (itSV->svEnd - c.bpWindowOffset);
	posEnd = (c.bpWindowOffset) ? (itSV->svEnd + c.bpWindowOffset) : (itSV->svEnd + 1);
      }
      for(int i=0;i<(posEnd-posStart);++i) {
	dataOut << idToName.find(itSV->id)->second << "\t" << ((bpOrder+1)/2) << "\t" << references[itSV->chr].RefName << "\t" << (i+posStart);
	// Iterate all samples
	for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	  // Get the sample name
	  std::string sampleName(c.files[file_c].stem().string());
	  TSampleSVPair sampleSVPair = std::make_pair(sampleName, bpOrder*itSV->id);
	  typename TCountMap::iterator countMapIt=normalCountMap.find(sampleSVPair);
	  typename TCountMap::iterator abCountMapIt=abnormalCountMap.find(sampleSVPair);
	  _outputQualities(dataOut, countMapIt->second[i], abCountMapIt->second[i]);
	}
	dataOut << std::endl;
      }
    }
  }
  
  // End
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
  return 0;
}


int main(int argc, char **argv) {
  Config c;

  // Define required options
  boost::program_options::options_description required("Required options");
  required.add_options()
    ("intervals,i", boost::program_options::value<boost::filesystem::path>(&c.int_file), "SV interval file")
    ;

  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("type,t", boost::program_options::value<std::string>(&c.svType)->default_value("DEL"), "SV analysis type (DEL, DUP, INV)")
    ("show-mapq,s", "use list of PE MAPQ instead of PE counts")
    ("quality-cut,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(0), "min. paired-end mapping quality")
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
    std::cout << "Usage: " << argv[0] << " [OPTIONS] <sample1.bam> <sample2.bam> ..." << std::endl;
    std::cout << visible_options << "\n"; 
    return -1; 
  }

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Do we need to include the qualities
  if (vm.count("show-mapq")) c.mapq=true;
  else c.mapq=false;

  // Run spanning coverage
  if (c.mapq) {
    if (c.svType == "DEL") return run(c, HitInterval<int32_t, uint16_t>(), std::vector<uint16_t>(), SVType<DeletionTag>());
    else if (c.svType == "DUP") return run(c, HitInterval<int32_t, uint16_t>(), std::vector<uint16_t>(), SVType<DuplicationTag>());
    else if (c.svType == "INV") return run(c, HitInterval<int32_t, uint16_t>(), std::vector<uint16_t>(), SVType<InversionTag>());
    else {
      std::cerr << "SV analysis type not supported by Delly: " << c.svType << std::endl;
      return -1;
    }
  } else {
    if (c.svType == "DEL") return run(c, HitInterval<int32_t, void>(), int(), SVType<DeletionTag>());
    else if (c.svType == "DUP") return run(c, HitInterval<int32_t, void>(), int(), SVType<DuplicationTag>());
    else if (c.svType == "INV") return run(c, HitInterval<int32_t, void>(), int(), SVType<InversionTag>());
    else {
      std::cerr << "SV analysis type not supported by Delly: " << c.svType << std::endl;
      return -1;
    }
  }
}
