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
#include <boost/unordered_map.hpp>
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
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/progress.hpp>

#include <htslib/vcf.h>
#include <htslib/sam.h>

#ifdef OPENMP
#include <omp.h>
#endif

#include "tags.h"
#include "coverage.h"
#include "version.h"
#include "util.h"

using namespace torali;

struct Config {
  unsigned int window_size;
  unsigned int window_offset;
  unsigned int window_num;
  uint16_t minGenoQual;
  bool bp_flag;
  bool avg_flag;
  bool inclCigar;
  boost::filesystem::path outfile;
  boost::filesystem::path int_file;
  std::vector<boost::filesystem::path> files;
};


inline int
run(Config const& c)
{
  // Create library objects
  typedef boost::unordered_map<std::string, LibraryInfo> TLibraryMap;
  typedef boost::unordered_map<std::string, TLibraryMap> TSampleLibrary;
  TSampleLibrary sampleLib;
  getLibraryParams(c.files, sampleLib, 0, 5);

  // Get references
  typedef std::vector<std::string> TRefNames;
  typedef std::vector<uint32_t> TRefLength;
  TRefNames refnames;
  TRefLength reflen;
  if (refnames.empty()) {
    samFile* samfile = sam_open(c.files[0].string().c_str(), "r");
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    for (int i = 0; i<hdr->n_targets; ++i) {
      refnames.push_back(hdr->target_name[i]);
      reflen.push_back(hdr->target_len[i]);
    }
    bam_hdr_destroy(hdr);
    sam_close(samfile);
  }

  // Read all SV intervals
  typedef std::vector<CovRecord> TSVs;
  TSVs svs;
  std::map<unsigned int, std::string> idToName;
  unsigned int intervalCount=1;
  if (boost::filesystem::exists(c.int_file) && boost::filesystem::is_regular_file(c.int_file) && boost::filesystem::file_size(c.int_file)) {
    typedef boost::unordered_map<std::string, unsigned int> TMapChr;
    TMapChr mapChr;
    for(unsigned int i = 0; i<refnames.size(); ++i) mapChr[ refnames[i] ] = i;
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
	      CovRecord sv;	  
	      sv.chr = mapChrIt->second;
	      sv.svStart = boost::lexical_cast<int32_t>(*tokIter++);
	      sv.svEnd = boost::lexical_cast<int32_t>(*tokIter++) + 1;
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
  } else {
    // Create artificial intervals
    for(int refIndex=0; refIndex < (int) refnames.size(); ++refIndex) {
      uint32_t pos = 0;
      unsigned int wSize = c.window_size;
      unsigned int wOffset = c.window_offset;
      if (c.window_num>0) {
	wSize=(reflen[refIndex] / c.window_num) + 1;
	wOffset=wSize;
      }
      while (pos < reflen[refIndex]) {
	uint32_t window_len = pos+wSize;
	if (window_len > reflen[refIndex]) window_len = reflen[refIndex];
	CovRecord sv;
	sv.chr = refIndex;
	sv.svStart = pos;
	sv.svEnd = window_len;
	std::stringstream s; 
	s << refnames[sv.chr] << ":" << sv.svStart << "-" << sv.svEnd;
	idToName.insert(std::make_pair(intervalCount, s.str()));
	sv.id = intervalCount++;
	svs.push_back(sv);
	pos += wOffset;
      }
    }
  }

  // Output data types
  typedef std::pair<std::string, int> TSampleSVPair;
  typedef std::pair<int, int> TBpRead;
  typedef std::map<TSampleSVPair, TBpRead> TCountMap;
  TCountMap countMap;

  // Annotate coverage
  if (c.inclCigar) annotateCoverage(c.files, c.minGenoQual, sampleLib, svs, countMap, BpLevelType<BpLevelCount>());
  else annotateCoverage(c.files, c.minGenoQual, sampleLib, svs, countMap, BpLevelType<NoBpLevelCount>());

  // Output library statistics
  std::cout << "Library statistics" << std::endl;
  TSampleLibrary::const_iterator sampleIt=sampleLib.begin();
  for(;sampleIt!=sampleLib.end();++sampleIt) {
    std::cout << "Sample: " << sampleIt->first << std::endl;
    TLibraryMap::const_iterator libIt=sampleIt->second.begin();
    for(;libIt!=sampleIt->second.end();++libIt) {
      std::cout << "RG: ID=" << libIt->first << ",Median=" << libIt->second.median << ",MAD=" << libIt->second.mad << ",Orientation=" << (int) libIt->second.defaultOrient << std::endl;
    }
  }

  // Output file
  boost::iostreams::filtering_ostream dataOut;
  dataOut.push(boost::iostreams::gzip_compressor());
  dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));

  // Print header
  dataOut << "#chr\tstart\tend\tid";
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    std::string sampleName(c.files[file_c].stem().string());
    dataOut << "\t";
    if (c.avg_flag) dataOut << sampleName << "_avgcov" << "\t";
    if (c.bp_flag) dataOut << sampleName << "_bpcount" << "\t";
    if ((c.bp_flag) || (c.avg_flag)) dataOut << sampleName << "_readcount";
    else dataOut << sampleName;
  }
  dataOut << std::endl;

  // Iterate all SVs
  TSVs::const_iterator itSV = svs.begin();
  TSVs::const_iterator itSVEnd = svs.end();
  for(;itSV!=itSVEnd;++itSV) {
    dataOut << refnames[itSV->chr] << "\t" << itSV->svStart << "\t" << itSV->svEnd << "\t" << idToName.find(itSV->id)->second;
    // Iterate all samples
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      // Get the sample name
      std::string sampleName(c.files[file_c].stem().string());
      TSampleSVPair sampleSVPair = std::make_pair(sampleName, itSV->id);
      TCountMap::iterator countMapIt=countMap.find(sampleSVPair);
      dataOut << "\t";
      if (c.avg_flag) dataOut << ( (countMapIt->second.first) / (double) (itSV->svEnd - itSV->svStart)) << "\t";
      if (c.bp_flag) dataOut << countMapIt->second.first << "\t";
      dataOut << countMapIt->second.second;
    }
    dataOut << std::endl;
  }

  // End
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
  return 0;
}


int main(int argc, char **argv) {
  Config c;

  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("avg-cov,a", "show average coverage")
    ("bp-count,b", "show base pair count")
    ("quality-cut,q", boost::program_options::value<uint16_t>(&c.minGenoQual)->default_value(0), "exclude all alignments with quality < q")
    ("outfile,f", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("cov.gz"), "coverage output file")
    ;

  // Define window options
  boost::program_options::options_description window("Window options");
  window.add_options()
    ("window-size,s", boost::program_options::value<unsigned int>(&c.window_size)->default_value(10000), "window size")
    ("window-offset,o", boost::program_options::value<unsigned int>(&c.window_offset)->default_value(10000), "window offset")
    ("window-num,n", boost::program_options::value<unsigned int>(&c.window_num)->default_value(0), "#windows per chr, used if #n>0")
    ;

  // Define interval options
  boost::program_options::options_description interval("Interval options");
  interval.add_options()
    ("interval-file,i", boost::program_options::value<boost::filesystem::path>(&c.int_file), "interval file")
    ;

  // Define hidden options
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
    ("license,l", "show license")
    ("warranty,w", "show warranty")
    ;
  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  // Set the visibility
  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(window).add(interval).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic).add(window).add(interval);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);


  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file"))) { 
    printTitle("Coverage calculation");
    if (vm.count("warranty")) {
      displayWarranty();
    } else if (vm.count("license")) {
      gplV3();
    } else {
      std::cout << "Usage: " << argv[0] << " [OPTIONS] <sample1.bam> <sample2.bam> ..." << std::endl;
      std::cout << visible_options << "\n"; 
    }
    return 1; 
  }
  if (vm.count("bp-count")) c.bp_flag = true;
  else c.bp_flag = false;
  if (vm.count("avg-cov")) c.avg_flag = true;
  else c.avg_flag = false;
  if ((c.bp_flag) || (c.avg_flag)) c.inclCigar = true;
  else c.inclCigar = false;

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;
 
  // Run coverage annotation
  return run(c);
}
