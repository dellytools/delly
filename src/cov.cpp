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
#include <boost/icl/split_interval_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/progress.hpp>

#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/sam.h>

#ifdef OPENMP
#include <omp.h>
#endif

#include "gc.h"
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
  unsigned short madCutoff;
  bool bp_flag;
  bool avg_flag;
  bool inclCigar;
  bool hasIntervalFile;
  bool hasInvariantRegions;
  boost::filesystem::path outfile;
  boost::filesystem::path invariant;
  boost::filesystem::path int_file;
  boost::filesystem::path genome;
  std::vector<boost::filesystem::path> files;
  std::vector<std::string> sampleName;
};


inline int
coverageRun(Config const& c)
{
  // Get references
  samFile* samfile = sam_open(c.files[0].string().c_str(), "r");
  bam_hdr_t* hdr = sam_hdr_read(samfile);

  // GC Correction
  typedef boost::unordered_map<std::string, LibraryInfo> TLibraryMap;
  typedef std::vector<TLibraryMap> TSampleLibrary;
  TSampleLibrary sampleLib(c.files.size());
  typedef boost::unordered_map<std::string, GCBias> TGCMap;
  typedef std::vector<TGCMap> TSampleGC;
  TSampleGC gcBias(c.files.size());
  if (c.hasInvariantRegions) {
    typedef boost::icl::interval_set<uint32_t> TChrIntervals;
    typedef TChrIntervals::interval_type TIVal;
    typedef std::vector<TChrIntervals> TRegionsGenome;
    TRegionsGenome invariantRegions;
    invariantRegions.resize(hdr->n_targets);
    std::ifstream file(c.invariant.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
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
	    invariantRegions[tid].insert(TIVal::right_open(start, end));
	  }
	}
      }
    }
    file.close();

    // Debug code
    //for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
    //for(typename TChrIntervals::const_iterator vRIt = invariantRegions[refIndex].begin(); vRIt != invariantRegions[refIndex].end(); ++vRIt) {
    //std::cerr << std::string(hdr->target_name[refIndex]) << "\t" << vRIt->lower() << "\t" << vRIt->upper() << std::endl;
    //}
    //}

    // Create library objects
    getLibraryParams(c, invariantRegions, sampleLib);

    // Estimate GC bias
    gcBiasPerRG(c, invariantRegions, sampleLib, gcBias);
  }
  
  // Read all SV intervals
  typedef std::vector<CovRecord> TSVs;
  TSVs svs;
  std::vector<std::string> idToName;
  unsigned int intervalCount = 0;
  if (c.hasIntervalFile) {
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
	  int32_t tid = bam_name2id(hdr, chrName.c_str());
	  if (tid >= 0) {
	    if (tokIter!=tokens.end()) {
	      CovRecord sv;	  
	      sv.chr = tid;
	      sv.svStart = boost::lexical_cast<int32_t>(*tokIter++);
	      sv.svEnd = boost::lexical_cast<int32_t>(*tokIter++);
	      std::string svName = *tokIter;
	      idToName.push_back(svName);
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
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      uint32_t pos = 0;
      unsigned int wSize = c.window_size;
      unsigned int wOffset = c.window_offset;
      if (c.window_num>0) {
	wSize=(hdr->target_len[refIndex] / c.window_num) + 1;
	wOffset=wSize;
      }
      while (pos < hdr->target_len[refIndex]) {
	uint32_t window_len = pos+wSize;
	if (window_len > hdr->target_len[refIndex]) window_len = hdr->target_len[refIndex];
	CovRecord sv;
	sv.chr = refIndex;
	sv.svStart = pos;
	sv.svEnd = window_len;
	std::string refname(hdr->target_name[refIndex]);
	refname += ":" + boost::lexical_cast<std::string>(sv.svStart) + "-" + boost::lexical_cast<std::string>(sv.svEnd);
	idToName.push_back(refname);
	sv.id = intervalCount++;
	svs.push_back(sv);
	pos += wOffset;
      }
    }
  }

  // Output data types
  typedef std::pair<int, int> TBpRead;
  typedef std::vector<TBpRead> TSVReadCount;
  typedef std::vector<TSVReadCount> TSampleSVReadCount;
  TSampleSVReadCount countMap;

  // Annotate coverage
  if (c.hasInvariantRegions) {
    if (c.inclCigar) annotateCoverageGCAware(c, gcBias, svs, countMap, BpLevelType<BpLevelCount>());
    else annotateCoverageGCAware(c, gcBias, svs, countMap, BpLevelType<NoBpLevelCount>());
  } else {
    if (c.inclCigar) annotateCoverage(c.files, c.minGenoQual, svs, countMap, BpLevelType<BpLevelCount>());
    else annotateCoverage(c.files, c.minGenoQual, svs, countMap, BpLevelType<NoBpLevelCount>());
  }




  
  // Output file
  boost::iostreams::filtering_ostream dataOut;
  dataOut.push(boost::iostreams::gzip_compressor());
  dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));

  // Print header
  dataOut << "chr\tstart\tend\tid";
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
    dataOut << hdr->target_name[itSV->chr] << "\t" << itSV->svStart << "\t" << itSV->svEnd << "\t" << idToName[itSV->id];
    // Iterate all samples
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      dataOut << "\t";
      if (c.avg_flag) dataOut << ( (countMap[file_c][itSV->id].first) / (double) (itSV->svEnd - itSV->svStart)) << "\t";
      if (c.bp_flag) dataOut << countMap[file_c][itSV->id].first << "\t";
      dataOut << countMap[file_c][itSV->id].second;
    }
    dataOut << std::endl;
  }

  // Clean-up
  bam_hdr_destroy(hdr);
  sam_close(samfile);

  // Output Statistics
  if (c.hasInvariantRegions) {
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      std::cout << "Sample: " << c.sampleName[file_c] << std::endl;
      for(TLibraryMap::const_iterator libIt = sampleLib[file_c].begin(); libIt != sampleLib[file_c].end(); ++libIt) {
	TGCMap::const_iterator gcIt = gcBias[file_c].find(libIt->first);
	std::cout << "RG: ID=" << libIt->first << ",ReadSize=" << libIt->second.rs << ",AvgDist=" << libIt->second.avgdist << ",EstCov=" << (double) libIt->second.rs / (double) libIt->second.avgdist << ",MappedAsPair=" << libIt->second.mappedAsPair << ",Median=" << libIt->second.median << ",MAD=" << libIt->second.mad << ",Layout=" << (int) libIt->second.defaultOrient << ",GC-Correction=[";
	for(uint32_t gcbin = 0; gcbin < gcIt->second.gc.size() - 1; ++gcbin) std::cout << gcIt->second.gc[gcbin] << ',';
	std::cout << gcIt->second.gc[gcIt->second.gc.size() - 1] << "]" << std::endl;
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
  c.madCutoff = 5;

  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("avg-cov,a", "show average coverage")
    ("bp-count,b", "show base pair count")
    ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
    ("quality-cut,q", boost::program_options::value<uint16_t>(&c.minGenoQual)->default_value(1), "exclude all alignments with quality < q")
    ("outfile,f", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("cov.gz"), "coverage output file")
    ;

  // Define window options
  boost::program_options::options_description window("Window options");
  window.add_options()
    ("window-size,s", boost::program_options::value<unsigned int>(&c.window_size)->default_value(10000), "window size")
    ("window-offset,o", boost::program_options::value<unsigned int>(&c.window_offset)->default_value(5000), "window offset")
    ("window-num,n", boost::program_options::value<unsigned int>(&c.window_num)->default_value(0), "#windows per chr, used if #n>0")
    ("interval-file,i", boost::program_options::value<boost::filesystem::path>(&c.int_file), "interval file, used if present")
    ;

  // Define interval options
  boost::program_options::options_description gc("GC content correction (slow)");
  gc.add_options()
    ("invariant,r", boost::program_options::value<boost::filesystem::path>(&c.invariant), "invariant regions")
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
  cmdline_options.add(generic).add(window).add(gc).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic).add(window).add(gc);
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
      std::cout << "Usage: " << argv[0] << " [OPTIONS] -g <ref.fa> <sample1.bam> <sample2.bam> ..." << std::endl;
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

  // Check reference
  if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
    std::cerr << "Reference file is missing! " << c.genome.string() << std::endl;
    return 1;
  } else {
    faidx_t* fai = fai_load(c.genome.string().c_str());
    if (fai == NULL) {
      if (fai_build(c.genome.string().c_str()) == -1) {
	std::cerr << "Fail to open genome fai index for " << c.genome.string() << std::endl;
	return 1;
      } else fai = fai_load(c.genome.string().c_str());
    }
    fai_destroy(fai);
  }

  // Check input files
  c.sampleName.resize(c.files.size());
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    if (!(boost::filesystem::exists(c.files[file_c]) && boost::filesystem::is_regular_file(c.files[file_c]) && boost::filesystem::file_size(c.files[file_c]))) {
      std::cerr << "Alignment file is missing: " << c.files[file_c].string() << std::endl;
      return 1;
    }
    samFile* samfile = sam_open(c.files[file_c].string().c_str(), "r");
    if (samfile == NULL) {
      std::cerr << "Fail to open file " << c.files[file_c].string() << std::endl;
      return 1;
    }
    hts_idx_t* idx = sam_index_load(samfile, c.files[file_c].string().c_str());
    if (idx == NULL) {
      std::cerr << "Fail to open index for " << c.files[file_c].string() << std::endl;
      return 1;
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    if (hdr == NULL) {
      std::cerr << "Fail to open header for " << c.files[file_c].string() << std::endl;
      return 1;
    }
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
      std::string tname(hdr->target_name[refIndex]);
      if (!faidx_has_seq(fai, tname.c_str())) {
	std::cerr << "BAM file chromosome " << hdr->target_name[refIndex] << " is NOT present in your reference file " << c.genome.string() << std::endl;
	return 1;
      }
    }
    fai_destroy(fai);
    std::string sampleName;
    if (!getSMTag(std::string(hdr->text), c.files[file_c].stem().string(), sampleName)) {
      std::cerr << "Only one sample (@RG:SM) is allowed per input BAM file " << c.files[file_c].string() << std::endl;
      return 1;
    } else c.sampleName[file_c] = sampleName;
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }

  // Check input intervals
  if (vm.count("interval-file")) {
    if (!(boost::filesystem::exists(c.int_file) && boost::filesystem::is_regular_file(c.int_file) && boost::filesystem::file_size(c.int_file))) {
      std::cerr << "Interval file is missing: " << c.int_file.string() << std::endl;
      return 1;
    }
    std::string oldChr;
    faidx_t* fai = fai_load(c.genome.string().c_str());
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
	  if (chrName.compare(oldChr) != 0) {
	    oldChr = chrName;
	    if (!faidx_has_seq(fai, chrName.c_str())) {
	      std::cerr << "Interval file chromosome " << chrName << " is NOT present in your reference file " << c.genome.string() << std::endl;
	      return 1;
	    }
	  }
	}
      }
      interval_file.close();
    }
    fai_destroy(fai);
    c.hasIntervalFile= true;
  } else c.hasIntervalFile = false;


  // Check invariant regions
  if (vm.count("invariant")) {
    if (!(boost::filesystem::exists(c.invariant) && boost::filesystem::is_regular_file(c.invariant) && boost::filesystem::file_size(c.invariant))) {
      std::cerr << "Invariant regions file is missing: " << c.invariant.string() << std::endl;
      return 1;
    }
    c.hasInvariantRegions = true;
  } else c.hasInvariantRegions = false;
  
  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;
 
  // Run coverage annotation
  return coverageRun(c);
}
