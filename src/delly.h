/*
============================================================================
DELLY: Structural variant discovery by integrated PE mapping and SR analysis
============================================================================
Copyright (C) 2012-2018 Tobias Rausch

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

#ifndef DELLY_H
#define DELLY_H


#include <iostream>
#include <fstream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/icl/split_interval_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/functional/hash.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/sam.h>

#include "version.h"
#include "util.h"
#include "bolog.h"
#include "tags.h"
#include "coverage.h"
#include "msa.h"
#include "split.h"
#include "shortpe.h"
#include "modvcf.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <stdio.h>

namespace torali
{

// Config arguments
struct Config {
  unsigned short minMapQual;
  unsigned short minTraQual;
  unsigned short minGenoQual;
  unsigned short madCutoff;
  int32_t minimumFlankSize;
  int32_t indelsize;
  uint32_t graphPruning;
  float flankQuality;
  bool indels;
  bool hasExcludeFile;
  bool hasVcfFile;
  bool isHaplotagged;
  bool dumpflag;
  bool svtcmd;
  std::set<int32_t> svtset;
  DnaScore<int> aliscore;
  boost::filesystem::path outfile;
  boost::filesystem::path vcffile;
  boost::filesystem::path genome;
  boost::filesystem::path exclude;
  boost::filesystem::path srpedump;
  std::vector<boost::filesystem::path> files;
  std::vector<std::string> sampleName;
};



template<typename TConfig, typename TSampleLib, typename TSVs, typename TCountMap, typename TSampleSVJunctionMap, typename TSpanningCoverage>
inline void
_annotateCoverage(TConfig& c, bam_hdr_t* hdr, TSampleLib& sampleLib, TSVs& svs, TCountMap& countMap, TSampleSVJunctionMap& juncMap, TSpanningCoverage& spanMap) 
{
  // Find Ns in the reference genome
  typedef boost::icl::interval_set<int> TNIntervals;
  typedef std::vector<TNIntervals> TNGenome;
  TNGenome ni( hdr->n_targets );

  // Find valid chromosomes
  typedef std::vector<bool> TValidChr;
  TValidChr validChr(hdr->n_targets, false);
  for (typename TSVs::const_iterator itSV = svs.begin(); itSV != svs.end(); ++itSV) {
    validChr[itSV->chr] = true;
    validChr[itSV->chr2] = true;
  }

  // Parse Ns
  faidx_t* fai = fai_load(c.genome.string().c_str());
  for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
    if (!validChr[refIndex]) continue;
    std::string tname(hdr->target_name[refIndex]);
    int32_t seqlen = -1;
    char* seq = NULL;
    seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
    bool nrun = false;
    int nstart = seqlen;
    for(int i=0; i<seqlen; ++i) {
      if ((seq[i] != 'n') && (seq[i] != 'N')) {
	if (nrun) {
	  ni[refIndex].add(boost::icl::discrete_interval<int>::right_open(nstart,i));
	  nrun = false;
	}
      } else {
	if (!nrun) {
	  nrun = true;
	  nstart = i;
	}
      }
    }
    if (nrun) ni[refIndex].add(boost::icl::discrete_interval<int>::right_open(nstart,seqlen));
    if (seq != NULL) free(seq);
  }
  fai_destroy(fai);
  
  // Add control regions
  typedef std::vector<CovRecord> TCovRecord;
  typedef std::vector<TCovRecord> TGenomicCovRecord;
  TGenomicCovRecord svc(hdr->n_targets, TCovRecord());
  uint32_t lastId = svs.size();
  typedef std::vector<int32_t> TSVSize;
  TSVSize svSize(lastId);
  for (typename TSVs::const_iterator itSV = svs.begin(); itSV != svs.end(); ++itSV) {
    int halfSize = (itSV->svEnd - itSV->svStart)/2;
    if ((_translocation(itSV->svt)) || (itSV->svt == 4)) halfSize = 500;

    // Left control region
    CovRecord sLeft;
    sLeft.id = lastId + itSV->id;
    sLeft.svStart = std::max(itSV->svStart - halfSize, 0);
    sLeft.svEnd = itSV->svStart;
    typename TNIntervals::const_iterator itO = ni[itSV->chr].find(boost::icl::discrete_interval<int>::right_open(sLeft.svStart, sLeft.svEnd));
    while (itO != ni[itSV->chr].end()) {
      sLeft.svStart = std::max(itO->lower() - halfSize, 0);
      sLeft.svEnd = itO->lower();
      itO = ni[itSV->chr].find(boost::icl::discrete_interval<int>::right_open(sLeft.svStart, sLeft.svEnd));
    }
    svc[itSV->chr].push_back(sLeft);

    // Actual SV
    CovRecord sMiddle;
    sMiddle.id = itSV->id;
    sMiddle.svStart = itSV->svStart;
    sMiddle.svEnd = itSV->svEnd;
    if ((_translocation(itSV->svt)) || (itSV->svt == 4)) {
      sMiddle.svStart = std::max(itSV->svStart - halfSize, 0);
      sMiddle.svEnd = itSV->svStart + halfSize;      
    }
    svSize[itSV->id] = (itSV->svEnd - itSV->svStart);
    if (_translocation(itSV->svt)) svSize[itSV->id] = 50000000;
    svc[itSV->chr].push_back(sMiddle);
    

    // Right control region
    CovRecord sRight;
    sRight.id = 2 * lastId + itSV->id;
    sRight.svStart = itSV->svEnd;
    sRight.svEnd = itSV->svEnd + halfSize;
    itO = ni[itSV->chr2].find(boost::icl::discrete_interval<int>::right_open(sRight.svStart, sRight.svEnd));
    while (itO != ni[itSV->chr2].end()) {
      sRight.svStart = itO->upper();
      sRight.svEnd = itO->upper() + halfSize;
      itO = ni[itSV->chr2].find(boost::icl::discrete_interval<int>::right_open(sRight.svStart, sRight.svEnd));
    }
    svc[itSV->chr2].push_back(sRight);
    //std::cerr << itSV->id << ':' << sLeft.svStart << '-' << sLeft.svEnd << ',' << itSV->svStart << '-' << itSV->svEnd << ',' << sRight.svStart << '-' << sRight.svEnd << std::endl;
  }
  
  typedef std::pair<int32_t, int32_t> TBpRead;
  typedef std::vector<TBpRead> TSVReadCount;
  typedef std::vector<TSVReadCount> TSampleSVReadCount;
  TSampleSVReadCount readCountMap;
  annotateCoverage(c, sampleLib, svc, readCountMap, svs, juncMap, spanMap);
  countMap.resize(c.files.size());
  for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
    countMap[file_c].resize(svs.size());
    for (uint32_t id = 0; id < svs.size(); ++id) {
      if ((c.indels) && (svSize[id] <= c.indelsize)) {
	countMap[file_c][id].rc = readCountMap[file_c][id].first;
	countMap[file_c][id].leftRC = readCountMap[file_c][id + lastId].first;
	countMap[file_c][id].rightRC = readCountMap[file_c][id + 2*lastId].first;
      } else {
	countMap[file_c][id].rc = readCountMap[file_c][id].second;
	countMap[file_c][id].leftRC = readCountMap[file_c][id + lastId].second;
	countMap[file_c][id].rightRC = readCountMap[file_c][id + 2*lastId].second;
      }
    }
  }
}

template<typename TConfigStruct>
inline int dellyRun(TConfigStruct& c) {
#ifdef PROFILE
  ProfilerStart("delly.prof");
#endif

  // Collect all promising structural variants
  typedef std::vector<StructuralVariantRecord> TVariants;
  TVariants svs;

  // Open header
  samFile* samfile = sam_open(c.files[0].string().c_str(), "r");
  bam_hdr_t* hdr = sam_hdr_read(samfile);

  // Exclude intervals
  typedef boost::icl::interval_set<uint32_t> TChrIntervals;
  typedef typename TChrIntervals::interval_type TIVal;
  typedef std::vector<TChrIntervals> TRegionsGenome;
  TRegionsGenome validRegions;
  validRegions.resize(hdr->n_targets);
  TRegionsGenome exclg;
  exclg.resize(hdr->n_targets);
  std::vector<bool> validChr;
  validChr.resize(hdr->n_targets, true);
  if (c.hasExcludeFile) {
    std::ifstream chrFile(c.exclude.string().c_str(), std::ifstream::in);
    if (chrFile.is_open()) {
      while (chrFile.good()) {
	std::string chrFromFile;
	getline(chrFile, chrFromFile);
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
	      exclg[tid].insert(TIVal::right_open(start, end));
	    } else validChr[tid] = false; // Exclude entire chromosome
	  }
	}
      }
      chrFile.close();
    }
  }
  // Create the valid regions
  for (int32_t i = 0; i<hdr->n_targets; ++i) {
    if (!validChr[i]) continue;
    uint32_t istart = 0;
    for(typename TChrIntervals::iterator it = exclg[i].begin(); it != exclg[i].end(); ++it) {
      if (istart + 1 < it->lower()) validRegions[i].insert(TIVal::right_open(istart, it->lower() - 1));
      istart = it->upper();
    }
    if (istart + 1 < hdr->target_len[i]) validRegions[i].insert(TIVal::right_open(istart, hdr->target_len[i]));
  }
  exclg.clear();

  // Debug code
  //for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
  //for(typename TChrIntervals::const_iterator vRIt = validRegions[refIndex].begin(); vRIt != validRegions[refIndex].end(); ++vRIt) {
  //std::cerr << std::string(hdr->target_name[refIndex]) << "\t" << vRIt->lower() << "\t" << vRIt->upper() << std::endl;
  //}
  //}
  
  // Create library objects
  typedef boost::unordered_map<std::string, LibraryInfo> TLibraryMap;
  typedef std::vector<TLibraryMap> TSampleLibrary;
  TSampleLibrary sampleLib(c.files.size());
  if (!getLibraryParams(c, validRegions, sampleLib)) {
    std::cerr << "Library parameters could not be estimated!" << std::endl;
    bam_hdr_destroy(hdr);
    sam_close(samfile);
    return 1;
  }

  // SV Discovery
  if (!c.hasVcfFile) shortPE(c, validRegions, svs, sampleLib);
  else vcfParse(c, hdr, svs);

  // Re-number SVs
  sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());    
  uint32_t cliqueCount = 0;
  for(typename TVariants::iterator svIt = svs.begin(); svIt != svs.end(); ++svIt) svIt->id = cliqueCount++;

  // Annotate junction reads
  typedef std::vector<JunctionCount> TSVJunctionMap;
  typedef std::vector<TSVJunctionMap> TSampleSVJunctionMap;
  TSampleSVJunctionMap junctionCountMap;

  // Annotate spanning coverage
  typedef std::vector<SpanningCount> TSVSpanningMap;
  typedef std::vector<TSVSpanningMap> TSampleSVSpanningMap;
  TSampleSVSpanningMap spanCountMap;

  // Annotate coverage
  typedef std::vector<ReadCount> TSVReadCount;
  typedef std::vector<TSVReadCount> TSampleSVReadCount;
  TSampleSVReadCount rcMap;

  // SV Genotyping
  if (!svs.empty()) _annotateCoverage(c, hdr, sampleLib, svs, rcMap, junctionCountMap, spanCountMap);
  
  // VCF output
  vcfOutput(c, svs, junctionCountMap, rcMap, spanCountMap);

  // Clean-up
  bam_hdr_destroy(hdr);
  sam_close(samfile);

  // Output library statistics
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Library statistics" << std::endl;
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    std::cout << "Sample: " << c.sampleName[file_c] << std::endl;
    for(TLibraryMap::const_iterator libIt = sampleLib[file_c].begin(); libIt != sampleLib[file_c].end(); ++libIt) {
      std::cout << "RG: ID=" << libIt->first << ",ReadSize=" << libIt->second.rs << ",Median=" << libIt->second.median << ",MAD=" << libIt->second.mad << ",UniqueDiscordantPairs=" << libIt->second.abnormal_pairs << std::endl;
    }
  }
  
#ifdef PROFILE
  ProfilerStop();
#endif
  
  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
  return 0;
}


int delly(int argc, char **argv) {
  Config c;
  c.isHaplotagged = false;

  // Define generic options
  std::string svtype;
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("svtype,t", boost::program_options::value<std::string>(&svtype)->default_value("ALL"), "SV type to compute [DEL, INS, DUP, INV, BND, ALL]")
    ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
    ("exclude,x", boost::program_options::value<boost::filesystem::path>(&c.exclude), "file with regions to exclude")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("sv.bcf"), "SV BCF output file")
    ;

  boost::program_options::options_description disc("Discovery options");
  disc.add_options()
    ("map-qual,q", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(1), "min. paired-end (PE) mapping quality")
    ("qual-tra,r", boost::program_options::value<unsigned short>(&c.minTraQual)->default_value(20), "min. PE quality for translocation")
    ("mad-cutoff,s", boost::program_options::value<unsigned short>(&c.madCutoff)->default_value(9), "insert size cutoff, median+s*MAD (deletions only)")
    ("noindels,n", "no small InDel calling")
    ;

  boost::program_options::options_description geno("Genotyping options");
  geno.add_options()
    ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input VCF/BCF file for genotyping")
    ("geno-qual,u", boost::program_options::value<unsigned short>(&c.minGenoQual)->default_value(5), "min. mapping quality for genotyping")
    ("dump,d", boost::program_options::value<boost::filesystem::path>(&c.srpedump), "gzipped output file for SV-reads (optional)")
    ;

  // Define hidden options
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
    ("pruning,j", boost::program_options::value<uint32_t>(&c.graphPruning)->default_value(1000), "PE graph pruning cutoff")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  // Set the visibility
  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(disc).add(geno).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic).add(disc).add(geno);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);


  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) { 
    std::cout << std::endl;
    std::cout << "Usage: delly " << argv[0] << " [OPTIONS] -g <ref.fa> <sample1.sort.bam> <sample2.sort.bam> ..." << std::endl;
    std::cout << visible_options << "\n";
    return 0;
  }

  // Only one SV type to compute?
  c.svtcmd = false;
  if (vm.count("svtype")) {
    c.svtcmd = true;
    if (svtype == "DEL") {
      c.svtset.insert(2);
    } else if (svtype == "INS") {
      c.svtset.insert(4);
    } else if (svtype == "DUP") {
      c.svtset.insert(3);
    } else if (svtype == "INV") {
      c.svtset.insert(0);
      c.svtset.insert(1);
    } else if (svtype == "BND") {
      c.svtset.insert(DELLY_SVT_TRANS + 0);
      c.svtset.insert(DELLY_SVT_TRANS + 1);
      c.svtset.insert(DELLY_SVT_TRANS + 2);
      c.svtset.insert(DELLY_SVT_TRANS + 3);
    } else {
      c.svtcmd = false;
    }
  }
  
  // Dump PE and SR support?
  if (vm.count("dump")) c.dumpflag = true;
  else c.dumpflag = false;

  // Check quality cuts
  if (c.minMapQual > c.minTraQual) c.minTraQual = c.minMapQual;
  
  // Check reference
  if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
    std::cerr << "Reference file is missing: " << c.genome.string() << std::endl;
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

  // Check exclude file
  if (vm.count("exclude")) {
    if (!(boost::filesystem::exists(c.exclude) && boost::filesystem::is_regular_file(c.exclude) && boost::filesystem::file_size(c.exclude))) {
      std::cerr << "Exclude file is missing: " << c.exclude.string() << std::endl;
      return 1;
    }
    c.hasExcludeFile = true;
  } else c.hasExcludeFile = false;

  // Check input VCF file
  if (vm.count("vcffile")) {
    if (!(boost::filesystem::exists(c.vcffile) && boost::filesystem::is_regular_file(c.vcffile) && boost::filesystem::file_size(c.vcffile))) {
      std::cerr << "Input VCF/BCF file is missing: " << c.vcffile.string() << std::endl;
      return 1;
    }
    htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
    if (ifile == NULL) {
      std::cerr << "Fail to open file " << c.vcffile.string() << std::endl;
      return 1;
    }
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);
    if (hdr == NULL) {
      std::cerr << "Fail to open index file " << c.vcffile.string() << std::endl;
      return 1;
    }
    bcf_hdr_destroy(hdr);
    bcf_close(ifile);
    c.hasVcfFile = true;
  } else c.hasVcfFile = false;

  // Check output directory
  try {
    boost::filesystem::path outdir;
    if (c.outfile.has_parent_path()) outdir = c.outfile.parent_path();
    else outdir = boost::filesystem::current_path();
    if (!boost::filesystem::exists(outdir)) {
      std::cerr << "Output directory does not exist: " << outdir << std::endl;
      return 1;
    } else {
      boost::filesystem::file_status s = boost::filesystem::status(outdir);
      boost::filesystem::ofstream file(c.outfile.string());
      file.close();
      if (!(boost::filesystem::exists(c.outfile) && boost::filesystem::is_regular_file(c.outfile))) {
	std::cerr << "Fail to open output file " << c.outfile.string() << std::endl;
	std::cerr << "Output directory permissions: " << s.permissions() << std::endl;
	return 1;
      } else {
	boost::filesystem::remove(c.outfile.string());
      }
    }
  } catch (boost::filesystem::filesystem_error const& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  
  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  std::cout << "delly ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Always ignore reads of mapping quality <5 for genotyping, otherwise het. is more likely!
  if (c.minGenoQual<5) c.minGenoQual=5;

  // Small InDels?
  if (vm.count("noindels")) c.indels = false;
  else c.indels = true;
  
  // Run main program
  c.aliscore = DnaScore<int>(5, -4, -10, -1);
  c.flankQuality = 0.95;
  c.minimumFlankSize = 13;
  c.indelsize = 500;
  return dellyRun(c);
}

}

#endif
