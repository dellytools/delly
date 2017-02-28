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
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/sam.h>

#include "version.h"
#include "util.h"
#include "bolog.h"
#include "tags.h"
#include "spanning.h"
#include "coverage.h"
#include "junction.h"
#include "msa.h"
#include "split.h"
#include "json.h"
#include "shortpe.h"

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
  DnaScore<int> aliscore;
  std::string svType;
  std::string format;
  boost::filesystem::path outfile;
  boost::filesystem::path vcffile;
  boost::filesystem::path genome;
  boost::filesystem::path exclude;
  std::vector<boost::filesystem::path> files;
  std::vector<std::string> sampleName;
};



template<typename TSize>
inline bool
_svSizeCheck(TSize const s, TSize const e, SVType<DeletionTag>) {
  return (( e - s ) >= 300);
}

template<typename TSize>
inline bool
_svSizeCheck(TSize const s, TSize const e, SVType<DuplicationTag>) {
  return (( e - s ) >= 100);
}

template<typename TSize>
inline bool
_svSizeCheck(TSize const s, TSize const e, SVType<InversionTag>) {
  return (( e - s ) >= 100);
}

template<typename TSize>
inline bool
_svSizeCheck(TSize const s, TSize const e, SVType<InsertionTag>) {
  return (( e - s ) >= 0);
}

template<typename TSize>
inline bool
_svSizeCheck(TSize const, TSize const, SVType<TranslocationTag>) {
  return true;
}

 
template<typename TConfig, typename TSVs, typename TCountMap, typename TTag>
inline void
_annotateCoverage(TConfig const& c, bam_hdr_t* hdr, TSVs& svs, TCountMap& countMap, SVType<TTag>) 
{
  // Find Ns in the reference genome
  typedef boost::icl::interval_set<int> TNIntervals;
  typedef std::vector<TNIntervals> TNGenome;
  TNGenome ni( hdr->n_targets );

  // Find valid chromosomes
  typedef std::vector<bool> TValidChr;
  TValidChr validChr(hdr->n_targets, false);
  for (typename TSVs::const_iterator itSV = svs.begin(); itSV != svs.end(); ++itSV) validChr[itSV->chr] = true;

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
  TCovRecord svc;
  uint32_t lastId = svs.size();
  typedef std::vector<int32_t> TSVSize;
  TSVSize svSize(lastId);
  for (typename TSVs::const_iterator itSV = svs.begin(); itSV != svs.end(); ++itSV) {
    int halfSize = (itSV->svEnd - itSV->svStart)/2;

    // Left control region
    CovRecord sLeft;
    sLeft.chr = itSV->chr;
    sLeft.id = lastId + itSV->id;
    sLeft.svStart = std::max(itSV->svStart - halfSize, 0);
    sLeft.svEnd = itSV->svStart;
    sLeft.peSupport = 0;
    typename TNIntervals::const_iterator itO = ni[sLeft.chr].find(boost::icl::discrete_interval<int>::right_open(sLeft.svStart, sLeft.svEnd));
    while (itO != ni[sLeft.chr].end()) {
      sLeft.svStart = std::max(itO->lower() - halfSize, 0);
      sLeft.svEnd = itO->lower();
      itO = ni[sLeft.chr].find(boost::icl::discrete_interval<int>::right_open(sLeft.svStart, sLeft.svEnd));
    }
    svc.push_back(sLeft);

    // Actual SV
    CovRecord sMiddle;
    sMiddle.chr = itSV->chr;
    sMiddle.id = itSV->id;
    sMiddle.svStart = itSV->svStart;
    sMiddle.svEnd = itSV->svEnd;
    sMiddle.peSupport = itSV->peSupport;
    svc.push_back(sMiddle);
    svSize[itSV->id] = (itSV->svEnd - itSV->svStart);

    // Right control region
    CovRecord sRight;
    sRight.chr = itSV->chr;
    sRight.id = 2 * lastId + itSV->id;
    sRight.svStart = itSV->svEnd;
    sRight.svEnd = itSV->svEnd + halfSize;
    sRight.peSupport = 0;
    itO = ni[sRight.chr].find(boost::icl::discrete_interval<int>::right_open(sRight.svStart, sRight.svEnd));
    while (itO != ni[sRight.chr].end()) {
      sRight.svStart = itO->upper();
      sRight.svEnd = itO->upper() + halfSize;
      itO = ni[sRight.chr].find(boost::icl::discrete_interval<int>::right_open(sRight.svStart, sRight.svEnd));
    }
    svc.push_back(sRight);
    //std::cerr << itSV->id << ':' << sLeft.svStart << '-' << sLeft.svEnd << ',' << itSV->svStart << '-' << itSV->svEnd << ',' << sRight.svStart << '-' << sRight.svEnd << std::endl;
  }
  
  typedef std::pair<int, int> TBpRead;
  typedef std::vector<TBpRead> TSVReadCount;
  typedef std::vector<TSVReadCount> TSampleSVReadCount;
  TSampleSVReadCount readCountMap;
  if (c.indels) annotateCoverage(c.files, c.minGenoQual, svc, readCountMap, BpLevelType<BpLevelCount>());
  else annotateCoverage(c.files, c.minGenoQual, svc, readCountMap, BpLevelType<NoBpLevelCount>());
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

template<typename TConfig, typename TSVs, typename TCountMap>
inline void
_annotateCoverage(TConfig const& c, bam_hdr_t*, TSVs& svs, TCountMap& countMap, SVType<TranslocationTag>) 
{
  countMap.resize(c.files.size());
  for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) countMap[file_c].resize(svs.size());
}

template<typename TConfig, typename TSVs, typename TCountMap>
inline void
_annotateCoverage(TConfig const& c, bam_hdr_t*, TSVs& svs, TCountMap& countMap, SVType<InsertionTag>) 
{
  countMap.resize(c.files.size());
  for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) countMap[file_c].resize(svs.size());
}





template<typename TSVType>
inline int dellyRun(Config& c, TSVType svType) {
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
  getLibraryParams(c, validRegions, sampleLib);

  // SV Discovery
  if (!c.hasVcfFile) shortPE(c, validRegions, svs, sampleLib, svType);
  else {
    // Read SV records from input file
    if (c.format == "json.gz") jsonParse(c, hdr, svs, svType);
    else vcfParse(c, hdr, svs, svType);
  }

  
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
  if (!svs.empty()) {
    annotateJunctionReads(c, svs, junctionCountMap, svType);
    annotateSpanningCoverage(c, sampleLib, svs, spanCountMap, svType);
    _annotateCoverage(c, hdr, svs, rcMap, svType);
  }
  
  // VCF output
  vcfOutput(c, svs, junctionCountMap, rcMap, spanCountMap, svType);

  // Clean-up
  bam_hdr_destroy(hdr);
  sam_close(samfile);

  // Output library statistics
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Library statistics" << std::endl;
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    std::cout << "Sample: " << c.sampleName[file_c] << std::endl;
    for(TLibraryMap::const_iterator libIt = sampleLib[file_c].begin(); libIt != sampleLib[file_c].end(); ++libIt) {
      std::cout << "RG: ID=" << libIt->first << ",ReadSize=" << libIt->second.rs << ",AvgDist=" << libIt->second.avgdist << ",EstCov=" << (double) libIt->second.rs / (double) libIt->second.avgdist << ",MappedAsPair=" << libIt->second.mappedAsPair << ",Median=" << libIt->second.median << ",MAD=" << libIt->second.mad << ",Layout=" << (int) libIt->second.defaultOrient << ",MaxSize=" << libIt->second.maxISizeCutoff << ",MinSize=" << libIt->second.minISizeCutoff << ",UniqueDiscordantPairs=" << libIt->second.abnormal_pairs << std::endl;
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
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("type,t", boost::program_options::value<std::string>(&c.svType)->default_value("DEL"), "SV type (DEL, DUP, INV, TRA, INS)")
    ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
    ("exclude,x", boost::program_options::value<boost::filesystem::path>(&c.exclude), "file with regions to exclude")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("sv.bcf"), "SV BCF output file")
    ;

  boost::program_options::options_description disc("Discovery options");
  disc.add_options()
    ("map-qual,q", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(1), "min. paired-end mapping quality")
    ("mad-cutoff,s", boost::program_options::value<unsigned short>(&c.madCutoff)->default_value(9), "insert size cutoff, median+s*MAD (deletions only)")
    ("noindels,n", "no small InDel calling")
    ;

  boost::program_options::options_description geno("Genotyping options");
  geno.add_options()
    ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input VCF/BCF file for re-genotyping")
    ("geno-qual,u", boost::program_options::value<unsigned short>(&c.minGenoQual)->default_value(5), "min. mapping quality for genotyping")
    ;

  // Define hidden options
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
    ("format,f", boost::program_options::value<std::string>(&c.format)->default_value("bcf"), "output format (bcf, json.gz)")
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
    if (c.format != "json.gz") {
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
    }
    c.hasVcfFile = true;
  } else c.hasVcfFile = false;

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  std::cout << "delly ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Always ignore reads of mapping quality <5 for genotyping, otherwise het. is more likely!
  if (c.minGenoQual<5) c.minGenoQual=5;

  // Small InDels?
  c.indels = false;
  if (!vm.count("noindels")) {
    if ((c.svType == "DEL") || (c.svType == "INS")) c.indels = true;
  }

  // Run main program
  c.aliscore = DnaScore<int>(5, -4, -10, -1);
  c.flankQuality = 0.95;
  c.minimumFlankSize = 13;
  c.indelsize = 500;
  if (c.svType == "DEL") return dellyRun(c, SVType<DeletionTag>());
  else if (c.svType == "DUP") return dellyRun(c, SVType<DuplicationTag>());
  else if (c.svType == "INV") return dellyRun(c, SVType<InversionTag>());
  else if (c.svType == "TRA") return dellyRun(c, SVType<TranslocationTag>());
  else if (c.svType == "INS") return dellyRun(c, SVType<InsertionTag>());
  else {
    std::cerr << "SV analysis type not supported by Delly: " << c.svType << std::endl;
    return 1;
  }
}

}

#endif
