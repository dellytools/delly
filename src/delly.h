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
    uint16_t minMapQual;
    uint16_t minTraQual;
    uint16_t minGenoQual;
    uint16_t madCutoff;
    int32_t nchr;
    int32_t minimumFlankSize;
    int32_t indelsize;
    uint32_t graphPruning;
    uint32_t minRefSep;
    uint32_t maxReadSep;
    uint32_t minClip;
    float flankQuality;
    bool hasExcludeFile;
    bool hasVcfFile;
    bool isHaplotagged;
    bool hasDumpFile;
    bool svtcmd;
    std::set<int32_t> svtset;
    DnaScore<int> aliscore;
    boost::filesystem::path outfile;
    boost::filesystem::path vcffile;
    boost::filesystem::path genome;
    boost::filesystem::path exclude;
    boost::filesystem::path dumpfile;
    std::vector<boost::filesystem::path> files;
    std::vector<std::string> sampleName;
  };


  template<typename TConfig, typename TRegionsGenome>
  inline int32_t
  _parseExcludeIntervals(TConfig const& c, bam_hdr_t* hdr, TRegionsGenome& validRegions) {
    typedef typename TRegionsGenome::value_type TChrIntervals;
    typedef typename TChrIntervals::interval_type TIVal;
    
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
		int32_t start = 0;
		try {
		  start = boost::lexical_cast<int32_t>(*tokIter++);
		} catch (boost::bad_lexical_cast&) {
		  std::cerr << "Exclude file needs to be in tab-delimited format: chr, start, end" << std::endl;
		  std::cerr << "Offending line: " << chrFromFile << std::endl;
		  return false;
		}
		if (tokIter!=tokens.end()) {
		  int32_t end = start + 1;
		  try {
		    end = boost::lexical_cast<int32_t>(*tokIter++);
		  } catch (boost::bad_lexical_cast&) {
		    std::cerr << "Exclude file needs to be in tab-delimited format: chr, start, end" << std::endl;
		    std::cerr << "Offending line: " << chrFromFile << std::endl;
		    return false;
		  }
		  if (start < end) {
		    exclg[tid].insert(TIVal::right_open(start, end));
		  } else {
		    std::cerr << "Exclude file needs to be in tab-delimited format (chr, start, end) and start < end." << std::endl;
		    std::cerr << "Offending line: " << chrFromFile << std::endl;
		    return false;
		  }
		} else {
		  std::cerr << "Exclude file needs to be in tab-delimited format: chr, start, end" << std::endl;
		  std::cerr << "Offending line: " << chrFromFile << std::endl;
		  return false;
		}
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
    return true;
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
    typedef std::vector<TChrIntervals> TRegionsGenome;
    TRegionsGenome validRegions;
    if (!_parseExcludeIntervals(c, hdr, validRegions)) {
      std::cerr << "Delly couldn't parse exclude intervals!" << std::endl;
      bam_hdr_destroy(hdr);
      sam_close(samfile);
      return 1;
    }
    
    // Debug code
    //for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
    //for(typename TChrIntervals::const_iterator vRIt = validRegions[refIndex].begin(); vRIt != validRegions[refIndex].end(); ++vRIt) {
    //std::cerr << std::string(hdr->target_name[refIndex]) << "\t" << vRIt->lower() << "\t" << vRIt->upper() << std::endl;
    //}
    //}
    
    // Create library objects
    typedef std::vector<LibraryInfo> TSampleLibrary;
    TSampleLibrary sampleLib(c.files.size(), LibraryInfo());
    getLibraryParams(c, validRegions, sampleLib);
    for(uint32_t i = 0; i<sampleLib.size(); ++i) {
      if (sampleLib[i].rs == 0) {
	std::cerr << "Sample has not enough data to estimate library parameters! File: " << c.files[i].string() << std::endl;
	bam_hdr_destroy(hdr);
	sam_close(samfile);
	return 1;
      }
    }
    
    // SV Discovery
    if (!c.hasVcfFile) {
      // Split-read SVs
      typedef std::vector<StructuralVariantRecord> TVariants;
      TVariants srSVs;
      
      // SR Store
      {
	typedef std::pair<int32_t, std::size_t> TPosRead;
	typedef boost::unordered_map<TPosRead, int32_t> TPosReadSV;
	typedef std::vector<TPosReadSV> TGenomicPosReadSV;
	TGenomicPosReadSV srStore(c.nchr, TPosReadSV());
	scanPEandSR(c, validRegions, svs, srSVs, srStore, sampleLib);
	
	// Assemble split-read calls
	assembleSplitReads(c, validRegions, srStore, srSVs);
      }

      // Sort and merge PE and SR calls
      mergeSort(svs, srSVs);
    } else vcfParse(c, hdr, svs);
    // Clean-up
    bam_hdr_destroy(hdr);
    sam_close(samfile);

    // Re-number SVs
    sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());    
    uint32_t cliqueCount = 0;
    for(typename TVariants::iterator svIt = svs.begin(); svIt != svs.end(); ++svIt, ++cliqueCount) svIt->id = cliqueCount;
    
    // Annotate junction reads
    typedef std::vector<JunctionCount> TSVJunctionMap;
    typedef std::vector<TSVJunctionMap> TSampleSVJunctionMap;
    TSampleSVJunctionMap jctMap;
    
    // Annotate spanning coverage
    typedef std::vector<SpanningCount> TSVSpanningMap;
    typedef std::vector<TSVSpanningMap> TSampleSVSpanningMap;
    TSampleSVSpanningMap spanMap;

    // Annotate coverage
    typedef std::vector<ReadCount> TSVReadCount;
    typedef std::vector<TSVReadCount> TSampleSVReadCount;
    TSampleSVReadCount rcMap;
    
    // SV Genotyping
    if (!svs.empty()) annotateCoverage(c, sampleLib, svs, rcMap, jctMap, spanMap);
    
    // VCF output
    vcfOutput(c, svs, jctMap, rcMap, spanMap);
    
    // Output library statistics
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Sample statistics" << std::endl;
    for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
      std::cout << "Sample:" << c.sampleName[file_c] << ",ReadSize=" << sampleLib[file_c].rs << ",Median=" << sampleLib[file_c].median << ",MAD=" << sampleLib[file_c].mad << ",UniqueDiscordantPairs=" << sampleLib[file_c].abnormal_pairs << std::endl;
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
      ("map-qual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(1), "min. paired-end (PE) mapping quality")
      ("qual-tra,r", boost::program_options::value<uint16_t>(&c.minTraQual)->default_value(20), "min. PE quality for translocation")
      ("mad-cutoff,s", boost::program_options::value<uint16_t>(&c.madCutoff)->default_value(9), "insert size cutoff, median+s*MAD (deletions only)")
      ("minclip,c", boost::program_options::value<uint32_t>(&c.minClip)->default_value(25), "min. clipping length")
      ("minrefsep,m", boost::program_options::value<uint32_t>(&c.minRefSep)->default_value(50), "min. reference separation")
      ("maxreadsep,n", boost::program_options::value<uint32_t>(&c.maxReadSep)->default_value(10), "max. read separation")
      ;
    
    boost::program_options::options_description geno("Genotyping options");
    geno.add_options()
      ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input VCF/BCF file for genotyping")
      ("geno-qual,u", boost::program_options::value<uint16_t>(&c.minGenoQual)->default_value(5), "min. mapping quality for genotyping")
      ("dump,d", boost::program_options::value<boost::filesystem::path>(&c.dumpfile), "gzipped output file for SV-reads (optional)")
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
    
    // SV types to compute?
    _svTypesToCompute(c, svtype, vm.count("svtype"));
    
    // Dump PE and SR support?
    if (vm.count("dump")) c.hasDumpFile = true;
    else c.hasDumpFile = false;
    
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
    c.nchr = 0;
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
      if (!c.nchr) c.nchr = hdr->n_targets;
      else {
	if (c.nchr != hdr->n_targets) {
	  std::cerr << "BAM files have different number of chromosomes!" << std::endl;
	  return 1;
	}
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
      std::string sampleName = "unknown";
      getSMTag(std::string(hdr->text), c.files[file_c].stem().string(), sampleName);
      c.sampleName[file_c] = sampleName;
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
    if (!_outfileValid(c.outfile)) return 1;
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "delly ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    // Always ignore reads of mapping quality <5 for genotyping, otherwise het. is more likely!
    if (c.minGenoQual<5) c.minGenoQual=5;
    
    // Run main program
    c.aliscore = DnaScore<int>(5, -4, -10, -1);
    c.flankQuality = 0.95;
    c.minimumFlankSize = 13;
    c.indelsize = 500;
    return dellyRun(c);
  }

}

#endif
