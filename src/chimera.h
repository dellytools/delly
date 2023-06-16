#ifndef CHIMERA_H
#define CHIMERA_H

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
#include <boost/icl/interval_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>


#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <stdio.h>

#include "edlib.h"
#include "tags.h"
#include "version.h"
#include "util.h"
#include "modvcf.h"

namespace torali
{

  struct ChimeraConfig {
    int32_t minSize;
    int32_t format;
    float divergence;
    boost::filesystem::path genome;
    boost::filesystem::path file;
  };

  inline int32_t     // -1: failure, 0: bam, 1: fasta, 2: fastq
  inputType(std::string const& path) {
    htsFile *hts_fp = hts_open(path.c_str(), "r");
    if (hts_fp == NULL) return -1;
    else {
      std::string ext = std::string(hts_format_file_extension(hts_get_format(hts_fp)));
      hts_close(hts_fp);
      if ((ext == "bam") || (ext == "cram")) return 0;
      else if (ext == "fa") return 1;
      else if (ext == "fq") return 2;
      else {
	std::cerr << ext << std::endl;
	return -1;
      }
    }
  }

  inline double
  computeChimeraScore(std::string const& sequence) {
    // Split into prefix and suffix
    int32_t halfPoint = sequence.size() / 2;
    std::string prefix = sequence.substr(0, halfPoint);
    std::string suffix = sequence.substr(halfPoint);
    reverseComplement(suffix);

    // Align
    EdlibAlignResult cigar = edlibAlign(prefix.c_str(), prefix.size(), suffix.c_str(), suffix.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
    double score = (double) cigar.editDistance / (double) halfPoint;
    //printAlignment(prefix, suffix, EDLIB_MODE_NW, cigar);
    edlibFreeAlignResult(cigar);
    
    return score;
  }

  inline double
  fastqInputChimera(ChimeraConfig& c) {                                                                                                     
    std::ifstream file(c.file.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    std::string header;
    std::string seq;
    uint64_t lcount = 0;
    uint64_t chimera = 0;
    uint64_t total = 0;
    while(std::getline(instream, gline)) {
      if (lcount % 4 == 0) header = gline;
      else if (lcount % 4 == 1) seq = gline;
      else if (lcount % 4 == 3) {
	double score = computeChimeraScore(seq);
	if (score < c.divergence) {
	  std::cout << header.substr(0, 60) << std::endl;
	  ++chimera;
	}
	++total;
      }
      ++lcount;
    }
    dataIn.pop();
    dataIn.pop();

    return (double) chimera / (double) total;
  }

  inline double
  bamInputChimera(ChimeraConfig& c) {
    // Open file handles
    samFile* samfile = sam_open(c.file.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.file.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    uint64_t chimera = 0;
    uint64_t total = 0;
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Processing... " << hdr->target_name[refIndex] << std::endl;

      // Iterate alignments
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;

	if (rec->core.l_qseq > c.minSize) {
	  // Load sequence
	  std::string sequence;
	  sequence.resize(rec->core.l_qseq, 'N');
	  uint8_t* seqptr = bam_get_seq(rec);
	  for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	  double score = computeChimeraScore(sequence);
	  if (score < c.divergence) {
	    std::cout << bam_get_qname(rec) << std::endl;
	    ++chimera;
	  }
	  ++total;
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);

    return (double) chimera / (double) total;
  }

	  
  inline int
  chimeraRun(ChimeraConfig& c) {
    double chimeraRate = 0;
    if (c.format == 0) chimeraRate = bamInputChimera(c);
    else chimeraRate = fastqInputChimera(c);

    // Chimera fraction
    std::cerr << "Chimera: " << chimeraRate << std::endl;
    
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    return 0;
  }


  int chimera(int argc, char **argv) {
    ChimeraConfig c;
    
    // Define generic options
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("minsize,m", boost::program_options::value<int32_t>(&c.minSize)->default_value(5000), "min. sequence size")
      ("divergence,d", boost::program_options::value<float>(&c.divergence)->default_value(0.05), "max. sequence divergence")
      ;
    
    // Define hidden options
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.file), "input BAM file")
      ;
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    // Set the visibility
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
      std::cerr << std::endl;
      std::cerr << "Usage: delly " << argv[0] << " [OPTIONS] -g <ref.fa> <input.bam>" << std::endl;
      std::cerr << visible_options << "\n";
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

    // Input format
    c.format = inputType(c.file.string());
    if (c.format == -1) {
      std::cerr << "Unknown input file format!" << std::endl;
      return 1;
    } else if (c.format == 0) {
      // Check BAM file
      if (!(boost::filesystem::exists(c.file) && boost::filesystem::is_regular_file(c.file) && boost::filesystem::file_size(c.file))) {
	std::cerr << "Alignment file is missing: " << c.file.string() << std::endl;
	return 1;
      }
      samFile* samfile = sam_open(c.file.string().c_str(), "r");
      if (samfile == NULL) {
	std::cerr << "Fail to open file " << c.file.string() << std::endl;
	return 1;
      }
      hts_idx_t* idx = sam_index_load(samfile, c.file.string().c_str());
      if (idx == NULL) {
	std::cerr << "Fail to open index for " << c.file.string() << std::endl;
	return 1;
      }
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.file.string() << std::endl;
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
      bam_hdr_destroy(hdr);
      hts_idx_destroy(idx);
      sam_close(samfile);
    }
      
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cerr << "delly ";
    for(int i=0; i<argc; ++i) { std::cerr << argv[i] << ' '; }
    std::cerr << std::endl;
    
    // Run comparison
    return chimeraRun(c);
  }
  
}

#endif
