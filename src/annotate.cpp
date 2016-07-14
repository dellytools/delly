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

#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>

#define BOOST_DISABLE_ASSERTS
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
#include <htslib/kseq.h>

#ifdef OPENMP
#include <omp.h>
#endif

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <stdio.h>

#include "tags.h"
#include "version.h"
#include "util.h"
#include "modvcf.h"
#include "junction.h"

using namespace torali;

struct ConfigAnnotate {
  int32_t maxlen;
  int32_t minimumFlankSize;
  float flankQuality;
  DnaScore<int> aliscore;
  std::string technology;
  std::string svType;
  boost::filesystem::path genome;
  boost::filesystem::path outfile;
  boost::filesystem::path infile;
};


template<typename TTag>
inline bool
_addRefAltAllele(SVType<TTag>) {
  return false;
}

inline bool
_addRefAltAllele(SVType<InsertionTag>) {
  return true;
}

inline bool
_addRefAltAllele(SVType<DeletionTag>) {
  return true;
}


template<typename TSVType>
inline int
runAnnotate(ConfigAnnotate const& c, TSVType svType)
{
    // Load bcf file
  htsFile* ifile = hts_open(c.infile.string().c_str(), "r");
  hts_idx_t* bcfidx = bcf_index_load(c.infile.string().c_str());
  bcf_hdr_t* hdr = bcf_hdr_read(ifile);

  // Open output VCF file
  htsFile *ofile = hts_open(c.outfile.string().c_str(), "wb");
  bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);
  bcf_hdr_remove(hdr_out, BCF_HL_INFO, "END");
  bcf_hdr_append(hdr_out, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">");
  bcf_hdr_remove(hdr_out, BCF_HL_INFO, "INSLEN");
  bcf_hdr_append(hdr_out, "##INFO=<ID=INSLEN,Number=1,Type=Integer,Description=\"Predicted length of the insertion\">");
  bcf_hdr_remove(hdr_out, BCF_HL_INFO, "SRQ");
  bcf_hdr_append(hdr_out, "##INFO=<ID=SRQ,Number=1,Type=Float,Description=\"Split-read consensus alignment quality\">");
  bcf_hdr_remove(hdr_out, BCF_HL_INFO, "CE");
  bcf_hdr_append(hdr_out, "##INFO=<ID=CE,Number=1,Type=Float,Description=\"Consensus sequence entropy\">");
  bcf_hdr_remove(hdr_out, BCF_HL_INFO, "HOMLEN");
  bcf_hdr_append(hdr_out, "##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description=\"Predicted microhomology length using a max. edit distance of 2\">");
  bcf_hdr_write(ofile, hdr_out);

  // VCF fields
  int32_t nsvend = 0;
  int32_t* svend = NULL;
  int32_t nsvt = 0;
  char* svt = NULL;
  int32_t ncons = 0;
  char* cons = NULL;
  
  // Get #sequences
  const char **seqnames = NULL;
  int nseq = 0;
  seqnames = bcf_hdr_seqnames(hdr, &nseq);
  if (seqnames != NULL) free(seqnames);

  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Annotating BCF file" << std::endl;
  boost::progress_display show_progress(nseq);

  // Parse genome
  kseq_t *seq;
  int l;
  gzFile fp = gzopen(c.genome.string().c_str(), "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    std::string seqname(seq->name.s);
    int32_t chrid = bcf_hdr_name2id(hdr, seqname.c_str());
    if (chrid < 0) continue;
    ++show_progress;
    hts_itr_t* itervcf = bcf_itr_queryi(bcfidx, chrid, 0, seq->seq.l);
    if (itervcf != NULL) {
      bcf1_t* rec = bcf_init1();
      while (true) {
	int32_t ret = bcf_itr_next(ifile, itervcf, rec);
	if (ret < 0) break;
	bcf_unpack(rec, BCF_UN_INFO);

	// Check SV type
	bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt);
	if ((svt != NULL) && (std::string(svt) != _addID(svType))) continue;

	// Check size
	bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend);
	int32_t svlen = 1;
	if (svend != NULL) svlen = *svend - rec->pos;

	// Check precise
	bool precise = false;
	if (bcf_get_info_flag(hdr, rec, "PRECISE", 0, 0) > 0) precise = true;

	// Try re-alignment
	bool useTags = true;
	if ((svlen <= c.maxlen) && (precise) && (_addRefAltAllele(svType))) {
	  useTags = false;
	  std::string chrName(bcf_hdr_id2name(hdr, rec->rid));
	    
	  // Get consensus sequence
	  bcf_get_info_string(hdr, rec, "CONSENSUS", &cons, &ncons);
	  std::string consensus = boost::to_upper_copy(std::string(cons));
	  int32_t consLen = consensus.size();

	  // Get reference sequence
	  int32_t regStart = std::max(rec->pos - consLen, 0);
	  int32_t regEnd = std::min((uint32_t) *svend + consLen, (uint32_t) seq->seq.l);
	  std::string svRefStr = boost::to_upper_copy(std::string(seq->seq.s + regStart, seq->seq.s + regEnd));

	  // Find breakpoint to reference
	  typedef boost::multi_array<char, 2> TAlign;
	  TAlign align;
	  if (!_consRefAlignment(consensus, svRefStr, align, svType)) useTags = true;
	  else {
	    AlignDescriptor ad;
	    if (!_findSplit(c, consensus, svRefStr, align, ad, svType)) useTags = true;
	    else {
	      std::string precChar = svRefStr.substr(ad.rStart - 1, 1);
	      std::string refPart = precChar;
	      if (ad.rEnd > ad.rStart + 1) refPart += svRefStr.substr(ad.rStart, (ad.rEnd - ad.rStart) - 1);
	      std::string altPart = precChar;
	      if (ad.cEnd > ad.cStart + 1) altPart += consensus.substr(ad.cStart, (ad.cEnd - ad.cStart) - 1);
	    
	      // Update BCF record
	      rec->pos = regStart + ad.rStart - 1;
	      std::string alleles;
	      alleles += refPart + "," + altPart;
	      bcf_update_alleles_str(hdr_out, rec, alleles.c_str());
	      int32_t tmpi = regStart + ad.rEnd;
	      bcf_update_info_int32(hdr_out, rec, "END", &tmpi, 1);
	      tmpi = ad.cEnd - ad.cStart - 1;
	      bcf_update_info_int32(hdr_out, rec, "INSLEN", &tmpi, 1);
	      tmpi = std::max(0, ad.homLeft + ad.homRight - 2);
	      bcf_update_info_int32(hdr_out, rec, "HOMLEN", &tmpi, 1);
	      float tmpf = ad.percId;
	      bcf_update_info_float(hdr_out, rec, "SRQ", &tmpf, 1);
	      tmpf = entropy(consensus);
	      bcf_update_info_float(hdr_out, rec, "CE", &tmpf, 1);
	    }
	  }
	}
	if (useTags) {
	  // Use tags
	  std::string alleles;
	  std::string refAllele = boost::to_upper_copy(std::string(seq->seq.s + rec->pos, seq->seq.s + rec->pos + 1));
	  alleles += refAllele + ",<" + _addID(svType) + ">";
	  bcf_update_alleles_str(hdr_out, rec, alleles.c_str());
	}
	bcf_write1(ofile, hdr_out, rec);
      }
      bcf_destroy(rec);
    }
    hts_itr_destroy(itervcf);
  }
  kseq_destroy(seq);
  gzclose(fp);

  // Clean-up
  if (svend != NULL) free(svend);
  if (svt != NULL) free(svt);
  if (cons != NULL) free(cons);
  
  // Close output VCF
  bcf_hdr_destroy(hdr_out);
  hts_close(ofile);

  // Build index
  bcf_index_build(c.outfile.string().c_str(), 14);

  // Close VCF
  bcf_hdr_destroy(hdr);
  hts_idx_destroy(bcfidx);
  bcf_close(ifile);

  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;

  return 0;
}


int main(int argc, char **argv) {
  ConfigAnnotate c;

  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("type,t", boost::program_options::value<std::string>(&c.svType)->default_value("DEL"), "SV type (DEL, DUP, INV, INS)")
    ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "Genomic reference file")
    ("maxlen,m", boost::program_options::value<int32_t>(&c.maxlen)->default_value(500), "max. SV size before tags (<DEL>) are used")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.bcf"), "output BCF file")
    ("technology,e", boost::program_options::value<std::string>(&c.technology)->default_value("illumina"), "technology (illumina, pacbio)")
    ;

  // Define hidden options
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("infile", boost::program_options::value<boost::filesystem::path>(&c.infile), "input file")
    ("license,l", "show license")
    ("warranty,w", "show warranty")
    ;
  boost::program_options::positional_options_description pos_args;
  pos_args.add("infile", -1);

  // Set the visibility
  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);


  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("infile"))) { 
    printTitle("SV annotation");
    if (vm.count("warranty")) {
      displayWarranty();
    } else if (vm.count("license")) {
      gplV3();
    } else {
      std::cout << "Usage: " << argv[0] << " [OPTIONS] <input.bcf>" << std::endl;
      std::cout << visible_options << "\n"; 
    }
    return 1; 
  }

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Check reference
  if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
    std::cerr << "Reference file is missing: " << c.genome.string() << std::endl;
    return 1;
  }

  // Check technology
  if (c.technology == "illumina") {
    c.aliscore = DnaScore<int>(5, -4, -10, -1);
    c.flankQuality = 0.95;
    c.minimumFlankSize = 13;
  } else if (c.technology == "pacbio") {
    c.aliscore = DnaScore<int>(5, -4, -2, -1);
    c.flankQuality = 0.8;
    c.minimumFlankSize = 13;
  }
  
  // Run SV annotation
  if (c.svType == "DEL") return runAnnotate(c, SVType<DeletionTag>());
  else if (c.svType == "INS") return runAnnotate(c, SVType<InsertionTag>());
  else if (c.svType == "DUP") return runAnnotate(c, SVType<DuplicationTag>());
  else if (c.svType == "INV") return runAnnotate(c, SVType<InversionTag>());
  else if (c.svType == "TRA") return runAnnotate(c, SVType<TranslocationTag>());
  else {
    std::cerr << "SV analysis type not supported " << c.svType << std::endl;
    return 1;
  }
}
