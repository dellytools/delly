#ifndef COMPVCF_H
#define COMPVCF_H

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

  struct MinimalSVRecord {
    int32_t tid;
    int32_t svStart;
    int32_t svEnd;
    int32_t svLen;
    int32_t svt;
    int32_t qual;
    std::string id;
    std::vector<int32_t> gt;
  };

  template<typename TSV>
  struct SortMinimalSVRecord : public std::binary_function<TSV, TSV, bool>
  {
    inline bool operator()(TSV const& sv1, TSV const& sv2) {
      return ((sv1.tid<sv2.tid) || ((sv1.tid==sv2.tid) && (sv1.svStart<sv2.svStart)) || ((sv1.tid==sv2.tid) && (sv1.svStart==sv2.svStart) && (sv1.svEnd<sv2.svEnd)));
    }
  };
  
  
  struct CompvcfConfig {
    bool filterForPass;
    int32_t qualthres;
    int32_t bpdiff;
    float sizeratio;
    boost::filesystem::path vcffile;
    boost::filesystem::path base;
  };

  inline bool
  _loadMinimalSVs(CompvcfConfig const& c, std::string const& filename, std::vector<MinimalSVRecord>& allsv) {
    bool success = true;
    std::set<std::string> allIds;
    
    // Load bcf file
    htsFile* ifile = hts_open(filename.c_str(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);

    // VCF fields
    int32_t nsvend = 0;
    int32_t* svend = NULL;
    int32_t svEndVal;
    int32_t nsvt = 0;
    char* svt = NULL;
    std::string svtVal;
    int32_t nct = 0;
    char* ct = NULL;
    std::string ctVal;
    int ngt = 0;
    int32_t* gt = NULL;
    int32_t inslenVal = 0;
    
    // Parse BCF
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Parsing VCF/BCF file" << std::endl;
    bcf1_t* rec = bcf_init1();
    while (bcf_read(ifile, hdr, rec) == 0) {
      bcf_unpack(rec, BCF_UN_INFO);

      // Check SV type
      if (_isKeyPresent(hdr, "SVTYPE")) {
	bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt);
	svtVal = std::string(svt);
	bcf_get_info_string(hdr, rec, "CT", &ct, &nct);
	ctVal = std::string(ct);
      } else {
	std::string refAllele = rec->d.allele[0];
	std::string altAllele = rec->d.allele[1];
	if (refAllele.size() > altAllele.size()) {
	  svtVal = "DEL";
	  ctVal = "3to5";
	  inslenVal = 0;
	} else {
	  svtVal = "INS";
	  ctVal = "NtoN";
	  inslenVal = altAllele.size() - refAllele.size();
	}
      }

      // Check size and PASS
      if (_isKeyPresent(hdr, "END")) {
	bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend);
	svEndVal = *svend;
      } else {
	std::string refAllele = rec->d.allele[0];
	std::string altAllele = rec->d.allele[1];
	if (refAllele.size() > altAllele.size()) svEndVal = rec->pos + (refAllele.size() - altAllele.size());
	else svEndVal = rec->pos + 1;
      }
      bool pass = true;
      if (c.filterForPass) pass = (bcf_has_filter(hdr, rec, const_cast<char*>("PASS"))==1);
      if ((rec->qual >= c.qualthres) && (pass)) {
	// Define SV event
	MinimalSVRecord sv;
	sv.tid = rec->rid;
	sv.svStart = rec->pos;
	sv.svEnd = sv.svStart + 1;
	sv.svLen = 1;
	sv.svt = _decodeOrientation(ctVal, svtVal);
	sv.qual = rec->qual;
	if (svtVal != "INS") {
	  sv.svLen = svEndVal - rec->pos;
	  sv.svEnd = svEndVal;
	} else {
	  sv.svLen = inslenVal;
	}

	// Check genotypes
	bcf_unpack(rec, BCF_UN_ALL);
	bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
	sv.gt.resize(bcf_hdr_nsamples(hdr), 0);
	for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
	  if ((bcf_gt_allele(gt[i*2]) != -1) && (bcf_gt_allele(gt[i*2 + 1]) != -1)) {
	    sv.gt[i] = bcf_gt_allele(gt[i*2]) + bcf_gt_allele(gt[i*2 + 1]);
	  }
	}
	sv.id = std::string(rec->d.id);
	if (allIds.find(sv.id) != allIds.end()) {
	  success=false;
	  std::cerr << "Error: Duplicate IDs " << sv.id << std::endl;
	} else {
	  std::cerr << sv.tid << ',' << sv.svStart << ',' << sv.svEnd << ',' << sv.id << ',' << sv.svLen << ',' << sv.svt << std::endl;
	  allIds.insert(sv.id);
	  allsv.push_back(sv);
	}
      }
    }
    bcf_destroy(rec);

    // Clean-up
    if (svend != NULL) free(svend);
    if (svt != NULL) free(svt);
    if (gt != NULL) free(gt);
    if (ct != NULL) free(ct);

    // Close VCF
    bcf_hdr_destroy(hdr);
    bcf_close(ifile);

    return success;    
  }

  inline int
  compvcfRun(CompvcfConfig const& c) {

    // Load SVs
    std::vector<MinimalSVRecord> basesv;
    if (!_loadMinimalSVs(c, c.base.string(), basesv)) return -1;

    std::vector<MinimalSVRecord> compsv;
    if (!_loadMinimalSVs(c, c.vcffile.string(), compsv)) return -1;

    std::cerr << basesv.size() << ',' << compsv.size() << std::endl;

    // Sort SVs
    //sort(allsv.begin(), allsv.end(), SortSVEvents<SVEvent>());

    // Mark duplicates
    //_markDuplicates(c, allsv);

    // Write non-duplicate SV sites
    //if (!_writeUniqueSVs(c, allsv)) return -1;

    return 0;
  }


  int compvcf(int argc, char **argv) {
    CompvcfConfig c;
    
    // Define generic options
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("base,b", boost::program_options::value<boost::filesystem::path>(&c.base), "base VCF/BCF file")
      ("quality,y", boost::program_options::value<int32_t>(&c.qualthres)->default_value(300), "min. SV site quality")
      ("bpdiff,b", boost::program_options::value<int32_t>(&c.bpdiff)->default_value(50), "max. SV breakpoint offset")
      ("sizeratio,s", boost::program_options::value<float>(&c.sizeratio)->default_value(0.8), "min. SV size ratio")
      ("pass,p", "Filter sites for PASS")
      ;
    
    // Define hidden options
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "comparison VCF/BCF file")
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
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cerr << std::endl;
      std::cerr << "Usage: delly " << argv[0] << " [OPTIONS] <input.bcf>" << std::endl;
      std::cerr << visible_options << "\n";
      return 0;
    }
    
    // Filter for PASS
    if (vm.count("pass")) c.filterForPass = true;
    else c.filterForPass = false;
        
    // Check input VCF file
    if (vm.count("input-file")) {
      if (!(boost::filesystem::exists(c.vcffile) && boost::filesystem::is_regular_file(c.vcffile) && boost::filesystem::file_size(c.vcffile))) {
	std::cerr << "Input VCF/BCF file is missing: " << c.vcffile.string() << std::endl;
	return 1;
      }
      htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
      if (ifile == NULL) {
	std::cerr << "Fail to open file " << c.vcffile.string() << std::endl;
      return 1;
      }
      hts_idx_t* bcfidx = NULL;
      tbx_t* tbx = NULL;
      if (hts_get_format(ifile)->format==vcf) tbx = tbx_index_load(c.vcffile.string().c_str());
      else bcfidx = bcf_index_load(c.vcffile.string().c_str());
      if ((bcfidx == NULL) && (tbx == NULL)) {
	std::cerr << "Fail to open index file for " << c.vcffile.string() << std::endl;
	return 1;
      }
      bcf_hdr_t* hdr = bcf_hdr_read(ifile);
      if (hdr == NULL) {
	std::cerr << "Fail to header for " << c.vcffile.string() << std::endl;
	return 1;
      }
      if (!(bcf_hdr_nsamples(hdr)>0)) {
	std::cerr << "BCF/VCF file has no sample genotypes!" << std::endl;
	return 1;
      }
      bcf_hdr_destroy(hdr);
      if (bcfidx) hts_idx_destroy(bcfidx);
      if (tbx) tbx_destroy(tbx);
      bcf_close(ifile);
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cerr << "delly ";
    for(int i=0; i<argc; ++i) { std::cerr << argv[i] << ' '; }
    std::cerr << std::endl;
    
    return compvcfRun(c);
  }
  
}

#endif
