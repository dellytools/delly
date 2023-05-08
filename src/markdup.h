#ifndef MARKDUP_H
#define MARKDUP_H

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

  struct SVEvent {
    bool duplicate;
    int32_t tid;
    int32_t svStart;
    int32_t svEnd;
    int32_t svLen;
    int32_t svt;
    std::string id;
    std::string consensus;
    std::vector<float> vaf;
    std::vector<float> gq;
    std::vector<int32_t> gt;
  };

  template<typename TSV>
  struct SortSVEvents : public std::binary_function<TSV, TSV, bool>
  {
    inline bool operator()(TSV const& sv1, TSV const& sv2) {
      return ((sv1.tid<sv2.tid) || ((sv1.tid==sv2.tid) && (sv1.svStart<sv2.svStart)) || ((sv1.tid==sv2.tid) && (sv1.svStart==sv2.svStart) && (sv1.svEnd<sv2.svEnd)));
    }
  };
  
  
  struct MarkdupConfig {
    bool filterForPass;
    bool softFilter;
    int32_t qualthres;
    int32_t bpdiff;
    float sizeratio;
    float divergence;
    boost::filesystem::path outfile;
    boost::filesystem::path vcffile;
  };

  inline bool
  _loadSVEvents(MarkdupConfig const& c, std::vector<SVEvent>& allsv) {
    bool success = true;
    std::set<std::string> allIds;
    
    // Load bcf file
    htsFile* ifile = hts_open(c.vcffile.string().c_str(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);

    // VCF fields
    int32_t nsvend = 0;
    int32_t* svend = NULL;
    int32_t nsvt = 0;
    char* svt = NULL;
    int32_t ninslen = 0;
    int32_t* inslen = NULL;
    int ngt = 0;
    int32_t* gt = NULL;
    int ngq = 0;
    int32_t* gq = NULL;
    float* gqf = NULL;
    int ndv = 0;
    int32_t* dv = NULL;
    int ndr = 0;
    int32_t* dr = NULL;
    int nrv = 0;
    int32_t* rv = NULL;
    int nrr = 0;
    int32_t* rr = NULL;
    int32_t ncons = 0;
    char* cons = NULL;
    
    // Parse BCF
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Parsing VCF/BCF file" << std::endl;
    bcf1_t* rec = bcf_init1();
    while (bcf_read(ifile, hdr, rec) == 0) {
      bcf_unpack(rec, BCF_UN_INFO);

      // Check SV type
      bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt);

      // Check size and PASS
      bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend);
      bool pass = true;
      if (c.filterForPass) pass = (bcf_has_filter(hdr, rec, const_cast<char*>("PASS"))==1);
      int32_t inslenVal = 0;
      if (bcf_get_info_int32(hdr, rec, "INSLEN", &inslen, &ninslen) > 0) inslenVal = *inslen;
      bool precise = false;
      if (bcf_get_info_flag(hdr, rec, "PRECISE", 0, 0) > 0) precise=true;
      if ((rec->qual >= c.qualthres) && (pass)) {
	// Define SV event
	SVEvent sv;
	sv.duplicate = false;
	sv.tid = rec->rid;
	sv.svStart = rec->pos;
	sv.svEnd = sv.svStart + 1;
	sv.svLen = 1;
	if (svend != NULL) {
	  sv.svLen = *svend - rec->pos;
	  sv.svEnd = *svend;
	}
	if (std::string(svt) == "INS") sv.svLen = inslenVal;
	if (precise) {
	  if (bcf_get_info_string(hdr, rec, "CONSENSUS", &cons, &ncons) > 0) {
	    sv.consensus = boost::to_upper_copy(std::string(cons));
	  }
	}
	
	// Check genotypes
	bcf_unpack(rec, BCF_UN_ALL);
	bool precise = false;
	if (bcf_get_info_flag(hdr, rec, "PRECISE", 0, 0) > 0) precise = true;
	bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
	if (_getFormatType(hdr, "GQ") == BCF_HT_INT) bcf_get_format_int32(hdr, rec, "GQ", &gq, &ngq);
	else if (_getFormatType(hdr, "GQ") == BCF_HT_REAL) bcf_get_format_float(hdr, rec, "GQ", &gqf, &ngq);
	bcf_get_format_int32(hdr, rec, "DV", &dv, &ndv);
	bcf_get_format_int32(hdr, rec, "DR", &dr, &ndr);
	bcf_get_format_int32(hdr, rec, "RV", &rv, &nrv);
	bcf_get_format_int32(hdr, rec, "RR", &rr, &nrr);

	sv.vaf.resize(bcf_hdr_nsamples(hdr), 0);
	sv.gq.resize(bcf_hdr_nsamples(hdr), 0);
	sv.gt.resize(bcf_hdr_nsamples(hdr), 0);
	for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
	  if ((bcf_gt_allele(gt[i*2]) != -1) && (bcf_gt_allele(gt[i*2 + 1]) != -1)) {
	    sv.gt[i] = bcf_gt_allele(gt[i*2]) + bcf_gt_allele(gt[i*2 + 1]);
	    if (_getFormatType(hdr, "GQ") == BCF_HT_INT) sv.gq[i] = gq[i];
	    else if (_getFormatType(hdr, "GQ") == BCF_HT_REAL) sv.gq[i] = gqf[i];
	    float rVar = 0;
	    if (!precise) rVar = (float) dv[i] / (float) (dr[i] + dv[i]);
	    else rVar = (float) rv[i] / (float) (rr[i] + rv[i]);
	    sv.vaf[i] = rVar;
	  }
	}
	sv.id = std::string(rec->d.id);
	if (allIds.find(sv.id) != allIds.end()) {
	  success=false;
	  std::cerr << "Error: Duplicate IDs " << sv.id << std::endl;
	} else {
	  allIds.insert(sv.id);
	  allsv.push_back(sv);
	}
      }
    }
    bcf_destroy(rec);

    // Clean-up
    if (svend != NULL) free(svend);
    if (svt != NULL) free(svt);
    if (inslen != NULL) free(inslen);
    if (gt != NULL) free(gt);
    if (gq != NULL) free(gq);
    if (gqf != NULL) free(gqf);
    if (dv != NULL) free(dv);
    if (dr != NULL) free(dr);
    if (rv != NULL) free(rv);
    if (rr != NULL) free(rr);
    if (cons != NULL) free(cons);

    // Close VCF
    bcf_hdr_destroy(hdr);
    bcf_close(ifile);

    return success;    
  }

  
  inline void
  _markDuplicates(MarkdupConfig const& c, std::vector<SVEvent>& allsv) {
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Mark duplicates" << std::endl;
    for(uint32_t i = 0; i < allsv.size(); ++i) {
      for(uint32_t j = i + 1; j < allsv.size(); ++j) {
	if (allsv[i].tid != allsv[j].tid) break;
	if (allsv[j].svStart - allsv[i].svStart > c.bpdiff) break;
	float sizerat = (float) allsv[i].svLen / (float) allsv[j].svLen;
	if (allsv[i].svLen > allsv[j].svLen) sizerat = (float) allsv[j].svLen / (float) allsv[i].svLen;
	if (sizerat < c.sizeratio) continue;
	std::string longc;
	std::string shortc;
	if (allsv[i].consensus.size() > allsv[j].consensus.size()) {
	  longc = allsv[i].consensus;
	  int32_t deslen = 0.8  * allsv[j].consensus.size();
	  int32_t offset = (allsv[j].consensus.size() - deslen)/2;
	  shortc = allsv[j].consensus.substr(offset, deslen);
	} else {
	  longc = allsv[j].consensus;
	  int32_t deslen = 0.8  * allsv[i].consensus.size();
	  int32_t offset = (allsv[i].consensus.size() - deslen)/2;
	  shortc = allsv[i].consensus.substr(offset, deslen);
	}
	EdlibAlignResult cigar = edlibAlign(shortc.c_str(), shortc.size(), longc.c_str(), longc.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
	//printAlignment(shortc, longc, EDLIB_MODE_HW, cigar);
	double score = (double) cigar.editDistance / (double) shortc.size();
	edlibFreeAlignResult(cigar);
	if (score > c.divergence) continue;

	// Duplicate
	std::cerr << allsv[i].tid << ',' << allsv[i].svStart << ',' << allsv[i].svEnd << ',' << allsv[i].id << ',' << allsv[i].svLen << std::endl;
	std::cerr << allsv[j].tid << ',' << allsv[j].svStart << ',' << allsv[j].svEnd << ',' << allsv[i].id << ',' << allsv[j].svLen << std::endl;
	std::cerr << std::endl;
      }
    }
  }
  

  inline int
  markdupRun(MarkdupConfig const& c) {

    // Load SVs
    std::vector<SVEvent> allsv;
    if (!_loadSVEvents(c, allsv)) return -1;

    // Sort SVs
    sort(allsv.begin(), allsv.end(), SortSVEvents<SVEvent>());

    // Mark duplicates
    _markDuplicates(c, allsv);
    
    /*

        // Load bcf file
    htsFile* ifile = hts_open(c.vcffile.string().c_str(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);

    // Open output VCF file
    std::string fmtout = "wb";
    if (c.outfile.string() == "-") fmtout = "w";
    htsFile *ofile = hts_open(c.outfile.string().c_str(), fmtout.c_str());
    bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);
    if (c.softFilter) {
      bcf_hdr_append(hdr_out, "##FILTER=<ID=Duplicate,Description=\"Marked duplicate.\">");
      bcf_hdr_append(hdr_out, "##FILTER=<ID=FAIL,Description=\"SV site fails quality check.\">");
    }
    if (bcf_hdr_write(ofile, hdr_out) != 0) std::cerr << "Error: Failed to write BCF header!" << std::endl;

    // VCF fields
    int32_t nsvend = 0;
    int32_t* svend = NULL;
    int32_t nsvt = 0;
    char* svt = NULL;
    int32_t ninslen = 0;
    int32_t* inslen = NULL;
    int ngt = 0;
    int32_t* gt = NULL;
    int ngq = 0;
    int32_t* gq = NULL;
    float* gqf = NULL;
    int ndv = 0;
    int32_t* dv = NULL;
    int ndr = 0;
    int32_t* dr = NULL;
    int nrv = 0;
    int32_t* rv = NULL;
    int nrr = 0;
    int32_t* rr = NULL;
    bool germline = false;

    // Parse BCF
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Parsing VCF/BCF file" << std::endl;
    bcf1_t* rec = bcf_init1();
    while (bcf_read(ifile, hdr, rec) == 0) {
      bcf_unpack(rec, BCF_UN_INFO);

      // Check SV type
      bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt);

      // Check size and PASS
      bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend);
      bool pass = true;
      if (c.filterForPass) pass = (bcf_has_filter(hdr, rec, const_cast<char*>("PASS"))==1);
      int32_t inslenVal = 0;
      if (bcf_get_info_int32(hdr, rec, "INSLEN", &inslen, &ninslen) > 0) inslenVal = *inslen;
      if ((rec->qual >= c.qualthres) && (pass)) {
	// Define SV event
	SVEvent sv;
	sv.svStart = rec->pos;
	sv.svEnd = sv.svStart + 1;
	sv.svLen = 1;
	if (svend != NULL) {
	  sv.svLen = *svend - rec->pos;
	  sv.svEnd = *svend;
	}
	if (std::string(svt) == "INS") sv.svLen = inslenVal;
	
	// Check genotypes
	bcf_unpack(rec, BCF_UN_ALL);
	bool precise = false;
	if (bcf_get_info_flag(hdr, rec, "PRECISE", 0, 0) > 0) precise = true;
	bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
	if (_getFormatType(hdr, "GQ") == BCF_HT_INT) bcf_get_format_int32(hdr, rec, "GQ", &gq, &ngq);
	else if (_getFormatType(hdr, "GQ") == BCF_HT_REAL) bcf_get_format_float(hdr, rec, "GQ", &gqf, &ngq);
	bcf_get_format_int32(hdr, rec, "DV", &dv, &ndv);
	bcf_get_format_int32(hdr, rec, "DR", &dr, &ndr);
	bcf_get_format_int32(hdr, rec, "RV", &rv, &nrv);
	bcf_get_format_int32(hdr, rec, "RR", &rr, &nrr);

	sv.vaf.resize(bcf_hdr_nsamples(hdr), 0);
	sv.gq.resize(bcf_hdr_nsamples(hdr), 0);
	sv.gt.resize(bcf_hdr_nsamples(hdr), 0);
	for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
	  if ((bcf_gt_allele(gt[i*2]) != -1) && (bcf_gt_allele(gt[i*2 + 1]) != -1)) {
	    sv.gt[i] = bcf_gt_allele(gt[i*2]) + bcf_gt_allele(gt[i*2 + 1]);
	    if (_getFormatType(hdr, "GQ") == BCF_HT_INT) sv.gq[i] = gq[i];
	    else if (_getFormatType(hdr, "GQ") == BCF_HT_REAL) sv.gq[i] = gqf[i];
	    float rVar = 0;
	    if (!precise) rVar = (float) dv[i] / (float) (dr[i] + dv[i]);
	    else rVar = (float) rv[i] / (float) (rr[i] + rv[i]);
	    sv.vaf[i] = rVar;
	  }
	}
	if (c.softFilter) {
	  int32_t tmpi = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "Duplicate");
	  bcf_update_filter(hdr_out, rec, &tmpi, 1);
	  bcf_write1(ofile, hdr_out, rec);
	} else bcf_write1(ofile, hdr_out, rec);
	
	std::cerr << bcf_hdr_id2name(hdr, rec->rid) << '\t' << (rec->pos + 1) << '\t' << *svend << '\t' << rec->d.id << '\t' << svlen << '\t' << ac[0] << '\t' << ac[1] << '\t' << std::string(svt) << '\t' << precise << std::endl;
      } else if (c.softFilter) {
	int32_t tmpi = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "FAIL");
	bcf_update_filter(hdr_out, rec, &tmpi, 1);
	bcf_write1(ofile, hdr_out, rec);
      }
    }
    bcf_destroy(rec);

    // Clean-up
    if (svend != NULL) free(svend);
    if (svt != NULL) free(svt);
    if (inslen != NULL) free(inslen);
    if (gt != NULL) free(gt);
    if (gq != NULL) free(gq);
    if (gqf != NULL) free(gqf);
    if (dv != NULL) free(dv);
    if (dr != NULL) free(dr);
    if (rv != NULL) free(rv);
    if (rr != NULL) free(rr);

    // Close output VCF
    bcf_hdr_destroy(hdr_out);
    hts_close(ofile);

    // Build index
    if (c.outfile.string() != "-") bcf_index_build(c.outfile.string().c_str(), 14);

    // Close VCF
    bcf_hdr_destroy(hdr);
    bcf_close(ifile);
    
    // End
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    */

    return 0;
  }


  int markdup(int argc, char **argv) {
    MarkdupConfig c;
    
    // Define generic options
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "Filtered SV BCF output file")
      ("quality,y", boost::program_options::value<int32_t>(&c.qualthres)->default_value(300), "min. SV site quality")
      ("bpdiff,b", boost::program_options::value<int32_t>(&c.bpdiff)->default_value(50), "max. SV breakpoint offset")
      ("sizeratio,s", boost::program_options::value<float>(&c.sizeratio)->default_value(0.8), "min. SV size ratio")
      ("divergence,d", boost::program_options::value<float>(&c.divergence)->default_value(0.1), "max. SV allele divergence")
      ("pass,p", "Filter sites for PASS")
      ("tag,t", "Tag duplicate marked sites in the FILTER column instead of removing them")
      ;
    
    // Define hidden options
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input file")
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
    
    // Soft filtering
    if (vm.count("tag")) c.softFilter = true;
    else c.softFilter = false;
    
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
    
    // Check outfile
    if (!vm.count("outfile")) c.outfile = "-";
    else {
      if (c.outfile.string() != "-") {
	if (!_outfileValid(c.outfile)) return 1;
      }
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cerr << "delly ";
    for(int i=0; i<argc; ++i) { std::cerr << argv[i] << ' '; }
    std::cerr << std::endl;
    
    return markdupRun(c);
  }
  
}

#endif
