#ifndef DPE_H
#define DPE_H

#include <iostream>
#include <fstream>
#include <boost/dynamic_bitset.hpp>
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

#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/sam.h>

#include "tags.h"
#include "coverage.h"
#include "version.h"
#include "util.h"
#include "modvcf.h"

namespace torali
{

  struct DoublePEConfig {
    int32_t wiggle;
    int32_t svsize;
    float carconc;
    boost::filesystem::path outfile;
    boost::filesystem::path infile;
  };

  struct SVCarrier {
    typedef boost::dynamic_bitset<> TBitSet;
    
    int32_t start;
    int32_t end;
    std::string id;
    TBitSet carrier;
    
    SVCarrier(int32_t s, int32_t e, std::string i, TBitSet c) : start(s), end(e), id(i), carrier(c) {}
  };

  struct DPERecord {
    int32_t start1;
    int32_t end1;
    int32_t start2;
    int32_t end2;
    float carconc;
    std::string id1;
    std::string id2;
    
    DPERecord(int32_t s1, int32_t e1, int32_t s2, int32_t e2, float cc, std::string i1, std::string i2) : start1(s1), end1(e1), start2(s2), end2(e2), carconc(cc), id1(i1), id2(i2) {}
  };
  
  
  inline int
  dpeRun(DoublePEConfig const& c)
  {
    
    // Open BCF file
    htsFile* ifile = bcf_open(c.infile.string().c_str(), "r");
    hts_idx_t* bcfidx = bcf_index_load(c.infile.string().c_str());
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);
    
    // Read BCF file
    int32_t nsvend = 0;
    int32_t* svend = NULL;
    int32_t nsvt = 0;
    char* svt = NULL;
    int32_t nchr2 = 0;
    char* chr2 = NULL;
    int32_t nct = 0;
    char* ct = NULL;
    int ngt = 0;
    int32_t* gt = NULL;
    
    // Get sequences
    int32_t nseq = 0;
    const char** seqnames = bcf_hdr_seqnames(hdr, &nseq);

    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Searching complex SVs" << std::endl;
    
    // Open output file
    std::string fmtout = "wb";
    if (c.outfile.string() == "-") fmtout = "w";
    htsFile *ofile = hts_open(c.outfile.string().c_str(), fmtout.c_str());
    bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);
    bcf_hdr_remove(hdr_out, BCF_HL_INFO, "LINKID");
    bcf_hdr_append(hdr_out, "##INFO=<ID=LINKID,Number=2,Type=String,Description=\"Linked paired-end IDs.\">");
    bcf_hdr_remove(hdr_out, BCF_HL_INFO, "REGION");
    bcf_hdr_append(hdr_out, "##INFO=<ID=REGION,Number=3,Type=String,Description=\"Entire spanning region of the complex SV as chr, start, end.\">");
    bcf_hdr_remove(hdr_out, BCF_HL_INFO, "REGION1");
    bcf_hdr_append(hdr_out, "##INFO=<ID=REGION1,Number=3,Type=String,Description=\"Sub-Region1 of the complex SV as chr, start, end.\">");
    bcf_hdr_remove(hdr_out, BCF_HL_INFO, "REGION2");
    bcf_hdr_append(hdr_out, "##INFO=<ID=REGION2,Number=3,Type=String,Description=\"Sub-Region2 of the complex SV as chr, start, end.\">");
    bcf_hdr_remove(hdr_out, BCF_HL_INFO, "REGION3");
    bcf_hdr_append(hdr_out, "##INFO=<ID=REGION3,Number=3,Type=String,Description=\"Sub-Region3 of the complex SV as chr, start, end.\">");
    bcf_hdr_remove(hdr_out, BCF_HL_INFO, "CARCONC");
    bcf_hdr_append(hdr_out, "##INFO=<ID=CARCONC,Number=1,Type=Float,Description=\"Carrier concordance of the linked paired-end calls.\">");
    if (bcf_hdr_write(ofile, hdr_out) != 0) std::cerr << "Error: Failed to write BCF header!" << std::endl;
    
    // Parse BCF
    for(int32_t refIndex = 0; refIndex < nseq; ++refIndex) {
      // Fetch SVs on this chromosome
      int32_t maxCTs = 5;
      typedef std::vector<SVCarrier> TSVCarrier;
      typedef std::vector<TSVCarrier> TCTs;
      TCTs cts(maxCTs);
      hts_itr_t* itervcf = bcf_itr_querys(bcfidx, hdr, bcf_hdr_id2name(hdr, refIndex));
      bcf1_t* rec = bcf_init();
      while (bcf_itr_next(ifile, itervcf, rec) >= 0) {
	// Fetch info
	bcf_unpack(rec, BCF_UN_ALL);
	bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
	bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend);
	bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt);
	std::string chr2Name("NA");
	if (bcf_get_info_string(hdr, rec, "CHR2", &chr2, &nchr2) > 0) chr2Name = std::string(chr2);
	uint8_t ict = 0;
	if (bcf_get_info_string(hdr, rec, "CT", &ct, &nct) > 0) ict = _decodeOrientation(std::string(ct));
	
	// Fetch carriers
	if ((*svend - rec->pos) < c.svsize) {
	  SVCarrier::TBitSet car(bcf_hdr_nsamples(hdr));
	  for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
	    if ((bcf_gt_allele(gt[i*2]) != -1) && (bcf_gt_allele(gt[i*2 + 1]) != -1)) {
	      int gt_type = bcf_gt_allele(gt[i*2]) + bcf_gt_allele(gt[i*2 + 1]);
	      if (gt_type > 0) car[i] = true;
	    }
	  }
	  cts[(int32_t) ict].push_back(SVCarrier(rec->pos, *svend, rec->d.id, car));
	}
      }
      bcf_destroy(rec);
      hts_itr_destroy(itervcf);
      
      // Process SVs
      typedef std::vector<DPERecord> Tdper;
      Tdper dper;
      typedef std::set<std::string> TSvIds;
      TSvIds svIds;
      for(int32_t i = 0; i<maxCTs; ++i) {
	if (!cts[i].empty()) {
	  for(int32_t j = i+1; j<maxCTs; ++j) {
	    if (!cts[j].empty()) {
	      // Compare these 2 CTs
	      typedef std::vector<int32_t> TJPCount;
	      TJPCount jpc(cts[j].size(), 0);
	      typedef std::map<std::pair<int32_t, int32_t>, float> TScoreMap;
	      TScoreMap scm;
	      for(int32_t ip = 0; ip < (int32_t) cts[i].size(); ++ip) {
		int32_t bestJP = -1;
		float bestCC = -1;
		for(int32_t jp = 0; jp < (int32_t) cts[j].size(); ++jp) {
		  if ((cts[i][ip].end < cts[j][jp].start) || (cts[j][jp].end < cts[i][ip].start)) continue;
		  if (((cts[i][ip].start - c.wiggle < cts[j][jp].start) && (cts[j][jp].start < cts[i][ip].end) && (cts[i][ip].end - c.wiggle < cts[j][jp].end)) || ((cts[j][jp].start - c.wiggle < cts[i][ip].start) && (cts[i][ip].start < cts[j][jp].end) && (cts[j][jp].end - c.wiggle < cts[i][ip].end))) {
		    int32_t common = (cts[i][ip].carrier & cts[j][jp].carrier).count();
		    int32_t all = (cts[i][ip].carrier | cts[j][jp].carrier).count();
		    float cc = 0;
		    if (all > 0) cc = (float) common / (float) all;
		    if ((cc >= c.carconc) && (cc > bestCC)) {
		      bestJP = jp;
		      bestCC = cc;
		    }
		  }
		}
		if (bestJP >= 0) {
		  scm.insert(std::make_pair(std::make_pair(ip, bestJP), bestCC));
		  ++jpc[bestJP];
		}
	      }
	      for(TScoreMap::iterator scmIt = scm.begin(); scmIt != scm.end(); ++scmIt) {
		int32_t ip = scmIt->first.first;
		int32_t jp = scmIt->first.second;
		float cc = scmIt->second;
		bool savePair = true;
		if (jpc[jp] > 1) {
		  for(TScoreMap::iterator scmSec = scm.begin(); scmSec != scm.end(); ++scmSec) {
		    if ((ip != scmSec->first.first) && (jp == scmSec->first.second)) {
		      if ((cc < scmSec->second) || ((cc == scmSec->second) && (ip > scmSec->first.first))) {
			savePair = false;
			break;
		      }
		    }
		  }
		}
		if (savePair) {
		  dper.push_back(DPERecord(cts[i][ip].start, cts[i][ip].end, cts[j][jp].start, cts[j][jp].end, cc, cts[i][ip].id, cts[j][jp].id));
		  if (!svIds.insert(cts[i][ip].id).second) std::cerr << "SV already exists!" << std::endl;
		  if (!svIds.insert(cts[j][jp].id).second) std::cerr << "SV already exists!" << std::endl;
		}
	      }
	    }
	  }
	}
      }
      
      
      hts_itr_t* ivcf = bcf_itr_querys(bcfidx, hdr, bcf_hdr_id2name(hdr, refIndex));
      bcf1_t* r = bcf_init();
      while (bcf_itr_next(ifile, ivcf, r) >= 0) {
	bcf_unpack(r, BCF_UN_ALL);
	bcf_get_info_int32(hdr, r, "END", &svend, &nsvend);
	std::string id = std::string(r->d.id);
	if (svIds.find(id) != svIds.end()) {
	  // Find matching DPERecord
	  for(int32_t i = 0; i < (int32_t) dper.size(); ++i) {
	    if (((dper[i].id1 == id) && (dper[i].start1 == r->pos) && (dper[i].end1 == *svend)) || ((dper[i].id2 == id) && (dper[i].start2 == r->pos) && (dper[i].end2 == *svend))) {
	      
	      std::string linkid = dper[i].id1 + "," + dper[i].id2;
	      _remove_info_tag(hdr_out, r, "LINKID");
	      bcf_update_info_string(hdr_out, r, "LINKID", linkid.c_str());
	      _remove_info_tag(hdr_out, r, "CARCONC");
	      bcf_update_info_float(hdr_out, r, "CARCONC", &dper[i].carconc, 1);
	      std::string reg = bcf_hdr_id2name(hdr, refIndex);
	      reg += "," + boost::lexical_cast<std::string>(std::min(dper[i].start1 + 1, dper[i].start2 + 1));
	      reg += "," + boost::lexical_cast<std::string>(std::max(dper[i].end1, dper[i].end2));
	      _remove_info_tag(hdr_out, r, "REGION");
	      bcf_update_info_string(hdr_out, r, "REGION", reg.c_str());
	      std::string reg1 = bcf_hdr_id2name(hdr, refIndex);
	      reg1 += "," + boost::lexical_cast<std::string>(std::min(dper[i].start1 + 1, dper[i].start2 + 1));
	      reg1 += "," + boost::lexical_cast<std::string>(std::max(dper[i].start1 + 1, dper[i].start2 + 1));
	      _remove_info_tag(hdr_out, r, "REGION1");
	      bcf_update_info_string(hdr_out, r, "REGION1", reg1.c_str());
	      std::string reg2 = bcf_hdr_id2name(hdr, refIndex);
	      reg2 += "," + boost::lexical_cast<std::string>(std::max(dper[i].start1 + 1, dper[i].start2 + 1));
	      reg2 += "," + boost::lexical_cast<std::string>(std::min(dper[i].end1, dper[i].end2));
	      _remove_info_tag(hdr_out, r, "REGION2");
	      bcf_update_info_string(hdr_out, r, "REGION2", reg2.c_str());
	      std::string reg3 = bcf_hdr_id2name(hdr, refIndex);
	      reg3 += "," + boost::lexical_cast<std::string>(std::min(dper[i].end1, dper[i].end2));
	      reg3 += "," + boost::lexical_cast<std::string>(std::max(dper[i].end1, dper[i].end2));
	      _remove_info_tag(hdr_out, r, "REGION3");
	      bcf_update_info_string(hdr_out, r, "REGION3", reg3.c_str());
	      bcf_write1(ofile, hdr_out, r);
	    }
	}
	}
      }
      bcf_destroy(r);
      hts_itr_destroy(ivcf);    
    }
    if (nseq) free(seqnames);
    
    // Close output BCF
    bcf_hdr_destroy(hdr_out);
    hts_close(ofile);
    
    // Build index
    if (c.outfile.string() != "-") bcf_index_build(c.outfile.string().c_str(), 14);
    
    // Clean-up
    if (svend != NULL) free(svend);
    if (svt != NULL) free(svt);
    if (chr2 != NULL) free(chr2);
    if (ct != NULL) free(ct);
    if (gt != NULL) free(gt);
    
    // BCF clean-up
    bcf_hdr_destroy(hdr);
    hts_idx_destroy(bcfidx);
    bcf_close(ifile);
    
    
    now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
    return 0;
  }
  
  
  int dpe(int argc, char **argv) {
    DoublePEConfig c;
    c.wiggle = 150;
    
    // Define generic options
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("svsize,s", boost::program_options::value<int32_t>(&c.svsize)->default_value(50000), "max. SV size")
      ("carconc,c", boost::program_options::value<float>(&c.carconc)->default_value(0.75), "min. carrier concordance")
      ("outfile,f", boost::program_options::value<boost::filesystem::path>(&c.outfile), "complex SV output file")
      ;
    
    // Define hidden options
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.infile), "input BCF file")
      ("license,l", "show license")
      ("warranty,w", "show warranty")
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
      printTitle("Complex SVs using double paired-end signatures");
      if (vm.count("warranty")) {
	displayWarranty();
      } else if (vm.count("license")) {
	bsd();
      } else {
	std::cerr << "Usage: " << argv[0] << " [OPTIONS] <deldup.bcf>" << std::endl;
	std::cerr << visible_options << "\n"; 
      }
      return 1; 
    }
    
    // Check input VCF file
    if (vm.count("input-file")) {
      if (!(boost::filesystem::exists(c.infile) && boost::filesystem::is_regular_file(c.infile) && boost::filesystem::file_size(c.infile))) {
	std::cerr << "Input BCF file is missing: " << c.infile.string() << std::endl;
	return 1;
      }
      htsFile* ifile = bcf_open(c.infile.string().c_str(), "r");
      if (ifile == NULL) {
	std::cerr << "Fail to open file " << c.infile.string() << std::endl;
      return 1;
      }
      hts_idx_t* bcfidx = bcf_index_load(c.infile.string().c_str());
      if (bcfidx == NULL) {
	std::cerr << "Fail to open index file for " << c.infile.string() << std::endl;
	return 1;
      }
      bcf_hdr_t* hdr = bcf_hdr_read(ifile);
      if (hdr == NULL) {
	std::cerr << "Fail to header for " << c.infile.string() << std::endl;
	return 1;
      }
      bcf_hdr_destroy(hdr);
      hts_idx_destroy(bcfidx);
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
    for(int i=0; i<argc; ++i) { std::cerr << argv[i] << ' '; }
    std::cerr << std::endl;
 
    // Run coverage annotation
    return dpeRun(c);
  }

}

#endif
