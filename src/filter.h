#ifndef FILTER_H
#define FILTER_H

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

#include "tags.h"
#include "version.h"
#include "util.h"
#include "modvcf.h"

namespace torali
{


struct FilterConfig {
  bool filterForPass;
  bool hasSampleFile;
  bool softFilter;
  int32_t minsize;
  int32_t maxsize;
  int32_t coverage;
  int32_t qualthres;
  float ratiogeno;
  float altaf;
  float controlcont;
  float gq;
  float rddel;
  float rddup;
  std::string filter;
  std::set<std::string> tumorSet;
  std::set<std::string> controlSet;
  boost::filesystem::path outfile;
  boost::filesystem::path samplefile;
  boost::filesystem::path vcffile;
};


template<typename TFilterConfig>
inline int
filterRun(TFilterConfig const& c) {

  // Load bcf file
  htsFile* ifile = hts_open(c.vcffile.string().c_str(), "r");
  bcf_hdr_t* hdr = bcf_hdr_read(ifile);

  // Open output VCF file
  std::string fmtout = "wb";
  if (c.outfile.string() == "-") fmtout = "w";
  htsFile *ofile = hts_open(c.outfile.string().c_str(), fmtout.c_str());
  bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);
  if (c.filter == "somatic") {
    bcf_hdr_remove(hdr_out, BCF_HL_INFO, "RDRATIO");
    bcf_hdr_append(hdr_out, "##INFO=<ID=RDRATIO,Number=1,Type=Float,Description=\"Read-depth ratio of tumor vs. normal.\">");
    bcf_hdr_remove(hdr_out, BCF_HL_INFO, "SOMATIC");
    bcf_hdr_append(hdr_out, "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic structural variant.\">");
    if (c.softFilter) {
      bcf_hdr_append(hdr_out, "##FILTER=<ID=FailDellyFilter,Description=\"Failed delly filter.\">");
      bcf_hdr_append(hdr_out, "##FILTER=<ID=FailSomatic,Description=\"Failed somatic filter (likely germline variant).\">");
    }
  } else if (c.filter == "germline") {
    bcf_hdr_remove(hdr_out, BCF_HL_INFO, "RDRATIO");
    bcf_hdr_append(hdr_out, "##INFO=<ID=RDRATIO,Number=1,Type=Float,Description=\"Read-depth ratio of SV carrier vs. non-carrier.\">");
    if (c.softFilter) {
      bcf_hdr_append(hdr_out, "##FILTER=<ID=FailDellyFilter,Description=\"Failed delly filter.\">");
      bcf_hdr_append(hdr_out, "##FILTER=<ID=FailGermline,Description=\"Failed germline filter.\">");
    }
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
  int nrc = 0;
  int32_t* rc = NULL;
  int nrcl = 0;
  int32_t* rcl = NULL;
  int nrcr = 0;
  int32_t* rcr = NULL;
  int ndv = 0;
  int32_t* dv = NULL;
  int ndr = 0;
  int32_t* dr = NULL;
  int nrv = 0;
  int32_t* rv = NULL;
  int nrr = 0;
  int32_t* rr = NULL;
  bool germline = false;
  if (c.filter == "germline") germline = true;

  // Parse BCF
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Filtering VCF/BCF file" << std::endl;
  bcf1_t* rec = bcf_init1();
  while (bcf_read(ifile, hdr, rec) == 0) {
    bcf_unpack(rec, BCF_UN_INFO);

    // Check SV type
    bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt);

    // Check size and PASS
    bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend);
    bool pass = true;
    if (c.filterForPass) pass = (bcf_has_filter(hdr, rec, const_cast<char*>("PASS"))==1);
    int32_t svlen = 1;
    if (svend != NULL) svlen = *svend - rec->pos;
    int32_t inslenVal = 0;
    if (bcf_get_info_int32(hdr, rec, "INSLEN", &inslen, &ninslen) > 0) inslenVal = *inslen;
    if ((rec->qual >= c.qualthres) && (pass) && ((std::string(svt) == "BND") || ((std::string(svt) == "INS") && (inslenVal >= c.minsize) && (inslenVal <= c.maxsize)) || ((std::string(svt) != "BND") && (std::string(svt) != "INS") && (svlen >= c.minsize) && (svlen <= c.maxsize)))) {
      // Check genotypes
      bcf_unpack(rec, BCF_UN_ALL);
      bool precise = false;
      if (bcf_get_info_flag(hdr, rec, "PRECISE", 0, 0) > 0) precise = true;
      bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
      if (_getFormatType(hdr, "GQ") == BCF_HT_INT) bcf_get_format_int32(hdr, rec, "GQ", &gq, &ngq);
      else if (_getFormatType(hdr, "GQ") == BCF_HT_REAL) bcf_get_format_float(hdr, rec, "GQ", &gqf, &ngq);
      bcf_get_format_int32(hdr, rec, "RC", &rc, &nrc);
      if (_isKeyPresent(hdr, "RCL")) bcf_get_format_int32(hdr, rec, "RCL", &rcl, &nrcl);
      if (_isKeyPresent(hdr, "RCR")) bcf_get_format_int32(hdr, rec, "RCR", &rcr, &nrcr);
      bcf_get_format_int32(hdr, rec, "DV", &dv, &ndv);
      bcf_get_format_int32(hdr, rec, "DR", &dr, &ndr);
      bcf_get_format_int32(hdr, rec, "RV", &rv, &nrv);
      bcf_get_format_int32(hdr, rec, "RR", &rr, &nrr);
      std::vector<float> rcraw;
      std::vector<float> rcControl;
      std::vector<float> rcTumor;
      std::vector<float> rcAlt;
      std::vector<float> rRefVar;
      std::vector<float> rAltVar;
      std::vector<float> gqRef;
      std::vector<float> gqAlt;
      uint32_t nCount = 0;
      uint32_t tCount = 0;
      uint32_t controlpass = 0;
      uint32_t tumorpass = 0;
      int32_t ac[2];
      ac[0] = 0;
      ac[1] = 0;
      for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
	if ((bcf_gt_allele(gt[i*2]) != -1) && (bcf_gt_allele(gt[i*2 + 1]) != -1)) {
	  int gt_type = bcf_gt_allele(gt[i*2]) + bcf_gt_allele(gt[i*2 + 1]);
	  ++ac[bcf_gt_allele(gt[i*2])];
	  ++ac[bcf_gt_allele(gt[i*2 + 1])];
	  if ((germline) || (c.controlSet.find(hdr->samples[i]) != c.controlSet.end())) {
	    // Control or population genomics
	    ++nCount;
	    if (gt_type == 0) {
	      rcraw.push_back(rc[i]);
	      if (_getFormatType(hdr, "GQ") == BCF_HT_INT) gqRef.push_back(gq[i]);
	      else if (_getFormatType(hdr, "GQ") == BCF_HT_REAL) gqRef.push_back(gqf[i]);
	      if ((rcl != NULL) && (rcr != NULL) && (rcl[i] + rcr[i] != 0)) rcControl.push_back((float) rc[i] / ((float) (rcl[i] + rcr[i])));
	      else rcControl.push_back(rc[i]);
	      float rVar = 0;
	      if (!precise) rVar = (float) dv[i] / (float) (dr[i] + dv[i]);
	      else rVar = (float) rv[i] / (float) (rr[i] + rv[i]);
	      rRefVar.push_back(rVar);
	      if (rVar <= c.controlcont) ++controlpass;
	    } else if ((germline) && (gt_type >= 1)) {
	      if (_getFormatType(hdr, "GQ") == BCF_HT_INT) gqAlt.push_back(gq[i]);
	      else if (_getFormatType(hdr, "GQ") == BCF_HT_REAL) gqAlt.push_back(gqf[i]);
	      if ((rcl != NULL) && (rcr != NULL) && (rcl[i] + rcr[i] != 0)) rcAlt.push_back((float) rc[i] / ((float) (rcl[i] + rcr[i])));
	      else rcAlt.push_back(rc[i]);
	      float rVar = 0;
	      if (!precise) rVar = (float) dv[i] / (float) (dr[i] + dv[i]);
	      else rVar = (float) rv[i] / (float) (rr[i] + rv[i]);
	      rAltVar.push_back(rVar);
	    }
	  } else if ((!germline) && (c.tumorSet.find(hdr->samples[i]) != c.tumorSet.end())) {
	    // Tumor
	    ++tCount;
	    if ((rcl != NULL) && (rcr != NULL) && (rcl[i] + rcr[i] != 0)) rcTumor.push_back((float) rc[i] / ((float) (rcl[i] + rcr[i])));
	    else rcTumor.push_back(rc[i]);
	    if (!precise) {
	      if ((((float) dv[i] / (float) (dr[i] + dv[i])) >= c.altaf) && (dr[i] + dv[i] >= c.coverage)) ++tumorpass;
	    } else {
	      if ((((float) rv[i] / (float) (rr[i] + rv[i])) >= c.altaf) && (rr[i] + rv[i] >= c.coverage)) ++tumorpass;
	    }
	  }
	}
      }
      if (c.filter == "somatic") {
	float genotypeRatio = (float) (nCount + tCount) / (float) (c.controlSet.size() + c.tumorSet.size());
	if ((controlpass) && (tumorpass) && (controlpass == nCount) && (genotypeRatio >= c.ratiogeno)) {
	  float rccontrolmed = 0;
	  getMedian(rcControl.begin(), rcControl.end(), rccontrolmed);
	  float rctumormed = 0;
	  getMedian(rcTumor.begin(), rcTumor.end(), rctumormed);
	  float rdRatio = 1;
	  if (rccontrolmed != 0) rdRatio = rctumormed/rccontrolmed;
	  _remove_info_tag(hdr_out, rec, "RDRATIO");
	  bcf_update_info_float(hdr_out, rec, "RDRATIO", &rdRatio, 1);
	  _remove_info_tag(hdr_out, rec, "SOMATIC");
	  bcf_update_info_flag(hdr_out, rec, "SOMATIC", NULL, 1);
	  bcf_write1(ofile, hdr_out, rec);
	} else if (c.softFilter) {
	  int32_t tmpi = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "FailSomatic");
	  bcf_update_filter(hdr_out, rec, &tmpi, 1);
	  bcf_write1(ofile, hdr_out, rec);
	} 
      } else if (c.filter == "germline") {
	float genotypeRatio = (float) (nCount + tCount) / (float) (bcf_hdr_nsamples(hdr));
	float rrefvarpercentile = 0;
	if (!rRefVar.empty()) getPercentile(rRefVar, 0.9, rrefvarpercentile);
	float raltvarmed = 0;
	if (!rAltVar.empty()) getMedian(rAltVar.begin(), rAltVar.end(), raltvarmed);
	float rccontrolmed = 0;
	if (!rcControl.empty()) getMedian(rcControl.begin(), rcControl.end(), rccontrolmed);
	float rcaltmed = 0;
	if (!rcAlt.empty()) getMedian(rcAlt.begin(), rcAlt.end(), rcaltmed);
	float rdRatio = 1;
	if (rccontrolmed != 0) rdRatio = rcaltmed/rccontrolmed;
	float gqaltmed = 0;
	if (!gqAlt.empty()) getMedian(gqAlt.begin(), gqAlt.end(), gqaltmed);
	float gqrefmed = 0;
	if (!gqRef.empty()) getMedian(gqRef.begin(), gqRef.end(), gqrefmed);
	float af = (float) ac[1] / (float) (ac[0] + ac[1]);
	
	//std::cerr << bcf_hdr_id2name(hdr, rec->rid) << '\t' << (rec->pos + 1) << '\t' << *svend << '\t' << rec->d.id << '\t' << svlen << '\t' << ac[1] << '\t' << af << '\t' << genotypeRatio << '\t' << std::string(svt) << '\t' << precise << '\t' << rrefvarpercentile << '\t' << raltvarmed << '\t' << gqrefmed << '\t' << gqaltmed << '\t' << rdRatio << std::endl;

	bool failgerm = false;
	if (!((af>0) && (gqaltmed >= c.gq) && (gqrefmed >= c.gq) && (raltvarmed >= c.altaf) && (genotypeRatio >= c.ratiogeno))) failgerm = true;
	if ((std::string(svt)=="DEL") && (rdRatio > c.rddel)) failgerm = true;
	if ((std::string(svt)=="DUP") && (rdRatio < c.rddup)) failgerm = true;
	if ((std::string(svt)!="DEL") && (std::string(svt)!="DUP") && (rrefvarpercentile > 0)) failgerm = true;
	if (!failgerm) {
	  _remove_info_tag(hdr_out, rec, "RDRATIO");
	  bcf_update_info_float(hdr_out, rec, "RDRATIO", &rdRatio, 1);
	  bcf_write1(ofile, hdr_out, rec);
	} else if (c.softFilter) {
	  int32_t tmpi = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "FailGermline");
	  bcf_update_filter(hdr_out, rec, &tmpi, 1);
	  bcf_write1(ofile, hdr_out, rec);
	}
      }
    } else if (c.softFilter) {
      int32_t tmpi = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "FailDellyFilter");
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
  if (rc != NULL) free(rc);
  if (rcl != NULL) free(rcl);
  if (rcr != NULL) free(rcr);
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
  now = boost::posix_time::second_clock::local_time();
  std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;

  return 0;
}


int filter(int argc, char **argv) {
  FilterConfig c;

  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("filter,f", boost::program_options::value<std::string>(&c.filter)->default_value("somatic"), "Filter mode (somatic, germline)")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "Filtered SV BCF output file")
    ("quality,y", boost::program_options::value<int32_t>(&c.qualthres)->default_value(300), "min. SV site quality")
    ("altaf,a", boost::program_options::value<float>(&c.altaf)->default_value(0.03), "min. fractional ALT support")
    ("minsize,m", boost::program_options::value<int32_t>(&c.minsize)->default_value(0), "min. SV size")
    ("maxsize,n", boost::program_options::value<int32_t>(&c.maxsize)->default_value(500000000), "max. SV size")
    ("ratiogeno,r", boost::program_options::value<float>(&c.ratiogeno)->default_value(0.75), "min. fraction of genotyped samples")
    ("pass,p", "Filter sites for PASS")
    ("tag,t", "Tag filtered sites in the FILTER column instead of removing them")
    ;

  // Define somatic options
  boost::program_options::options_description somatic("Somatic options");
  somatic.add_options()
    ("samples,s", boost::program_options::value<boost::filesystem::path>(&c.samplefile), "Two-column sample file listing sample name and tumor or control")
    ("coverage,v", boost::program_options::value<int32_t>(&c.coverage)->default_value(10), "min. coverage in tumor")
    ("controlcontamination,c", boost::program_options::value<float>(&c.controlcont)->default_value(0.0), "max. fractional ALT support in control")
    ;

  // Define germline options
  boost::program_options::options_description germline("Germline options");
  germline.add_options()
    ("gq,q", boost::program_options::value<float>(&c.gq)->default_value(15), "min. median GQ for carriers and non-carriers")
    ("rddel,e", boost::program_options::value<float>(&c.rddel)->default_value(0.8), "max. read-depth ratio of carrier vs. non-carrier for a deletion")
    ("rddup,u", boost::program_options::value<float>(&c.rddup)->default_value(1.2), "min. read-depth ratio of carrier vs. non-carrier for a duplication")
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
  cmdline_options.add(generic).add(somatic).add(germline).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic).add(somatic).add(germline);
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

  // Population Genomics
  if (c.filter == "germline") c.controlcont = 1.0;

  // Check sample file
  std::set<std::string> tSet;
  std::set<std::string> cSet;
  if (c.filter == "somatic") {
    c.hasSampleFile = true;
    if (!(boost::filesystem::exists(c.samplefile) && boost::filesystem::is_regular_file(c.samplefile) && boost::filesystem::file_size(c.samplefile))) {
      std::cerr << "Sample file is missing " << c.samplefile.string() << std::endl;
      return 1;
    } else {
      // Get samples
      std::ifstream sampleFile(c.samplefile.string().c_str(), std::ifstream::in);
      if (sampleFile.is_open()) {
	while (sampleFile.good()) {
	  std::string sampleFromFile;
	  getline(sampleFile, sampleFromFile);
	  typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	  boost::char_separator<char> sep(",\t ");
	  Tokenizer tokens(sampleFromFile, sep);
	  Tokenizer::iterator tokIter = tokens.begin();
	  if (tokIter != tokens.end()) {
	    std::string sample = *tokIter++;
	    if (tokIter != tokens.end()) {
	      std::string type = *tokIter;
	      if (type == "control") cSet.insert(sample);
	      else if (type == "tumor") tSet.insert(sample);
	      else {
		std::cerr << "Sample type for " << sample << " is neither tumor nor control" << std::endl;
		return 1;
	      }
	    }
	  }
	}
	sampleFile.close();
      }
      if (tSet.empty()) {
	std::cerr << "No tumor samples specified." << std::endl;
	return 1;
      }
      if (cSet.empty()) {
	std::cerr << "No control samples specified." << std::endl;
	return 1;
      }
      std::vector<std::string> intersection;
      std::set_intersection(cSet.begin(), cSet.end(), tSet.begin(), tSet.end(), std::back_inserter(intersection));
      if (!intersection.empty()) {
	std::cerr << "Sample " << intersection[0] << " is both a tumor and control sample." << std::endl;
	return 1;
      }
    }
  } else c.hasSampleFile = false;

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
    // Check sample names
    if (c.filter == "somatic") {
      for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
	if (tSet.find(hdr->samples[i]) != tSet.end()) c.tumorSet.insert(hdr->samples[i]);
	else if (cSet.find(hdr->samples[i]) != cSet.end()) c.controlSet.insert(hdr->samples[i]);
	else std::cerr << "Warning: Sample " << hdr->samples[i] << " is missing in sample file." << std::endl;
      }
      if (c.tumorSet.empty()) {
	std::cerr << "No tumor samples specified." << std::endl;
	return 1;
      }
      if (c.controlSet.empty()) {
	std::cerr << "No control samples specified." << std::endl;
	return 1;
      }
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

  return filterRun(c);
}

}

#endif
