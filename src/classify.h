#ifndef CLASSIFY_H
#define CLASSIFY_H

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
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

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


struct ClassifyConfig {
  bool filterForPass;
  bool hasSampleFile;
  int32_t minsize;
  int32_t maxsize;
  int32_t qual;
  uint16_t ploidy;
  float pgerm;
  float maxsd;
  float cn_offset;
  std::string filter;
  std::set<std::string> tumorSet;
  std::set<std::string> controlSet;
  boost::filesystem::path outfile;
  boost::filesystem::path samplefile;
  boost::filesystem::path vcffile;
};


template<typename TClassifyConfig>
inline int
classifyRun(TClassifyConfig const& c) {

  // Load bcf file
  htsFile* ifile = hts_open(c.vcffile.string().c_str(), "r");
  bcf_hdr_t* hdr = bcf_hdr_read(ifile);

  // Open output VCF file
  std::string fmtout = "wb";
  if (c.outfile.string() == "-") fmtout = "w";
  htsFile *ofile = hts_open(c.outfile.string().c_str(), fmtout.c_str());
  bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);
  if (c.filter == "somatic") {
    bcf_hdr_remove(hdr_out, BCF_HL_INFO, "SOMATIC");
    bcf_hdr_append(hdr_out, "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic copy-number variant.\">");
    bcf_hdr_remove(hdr_out, BCF_HL_INFO, "PGERM");
    bcf_hdr_append(hdr_out, "##INFO=<ID=PGERM,Number=1,Type=Float,Description=\"Probability of being germline.\">");
    bcf_hdr_remove(hdr_out, BCF_HL_INFO, "CNDIFF");
    bcf_hdr_append(hdr_out, "##INFO=<ID=CNDIFF,Number=1,Type=Float,Description=\"Absolute tumor-normal CN difference.\">");
  } else {
    bcf_hdr_remove(hdr_out, BCF_HL_INFO, "CNSHIFT");
    bcf_hdr_append(hdr_out, "##INFO=<ID=CNSHIFT,Number=1,Type=Float,Description=\"Estimated CN shift.\">");
    bcf_hdr_remove(hdr_out, BCF_HL_INFO, "CNSD");
    bcf_hdr_append(hdr_out, "##INFO=<ID=CNSD,Number=1,Type=Float,Description=\"Estimated CN standard deviation.\">");
  }
  if (bcf_hdr_write(ofile, hdr_out) != 0) std::cerr << "Error: Failed to write BCF header!" << std::endl;

  // VCF fields
  int32_t nsvend = 0;
  int32_t* svend = NULL;
  int32_t nsvt = 0;
  char* svt = NULL;
  int ngqval = 0;
  int32_t* gqval = NULL;
  int ncnval = 0;
  int32_t* cnval = NULL;
  int ncnl = 0;
  float* cnl = NULL;
  int nrdcn = 0;
  float* rdcn = NULL;
  int nrdsd = 0;
  float* rdsd = NULL;
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
    if (std::string(svt) != "CNV") continue;

    // Check PASS
    bool pass = true;
    if (c.filterForPass) pass = (bcf_has_filter(hdr, rec, const_cast<char*>("PASS"))==1);
    if (!pass) continue;

    // Check size
    int32_t svStart= rec->pos - 1;
    bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend);
    int32_t svEnd = *svend;
    if (svStart > svEnd) continue;
    int32_t svlen = svEnd - svStart;
    if ((svlen < c.minsize) || (svlen > c.maxsize)) continue;

    // Check copy-number
    bcf_unpack(rec, BCF_UN_ALL);
    bcf_get_format_int32(hdr, rec, "GQ", &gqval, &ngqval);
    bcf_get_format_int32(hdr, rec, "CN", &cnval, &ncnval);
    bcf_get_format_float(hdr, rec, "CNL", &cnl, &ncnl);
    bcf_get_format_float(hdr, rec, "RDCN", &rdcn, &nrdcn);
    bcf_get_format_float(hdr, rec, "RDSD", &rdsd, &nrdsd);

    typedef std::pair<float, float> TCnSd;
    typedef std::vector<TCnSd> TSampleDist;
    TSampleDist control;
    TSampleDist tumor;
    bool invalidCNV = false;
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
      if ((!std::isfinite(rdcn[i])) || (rdcn[i] == -1)) {
	invalidCNV = true;
	break;
      }
      if ((germline) || (c.controlSet.find(hdr->samples[i]) != c.controlSet.end())) {
	// Control or population genomics
	control.push_back(std::make_pair(rdcn[i], rdsd[i]));
      } else if ((!germline) && (c.tumorSet.find(hdr->samples[i]) != c.tumorSet.end())) {
	// Tumor
	tumor.push_back(std::make_pair(rdcn[i], rdsd[i]));
      }
    }
    if (invalidCNV) continue;

    // Classify
    if (!germline) {
      // Somatic mode
      double bestCnOffset = 0;
      bool somaticcnv = false;
      double lowestp = 1;
      for(uint32_t i = 0; i < tumor.size(); ++i) {
	bool germcnv = false;
	double highestprob = 0;
	double tcnoffset = -1;
	for(uint32_t k = 0; k < control.size(); ++k) {
	  boost::math::normal s1(control[k].first, control[k].second);
	  double prob1 = boost::math::pdf(s1, tumor[i].first);
	  boost::math::normal s2(tumor[i].first, tumor[i].second);
	  double prob2 = boost::math::pdf(s2, control[k].first);
	  double prob = std::max(prob1, prob2);
	  if (prob > c.pgerm) germcnv = true;
	  else {
	    // Among all controls, take highest p-value (most likely germline CNV)
	    if (prob > highestprob) highestprob = prob;
	  }
	  double cndiff = std::abs(tumor[i].first - control[k].first);
	  if (cndiff < c.cn_offset) germcnv = true;
	  else {
	    // Among all controls, take smallest CN difference
	    if ((tcnoffset == -1) || (cndiff < tcnoffset)) tcnoffset = cndiff;
	  }
	}
	// Among all tumors take best CN difference and lowest p-value
	if (!germcnv) {
	  somaticcnv = true;
	  if ((highestprob < lowestp) && (tcnoffset > bestCnOffset)) {
	    lowestp = highestprob;
	    bestCnOffset = tcnoffset;
	  }
	}
      }
      if (!somaticcnv) continue;
      _remove_info_tag(hdr_out, rec, "SOMATIC");
      bcf_update_info_flag(hdr_out, rec, "SOMATIC", NULL, 1);
      float pgerm = (float) lowestp;
      _remove_info_tag(hdr_out, rec, "PGERM");
      bcf_update_info_float(hdr_out, rec, "PGERM", &pgerm, 1);
      float cndiv = (float) bestCnOffset;
      _remove_info_tag(hdr_out, rec, "CNDIFF");
      bcf_update_info_float(hdr_out, rec, "CNDIFF", &cndiv, 1);
    } else {
      // Correct CN shift
      int32_t cnmain = 0;
      {
	std::vector<int32_t> cncount(MAX_CN, 0);
	{
	  bool validsite = true;
	  boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance> > acc;
	  for(uint32_t k = 0; k < control.size(); ++k) {
	    if ((boost::math::isinf(control[k].first)) || (boost::math::isnan(control[k].first))) validsite = false;
	    else acc(boost::math::round(control[k].first) - control[k].first);
	  }
	  if (!validsite) continue;
	  double cnshift = boost::accumulators::mean(acc);
	  float cnshiftval = cnshift;
	  _remove_info_tag(hdr_out, rec, "CNSHIFT");
	  bcf_update_info_float(hdr_out, rec, "CNSHIFT", &cnshiftval, 1);
	  for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
	    rdcn[i] += cnshift;
	    cnval[i] = boost::math::round(rdcn[i]);
	    if ((cnval[i] >= 0) && (cnval[i] < MAX_CN)) ++cncount[cnval[i]];
	  }
	}
	
	// Find max CN
	for(uint32_t i = 1; i < MAX_CN; ++i) {
	  if (cncount[i] > cncount[cnmain]) cnmain = i;
	}
      }

      // Calculate SD
      boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance> > accLocal;
      for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
	if (cnval[i] == cnmain) accLocal(rdcn[i]);
      }
      double sd = sqrt(boost::accumulators::variance(accLocal));
      if (sd < 0.025) sd = 0.025;
      float cnsdval = sd;
      _remove_info_tag(hdr_out, rec, "CNSD");
      bcf_update_info_float(hdr_out, rec, "CNSD", &cnsdval, 1);
      if (cnsdval > c.maxsd) continue;
      
      // Re-compute CNLs
      std::vector<std::string> ftarr(bcf_hdr_nsamples(hdr));
      int32_t altqual = 0;
      int32_t altcount = 0;
      for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
	int32_t qval = _computeCNLs(c, rdcn[i], sd, cnl, gqval, i);
	if (cnval[i] != c.ploidy) {
	  altqual += qval;
	  ++altcount;
	}
	if (gqval[i] < 15) ftarr[i] = "LowQual";
	else ftarr[i] = "PASS";
      }
      if (altcount == 0) continue;
      altqual /= altcount;
      if (altqual < c.qual) continue;
      if (altqual > 10000) altqual = 10000;
      
      // Update QUAL and FILTER
      rec->qual = altqual;
      int32_t tmpi = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "PASS");
      if (rec->qual < 15) tmpi = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "LowQual");
      bcf_update_filter(hdr_out, rec, &tmpi, 1);

      // Update GT fields
      std::vector<const char*> strp(bcf_hdr_nsamples(hdr));
      std::transform(ftarr.begin(), ftarr.end(), strp.begin(), cstyle_str());
      bcf_update_format_int32(hdr_out, rec, "CN", cnval, bcf_hdr_nsamples(hdr));
      bcf_update_format_float(hdr_out, rec, "CNL",  cnl, bcf_hdr_nsamples(hdr) * MAX_CN);
      bcf_update_format_int32(hdr_out, rec, "GQ", gqval, bcf_hdr_nsamples(hdr));
      bcf_update_format_string(hdr_out, rec, "FT", &strp[0], bcf_hdr_nsamples(hdr));
      bcf_update_format_float(hdr_out, rec, "RDCN",  rdcn, bcf_hdr_nsamples(hdr));
    }
    bcf_write1(ofile, hdr_out, rec);
  }
  bcf_destroy(rec);

  // Clean-up
  if (svend != NULL) free(svend);
  if (svt != NULL) free(svt);
  if (gqval != NULL) free(gqval);
  if (cnval != NULL) free(cnval);
  if (cnl != NULL) free(cnl);
  if (rdcn != NULL) free(rdcn);
  if (rdsd != NULL) free(rdsd);

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


int classify(int argc, char **argv) {
  ClassifyConfig c;

  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("filter,f", boost::program_options::value<std::string>(&c.filter)->default_value("somatic"), "Filter mode (somatic, germline)")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "Filtered CNV BCF output file")
    ("minsize,m", boost::program_options::value<int32_t>(&c.minsize)->default_value(1000), "min. CNV size")
    ("maxsize,n", boost::program_options::value<int32_t>(&c.maxsize)->default_value(500000000), "max. CNV size")
    ("pass,p", "Filter sites for PASS")
    ;

  // Define somatic options
  boost::program_options::options_description somatic("Somatic options");
  somatic.add_options()
    ("samples,s", boost::program_options::value<boost::filesystem::path>(&c.samplefile), "Two-column sample file listing sample name and tumor or control")
    ("pgerm,e", boost::program_options::value<float>(&c.pgerm)->default_value(0.001), "probability germline")
    ("cn-offset,t", boost::program_options::value<float>(&c.cn_offset)->default_value(0.2), "min. CN offset")
    ;

  // Define germline options
  boost::program_options::options_description germline("Germline options");
  germline.add_options()
    ("ploidy,y", boost::program_options::value<uint16_t>(&c.ploidy)->default_value(2), "baseline ploidy")
    ("qual,q", boost::program_options::value<int32_t>(&c.qual)->default_value(50), "min. site quality")
    ("maxsd,x", boost::program_options::value<float>(&c.maxsd)->default_value(0.15), "max. population SD")
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

  return classifyRun(c);
}

}

#endif
