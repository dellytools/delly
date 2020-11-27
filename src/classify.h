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
#include <boost/progress.hpp>

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
  float gq;
  float pgerm;
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
  htsFile *ofile = hts_open(c.outfile.string().c_str(), "wb");
  bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);
  if (c.filter == "somatic") {
    bcf_hdr_remove(hdr_out, BCF_HL_INFO, "SOMATIC");
    bcf_hdr_append(hdr_out, "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic copy-number variant.\">");
    bcf_hdr_remove(hdr_out, BCF_HL_INFO, "PGERM");
    bcf_hdr_append(hdr_out, "##INFO=<ID=PGERM,Number=1,Type=Float,Description=\"Probability of being germline.\">");
  }
  if (bcf_hdr_write(ofile, hdr_out) != 0) std::cerr << "Error: Failed to write BCF header!" << std::endl;

  // VCF fields
  int32_t nsvend = 0;
  int32_t* svend = NULL;
  int32_t nsvt = 0;
  char* svt = NULL;
  int ngq = 0;
  int32_t* gq = NULL;
  int ncn = 0;
  int32_t* cn = NULL;
  int nrdcn = 0;
  float* rdcn = NULL;
  int nrdsd = 0;
  float* rdsd = NULL;
  bool germline = false;
  if (c.filter == "germline") germline = true;

  // Parse BCF
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Filtering VCF/BCF file" << std::endl;
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
    bcf_get_format_int32(hdr, rec, "GQ", &gq, &ngq);
    bcf_get_format_int32(hdr, rec, "CN", &cn, &ncn);
    bcf_get_format_float(hdr, rec, "RDCN", &rdcn, &nrdcn);
    bcf_get_format_float(hdr, rec, "RDSD", &rdsd, &nrdsd);

    typedef std::pair<float, float> TCnSd;
    typedef std::vector<TCnSd> TSampleDist;
    TSampleDist control;
    TSampleDist tumor;
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
      if ((germline) || (c.controlSet.find(hdr->samples[i]) != c.controlSet.end())) {
	// Control or population genomics
	control.push_back(std::make_pair(rdcn[i], rdsd[i]));
      } else if ((!germline) && (c.tumorSet.find(hdr->samples[i]) != c.tumorSet.end())) {
	// Tumor
	tumor.push_back(std::make_pair(rdcn[i], rdsd[i]));
      }
    }

    // Classify
    if (!germline) {
      // Somatic mode
      bool somaticcnv = false;
      double worstp = 0;
      for(uint32_t i = 0; i < tumor.size(); ++i) {
	bool germcnv = false;
	double highestprob = 0;
	for(uint32_t k = 0; k < control.size(); ++k) {
	  boost::math::normal s(control[k].first, control[k].second);
	  double prob = boost::math::pdf(s, tumor[i].first);
	  if (prob > c.pgerm) germcnv = true;
	  else {
	    if (prob > highestprob) highestprob = prob;
	  }
	}
	if (!germcnv) {
	  somaticcnv = true;
	  if (highestprob > worstp) worstp = highestprob;
	}
      }
      if (!somaticcnv) continue;
      _remove_info_tag(hdr_out, rec, "SOMATIC");
      bcf_update_info_flag(hdr_out, rec, "SOMATIC", NULL, 1);
      float pgerm = (float) worstp;
      _remove_info_tag(hdr_out, rec, "PGERM");
      bcf_update_info_float(hdr_out, rec, "PGERM", &pgerm, 1);
    } else {
      /*
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
      if ((af>0) && (gqaltmed >= c.gq) && (gqrefmed >= c.gq) && (raltvarmed >= c.altaf) && (genotypeRatio >= c.ratiogeno)) {
      if ((std::string(svt)=="DEL") && (rdRatio > c.rddel)) continue;
      if ((std::string(svt)=="DUP") && (rdRatio < c.rddup)) continue;
      if ((std::string(svt)!="DEL") && (std::string(svt)!="DUP") && (rrefvarpercentile > 0)) continue;
      _remove_info_tag(hdr_out, rec, "RDRATIO");
      bcf_update_info_float(hdr_out, rec, "RDRATIO", &rdRatio, 1);
      */
    }
    bcf_write1(ofile, hdr_out, rec);
  }
  bcf_destroy(rec);

  // Clean-up
  if (svend != NULL) free(svend);
  if (svt != NULL) free(svt);
  if (gq != NULL) free(gq);
  if (cn != NULL) free(cn);
  if (rdcn != NULL) free(rdcn);
  if (rdsd != NULL) free(rdsd);

  // Close output VCF
  bcf_hdr_destroy(hdr_out);
  hts_close(ofile);

  // Build index
  bcf_index_build(c.outfile.string().c_str(), 14);

  // Close VCF
  bcf_hdr_destroy(hdr);
  bcf_close(ifile);

  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;

  return 0;
}


int classify(int argc, char **argv) {
  ClassifyConfig c;

  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("filter,f", boost::program_options::value<std::string>(&c.filter)->default_value("somatic"), "Filter mode (somatic, germline)")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("cnv.bcf"), "Filtered CNV BCF output file")
    ("minsize,m", boost::program_options::value<int32_t>(&c.minsize)->default_value(1000), "min. CNV size")
    ("maxsize,n", boost::program_options::value<int32_t>(&c.maxsize)->default_value(500000000), "max. CNV size")
    ("pass,p", "Filter sites for PASS")
    ;

  // Define somatic options
  boost::program_options::options_description somatic("Somatic options");
  somatic.add_options()
    ("samples,s", boost::program_options::value<boost::filesystem::path>(&c.samplefile), "Two-column sample file listing sample name and tumor or control")
    ("pgerm,e", boost::program_options::value<float>(&c.pgerm)->default_value(0.05), "probability germline")
    ;

  // Define germline options
  boost::program_options::options_description germline("Germline options");
  germline.add_options()
    ("gq,q", boost::program_options::value<float>(&c.gq)->default_value(15), "min. median GQ for carriers and non-carriers")
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
    std::cout << std::endl;
    std::cout << "Usage: delly " << argv[0] << " [OPTIONS] <input.bcf>" << std::endl;
    std::cout << visible_options << "\n";
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

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  std::cout << "delly ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  return classifyRun(c);
}

}

#endif
