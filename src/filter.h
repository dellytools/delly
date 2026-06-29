#ifndef FILTER_H
#define FILTER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstdlib>
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
#include "popgen.h"

namespace torali
{
  
  
  struct FilterConfig {
    bool filterForPass;
    bool hasSampleFile;
    bool softFilter;
    bool noRefine;
    bool noCollapse;
    int32_t minsize;
    int32_t maxsize;
    int32_t coverage;
    int32_t qualthres;
    int32_t rdist;
    int32_t rminshared;
    uint32_t maxiter;
    double epsilon;
    float ratiogeno;
    float altaf;
    float controlcont;
    float gq;
    float genogq;
    float hwe;
    float exhet;
    float minrsq;
    float rsize;
    float rcorr;
    float rddel;
    float rddup;
    std::string filter;
    std::set<std::string> tumorSet;
    std::set<std::string> controlSet;
    boost::filesystem::path outfile;
    boost::filesystem::path samplefile;
    boost::filesystem::path vcffile;
  };
  
  
  // One buffered record for the redundant SV search
  struct RedRec {
    bcf1_t* rec;
    std::vector<int8_t> dos;   // GTs (0/1/2, -1 missing)
    std::string svtype;
    int32_t spos;
    int32_t epos;
    int32_t len;
    int32_t ac;
    int32_t ncalled;
    float qual;
    bool precise;
    bool eligible;
    bool redundant;
  };
  
  // Redundant record sort
  inline bool
  _redBetter(RedRec const& a, RedRec const& b) {
    if (a.precise != b.precise) return a.precise;
    if (a.ncalled != b.ncalled) return (a.ncalled > b.ncalled);
    if (a.ac != b.ac) return (a.ac > b.ac);
    return (a.qual > b.qual);
  }

  // Proximity of two SV sites
  inline bool
  _redProximal(RedRec const& a, RedRec const& b, int32_t const rdist, float const rsize) {
    if (a.svtype != b.svtype) return false;
    if (a.svtype == "BND") return false; // ignore BND
    if (std::abs(a.spos - b.spos) > rdist) return false;
    if (a.svtype != "INS") {
      if (std::abs(a.epos - b.epos) > rdist) return false;
    }
    int32_t mn = std::min(a.len, b.len);
    int32_t mx = std::max(a.len, b.len);
    if ((mx > 0) && ((float) mn / (float) mx < rsize)) return false;
    return true;
  }

  // remove redundant SVs
  inline void
  _flushRedundancy(std::vector<RedRec>& win, int32_t const flushBelow, bool const flushAll, htsFile* ofile, bcf_hdr_t* hdr_out, bool const softFilter, int32_t const redId) {
    std::vector<RedRec> keep;
    keep.reserve(win.size());
    for(std::vector<RedRec>::iterator it = win.begin(); it != win.end(); ++it) {
      if ((!flushAll) && (it->spos >= flushBelow)) {
	keep.push_back(*it);
	continue;
      }
      bool drop = false;
      if ((it->eligible) && (it->redundant)) {
	if (softFilter) bcf_update_filter(hdr_out, it->rec, const_cast<int32_t*>(&redId), 1);
	else drop = true;
      }
      if (!drop) bcf_write1(ofile, hdr_out, it->rec);
      bcf_destroy(it->rec);
    }
    win.swap(keep);
  }

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
      const char* popTags[] = {"AFmle","ACmle","GFmle","FIC","RSQ","HWEpval"};
      for(int t=0; t<6; ++t) bcf_hdr_remove(hdr_out, BCF_HL_INFO, popTags[t]);
      bcf_hdr_append(hdr_out, "##INFO=<ID=AFmle,Number=1,Type=Float,Description=\"AF estimated from genotype likelihoods.\">");
      bcf_hdr_append(hdr_out, "##INFO=<ID=ACmle,Number=1,Type=Integer,Description=\"Allele count estimated from genotype likelihoods.\">");
      bcf_hdr_append(hdr_out, "##INFO=<ID=GFmle,Number=G,Type=Float,Description=\"GT frequencies estimated from GLs.\">");
      bcf_hdr_append(hdr_out, "##INFO=<ID=FIC,Number=1,Type=Float,Description=\"Inbreeding coefficient.\">");
      bcf_hdr_append(hdr_out, "##INFO=<ID=RSQ,Number=1,Type=Float,Description=\"Imputation quality R^2.\">");
      bcf_hdr_append(hdr_out, "##INFO=<ID=HWEpval,Number=1,Type=Float,Description=\"HWE likelihood-ratio test p-value.\">");
      if (c.softFilter) {
	bcf_hdr_append(hdr_out, "##FILTER=<ID=RedundantSV,Description=\"Redundant SV site.\">");
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
    int npl = 0;
    int32_t* pl = NULL;
    bool germline = false;
    if (c.filter == "germline") germline = true;

    // Redundant-SV collapse window (germline only)
    bool collapse = ((germline) && (!c.noCollapse));
    std::vector<RedRec> redWin;
    int32_t redChrom = -1;
    int32_t redId = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "RedundantSV");

    // Parse BCF
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Filtering VCF/BCF file" << std::endl;
    bcf1_t* rec = bcf_init1();
    while (bcf_read(ifile, hdr, rec) == 0) {
      bcf_unpack(rec, BCF_UN_INFO);

      // Advance the collapse window: flush records that can no longer have a redundant partner
      if (collapse) {
	if (rec->rid != redChrom) { _flushRedundancy(redWin, 0, true, ofile, hdr_out, c.softFilter, redId); redChrom = rec->rid; }
	else _flushRedundancy(redWin, rec->pos - c.rdist, false, ofile, hdr_out, c.softFilter, redId);
      }

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
	
	// Population refinement from GLs
	bool refined = false;
	double hwepvalStore = 1;
	double ficStore = 0;
	double rsqStore = 0;
	if ((germline) && (!c.noRefine) && (rec->n_allele == 2) && (bcf_get_format_int32(hdr, rec, "PL", &pl, &npl) > 0)) {
	  int32_t nsmpl = bcf_hdr_nsamples(hdr);
	  int32_t plStride = npl / nsmpl;
	  if (plStride >= 3) {
	    typedef std::vector<double> TGLs;
	    std::vector<TGLs> glVector;
	    glVector.reserve(nsmpl);
	    for (int i = 0; i < nsmpl; ++i) {
	      if ((bcf_gt_allele(gt[i*2]) == -1) || (bcf_gt_allele(gt[i*2 + 1]) == -1)) continue;
	      int32_t* plp = pl + i * plStride;
	      if ((plp[0] == bcf_int32_missing) || (plp[0] == bcf_int32_vector_end)) continue;
	      TGLs glTriple(3);
	      for(int k = 0; k < 3; ++k) glTriple[k] = std::pow((double) 10.0, (double) -plp[k] / 10.0);
	      glVector.push_back(glTriple);
	    }
	    if (!glVector.empty()) {
	      refined = true;
	      double hweAF[2];
	      hweAF[0] = 0.5;
	      hweAF[1] = 0.5;
	      _estBiallelicAF(c, glVector, hweAF);
	      double mleGTFreq[3];
	      mleGTFreq[0] = 0;
	      mleGTFreq[1] = 0;
	      mleGTFreq[2] = 0;
	      _estBiallelicGTFreq(c, glVector, mleGTFreq);
	      double Fic = 0;
	      _estBiallelicFIC(glVector, hweAF, Fic);
	      double rsq = 0;
	      _estBiallelicRSQ(glVector, hweAF, rsq);
	      double pval = 1;
	      _estBiallelicHWE_LRT(glVector, hweAF, mleGTFreq, pval);
	      hwepvalStore = pval;
	      ficStore = Fic;
	      rsqStore = rsq;

	      // Posterior GQ (only GQ and missingness change)
	      bool gqIsInt = (_getFormatType(hdr, "GQ") == BCF_HT_INT);
	      for (int i = 0; i < nsmpl; ++i) {
		if ((bcf_gt_allele(gt[i*2]) == -1) || (bcf_gt_allele(gt[i*2 + 1]) == -1)) continue;
		int32_t* plp = pl + i * plStride;
		if ((plp[0] == bcf_int32_missing) || (plp[0] == bcf_int32_vector_end)) continue;
		double pp[3];
		int bestIdx = 0;
		for(int k = 0; k < 3; ++k) {
		  pp[k] = mleGTFreq[k] * std::pow((double) 10.0, (double) -plp[k] / 10.0);
		  if (plp[k] < plp[bestIdx]) bestIdx = k;
		}
		double sumPP = pp[0] + pp[1] + pp[2];
		double sampleGq = 0;
		if (sumPP > 0) sampleGq = (double) -10.0 * std::log10((double) 1.0 - pp[bestIdx] / sumPP);
		if (sampleGq > 99) sampleGq = 99;
		if (sampleGq < 0) sampleGq = 0;
		if (sampleGq < c.genogq) {
		  gt[i*2] = bcf_gt_missing;
		  gt[i*2 + 1] = bcf_gt_missing;
		}
		if (gqIsInt) gq[i] = (int32_t) (sampleGq + 0.5);
		else gqf[i] = (float) sampleGq;
	      }
	      bcf_update_genotypes(hdr_out, rec, gt, nsmpl * 2);
	      if (gqIsInt) bcf_update_format_int32(hdr_out, rec, "GQ", gq, nsmpl);
	      else bcf_update_format_float(hdr_out, rec, "GQ", gqf, nsmpl);

	      // Annotations
	      float afmle = (float) hweAF[1];
	      int32_t acmle = (int32_t) boost::math::iround(hweAF[1] * 2.0 * (double) glVector.size());
	      float gfmle[3]; gfmle[0] = (float) mleGTFreq[0]; gfmle[1] = (float) mleGTFreq[1]; gfmle[2] = (float) mleGTFreq[2];
	      float ficv = (float) Fic;
	      float rsqv = (float) rsq;
	      float hwev = (float) pval;
	      _remove_info_tag(hdr_out, rec, "AFmle");
	      bcf_update_info_float(hdr_out, rec, "AFmle", &afmle, 1);
	      _remove_info_tag(hdr_out, rec, "ACmle");
	      bcf_update_info_int32(hdr_out, rec, "ACmle", &acmle, 1);
	      _remove_info_tag(hdr_out, rec, "GFmle");
	      bcf_update_info_float(hdr_out, rec, "GFmle", &gfmle, 3);
	      _remove_info_tag(hdr_out, rec, "FIC");
	      bcf_update_info_float(hdr_out, rec, "FIC", &ficv, 1);
	      _remove_info_tag(hdr_out, rec, "RSQ");
	      bcf_update_info_float(hdr_out, rec, "RSQ", &rsqv, 1);
	      _remove_info_tag(hdr_out, rec, "HWEpval");
	      bcf_update_info_float(hdr_out, rec, "HWEpval", &hwev, 1);
	    }
	  }
	}
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
	  if ((refined) && (c.hwe > 0) && (ficStore < 0) && (hwepvalStore < c.hwe)) failgerm = true;
	  if ((refined) && (c.exhet > 0) && (ficStore < -c.exhet)) failgerm = true;
	  if ((refined) && (c.minrsq > 0) && (rsqStore < c.minrsq)) failgerm = true;
	  if (!failgerm) {
	    _remove_info_tag(hdr_out, rec, "RDRATIO");
	    bcf_update_info_float(hdr_out, rec, "RDRATIO", &rdRatio, 1);
	    if (collapse) {
	      // Enqueue as a collapse-eligible record
	      RedRec rr;
	      rr.rec = bcf_dup(rec);
	      rr.svtype = std::string(svt);
	      rr.spos = rec->pos;
	      rr.epos = (svend != NULL) ? (*svend) : rec->pos;
	      rr.len = (rr.svtype == "INS") ? inslenVal : std::abs(svlen);
	      rr.qual = rec->qual;
	      rr.precise = precise;
	      rr.eligible = true;
	      rr.redundant = false;
	      int32_t nsmpl = bcf_hdr_nsamples(hdr);
	      rr.dos.assign(nsmpl, (int8_t) -1);
	      rr.ac = 0; rr.ncalled = 0;
	      for (int i = 0; i < nsmpl; ++i) {
		int a0 = bcf_gt_allele(gt[i*2]);
		int a1 = bcf_gt_allele(gt[i*2 + 1]);
		if ((a0 >= 0) && (a1 >= 0)) {
		  int8_t d = (int8_t) ((a0 > 0 ? 1 : 0) + (a1 > 0 ? 1 : 0));
		  rr.dos[i] = d;
		  ++rr.ncalled;
		  rr.ac += d;
		}
	      }
	      // Collapse
	      for(std::vector<RedRec>::iterator wit = redWin.begin(); wit != redWin.end(); ++wit) {
		if ((!wit->eligible) || (wit->redundant)) continue;
		if (_redProximal(rr, *wit, c.rdist, c.rsize)) {
		  double r2 = _dosageR2(rr.dos, wit->dos, c.rminshared);
		  if (r2 >= (double) c.rcorr) {
		    if (_redBetter(rr, *wit)) wit->redundant = true;
		    else {
		      rr.redundant = true;
		      break;
		    }
		  }
		}
	      }
	      redWin.push_back(rr);
	    } else bcf_write1(ofile, hdr_out, rec);
	  } else if (c.softFilter) {
	    int32_t tmpi = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "FailGermline");
	    bcf_update_filter(hdr_out, rec, &tmpi, 1);
	    if (collapse) {
	      RedRec rr;
	      rr.rec = bcf_dup(rec);
	      rr.svtype = std::string(svt);
	      rr.spos = rec->pos;
	      rr.epos = (svend != NULL) ? (*svend) : rec->pos;
	      rr.len = 0; rr.qual = rec->qual; rr.precise = precise;
	      rr.eligible = false; rr.redundant = false; rr.ac = 0; rr.ncalled = 0;
	      redWin.push_back(rr);
	    } else bcf_write1(ofile, hdr_out, rec);
	  }
	}
      } else if (c.softFilter) {
	int32_t tmpi = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "FailDellyFilter");
	bcf_update_filter(hdr_out, rec, &tmpi, 1);
	if (collapse) {
	  RedRec rr;
	  rr.rec = bcf_dup(rec);
	  rr.svtype = std::string(svt);
	  rr.spos = rec->pos;
	  rr.epos = (svend != NULL) ? (*svend) : rec->pos;
	  rr.len = 0; rr.qual = rec->qual; rr.precise = false;
	  rr.eligible = false; rr.redundant = false; rr.ac = 0; rr.ncalled = 0;
	  redWin.push_back(rr);
	} else bcf_write1(ofile, hdr_out, rec);
      }
    }
    // Flush any records still buffered in the collapse window
    if (collapse) _flushRedundancy(redWin, 0, true, ofile, hdr_out, c.softFilter, redId);
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
    if (pl != NULL) free(pl);
    
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
      ("genogq,j", boost::program_options::value<float>(&c.genogq)->default_value(15), "set genotypes below this posterior GQ to missing")
      ("hwe,w", boost::program_options::value<float>(&c.hwe)->default_value(0.000001), "min. HWE p-value for excess-heterozygosity (one-sided; set 0 to disable)")
      ("exhet", boost::program_options::value<float>(&c.exhet)->default_value(0), "max. excess heterozygosity (fail if FIC < -exhet; 0 disables)")
      ("minrsq", boost::program_options::value<float>(&c.minrsq)->default_value(0), "min. imputation R-squared RSQ (0 disables)")
      ("no-refine", boost::program_options::bool_switch(&c.noRefine), "disable population refinement")
      ("no-collapse", boost::program_options::bool_switch(&c.noCollapse), "disable redundant-SV collapse")
      ("rdist", boost::program_options::value<int32_t>(&c.rdist)->default_value(250), "max. breakpoint distance for redundant SVs")
      ("rsize", boost::program_options::value<float>(&c.rsize)->default_value(0.8), "min. size ratio for redundant SVs")
      ("rcorr", boost::program_options::value<float>(&c.rcorr)->default_value(0.8), "min. genotype r-squared for redundant SVs")
      ("rminshared", boost::program_options::value<int32_t>(&c.rminshared)->default_value(20), "min. samples to assess GT concordance")
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
    
    // EM parameters for population refinement
    c.epsilon = 1e-20;
    c.maxiter = 1000;
    
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
