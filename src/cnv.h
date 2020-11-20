#ifndef CNV_H
#define CNV_H

#include <boost/filesystem.hpp>
#include <boost/multi_array.hpp>
#include <boost/progress.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/sam.h>


namespace torali
{

  struct CNV {
    int32_t chr;
    int32_t start;
    int32_t end;
    int32_t ciposlow;
    int32_t ciposhigh;
    int32_t ciendlow;
    int32_t ciendhigh;
    double cn;
    double zscore;
    double mappable;
    
    CNV(int32_t const c, int32_t const s, int32_t const e, int32_t const cil, int32_t const cih, int32_t const cel, int32_t ceh, double const estcn, double const zs, double const mp) : chr(c), start(s), end(e), ciposlow(cil), ciposhigh(cih), ciendlow(cel), ciendhigh(ceh), cn(estcn), zscore(zs), mappable(mp) {}
  };

  template<typename TCNV>
  struct SortCNVs : public std::binary_function<TCNV, TCNV, bool>
  {
    inline bool operator()(TCNV const& sv1, TCNV const& sv2) {
      return ((sv1.chr<sv2.chr) || ((sv1.chr==sv2.chr) && (sv1.start<sv2.start)) || ((sv1.chr==sv2.chr) && (sv1.start==sv2.start) && (sv1.end<sv2.end)) || ((sv1.chr==sv2.chr) && (sv1.start==sv2.start) && (sv1.end==sv2.end) && (sv1.cn < sv2.cn)));
    }
  };


  struct BpCNV {
    int32_t start;
    int32_t end;
    double zscore;

    BpCNV(int32_t const s, int32_t const e, double const z) : start(s), end(e), zscore(z) {}
  };
  
  template<typename TConfig, typename TGcBias, typename TCoverage>
  inline void
  callCNVs(TConfig const& c, std::pair<uint32_t, uint32_t> const& gcbound, std::vector<uint16_t> const& gcContent, std::vector<uint16_t> const& uniqContent, TGcBias const& gcbias, TCoverage const& cov, bam_hdr_t const* hdr, int32_t const refIndex, std::vector<CNV>& chrcnv) {

    // Parameters
    int32_t smallestWin = 100;
    int32_t biggestWin = 15000;
    uint32_t chain = 10;

    // Find breakpoints
    std::vector<BpCNV> bpmax;
    if (bpmax.empty()) {
      // Scanning window sizes
      std::vector<int32_t> winsize;
      int32_t wsize = smallestWin;
      while (wsize < biggestWin) {
	winsize.push_back(wsize);
	wsize *= 2;
      }

      // Iterate window sizes
      typedef int32_t TCnVal;
      typedef std::vector<TCnVal> TCN;
      typedef std::vector<int32_t> TChrPos;
      std::vector<BpCNV> bpvec;
      for(uint32_t idx = 0; idx < winsize.size(); ++idx) {
	uint32_t idxOffset = winsize[idx] / winsize[0];
	//std::cerr << idx << ',' << winsize[idx] << ',' << idxOffset << ',' << bpvec.size() << ',' << hdr->target_len[refIndex] << std::endl;
	TCN cnvec;
	TChrPos wpos;
	uint32_t wstart = 0;
	while(wstart < hdr->target_len[refIndex]) {
	  double covsum = 0;
	  double expcov = 0;
	  int32_t winlen = 0;
	  uint32_t pos = wstart;
	  while ((winlen < winsize[idx]) && (pos < hdr->target_len[refIndex])) {
	    if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	      covsum += cov[pos];
	      expcov += gcbias[gcContent[pos]].coverage;
	      ++winlen;
	    }
	    ++pos;
	  }
	  if (winlen == winsize[idx]) {
	    // Full window
	    cnvec.push_back((int32_t) boost::math::round(c.ploidy * covsum / expcov * 100.0));
	    wpos.push_back(wstart);
	  }
	  wstart = pos;
	}
	
	// Identify breakpoints
	TCN pre(chain, -1);
	TCN suc(chain, -1);
	TChrPos prep(chain, 0);
	TChrPos sucp(chain, 0);
	uint32_t idxbp = 0;
	for(uint32_t k = 0; k < cnvec.size(); ++k) {
	  if (k < chain) {
	    pre[k % chain] = cnvec[k];
	    prep[k % chain] = wpos[k];
	    if (k + 1 < cnvec.size()) {
	      if (idx == 0 ) bpvec.push_back(BpCNV(wpos[k], wpos[k+1], 0));
	      else idxbp += idxOffset;
	    }
	  } else if (k < 2 * chain) {
	    suc[k % chain] = cnvec[k];
	    sucp[k % chain] = wpos[k];
	  } else {
	    // Midpoint
	    TCnVal val = suc[k%chain];	  
	    int32_t pos = sucp[k%chain];
	    int32_t posNext = sucp[(k+1)%chain];
	    suc[k%chain] = cnvec[k];
	    sucp[k%chain] = wpos[k];
	    
	    // Debug
	    //for(uint32_t m = 0; m < pre.size(); ++m) std::cerr << prep[m] << '\t' << pre[m] << std::endl;
	    //std::cerr << "M:" << pos << '\t' << val << std::endl;
	    //for(uint32_t m = 0; m < suc.size(); ++m) std::cerr << sucp[m] << '\t' << suc[m] << std::endl;
	    
	    // Any shift in CN?
	    boost::accumulators::accumulator_set<TCnVal, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> accpre;
	    boost::accumulators::accumulator_set<TCnVal, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> accsuc;
	    for(uint32_t m = 0; m < pre.size(); ++m) accpre(pre[m]);
	    for(uint32_t m = 0; m < suc.size(); ++m) accsuc(suc[m]);
	    double diff = std::abs(boost::accumulators::mean(accsuc) - boost::accumulators::mean(accpre));
	    // Breakpoint candidate
	    double zscore = 0;
	    if ((diff > c.stringency * sqrt(boost::accumulators::variance(accpre))) && (diff > c.stringency * sqrt(boost::accumulators::variance(accsuc)))) {
	      zscore = diff / std::max(sqrt(boost::accumulators::variance(accpre)), sqrt(boost::accumulators::variance(accsuc)));
	    }
	    if (idx == 0) bpvec.push_back(BpCNV(pos, posNext, zscore));
	    else {
	      for(uint32_t sub = idxbp; sub < idxbp + idxOffset; ++sub) bpvec[sub].zscore += zscore;
	      idxbp += idxOffset;
	    }
	    pre[k%chain] = val;
	    prep[k%chain] = pos;
	  }
	}
      }
    
      // Local maxima
      if (bpvec.size()) {
	int32_t pos = bpvec[0].start;
	int32_t posNext = bpvec[0].end;
	double bestDiff = bpvec[0].zscore;
	for(uint32_t n = 1; n < bpvec.size(); ++n) {
	  //std::cerr << "B:" << bpvec[n].start << '-' << bpvec[n].end << ':' << bpvec[n].zscore << std::endl;
	  if (bpvec[n].zscore == 0) {
	    if (bestDiff != 0) {
	      //std::cerr << "M:" << pos << '-' << posNext << ':' << bestDiff << std::endl;
	      bpmax.push_back(BpCNV(pos, posNext, bestDiff));
	      pos = bpvec[n].start;
	      posNext = bpvec[n].end;
	      bestDiff = bpvec[n].zscore;
	    }
	  } else {
	    if (bpvec[n].zscore > bestDiff) {
	      // Replace local max
	      pos = bpvec[n].start;
	      posNext = bpvec[n].end;
	      bestDiff = bpvec[n].zscore;
	    } else if (bpvec[n].zscore == bestDiff) {
	      // Extend local max
	      posNext = bpvec[n].end;
	    }
	  }
	}
      }
    }

    // Breakpoints
    for(uint32_t n = 0; n <= bpmax.size(); ++n) {
      int32_t cil = 0;
      int32_t cih = 0;
      double zscorepre = 0;
      if (n > 0) {
	cil = bpmax[n-1].start;
	cih = bpmax[n-1].end;
	zscorepre = bpmax[n-1].zscore;
      }
      int32_t cel = hdr->target_len[refIndex] - 1;
      int32_t ceh = hdr->target_len[refIndex] - 1;
      double zscoresuc = 0;
      if (n < bpmax.size()) {
	cel = bpmax[n].start;
	ceh = bpmax[n].end;
	zscoresuc = bpmax[n].zscore;
      }
      int32_t cnvstart = (int32_t) ((cil + cih)/2);
      int32_t cnvend = (int32_t) ((cel + ceh)/2);
      double covsum = 0;
      double expcov = 0;
      double obsexp = 0;
      int32_t winlen = 0;
      int32_t pos = cnvstart;
      while((pos < cnvend) && (pos < (int32_t) hdr->target_len[refIndex])) {
	if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	  covsum += cov[pos];
	  obsexp += gcbias[gcContent[pos]].obsexp;
	  expcov += gcbias[gcContent[pos]].coverage;
	  ++winlen;
	}
	++pos;
      }
      obsexp /= (double) winlen;
      double cn = c.ploidy * covsum / expcov;
      double zsc = (zscorepre + zscoresuc) / 2.0;
      double mp = (double) winlen / (double) (cnvend - cnvstart);
      chrcnv.push_back(CNV(refIndex, cnvstart, cnvend, cil, cih, cel, ceh, cn, zsc, mp));
      //std::cerr << hdr->target_name[refIndex] << '\t' << cnvstart << '\t' << cnvend << '\t' << '(' << cil << ',' << cih << ')' << '\t' << '(' << cel << ',' << ceh << ')' << '\t' << cn << '\t' << zsc << '\t' << mp << std::endl;
    }    
  }

  template<typename TConfig>
  inline void
  mergeCNVs(TConfig const& c, std::vector<CNV>& chrcnv, std::vector<CNV>& cnvs) {
    // Merge neighboring segments if too similar
    bool merged = true;
    std::vector<CNV> newcnv;
    while(merged) {
      int32_t k = -1;
      for(int32_t i = 0; i < (int32_t) chrcnv.size(); ++i) {
	if (i <= k) continue;
	k = i;
	for(int32_t j = i + 1; j < (int32_t) chrcnv.size(); ++j) {
	  bool allValid = true;
	  for(int32_t pre = i; pre < j; ++pre) {
	    double diff = std::abs(chrcnv[pre].cn - chrcnv[j].cn);
	    if (diff >= c.cn_offset) {
	      allValid = false;
	      break;
	    }
	  }
	  if (allValid) k = j;
	  else break;
	}
	if (k > i) {
	  // Merge
	  double cn = (chrcnv[i].cn + chrcnv[k].cn) / 2.0;
	  double zsc = (chrcnv[i].zscore + chrcnv[k].zscore) / 2.0;
	  double mp = (chrcnv[i].mappable + chrcnv[k].mappable) / 2.0;
	  newcnv.push_back(CNV(chrcnv[i].chr, chrcnv[i].start, chrcnv[k].end, chrcnv[i].ciposlow, chrcnv[i].ciposhigh, chrcnv[k].ciendlow, chrcnv[k].ciendhigh, cn, zsc, mp));	  
	} else {
	  newcnv.push_back(chrcnv[i]);
	}
      }
      if (newcnv.size() == chrcnv.size()) merged = false;
      else {
	chrcnv = newcnv;
	newcnv.clear();
      }
    }

    // Insert into global CNV vector
    for(uint32_t i = 0; i < chrcnv.size(); ++i) {
      cnvs.push_back(chrcnv[i]);
      //std::cerr << chrcnv[i].chr << '\t' << chrcnv[i].start << '\t' << chrcnv[i].end << '\t' << '(' << chrcnv[i].ciposlow << ',' << chrcnv[i].ciposhigh << ')' << '\t' << '(' << chrcnv[i].ciendlow << ',' << chrcnv[i].ciendhigh << ')' << '\t' << chrcnv[i].cn << '\t' << chrcnv[i].zscore << '\t' << chrcnv[i].mappable << std::endl;
    }
  }


  template<typename TConfig>
  inline void
  cnvVCF(TConfig const& c, std::vector<CNV> const& cnvs) {
    // BoLog class
    BoLog<double> bl;

    // Open one bam file header
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* bamhd = sam_hdr_read(samfile);

    // Output all copy-number variants
    htsFile *fp = hts_open(c.cnvfile.string().c_str(), "wb");
    bcf_hdr_t *hdr = bcf_hdr_init("w");

    // Print vcf header
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    boost::gregorian::date today = now.date();
    std::string datestr("##fileDate=");
    datestr += boost::gregorian::to_iso_string(today);
    bcf_hdr_append(hdr, datestr.c_str());
    bcf_hdr_append(hdr, "##ALT=<ID=MCNV,Description=\"multi-allelic copy-number variants\">");
    bcf_hdr_append(hdr, "##FILTER=<ID=LowQual,Description=\"Poor quality and insufficient number of PEs and SRs.\">");
    bcf_hdr_append(hdr, "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END\">");
    bcf_hdr_append(hdr, "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS\">");
    bcf_hdr_append(hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the copy-number variant\">");
    bcf_hdr_append(hdr, "##INFO=<ID=MP,Number=1,Type=Float,Description=\"Mappable fraction of CNV\">");
    bcf_hdr_append(hdr, "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect CNV\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Log10-scaled genotype likelihoods for RR,RA,AA genotypes\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Per-sample genotype filter\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Integer copy-number\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=RDCN,Number=1,Type=Float,Description=\"Read-depth based copy-number estimate\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=BP,Number=1,Type=Float,Description=\"Z-score breakpoint RD shift\">");

    // Add reference
    std::string refloc("##reference=");
    refloc += c.genome.string();
    bcf_hdr_append(hdr, refloc.c_str());
    for (int i = 0; i<bamhd->n_targets; ++i) {
      std::string refname("##contig=<ID=");
      refname += std::string(bamhd->target_name[i]) + ",length=" + boost::lexical_cast<std::string>(bamhd->target_len[i]) + ">";
      bcf_hdr_append(hdr, refname.c_str());
    }
    // Add samples
    bcf_hdr_add_sample(hdr, c.sampleName.c_str());
    bcf_hdr_add_sample(hdr, NULL);
    if (bcf_hdr_write(fp, hdr) != 0) std::cerr << "Error: Failed to write BCF header!" << std::endl;

    if (!cnvs.empty()) {
      // Genotype arrays
      int32_t *gts = (int*) malloc(bcf_hdr_nsamples(hdr) * 2 * sizeof(int));
      float *gls = (float*) malloc(bcf_hdr_nsamples(hdr) * 3 * sizeof(float));
      int32_t *gqval = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
      int32_t *cnval = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
      float *cnrdval = (float*) malloc(bcf_hdr_nsamples(hdr) * sizeof(float));
      float *zscoreval = (float*) malloc(bcf_hdr_nsamples(hdr) * sizeof(float));

      std::vector<std::string> ftarr;
      ftarr.resize(bcf_hdr_nsamples(hdr));
    
      // Iterate all structural variants
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Genotyping" << std::endl;
      boost::progress_display show_progress( cnvs.size() );
      bcf1_t *rec = bcf_init();
      for(uint32_t i = 0; i < cnvs.size(); ++i) {
	++show_progress;

	// Integer copy-number
	int32_t absCN = (int32_t) boost::math::round(cnvs[i].cn);
	if (absCN == c.ploidy) continue;
      
	// Output main vcf fields
	rec->rid = bcf_hdr_name2id(hdr, bamhd->target_name[cnvs[i].chr]);
	int32_t svStartPos = cnvs[i].start + 1;
	int32_t svEndPos = cnvs[i].end;
	if (svEndPos >= (int32_t) bamhd->target_len[cnvs[i].chr]) svEndPos = bamhd->target_len[cnvs[i].chr] - 1;
	rec->pos = svStartPos;
	std::string id("MCNV");
	std::string padNumber = boost::lexical_cast<std::string>(i+1);
	padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
	id += padNumber;
	bcf_update_id(hdr, rec, id.c_str());
	std::string svtype = "MCNV";
	std::string alleles = "N,<" + svtype + ">";
	bcf_update_alleles_str(hdr, rec, alleles.c_str());
	int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
	bcf_update_filter(hdr, rec, &tmpi, 1);
      
	// Add INFO fields
	bcf_update_info_flag(hdr, rec, "IMPRECISE", NULL, 1);

	bcf_update_info_string(hdr, rec, "SVTYPE", svtype.c_str());
	std::string dellyVersion("EMBL.DELLYv");
	dellyVersion += dellyVersionNumber;
	bcf_update_info_string(hdr,rec, "SVMETHOD", dellyVersion.c_str());
	tmpi = svEndPos;
	bcf_update_info_int32(hdr, rec, "END", &tmpi, 1);
	int32_t ciend[2];
	ciend[0] = cnvs[i].ciendlow - svEndPos;
	ciend[1] = cnvs[i].ciendhigh - svEndPos;
	int32_t cipos[2];
	cipos[0] = cnvs[i].ciposlow - svStartPos;
	cipos[1] = cnvs[i].ciposhigh - svStartPos;
	bcf_update_info_int32(hdr, rec, "CIPOS", cipos, 2);
	bcf_update_info_int32(hdr, rec, "CIEND", ciend, 2);
	float tmpf = cnvs[i].mappable;
	bcf_update_info_float(hdr, rec, "MP", &tmpf, 1);

	// Genotyping
	gts[0] = bcf_gt_missing;
	gts[1] = bcf_gt_missing;
	gls[0] = 0;
	gls[1] = 0;
	gls[2] = 0;
	gqval[0] = 0;
	if (gqval[0] < 15) ftarr[0] = "LowQual";
	else ftarr[0] = "PASS";
	cnval[0] = absCN;
	cnrdval[0] = cnvs[i].cn;
	zscoreval[0] = cnvs[i].zscore;
	

	
	bcf_update_genotypes(hdr, rec, gts, bcf_hdr_nsamples(hdr) * 2);
	bcf_update_format_float(hdr, rec, "GL",  gls, bcf_hdr_nsamples(hdr) * 3);
	bcf_update_format_int32(hdr, rec, "GQ", gqval, bcf_hdr_nsamples(hdr));
	std::vector<const char*> strp(bcf_hdr_nsamples(hdr));
	std::transform(ftarr.begin(), ftarr.end(), strp.begin(), cstyle_str());
	bcf_update_format_string(hdr, rec, "FT", &strp[0], bcf_hdr_nsamples(hdr));
	bcf_update_format_int32(hdr, rec, "CN", cnval, bcf_hdr_nsamples(hdr));
	bcf_update_format_float(hdr, rec, "RDCN",  cnrdval, bcf_hdr_nsamples(hdr));
	bcf_update_format_float(hdr, rec, "BP",  zscoreval, bcf_hdr_nsamples(hdr));
	bcf_write1(fp, hdr, rec);
	bcf_clear1(rec);
      }
      bcf_destroy1(rec);
    
      // Clean-up
      free(gts);
      free(gls);
      free(gqval);
      free(cnval);
      free(cnrdval);
      free(zscoreval);
    }

    // Close BAM file
    bam_hdr_destroy(bamhd);
    sam_close(samfile);
    
    // Close VCF file
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    
    // Build index
    bcf_index_build(c.outfile.string().c_str(), 14);
  }
 

}

#endif
