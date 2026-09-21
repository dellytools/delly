#ifndef PLOIDY_H
#define PLOIDY_H

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>

namespace torali
{

  // 0 = unknown (diploid GTs), 1 = male, 2 = female
  struct SexModel {
    typedef std::vector<std::pair<int32_t, int32_t> > TRegions;

    int32_t xTid;
    int32_t yTid;
    TRegions parX;
    TRegions parY;
    std::vector<uint8_t> sex;

    SexModel() : xTid(-1), yTid(-1) {}
  };

  inline bool
  _inRegions(SexModel::TRegions const& reg, int32_t const pos) {
    for(SexModel::TRegions::const_iterator it = reg.begin(); it != reg.end(); ++it) {
      if ((pos >= it->first) && (pos < it->second)) return true;
    }
    return false;
  }

  // Ploidy of a sample at a genomic position
  inline uint8_t
  _ploidy(SexModel const& sm, uint8_t const sex, int32_t const tid, int32_t const pos) {
    if (sex == 0) return 2;
    if (tid == sm.xTid) {
      if ((sex == 1) && (!_inRegions(sm.parX, pos))) return 1;
      return 2;
    }
    if (tid == sm.yTid) {
      if (sex == 2) return 0;
      if (_inRegions(sm.parY, pos)) return 2;
      return 1;
    }
    return 2;
  }

  // PAR regions
  inline bool
  _parRegions(uint32_t const xlen, SexModel& sm) {
    sm.parX.clear();
    sm.parY.clear();
    if (xlen == 156040895) {
      // GRCh38
      sm.parX.push_back(std::make_pair(10000, 2781479));
      sm.parX.push_back(std::make_pair(155701382, 156030895));
      sm.parY.push_back(std::make_pair(10000, 2781479));
      sm.parY.push_back(std::make_pair(56887902, 57217415));
    } else if (xlen == 155270560) {
      // GRCh37
      sm.parX.push_back(std::make_pair(60000, 2699520));
      sm.parX.push_back(std::make_pair(154931043, 155260560));
      sm.parY.push_back(std::make_pair(10000, 2649520));
      sm.parY.push_back(std::make_pair(59034049, 59363566));
    } else if (xlen == 154259566) {
      // T2T-CHM13
      sm.parX.push_back(std::make_pair(0, 2394410));
      sm.parX.push_back(std::make_pair(153925834, 154259566));
      sm.parY.push_back(std::make_pair(0, 2458320));
      sm.parY.push_back(std::make_pair(62122809, 62460029));
    } else return false;
    return true;
  }

  // Any sex chromosomes (BAM)
  inline bool
  _sexChromosomes(bam_hdr_t const* hdr, SexModel& sm) {
    sm.xTid = -1;
    sm.yTid = -1;
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      std::string tname(hdr->target_name[refIndex]);
      if ((tname == "chrX") || (tname == "X")) sm.xTid = refIndex;
      else if ((tname == "chrY") || (tname == "Y")) sm.yTid = refIndex;
    }
    if (sm.xTid == -1) return false;
    if (!_parRegions(hdr->target_len[sm.xTid], sm)) {
      // Unknown genome build
      sm.xTid = -1;
      sm.yTid = -1;
      return false;
    }
    return true;
  }

  // Any sex chromosomes (BCF)
  inline bool
  _sexChromosomes(bcf_hdr_t const* hdr, SexModel& sm) {
    sm.xTid = -1;
    sm.yTid = -1;
    int32_t nseq = 0;
    const char** seqnames = bcf_hdr_seqnames(hdr, &nseq);
    uint32_t xlen = 0;
    for(int32_t i = 0; i < nseq; ++i) {
      std::string tname(seqnames[i]);
      int32_t rid = bcf_hdr_name2id(hdr, seqnames[i]);
      if ((tname == "chrX") || (tname == "X")) {
	sm.xTid = rid;
	bcf_hrec_t* hrec = bcf_hdr_get_hrec(hdr, BCF_HL_CTG, "ID", seqnames[i], NULL);
	if (hrec != NULL) {
	  int32_t k = bcf_hrec_find_key(hrec, "length");
	  if (k >= 0) xlen = (uint32_t) std::atoll(hrec->vals[k]);
	}
      } else if ((tname == "chrY") || (tname == "Y")) sm.yTid = rid;
    }
    free(seqnames);
    if (sm.xTid == -1) return false;
    if (!_parRegions(xlen, sm)) {
      sm.xTid = -1;
      sm.yTid = -1;
      return false;
    }
    return true;
  }

  // Infer sex
  inline uint8_t
  _inferSexFromReads(samFile* samfile, hts_idx_t* idx, bam_hdr_t const* hdr, std::string const& genome, SexModel const& sm, double& xRatio, double& yRatio) {
    xRatio = -1;
    yRatio = -1;
    if (sm.xTid == -1) return 0;
    hts_set_fai_filename(samfile, genome.c_str());
    hts_set_opt(samfile, CRAM_OPT_REQUIRED_FIELDS, SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ);
    int32_t const nwin = 20;
    int32_t const wlen = 200000;
    typedef std::vector<std::pair<int32_t, std::pair<int32_t, int32_t> > > TIntervals;
    TIntervals autoIv;
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      std::string tname(hdr->target_name[refIndex]);
      if (tname.substr(0, 3) == "chr") tname = tname.substr(3);
      bool isNumber = (!tname.empty());
      for(std::size_t i = 0; i < tname.size(); ++i) if (!std::isdigit((unsigned char) tname[i])) isNumber = false;
      if ((isNumber) && (hdr->target_len[refIndex] > (uint32_t) (2 * nwin * wlen))) autoIv.push_back(std::make_pair(refIndex, std::make_pair(wlen, (int32_t) hdr->target_len[refIndex] - wlen)));
    }
    if (autoIv.empty()) return 0;
    int32_t xBeg = wlen;
    int32_t xEnd = hdr->target_len[sm.xTid] - wlen;
    if (!sm.parX.empty()) {
      xBeg = sm.parX.front().second + wlen;
      xEnd = sm.parX.back().first - wlen;
    }
    int32_t yBeg = -1;
    int32_t yEnd = -1;
    if (sm.yTid != -1) {
      // Exclude heterochromatic Yq
      yBeg = (sm.parY.empty()) ? wlen : (sm.parY.front().second + wlen);
      yEnd = (int32_t) (0.45 * hdr->target_len[sm.yTid]);
    }
    // Read counts per window
    std::vector<double> cnt[3];   // autosomes, X, Y
    for(int32_t cls = 0; cls < 3; ++cls) {
      TIntervals iv;
      if (cls == 0) iv = autoIv;
      else if (cls == 1) iv.push_back(std::make_pair(sm.xTid, std::make_pair(xBeg, xEnd)));
      else {
	if (sm.yTid == -1) continue;
	iv.push_back(std::make_pair(sm.yTid, std::make_pair(yBeg, yEnd)));
      }
      for(TIntervals::const_iterator it = iv.begin(); it != iv.end(); ++it) {
	int32_t span = it->second.second - it->second.first;
	if (span < nwin * wlen) continue;
	int32_t step = span / nwin;
	for(int32_t w = 0; w < nwin; ++w) {
	  int32_t beg = it->second.first + w * step;
	  hts_itr_t* iter = sam_itr_queryi(idx, it->first, beg, beg + wlen);
	  if (iter == NULL) continue;
	  bam1_t* rec = bam_init1();
	  uint32_t c = 0;
	  while (sam_itr_next(samfile, iter, rec) >= 0) {
	    if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	    if (rec->core.qual < 20) continue;
	    ++c;
	  }
	  bam_destroy1(rec);
	  hts_itr_destroy(iter);
	  cnt[cls].push_back(c);
	}
      }
    }
    if ((cnt[0].size() < 10) || (cnt[1].size() < 5)) return 0;
    double med[3] = {0, 0, 0};
    for(int32_t cls = 0; cls < 3; ++cls) {
      if (cnt[cls].empty()) continue;
      std::sort(cnt[cls].begin(), cnt[cls].end());
      med[cls] = cnt[cls][cnt[cls].size() / 2];
    }
    if (med[0] <= 0) return 0;
    xRatio = med[1] / med[0];
    yRatio = (cnt[2].empty()) ? -1 : (med[2] / med[0]);
    if (xRatio < 0.75) {
      if ((yRatio < 0) || (yRatio > 0.1)) return 1;
    } else {
      if ((yRatio < 0) || (yRatio < 0.15)) return 2;
    }
    return 0;
  }

  // Infer sample sex
  inline void
  _inferSexFromCallset(std::string const& vcffile, bcf_hdr_t const* hdr, SexModel& sm) {
    int32_t nsmpl = bcf_hdr_nsamples(hdr);
    sm.sex.assign(nsmpl, 0);
    if (sm.xTid == -1) return;
    bcf_srs_t* sr = bcf_sr_init();
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
    std::string regions = std::string(bcf_hdr_id2name(hdr, sm.xTid));
    if (sm.yTid != -1) regions += "," + std::string(bcf_hdr_id2name(hdr, sm.yTid));
    if (bcf_sr_set_regions(sr, regions.c_str(), 0) != 0) {
      bcf_sr_destroy(sr);
      return;
    }
    if (!bcf_sr_add_reader(sr, vcffile.c_str())) {
      bcf_sr_destroy(sr);
      return;
    }
    std::vector<uint32_t> xCalled(nsmpl, 0);
    std::vector<uint32_t> xHet(nsmpl, 0);
    std::vector<uint32_t> ySites(nsmpl, 0);
    std::vector<uint32_t> yCalled(nsmpl, 0);
    int32_t ngt = 0;
    int32_t* gt = NULL;
    while (bcf_sr_next_line(sr)) {
      bcf1_t* rec = bcf_sr_get_line(sr, 0);
      if (rec == NULL) continue;
      bool onX = (rec->rid == sm.xTid);
      if ((onX) && (_inRegions(sm.parX, rec->pos))) continue;
      if ((!onX) && (rec->rid != sm.yTid)) continue;
      if (bcf_get_genotypes(sr->readers[0].header, rec, &gt, &ngt) < 1) continue;
      int32_t stride = ngt / nsmpl;
      for(int32_t i = 0; i < nsmpl; ++i) {
	int32_t* g = gt + i * stride;
	if (g[0] == bcf_gt_missing) {
	  if (!onX) ++ySites[i];
	  continue;
	}
	int32_t a0 = bcf_gt_allele(g[0]);
	int32_t a1 = ((stride > 1) && (g[1] != bcf_int32_vector_end) && (g[1] != bcf_gt_missing)) ? bcf_gt_allele(g[1]) : a0;
	if (onX) {
	  ++xCalled[i];
	  if (a0 != a1) ++xHet[i];
	} else {
	  ++ySites[i];
	  ++yCalled[i];
	}
      }
    }
    if (gt != NULL) free(gt);
    bcf_sr_destroy(sr);

    // chrY call rate
    std::vector<double> hetX(nsmpl, -1);
    for(int32_t i = 0; i < nsmpl; ++i) {
      if (xCalled[i] >= 100) hetX[i] = (double) xHet[i] / (double) xCalled[i];
    }
    double hetSplit = -1;
    {
      std::vector<double> vals;
      for(int32_t i = 0; i < nsmpl; ++i) if (hetX[i] >= 0) vals.push_back(hetX[i]);
      if (vals.size() >= 6) {
	std::sort(vals.begin(), vals.end());
	double total = 0;
	for(std::size_t i = 0; i < vals.size(); ++i) total += vals[i];
	double best = 0;
	double cum = 0;
	for(std::size_t k = 1; k < vals.size(); ++k) {
	  cum += vals[k - 1];
	  double m1 = cum / k;
	  double m2 = (total - cum) / (vals.size() - k);
	  double between = (double) k * (double) (vals.size() - k) * (m1 - m2) * (m1 - m2);
	  if ((between > best) && (k >= 3) && (vals.size() - k >= 3) && (m2 > 1.5 * m1)) {
	    best = between;
	    hetSplit = (vals[k - 1] + vals[k]) / 2.0;
	  }
	}
      }
    }
    for(int32_t i = 0; i < nsmpl; ++i) {
      double calledY = (ySites[i] >= 20) ? ((double) yCalled[i] / (double) ySites[i]) : -1;
      bool hetMale = ((hetSplit > 0) && (hetX[i] >= 0) && (hetX[i] < hetSplit));
      bool hetFemale = ((hetSplit > 0) && (hetX[i] >= 0) && (hetX[i] >= hetSplit));
      if (calledY >= 0) {
	if ((calledY >= 0.5) && (!hetFemale)) sm.sex[i] = 1;
	else if ((calledY < 0.5) && (!hetMale)) sm.sex[i] = 2;
      } else {
	if (hetMale) sm.sex[i] = 1;
	else if (hetFemale) sm.sex[i] = 2;
      }
    }
  }

  template<typename TConfig>
  inline bool
  _parseSex(TConfig& c, std::string const& sexArg, SexModel& sm, bool& autoSex) {
    autoSex = false;
    sm.sex.assign(c.sampleName.size(), 0);
    if (sexArg == "auto") autoSex = true;
    else if (sexArg == "none") return true;
    else if (sexArg == "male") sm.sex.assign(c.sampleName.size(), 1);
    else if (sexArg == "female") sm.sex.assign(c.sampleName.size(), 2);
    else {
      // Parse sex file
      std::ifstream sexFile(sexArg.c_str());
      if (!sexFile.is_open()) {
	std::cerr << "Sex file cannot be opened: " << sexArg << std::endl;
	return false;
      }
      std::map<std::string, uint8_t> sexMap;
      std::string line;
      while (std::getline(sexFile, line)) {
	std::vector<std::string> tokens;
	boost::split(tokens, line, boost::is_any_of("\t ,"), boost::token_compress_on);
	if (tokens.size() < 2) continue;
	std::string s = boost::to_lower_copy(tokens[1]);
	if ((s == "male") || (s == "m") || (s == "1")) sexMap[tokens[0]] = 1;
	else if ((s == "female") || (s == "f") || (s == "2")) sexMap[tokens[0]] = 2;
      }
      for(std::size_t i = 0; i < c.sampleName.size(); ++i) {
	if (sexMap.find(c.sampleName[i]) != sexMap.end()) sm.sex[i] = sexMap[c.sampleName[i]];
      }
    }
    return true;
  }

  inline std::string
  _sexName(uint8_t const sex) {
    if (sex == 1) return "male";
    if (sex == 2) return "female";
    return "unknown";
  }

}

#endif
