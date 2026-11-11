#ifndef BAF_H
#define BAF_H

#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

#include <htslib/sam.h>
#include <htslib/vcf.h>

namespace torali
{

  // Common SNV sites
  struct SnvSite {
    int32_t pos;
    uint16_t ref;
    uint16_t alt;
    uint8_t refBase;
    uint8_t altBase;

    SnvSite(int32_t const p, uint8_t const r, uint8_t const a) : pos(p), ref(0), alt(0), refBase(r), altBase(a) {}
  };

  // Allele counts
  struct AlleleCount {
    int32_t pos;
    uint16_t ref;
    uint16_t alt;

    AlleleCount(int32_t const p, uint16_t const r, uint16_t const a) : pos(p), ref(r), alt(a) {}
  };

  // BAF
  struct BinBaf {
    uint32_t nsites;
    uint32_t nhet;
    double baf;

    BinBaf() : nsites(0), nhet(0), baf(-1) {}
  };

  // Common SNV sites
  struct SnvSitesFile {
    htsFile* fp;
    bcf_hdr_t* hdr;
    hts_idx_t* idx;

    SnvSitesFile() : fp(NULL), hdr(NULL), idx(NULL) {}

    bool open(std::string const& path) {
      fp = bcf_open(path.c_str(), "r");
      if (fp == NULL) return false;
      if (hts_get_format(fp)->format != bcf) return false;
      hdr = bcf_hdr_read(fp);
      if (hdr == NULL) return false;
      idx = bcf_index_load(path.c_str());
      return (idx != NULL);
    }

    void close() {
      if (idx != NULL) hts_idx_destroy(idx);
      if (hdr != NULL) bcf_hdr_destroy(hdr);
      if (fp != NULL) bcf_close(fp);
      idx = NULL;
      hdr = NULL;
      fp = NULL;
    }
  };

  inline bool
  _acgt(uint8_t const b) {
    return ((b == 1) || (b == 2) || (b == 4) || (b == 8));
  }

  // Biallelic SNVs
  inline void
  loadSnvSites(SnvSitesFile const& sf, std::string const& chrName, std::vector<SnvSite>& sites) {
    sites.clear();
    int32_t rid = bcf_hdr_name2id(sf.hdr, chrName.c_str());
    if (rid < 0) {
      // Toggle chr prefix
      std::string altName = (chrName.substr(0, 3) == "chr") ? chrName.substr(3) : ("chr" + chrName);
      rid = bcf_hdr_name2id(sf.hdr, altName.c_str());
    }
    if (rid < 0) return;
    hts_itr_t* itr = bcf_itr_queryi(sf.idx, rid, 0, HTS_POS_MAX);
    if (itr == NULL) return;
    bcf1_t* rec = bcf_init();
    int64_t lastPos = -1;
    while (bcf_itr_next(sf.fp, itr, rec) >= 0) {
      if (rec->pos == lastPos) {
	if ((!sites.empty()) && (sites.back().pos == rec->pos)) sites.pop_back();
	continue;
      }
      lastPos = rec->pos;
      if ((rec->n_allele != 2) || (rec->rlen != 1)) continue;
      bcf_unpack(rec, BCF_UN_STR);
      if ((rec->d.allele[0][1] != '\0') || (rec->d.allele[1][1] != '\0')) continue;
      uint8_t r = seq_nt16_table[(uint8_t) rec->d.allele[0][0]];
      uint8_t a = seq_nt16_table[(uint8_t) rec->d.allele[1][0]];
      if ((!_acgt(r)) || (!_acgt(a)) || (r == a)) continue;
      sites.push_back(SnvSite(rec->pos, r, a));
    }
    bcf_destroy(rec);
    hts_itr_destroy(itr);
  }

  // Count REF and ALT
  inline void
  addAlleleCounts(bam1_t const* rec, std::vector<SnvSite>& sites, uint8_t const minBaseQual) {
    if ((sites.empty()) || (rec->core.l_qseq == 0)) return;
    uint16_t const maxCount = std::numeric_limits<uint16_t>::max();
    uint32_t const* cigar = bam_get_cigar(rec);
    uint8_t const* seq = bam_get_seq(rec);
    uint8_t const* qual = bam_get_qual(rec);
    bool const hasQual = (qual[0] != 0xff);
    int32_t rp = rec->core.pos;
    int32_t qp = 0;
    std::vector<SnvSite>::iterator it = std::lower_bound(sites.begin(), sites.end(), rp, [](SnvSite const& s, int32_t const p) { return s.pos < p; });
    for(uint32_t i = 0; ((i < rec->core.n_cigar) && (it != sites.end())); ++i) {
      uint32_t const op = bam_cigar_op(cigar[i]);
      int32_t const ol = (int32_t) bam_cigar_oplen(cigar[i]);
      if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF)) {
	while ((it != sites.end()) && (it->pos < rp + ol)) {
	  int32_t const q = qp + (it->pos - rp);
	  if ((!hasQual) || (qual[q] >= minBaseQual)) {
	    uint8_t const b = bam_seqi(seq, q);
	    if (b == it->refBase) {
	      if (it->ref < maxCount) ++it->ref;
	    }
	    else if (b == it->altBase) {
	      if (it->alt < maxCount) ++it->alt;
	    }
	  }
	  ++it;
	}
	rp += ol;
	qp += ol;
      } else if ((op == BAM_CINS) || (op == BAM_CSOFT_CLIP)) qp += ol;
      else if ((op == BAM_CDEL) || (op == BAM_CREF_SKIP)) {
	rp += ol;
	while ((it != sites.end()) && (it->pos < rp)) ++it;
      }
    }
  }

  // Keep the covered sites
  inline void
  coveredSites(std::vector<SnvSite> const& sites, uint32_t const minDepth, std::vector<AlleleCount>& cov) {
    cov.clear();
    for(uint32_t i = 0; i < sites.size(); ++i) {
      if ((uint32_t) sites[i].ref + (uint32_t) sites[i].alt >= minDepth) cov.push_back(AlleleCount(sites[i].pos, sites[i].ref, sites[i].alt));
    }
  }

  // Het. sites
  inline bool
  _hetSite(uint32_t const ref, uint32_t const alt, double const minFrac) {
    double const d = (double) ref + (double) alt;
    double const f = (double) alt / d;
    return ((f >= minFrac) && (f <= 1.0 - minFrac) && (ref >= 2) && (alt >= 2));
  }

  // Het. site filter
  struct HetBand {
    uint32_t maxDepth;
    double minFrac;
    std::vector<double> grid;
    std::vector<std::vector<double> > lp;

    HetBand(double const mf) : maxDepth(500), minFrac(mf) {
      for(double b = mf; b <= 0.5 + 1e-9; b += 0.005) grid.push_back(b);
      lp.resize(maxDepth + 1, std::vector<double>(grid.size(), 0));
      std::vector<double> lchoose(maxDepth + 1);
      for(uint32_t d = 1; d <= maxDepth; ++d) {
	lchoose[0] = 0;
	for(uint32_t a = 1; a <= d; ++a) lchoose[a] = lchoose[a - 1] + std::log((double) (d - a + 1)) - std::log((double) a);
	for(uint32_t k = 0; k < grid.size(); ++k) {
	  double const lb = std::log(grid[k]);
	  double const l1b = std::log(1.0 - grid[k]);
	  double sum = 0;
	  for(uint32_t a = 0; a <= d; ++a) {
	    if (!_hetSite(d - a, a, minFrac)) continue;
	    sum += std::exp(lchoose[a] + a * lb + (d - a) * l1b);
	  }
	  lp[d][k] = (sum > 0) ? std::log(sum) : -1000.0;
	}
      }
    }

    double logProb(uint32_t const depth, uint32_t const k) const {
      return lp[std::min(depth, maxDepth)][k];
    }
  };

  // BAF
  inline BinBaf
  binBaf(std::vector<AlleleCount> const& cov, int32_t const start, int32_t const end, uint32_t const minHet, HetBand const& hb) {
    BinBaf bb;
    std::vector<AlleleCount>::const_iterator itS = std::lower_bound(cov.begin(), cov.end(), start, [](AlleleCount const& s, int32_t const p) { return s.pos < p; });
    std::vector<AlleleCount>::const_iterator itE = std::lower_bound(cov.begin(), cov.end(), end, [](AlleleCount const& s, int32_t const p) { return s.pos < p; });
    std::vector<AlleleCount const*> het;
    for(std::vector<AlleleCount>::const_iterator it = itS; it != itE; ++it) {
      ++bb.nsites;
      if (_hetSite(it->ref, it->alt, hb.minFrac)) het.push_back(&(*it));
    }
    bb.nhet = het.size();
    if (bb.nhet < minHet) return bb;

    // Grid
    double bestLL = -std::numeric_limits<double>::max();
    for(uint32_t k = 0; k < hb.grid.size(); ++k) {
      double const lb = std::log(hb.grid[k]);
      double const l1b = std::log(1.0 - hb.grid[k]);
      double ll = 0;
      for(uint32_t i = 0; i < het.size(); ++i) {
	double const a = het[i]->alt;
	double const r = het[i]->ref;
	double const x = a * lb + r * l1b;
	double const y = a * l1b + r * lb;
	ll += std::max(x, y) + std::log1p(std::exp(std::min(x, y) - std::max(x, y))) - hb.logProb(het[i]->ref + het[i]->alt, k);
      }
      if (ll > bestLL) {
	bestLL = ll;
	bb.baf = hb.grid[k];
      }
    }
    return bb;
  }

}

#endif
