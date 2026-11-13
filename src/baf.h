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

  // Allele-count model
  struct BafModel {
    uint32_t maxDepth;
    double minFrac;
    double outlier;
    double floorBaf;
    double err;
    double hetPrior;
    double rhoHet;
    double rhoHom;
    std::vector<double> grid;
    std::vector<double> lgFac;
    std::vector< std::vector<double> > lgHetA;
    std::vector< std::vector<double> > lgHetB;
    std::vector<double> lgHetS;
    std::vector<double> lgHetC;
    std::vector<double> lgHomA;
    std::vector<double> lgHomB;
    std::vector<double> lgHomS;
    double lgHomC;

    BafModel(double const mf) : maxDepth(500), minFrac(mf), outlier(0.05), floorBaf(0.02), err(0.001), hetPrior(0.1), rhoHet(0.02), rhoHom(0.02), lgHomC(0) {
      for(double b = mf; b <= 0.5 + 1e-9; b += 0.005) grid.push_back(b);
      lgFac.resize(maxDepth + 1);
      for(uint32_t i = 0; i <= maxDepth; ++i) lgFac[i] = std::lgamma((double) i + 1.0);
      setHet(rhoHet);
      setHom(rhoHom);
    }

    // Beta-binomial tables
    void _tables(double const p, double const rho, std::vector<double>& lgA, std::vector<double>& lgB, std::vector<double>& lgS, double& lgC) const {
      double const s = (1.0 - rho) / rho;
      double const al = p * s;
      double const be = (1.0 - p) * s;
      lgA.resize(maxDepth + 1);
      lgB.resize(maxDepth + 1);
      lgS.resize(maxDepth + 1);
      for(uint32_t i = 0; i <= maxDepth; ++i) {
	lgA[i] = std::lgamma((double) i + al);
	lgB[i] = std::lgamma((double) i + be);
	lgS[i] = std::lgamma((double) i + s);
      }
      lgC = std::lgamma(al) + std::lgamma(be) - std::lgamma(s);
    }

    void setHet(double const rho) {
      rhoHet = rho;
      lgHetA.resize(grid.size());
      lgHetB.resize(grid.size());
      lgHetC.resize(grid.size());
      for(uint32_t k = 0; k < grid.size(); ++k) _tables(grid[k], rho, lgHetA[k], lgHetB[k], lgHetS, lgHetC[k]);
    }

    void setHom(double const rho) {
      rhoHom = rho;
      _tables(err, rho, lgHomA, lgHomB, lgHomS, lgHomC);
    }

    // Scale counts to the table size
    inline void counts(uint32_t const ref, uint32_t const alt, uint32_t& a, uint32_t& d) const {
      d = ref + alt;
      a = alt;
      if (d > maxDepth) {
	a = (uint32_t) ((double) alt * (double) maxDepth / (double) d + 0.5);
	d = maxDepth;
      }
    }

    inline double lchoose(uint32_t const a, uint32_t const d) const { return lgFac[d] - lgFac[a] - lgFac[d - a]; }
    inline double lhet(uint32_t const a, uint32_t const d, uint32_t const k) const { return lgHetA[k][a] + lgHetB[k][d - a] - lgHetS[d] - lgHetC[k]; }
    inline double lhom(uint32_t const a, uint32_t const d) const { return lgHomA[a] + lgHomB[d - a] - lgHomS[d] - lgHomC; }
  };

  inline double
  _logSumExp(double const x, double const y) {
    double const m = std::max(x, y);
    return m + std::log1p(std::exp(std::min(x, y) - m));
  }

  // BAF
  inline BinBaf
  binBaf(std::vector<AlleleCount> const& cov, int32_t const start, int32_t const end, uint32_t const minHet, BafModel const& bm, double& loglik) {
    BinBaf bb;
    loglik = 0;
    std::vector<AlleleCount>::const_iterator itS = std::lower_bound(cov.begin(), cov.end(), start, [](AlleleCount const& s, int32_t const p) { return s.pos < p; });
    std::vector<AlleleCount>::const_iterator itE = std::lower_bound(cov.begin(), cov.end(), end, [](AlleleCount const& s, int32_t const p) { return s.pos < p; });
    bb.nsites = itE - itS;
    if (bb.nsites == 0) return bb;

    // Site likelihoods
    std::vector<uint32_t> a(bb.nsites);
    std::vector<uint32_t> d(bb.nsites);
    std::vector<double> lhom(bb.nsites);
    std::vector<double> lout(bb.nsites);
    std::vector<double> lhet(bb.nsites);
    std::vector<double> post(bb.nsites);
    double const log2 = std::log(2.0);
    uint32_t i = 0;
    for(std::vector<AlleleCount>::const_iterator it = itS; it != itE; ++it, ++i) {
      bm.counts(it->ref, it->alt, a[i], d[i]);
      lhom[i] = _logSumExp(bm.lhom(a[i], d[i]), bm.lhom(d[i] - a[i], d[i])) - log2;
      lout[i] = -std::log((double) d[i] + 1.0) - bm.lchoose(a[i], d[i]);
    }

    // Grid search
    double const lOut = std::log(bm.outlier);
    double bestLL = -std::numeric_limits<double>::max();
    double bestHet = 0;
    for(uint32_t k = 0; k < bm.grid.size(); ++k) {
      for(i = 0; i < bb.nsites; ++i) lhet[i] = _logSumExp(bm.lhet(a[i], d[i], k), bm.lhet(d[i] - a[i], d[i], k)) - log2;
      double pi = bm.hetPrior;
      double ll = 0;
      for(uint32_t iter = 0; iter < 3; ++iter) {
	double const lHet = std::log(pi);
	double const lHom = std::log(std::max(1.0 - pi - bm.outlier, 1e-3));
	ll = 0;
	double sum = 0;
	for(i = 0; i < bb.nsites; ++i) {
	  double const lp = _logSumExp(_logSumExp(lHet + lhet[i], lHom + lhom[i]), lOut + lout[i]);
	  post[i] = std::exp(lHet + lhet[i] - lp);
	  sum += post[i];
	  ll += lp;
	}
	pi = std::min(std::max(sum / (double) bb.nsites, 1e-3), 1.0 - bm.outlier - 1e-3);
      }
      if (ll > bestLL) {
	bestLL = ll;
	bb.baf = bm.grid[k];
	bestHet = pi * (double) bb.nsites;
      }
    }
    loglik = bestLL;
    bb.nhet = (uint32_t) (bestHet + 0.5);
    if ((bb.nhet < minHet) || (bb.baf < bm.minFrac + bm.floorBaf)) bb.baf = -1;   // too few het. sites or unreliable BAF
    return bb;
  }

  // Calibrate
  struct BafRegion {
    int32_t refIndex;
    int32_t start;
    int32_t end;

    BafRegion(int32_t const ri, int32_t const s, int32_t const e) : refIndex(ri), start(s), end(e) {}
  };

  inline void
  calibrateBafModel(std::vector< std::vector<AlleleCount> > const& snv, std::vector<BafRegion> const& regions, uint32_t const minHet, BafModel& bm) {
    // Error rate and het. fraction from clean hom and het sites
    double minor = 0;
    double depth = 0;
    uint64_t nhet = 0;
    uint64_t nsites = 0;
    for(uint32_t ri = 0; ri < snv.size(); ++ri) {
      for(uint32_t i = 0; i < snv[ri].size(); ++i) {
	double const d = (double) snv[ri][i].ref + (double) snv[ri][i].alt;
	double const f = (double) snv[ri][i].alt / d;
	++nsites;
	if ((f < bm.minFrac) || (f > 1.0 - bm.minFrac)) {
	  minor += std::min(snv[ri][i].ref, snv[ri][i].alt);
	  depth += d;
	} else ++nhet;
      }
    }
    if (depth > 0) bm.err = std::max(0.001, minor / depth);
    if (nsites > 0) bm.hetPrior = std::min(std::max((double) nhet / (double) nsites, 0.01), 0.9);

    // Over-dispersion
    double const rhos[] = {0.002, 0.005, 0.01, 0.02, 0.04, 0.08, 0.15};
    uint32_t const nrho = sizeof(rhos) / sizeof(rhos[0]);
    for(uint32_t comp = 0; comp < 2; ++comp) {
      double best = -std::numeric_limits<double>::max();
      double bestRho = (comp == 0) ? bm.rhoHom : bm.rhoHet;
      for(uint32_t r = 0; r < nrho; ++r) {
	if (comp == 0) bm.setHom(rhos[r]);
	else bm.setHet(rhos[r]);
	double sum = 0;
	for(uint32_t i = 0; i < regions.size(); ++i) {
	  double ll = 0;
	  binBaf(snv[regions[i].refIndex], regions[i].start, regions[i].end, minHet, bm, ll);
	  sum += ll;
	}
	if (sum > best) {
	  best = sum;
	  bestRho = rhos[r];
	}
      }
      if (comp == 0) bm.setHom(bestRho);
      else bm.setHet(bestRho);
    }
  }

}

#endif
