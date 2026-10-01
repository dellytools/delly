#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include <htslib/sam.h>

#include "util.h"
#include "methyl.h"
#include "ploidy.h"

namespace torali
{

  // Spanning TR read
  struct TrRead {
    int32_t len;
    int32_t ps;
    uint8_t hp;
  };

  // TR noise
  struct TrNoise {
    double mu[4];
    double b[4];
    double eps[4];
    std::vector<int32_t> dev[4];

    TrNoise() {
      for(int32_t i = 0; i < 4; ++i) { mu[i] = 0; b[i] = 8; eps[i] = 0.03; }
    }
  };

  inline int32_t
  _editDistanceNW(std::string const& query, std::string const& target) {
    EdlibAlignResult align = edlibAlign(query.c_str(), query.size(), target.c_str(), target.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
    // Debug: requires EDLIB_TASK_PATH otherwise EDLIB_TASK_DISTANCE
    //printAlignment(query, target, EDLIB_MODE_NW, align);
    int32_t ed = align.editDistance;
    edlibFreeAlignResult(align);
    return ed;
  }
  
  inline int32_t
  _readStart(bam1_t* rec) {
    uint32_t rp = rec->core.pos;
    const uint32_t* cigar = bam_get_cigar(rec);
    if (rec->core.n_cigar) {
      if ((bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP)) {
	if (rp > bam_cigar_oplen(cigar[0])) rp -= bam_cigar_oplen(cigar[0]);
	else rp = 0;
      }
    }
    return rp;
  }

  inline int32_t
  _readEnd(bam1_t* rec) {
    uint32_t rp = rec->core.pos;
    const uint32_t* cigar = bam_get_cigar(rec);
    if (rec->core.n_cigar) {
      for (uint32_t i = 0; i < rec->core.n_cigar; ++i) {
	if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF) || (bam_cigar_op(cigar[i]) == BAM_CDEL) || (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP)) rp += bam_cigar_oplen(cigar[i]);
      }
      if ((bam_cigar_op(cigar[rec->core.n_cigar - 1]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[rec->core.n_cigar - 1]) == BAM_CHARD_CLIP)) {
	rp += bam_cigar_oplen(cigar[rec->core.n_cigar - 1]);
      }
    }
    return rp;
  }

  inline int32_t
  _findSeqBp(bam1_t* rec, uint32_t const pos) {
    uint32_t rp = rec->core.pos; // reference pointer
    uint32_t sp = 0; // sequence pointer

    // Parse the CIGAR
    const uint32_t* cigar = bam_get_cigar(rec);
    if (rec->core.n_cigar) {
      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	  for(uint32_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k, ++rp, ++sp) {
	    if (rp >= pos) return sp;
	  }
	} else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	  rp += bam_cigar_oplen(cigar[i]);
	  if (rp >= pos) return sp;
	} else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	  sp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	  rp += bam_cigar_oplen(cigar[i]);
	  if (rp >= pos) return sp;
	} else if ((bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) {
	  sp += bam_cigar_oplen(cigar[i]);
	} else {
	  std::cerr << "Unknown Cigar options" << std::endl;
	}
      }
      if ((bam_cigar_op(cigar[rec->core.n_cigar - 1]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[rec->core.n_cigar - 1]) == BAM_CHARD_CLIP)) {
	return sp - bam_cigar_oplen(cigar[rec->core.n_cigar - 1]);
      }
    }
    return -1;
  }

  inline int32_t
  _trBin(int32_t const tractLen) {
    if (tractLen < 100) return 0;
    if (tractLen < 300) return 1;
    if (tractLen < 1000) return 2;
    return 3;
  }

  inline bool
  _trAlleleLength(StructuralVariantRecord const& sv, int32_t& d) {
    std::size_t comma = sv.alleles.find(',');
    if ((comma != std::string::npos) && (sv.alleles.find('<') == std::string::npos)) d = (int32_t) (sv.alleles.size() - comma - 1) - (int32_t) comma;
    else if (sv.svt == 4) d = sv.insLen;
    else if (sv.svt == 2) d = -(sv.svEnd - sv.svStart);
    else return false;
    return true;
  }

  // Calibrate the read noise
  template<typename TGroups, typename TGroupReads>
  inline void
  _trCalibrate(std::vector<StructuralVariantRecord> const& svs, TGroups const& alleleGroup, TGroupReads const& groupReads, std::vector<int32_t> const& grpWs, std::vector<int32_t> const& grpWe, int32_t const refIndex, TrNoise& nm) {
    for(uint32_t g = 0; g < alleleGroup.size(); ++g) {
      if ((alleleGroup[g].empty()) || (grpWs[g] < 0)) continue;
      StructuralVariantRecord const& first = svs[alleleGroup[g][0]];
      if (first.chr != refIndex) continue;
      if (groupReads[g].size() < 8) continue;
      uint32_t nRef = 0;
      for(uint32_t i = 0; i < groupReads[g].size(); ++i) if (std::abs(groupReads[g][i].len) <= 40) ++nRef;
      if (nRef < 0.8 * groupReads[g].size()) continue;
      int32_t bin = _trBin(grpWe[g] - grpWs[g]);
      if (nm.dev[bin].size() >= 100000) continue;
      for(uint32_t i = 0; i < groupReads[g].size(); ++i) nm.dev[bin].push_back(groupReads[g][i].len);
    }
    for(int32_t bin = 0; bin < 4; ++bin) {
      if (nm.dev[bin].size() < 200) continue;
      std::vector<int32_t> dev(nm.dev[bin]);
      std::sort(dev.begin(), dev.end());
      double mu = dev[dev.size() / 2];
      double sum = 0;
      uint32_t n = 0;
      uint32_t out = 0;
      for(uint32_t i = 0; i < dev.size(); ++i) {
	double ad = std::abs((double) dev[i] - mu);
	if (ad <= 50) { sum += ad; ++n; }
	else ++out;
      }
      if (n == 0) continue;
      nm.mu[bin] = mu;
      nm.b[bin] = std::max(1.0, sum / (double) n);
      nm.eps[bin] = std::max(0.005, (double) out / (double) dev.size());
    }
  }

  inline double
  _trLogSum(double const a, double const b) {
    if (a > b) return a + std::log10(1.0 + std::pow(10.0, b - a));
    return b + std::log10(1.0 + std::pow(10.0, a - b));
  }

  // GT for tandem repeats
  template<typename TMembers, typename TReads, typename TJctVector>
  inline void
  _trLocusGenotype(std::vector<StructuralVariantRecord> const& svs, TMembers const& members, TReads const& reads, TrNoise const& nm, int32_t const ws, int32_t const we, int32_t const flank, uint8_t const ploidy, TJctVector& jct) {
    if ((ploidy == 0) || (reads.size() < 3)) return;
    int32_t const tractLen = we - ws;
    int32_t bin = _trBin(tractLen);
    double const mu = nm.mu[bin];
    double const b = nm.b[bin];
    double const eps = nm.eps[bin];
    double const W = 2000;
    double const lOther = std::log10(1.0 / W);

    // Alleles (0 = REF)
    std::vector<int32_t> alen(1, 0);
    std::vector<int32_t> amember(1, -1);
    for(uint32_t k = 0; k < members.size(); ++k) {
      int32_t d = 0;
      if (!_trAlleleLength(svs[members[k]], d)) continue;
      StructuralVariantRecord const& sv = svs[members[k]];
      int32_t svE = (sv.svt == 4) ? sv.svStart : sv.svEnd;
      if ((sv.svStart < ws - flank) || (svE > we + flank)) continue;
      alen.push_back(d);
      amember.push_back(members[k]);
    }
    if (alen.size() < 2) return;

    // Candidate alleles
    std::vector<int32_t> cand(1, 0);
    for(uint32_t a = 1; a < alen.size(); ++a) {
      for(uint32_t r = 0; r < reads.size(); ++r) {
	if (std::abs((double) reads[r].len - (double) alen[a] - mu) <= 5 * b) { cand.push_back(a); break; }
      }
    }
    int32_t const other = (int32_t) alen.size();
    cand.push_back(other);
    int32_t nc = (int32_t) cand.size();

    // Read likelihoods
    std::vector<std::vector<double> > lp(reads.size(), std::vector<double>(nc, lOther));
    for(int32_t ci = 0; ci < nc - 1; ++ci) {
      int32_t abin = _trBin(std::max(0, tractLen + alen[cand[ci]]));
      for(uint32_t r = 0; r < reads.size(); ++r) {
	double d = std::abs((double) reads[r].len - (double) alen[cand[ci]] - nm.mu[abin]);
	lp[r][ci] = std::log10((1.0 - nm.eps[abin]) * std::exp(-d / nm.b[abin]) / (2.0 * nm.b[abin]) + nm.eps[abin] / W);
      }
    }

    // Alleles too close
    double const resol = std::max(3.0, 2.83 * b / std::sqrt(std::max(1.0, (double) reads.size() / (double) ploidy)));
    std::vector<int32_t> lenClass(nc, 0);
    {
      std::vector<int32_t> ord;
      for(int32_t k = 1; k < nc - 1; ++k) ord.push_back(k);
      std::sort(ord.begin(), ord.end(), [&](int32_t x, int32_t y) { return (alen[cand[x]] < alen[cand[y]]) || ((alen[cand[x]] == alen[cand[y]]) && (x < y)); });
      int32_t cid = 0;
      int32_t anchor = 0;
      for(uint32_t i = 0; i < ord.size(); ++i) {
	if ((i == 0) || ((double) (alen[cand[ord[i]]] - anchor) > resol)) {
	  ++cid;
	  anchor = alen[cand[ord[i]]];
	}
	lenClass[ord[i]] = cid;
      }
      lenClass[nc - 1] = cid + 1;
    }

    // Genotype likelihoods
    int32_t bestX = 0;
    int32_t bestY = 0;
    double bestL = -std::numeric_limits<double>::max();
    std::vector<std::vector<double> > gl(nc, std::vector<double>(nc, 0));
    auto validGt = [&](int32_t x, int32_t y) {
      if (x == y) return true;
      if (ploidy == 1) return false;
      return (lenClass[x] != lenClass[y]);
    };
    for(int32_t x = 0; x < nc; ++x) {
      for(int32_t y = x; y < nc; ++y) {
	if (!validGt(x, y)) continue;
	double l = 0;
	if (x == y) for(uint32_t r = 0; r < reads.size(); ++r) l += lp[r][x];
	else for(uint32_t r = 0; r < reads.size(); ++r) l += std::log10(0.5) + _trLogSum(lp[r][x], lp[r][y]);
	gl[x][y] = l;
	if (l > bestL) { bestL = l; bestX = x; bestY = y; }
      }
    }

    // Alleles
    auto classRep = [&](int32_t k) {
      if ((k == 0) || (k == nc - 1)) return k;
      int32_t rep = k;
      for(int32_t j = 1; j < nc - 1; ++j) {
	if ((lenClass[j] == lenClass[k]) && (jct[amember[cand[j]]].alt.size() > jct[amember[cand[rep]]].alt.size())) rep = j;
      }
      return rep;
    };
    bestX = classRep(bestX);
    bestY = classRep(bestY);

    // Allele records
    double const small = -std::numeric_limits<double>::max();
    double total = small;
    for(int32_t x = 0; x < nc; ++x) {
      for(int32_t y = x; y < nc; ++y) {
	if (!validGt(x, y)) continue;
	total = (total == small) ? gl[x][y] : _trLogSum(total, gl[x][y]);
      }
    }
    double const lFar = std::log10(eps / W);
    // GT confidence
    auto classMargin = [&](int32_t k) {
      double g[3] = { small, small, small };
      for(int32_t x = 0; x < nc; ++x) {
	for(int32_t y = x; y < nc; ++y) {
	  if (!validGt(x, y)) continue;
	  int32_t copies = ((lenClass[x] == lenClass[k]) ? 1 : 0) + ((lenClass[y] == lenClass[k]) ? 1 : 0);
	  if ((ploidy == 1) && (copies > 0)) copies = 2;
	  g[copies] = (g[copies] == small) ? gl[x][y] : _trLogSum(g[copies], gl[x][y]);
	}
      }
      int32_t called = ((lenClass[bestX] == lenClass[k]) ? 1 : 0) + ((lenClass[bestY] == lenClass[k]) ? 1 : 0);
      if ((ploidy == 1) && (called > 0)) called = 2;
      double other = small;
      for(int32_t d = 0; d < 3; ++d) if ((d != called) && (g[d] != small) && (g[d] > other)) other = g[d];
      if ((g[called] == small) || (other == small)) return (double) -SMALLEST_GL;
      return std::max(0.0, g[called] - other);
    };
    // By length
    std::vector<int32_t> cls(nc);
    for(int32_t k = 0; k < nc; ++k) {
      cls[k] = k;
      if ((k != bestX) && (lenClass[k] == lenClass[bestX])) cls[k] = bestX;
      else if ((k != bestY) && (lenClass[k] == lenClass[bestY])) cls[k] = bestY;
    }
    // Alleles without a read
    double const farG2 = (double) reads.size() * lFar;
    double farG1 = small;
    if (ploidy == 2) {
      for(int32_t y = 0; y < nc; ++y) {
	double l = 0;
	for(uint32_t r = 0; r < reads.size(); ++r) l += std::log10(0.5) + _trLogSum(lFar, lp[r][y]);
	farG1 = (farG1 == small) ? l : _trLogSum(farG1, l);
      }
    }
    for(uint32_t a = 1; a < alen.size(); ++a) {
      int32_t ci = -1;
      for(int32_t k = 0; k < nc; ++k) if (cand[k] == (int32_t) a) ci = k;
      double g0 = small;
      double g1 = small;
      double g2 = small;
      if (ci >= 0) {
	bool called = ((ci == bestX) || (ci == bestY));
	for(int32_t x = 0; x < nc; ++x) {
	  for(int32_t y = x; y < nc; ++y) {
	    if (!validGt(x, y)) continue;
	    int32_t copies = ((x == ci) ? 1 : 0) + ((y == ci) ? 1 : 0);
	    if (called) copies = ((cls[x] == ci) ? 1 : 0) + ((cls[y] == ci) ? 1 : 0);
	    if (copies == 0) g0 = (g0 == small) ? gl[x][y] : _trLogSum(g0, gl[x][y]);
	    else if (copies == 1) g1 = (g1 == small) ? gl[x][y] : _trLogSum(g1, gl[x][y]);
	    else g2 = (g2 == small) ? gl[x][y] : _trLogSum(g2, gl[x][y]);
	  }
	}
      } else {
	// No read is close to this allele
	g0 = total;
	g1 = farG1;
	g2 = farG2;
      }
      if (ploidy == 1) { g1 = small; if (g2 == small) g2 = SMALLEST_GL + g0; }
      double mx = std::max(g0, std::max(g1, g2));
      double v[3] = { g0, g1, g2 };
      for(int32_t k = 0; k < 3; ++k) {
	v[k] = (v[k] == small) ? (double) SMALLEST_GL : (v[k] - mx);
	if (v[k] < SMALLEST_GL) v[k] = SMALLEST_GL;
      }
      int32_t copies = 0;
      if (ci >= 0) copies = ((bestX == ci) ? 1 : 0) + ((bestY == ci) ? 1 : 0);
      if ((ploidy == 1) && (copies > 0)) copies = 2;
      if ((ci >= 0) && (cls[ci] != ci)) {
	double m = classMargin(cls[ci]);
	v[0] = 0;
	v[1] = std::max((double) SMALLEST_GL, -m);
	v[2] = std::max((double) SMALLEST_GL, -2.0 * m);
	if (ploidy == 1) { v[2] = v[1]; v[1] = SMALLEST_GL; }
      }
      double sum = std::pow(10.0, v[0]) + std::pow(10.0, v[1]) + std::pow(10.0, v[2]);
      double post = std::pow(10.0, v[copies]) / sum;
      int32_t gq = (post >= 1.0) ? 10000 : (int32_t) std::round(-10.0 * std::log10(1.0 - post));
      if (gq > 10000) gq = 10000;
      if (gq < 0) gq = 0;

      typename TJctVector::value_type& jc = jct[amember[a]];
      jc.gl[0] = (float) v[0];
      jc.gl[1] = (float) v[1];
      jc.gl[2] = (float) v[2];
      jc.locusGt = (int8_t) copies;
      jc.locusGq = gq;
      jc.ref.clear(); jc.alt.clear();
      jc.hp1ref.clear(); jc.hp1alt.clear(); jc.hp2ref.clear(); jc.hp2alt.clear();
      jc.ps = -1;
      for(uint32_t r = 0; r < reads.size(); ++r) {
	// Assign to the best matching allele
	int32_t assigned = (lp[r][bestX] >= lp[r][bestY]) ? bestX : bestY;
	bool isAlt = ((ci >= 0) && (assigned == ci));
	uint8_t qual = 20;
	if (isAlt) {
	  jc.alt.push_back(qual);
	  if (reads[r].hp == 1) jc.hp1alt.push_back(qual);
	  else if (reads[r].hp == 2) jc.hp2alt.push_back(qual);
	  if ((reads[r].hp > 0) && (reads[r].ps >= 0) && (jc.ps < 0)) jc.ps = reads[r].ps;
	} else {
	  jc.ref.push_back(qual);
	  if (reads[r].hp == 1) jc.hp1ref.push_back(qual);
	  else if (reads[r].hp == 2) jc.hp2ref.push_back(qual);
	}
      }
    }
  }

  template<typename TConfig, typename TJunctionMap, typename TReadCountMap, typename TMethylMap>
  inline void
  genotypeLR(TConfig& c, std::vector<StructuralVariantRecord>& svs, TJunctionMap& jctMap, TReadCountMap& covMap, TMethylMap& methylMap) {
    typedef std::vector<StructuralVariantRecord> TSVs;
    if (svs.empty()) return;

    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    typedef std::vector<bam_hdr_t*> THeader;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    THeader hdr(c.files.size());
    for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
      hdr[file_c] = sam_hdr_read(samfile[file_c]);
    }

    // Count aligned reads per SV
    typedef std::vector<uint32_t> TSVReadCount;
    typedef std::vector<TSVReadCount> TSVFileReadCount;
    TSVFileReadCount readSV(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      readSV[file_c].resize(svs.size(), 0);
    }

    // Methylation
    typedef std::vector<MethylAccum> TSVMethylAccum;
    typedef std::vector<TSVMethylAccum> TFileMethylAccum;
    TFileMethylAccum methylAccum(c.files.size());
    for (unsigned int file_c = 0; file_c < c.files.size(); ++file_c) methylAccum[file_c].resize(svs.size());

    // Dump file
    boost::iostreams::filtering_ostream dumpOut;
    if (c.hasDumpFile) {
      dumpOut.push(boost::iostreams::gzip_compressor());
      dumpOut.push(boost::iostreams::file_sink(c.dumpfile.string(), std::ios_base::out | std::ios_base::binary));
      dumpOut << "#svid\tbam\tqname\tchr\tpos\tmapq\ttype" << std::endl;
    }

    // Genotype SVs
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "SV annotation" << std::endl;
    
    // Multi-allelic locus
    typedef std::vector<int32_t> TSvIds;
    std::vector<TSvIds> alleleGroup;
    std::vector<int32_t> groupOf(svs.size(), -1);
    {
      std::map<int32_t, int32_t> aidIdx;
      for(uint32_t i = 0; i < svs.size(); ++i) {
	if ((svs[i].alleleid < 0) || (svs[i].nallele < 2)) continue;
	std::map<int32_t, int32_t>::iterator it = aidIdx.find(svs[i].alleleid);
	if (it == aidIdx.end()) {
	  it = aidIdx.insert(std::make_pair(svs[i].alleleid, (int32_t) alleleGroup.size())).first;
	  alleleGroup.push_back(TSvIds());
	}
	alleleGroup[it->second].push_back(svs[i].id);
      }
      for(uint32_t g = 0; g < alleleGroup.size(); ++g) {
	if (alleleGroup[g].size() < 2) continue;
	for(uint32_t k = 0; k < alleleGroup[g].size(); ++k) {
	  groupOf[alleleGroup[g][k]] = g;
	  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) jctMap[file_c][alleleGroup[g][k]].joint = true;
	}
      }
    }

    // Tandem repeats
    int32_t const trFlank = 200;
    std::vector<TSvIds> trUnit;
    std::vector<int32_t> unitWs;
    std::vector<int32_t> unitWe;
    std::vector<int32_t> unitOf(svs.size(), -1);
    for(uint32_t g = 0; g < alleleGroup.size(); ++g) {
      if (alleleGroup[g].size() < 2) continue;
      typedef std::pair<int32_t, int32_t> TTract;
      std::vector<TTract> tracts;
      bool indel = true;
      for(uint32_t k = 0; k < alleleGroup[g].size(); ++k) {
	StructuralVariantRecord const& sv = svs[alleleGroup[g][k]];
	if ((sv.svt != 2) && (sv.svt != 4)) { indel = false; break; }
	if (sv.trStart >= 0) tracts.push_back(std::make_pair(sv.trStart, sv.trEnd));
      }
      if ((!indel) || (tracts.empty())) continue;
      std::sort(tracts.begin(), tracts.end());
      std::vector<TTract> win;
      for(uint32_t i = 0; i < tracts.size(); ++i) {
	if ((!win.empty()) && (tracts[i].first <= win.back().second + trFlank)) win.back().second = std::max(win.back().second, tracts[i].second);
	else win.push_back(tracts[i]);
      }
      for(uint32_t w = 0; w < win.size(); ++w) {
	if (win[w].second - win[w].first > 10000) continue;
	TSvIds mem;
	for(uint32_t k = 0; k < alleleGroup[g].size(); ++k) {
	  StructuralVariantRecord const& sv = svs[alleleGroup[g][k]];
	  int32_t svE = (sv.svt == 4) ? sv.svStart : sv.svEnd;
	  if ((unitOf[sv.id] < 0) && (sv.svStart >= win[w].first - trFlank) && (svE <= win[w].second + trFlank)) mem.push_back(sv.id);
	}
	if (mem.size() < 2) continue;
	for(uint32_t k = 0; k < mem.size(); ++k) unitOf[mem[k]] = (int32_t) trUnit.size();
	trUnit.push_back(mem);
	unitWs.push_back(win[w].first);
	unitWe.push_back(win[w].second);
      }
    }

    // Tandem repeats
    typedef std::vector<TrRead> TTrReads;
    std::vector<std::vector<TTrReads> > trReads(c.files.size(), std::vector<TTrReads>(trUnit.size()));
    std::vector<TrNoise> trNoise(c.files.size());

    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < (int32_t) hdr[0]->n_targets; ++refIndex) {
      // Fetch breakpoints
      typedef std::multimap<int32_t, int32_t> TBreakpointMap;
      TBreakpointMap bpMap;
      for(typename TSVs::iterator itSV = svs.begin(); itSV != svs.end(); ++itSV) {
	if (itSV->chr == refIndex) bpMap.insert(std::make_pair(itSV->svStart, itSV->id));
	if (itSV->chr2 == refIndex) bpMap.insert(std::make_pair(itSV->svEnd, itSV->id));
      }
      if (bpMap.empty()) continue;
      
      // Load sequence
      int32_t seqlen = -1;
      std::string tname(hdr[0]->target_name[refIndex]);
      char* seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr[0]->target_len[refIndex], &seqlen);

      // Take care of symbolic ALTs and SV annotation
      for(typename TSVs::iterator itSV = svs.begin(); itSV != svs.end(); ++itSV) {
	if ((itSV->chr == refIndex) && (itSV->alleles.empty())) itSV->alleles = _addAlleles(_refAnchor(seq, itSV->svStart, hdr[0]->target_len[refIndex]), std::string(hdr[0]->target_name[itSV->chr2]), *itSV, itSV->svt);

	// Annotate SVs
	if ((itSV->chr == refIndex) && (!_translocation(itSV->svt))) {
	  annotateSV(c, hdr[0], seq, *itSV);
	}
      }
      
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	
	// Coverage track
	typedef uint16_t TMaxCoverage;
	uint32_t maxCoverage = std::numeric_limits<TMaxCoverage>::max();
	typedef std::vector<TMaxCoverage> TBpCoverage;
	TBpCoverage covBases(hdr[file_c]->target_len[refIndex], 0);
	
	// Parse BAM
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr[file_c]->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  // Coverage track
	  if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
	  if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;

	  // Annotate coverage
	  {
	    uint32_t rp = rec->core.pos; // reference pointer
	    const uint32_t* cigar = bam_get_cigar(rec);
	    for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
		for(uint32_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k, ++rp) {
		  if ((rp < hdr[file_c]->target_len[refIndex]) && (covBases[rp] < maxCoverage - 1)) ++covBases[rp];
		}
	      } else if ((bam_cigar_op(cigar[i]) == BAM_CDEL) || (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP)) {
		rp += bam_cigar_oplen(cigar[i]);
	      }
	    }
	  }
	  
	  // Only primary alignments for genotyping (full sequence)
	  if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) continue;
	  if (rec->core.l_qseq < 2 * c.minimumFlankSize) continue;
	  
	  // Overlaps any SV breakpoints?
	  typedef std::set<int32_t> TSVSet;
	  TSVSet process;
	  int32_t rStart = _readStart(rec) + c.minimumFlankSize;
	  int32_t rEnd = _readEnd(rec);
	  if (rEnd > c.minimumFlankSize) {
	    rEnd -= c.minimumFlankSize;
	    if (rStart < rEnd) {
	      TBreakpointMap::const_iterator itBegin = bpMap.lower_bound(rStart);
	      TBreakpointMap::const_iterator itEnd = bpMap.upper_bound(rEnd);
	      for(;((itBegin != itEnd) && (itBegin != bpMap.end())); ++itBegin) process.insert(itBegin->second);
	    }
	  }

	  // Genotype SVs
	  // Read HP (haplotype) and PS (phase set) tags once per read
	  uint8_t hp = 0;
	  int32_t ps = -1;
	  {
	    uint8_t* hpTag = bam_aux_get(rec, "HP");
	    if (hpTag) hp = (uint8_t) bam_aux2i(hpTag);
	    uint8_t* psTag = bam_aux_get(rec, "PS");
	    if (psTag) ps = bam_aux2i(psTag);
	  }
	  std::string sequence;
	  std::vector<int8_t> methCall;
	  bool methCallBuilt = false;
	  bool hasMethyl = false;
	  // Edit distance scoring
	  auto scoreSV = [&](int32_t svid, std::vector<int32_t>& candidates, int32_t& refedsum, int32_t& altedsum, int32_t& nInform) {
	    refedsum = 0;
	    altedsum = 0;
	    nInform = 0;
	    // Which SV breakpoint does the read overlap
	    if ((svs[svid].chr == refIndex) && (svs[svid].svStart >= rStart) && (svs[svid].svStart <= rEnd)) candidates.push_back(svs[svid].svStart);
	    if ((svs[svid].chr2 == refIndex) && (svs[svid].svEnd >= rStart) && (svs[svid].svEnd <= rEnd)) candidates.push_back(svs[svid].svEnd);
	    if (candidates.empty()) return;

	    // Genotype breakpoints
	    std::vector<double> scoreR(candidates.size());
	    std::vector<double> scoreA(candidates.size());
	    for(uint32_t i = 0; i < candidates.size(); ++i) {
	      int32_t pos = candidates[i];
	      int32_t spBp = _findSeqBp(rec, pos);
	      int32_t consBp = svs[svid].consBp;
	      if (pos == svs[svid].svEnd) consBp += svs[svid].insLen;

	      // Find flanking sequence offsets
	      int32_t rStartOffset = pos - std::max(0, pos - spBp);
	      int32_t rEndOffset = std::min(pos + rec->core.l_qseq - spBp, (int32_t) hdr[file_c]->target_len[refIndex]) - pos;
	      int32_t cStartOffset = consBp - std::max(0, consBp - spBp);
	      int32_t cEndOffset = std::min(consBp + rec->core.l_qseq - spBp, (int32_t) svs[svid].consensus.size()) - consBp;
	      // Breakpoint should be in the middle so flanking sequences do not bias edit distance
	      int32_t offset = std::min(std::min(rStartOffset, cStartOffset), std::min(rEndOffset, cEndOffset));
	      if (offset < c.minimumFlankSize) continue;
	      if (!_translocation(svs[svid].svt) && (2 * offset < c.minConsWindow)) continue;

	      // Load sequence
	      if (sequence.empty()) {
		sequence.resize(rec->core.l_qseq);
		const uint8_t* seqptr = bam_get_seq(rec);
		for (int ik = 0; ik < rec->core.l_qseq; ++ik) sequence[ik] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, ik)];
	      }
	      std::string ref = boost::to_upper_copy(std::string(seq + pos - offset, seq + pos + offset));
	      std::string alt = svs[svid].consensus.substr(consBp - offset, 2 * offset);
	      std::string probe = sequence.substr(spBp - offset, 2 * offset);

	      // Edit distances
	      int32_t refScore = _editDistanceNW(ref, probe);
	      if ( ((svs[svid].svt == 0) && (pos == svs[svid].svEnd)) ||
		   ((svs[svid].svt == 1) && (pos == svs[svid].svStart)) ||
		   ((svs[svid].svt == 5) && (pos == svs[svid].svEnd)) ||
		   ((svs[svid].svt == 6) && (pos == svs[svid].svStart))
		   ) {
		reverseComplement(probe);
	      }
	      int32_t altScore = _editDistanceNW(alt, probe);

	      scoreA[i] = (1.0 - c.flankQuality) * alt.size();
	      scoreR[i] = (1.0 - c.flankQuality) * ref.size();
	      scoreA[i] = scoreA[i] / (double) (altScore + 1);
	      scoreR[i] = scoreR[i] / (double) (refScore + 1);
	      // Only breakpoints where the read matches at least one allele are informative
	      if ((scoreR[i] > 0.6) || (scoreA[i] > 0.6)) {
		refedsum += refScore;
		altedsum += altScore;
		++nInform;
	      }
	    }
	  };

	  // Read quality
	  auto edQual = [&](int32_t adelta) {
	    double w = std::log10((double) c.flankQuality / (double) (1.0 - c.flankQuality));
	    double ex = (double) adelta * w;
	    if (ex > 4.0) ex = 4.0;
	    uint32_t mq = (uint32_t) (10.0 * std::log10(1.0 + std::pow(10.0, ex)));
	    if (mq > (uint32_t) c.genoCap) mq = (uint32_t) c.genoCap;
	    return (uint8_t) mq;
	  };

	  // Assign the read to REF or ALT
	  auto assignRead = [&](int32_t svid, bool isAlt, uint8_t qual, std::vector<int32_t> const& candidates) {
	    if (!methCallBuilt) {
	      methCallBuilt = true;
	      hasMethyl = buildMethylCalls(rec, (uint8_t)c.methylProb, methCall);
	    }
	    if (hasMethyl) accumulateMethyl(c, rec, methCall, svs[svid], refIndex, (int32_t)hdr[file_c]->target_len[refIndex], isAlt, candidates, methylAccum[file_c][svid], sequence);
	    if (!isAlt) {
	      // REF-supporting read
	      jctMap[file_c][svid].ref.push_back(qual);
	      if (hp == 1) jctMap[file_c][svid].hp1ref.push_back(qual);
	      else if (hp == 2) jctMap[file_c][svid].hp2ref.push_back(qual);
	    } else {
	      // ALT-supporting read
	      if (c.hasDumpFile) {
		std::string svidStr(_addID(svs[svid].svt));
		std::string padNumber = boost::lexical_cast<std::string>(svid);
		padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
		svidStr += padNumber;
		dumpOut << svidStr << "\t" << c.files[file_c].string() << "\t" << bam_get_qname(rec) << "\t" << hdr[file_c]->target_name[rec->core.tid] << "\t" << rec->core.pos << "\t" << (int32_t) rec->core.qual << "\tSR" << std::endl;
	      }
	      jctMap[file_c][svid].alt.push_back(qual);
	      if (hp == 1) jctMap[file_c][svid].hp1alt.push_back(qual);
	      else if (hp == 2) jctMap[file_c][svid].hp2alt.push_back(qual);
	      if ((hp > 0) && (ps >= 0) && (jctMap[file_c][svid].ps < 0)) jctMap[file_c][svid].ps = ps;
	    }
	  };

	  std::set<int32_t> groupDone;

	  std::set<int32_t> trSeen;
	  for(typename TSVSet::const_iterator it = process.begin(); it != process.end(); ++it) {
	    int32_t svid = *it;
	    // Tandem repeat locus
	    if ((unitOf[svid] >= 0) && (trSeen.insert(unitOf[svid]).second)) {
	      int32_t u = unitOf[svid];
	      int32_t ws = unitWs[u] - trFlank;
	      int32_t we = unitWe[u] + trFlank;
	      if ((ws >= 0) && (we <= (int32_t) hdr[file_c]->target_len[refIndex]) && ((int32_t) rec->core.pos <= ws) && ((int32_t) bam_endpos(rec) >= we) && (trReads[file_c][u].size() < c.maxGenoReadCount)) {
		int32_t qa = _findSeqBp(rec, ws);
		int32_t qb = _findSeqBp(rec, we);
		if ((qa >= 0) && (qb >= 0)) {
		  TrRead tr;
		  tr.len = (qb - qa) - (we - ws);
		  tr.ps = ps;
		  tr.hp = hp;
		  trReads[file_c][u].push_back(tr);
		}
	      }
	    }
	    if (groupOf[svid] < 0) {
	      // Bi-allelic site
	      if ((jctMap[file_c][svid].ref.size() + jctMap[file_c][svid].alt.size()) >= c.maxGenoReadCount) continue;
	      // Enough candidates?
	      if (readSV[file_c][svid] >= c.maxGenoReadCount) continue;
	      ++readSV[file_c][svid];
	      std::vector<int32_t> candidates;
	      int32_t refedsum = 0;
	      int32_t altedsum = 0;
	      int32_t nInform = 0;
	      scoreSV(svid, candidates, refedsum, altedsum, nInform);
	      if (nInform == 0) continue;
	      // Use edit distance instead of mapq
	      int32_t delta = refedsum - altedsum;
	      int32_t adelta = (delta < 0) ? -delta : delta;
	      assignRead(svid, (delta > 0), edQual(adelta), candidates);
	    } else {
	      // Multi-allelic locus
	      int32_t g = groupOf[svid];
	      if (groupDone.find(g) != groupDone.end()) continue;
	      groupDone.insert(g);
	      std::vector<int32_t> members;
	      std::vector<int32_t> deltas;
	      std::vector<std::vector<int32_t> > memberCand;
	      for(uint32_t k = 0; k < alleleGroup[g].size(); ++k) {
		int32_t m = alleleGroup[g][k];
		if (process.find(m) == process.end()) continue;
		if ((jctMap[file_c][m].ref.size() + jctMap[file_c][m].alt.size()) >= c.maxGenoReadCount) continue;
		if (readSV[file_c][m] >= c.maxGenoReadCount) continue;
		++readSV[file_c][m];
		std::vector<int32_t> candidates;
		int32_t refedsum = 0;
		int32_t altedsum = 0;
		int32_t nInform = 0;
		scoreSV(m, candidates, refedsum, altedsum, nInform);
		if (nInform == 0) continue;
		members.push_back(m);
		deltas.push_back(refedsum - altedsum);
		memberCand.push_back(candidates);
	      }
	      for(uint32_t j = 0; j < members.size(); ++j) {
		int32_t bestOther = 0;
		for(uint32_t k = 0; k < members.size(); ++k) {
		  if ((k != j) && (deltas[k] > bestOther)) bestOther = deltas[k];
		}
		int32_t eff = bestOther - deltas[j];
		if (eff == 0) continue;
		assignRead(members[j], (eff < 0), edQual((eff < 0) ? -eff : eff), memberCand[j]);
	      }
	    }
	  }
	}
	// Clean-up
	bam_destroy1(rec);
	hts_itr_destroy(iter);

	// Tandem repeat loci
	{
	  uint8_t sex = (file_c < c.sexModel.sex.size()) ? c.sexModel.sex[file_c] : 0;
	  _trCalibrate(svs, trUnit, trReads[file_c], unitWs, unitWe, refIndex, trNoise[file_c]);
	  for(uint32_t u = 0; u < trUnit.size(); ++u) {
	    StructuralVariantRecord const& first = svs[trUnit[u][0]];
	    if (first.chr != refIndex) continue;
	    uint8_t ploidy = _svPloidy(c.sexModel, sex, first.chr, unitWs[u], first.chr, unitWe[u]);
	    _trLocusGenotype(svs, trUnit[u], trReads[file_c][u], trNoise[file_c], unitWs[u], unitWe[u], trFlank, ploidy, jctMap[file_c]);
	    TTrReads().swap(trReads[file_c][u]);
	  }
	}

	// Coverage annotation
	for(uint32_t i = 0; i < svs.size(); ++i) {
	  if (svs[i].chr == refIndex) {
	    int32_t halfSize = (svs[i].svEnd - svs[i].svStart)/2;
	    if ((_translocation(svs[i].svt)) || (svs[i].svt == 4)) halfSize = 500;

	    // Left region
	    int32_t lstart = std::max(svs[i].svStart - halfSize, 0);
	    int32_t lend = svs[i].svStart;
	    int32_t covbase = 0;
	    for(uint32_t k = lstart; ((k < (uint32_t) lend) && (k < hdr[file_c]->target_len[refIndex])); ++k) covbase += covBases[k];
	    covMap[file_c][svs[i].id].leftRC = covbase;
	  
	    // Actual SV
	    covbase = 0;
	    int32_t mstart = svs[i].svStart;
	    int32_t mend = svs[i].svEnd;
	    if ((_translocation(svs[i].svt)) || (svs[i].svt == 4)) {
	      mstart = std::max(svs[i].svStart - halfSize, 0);
	      mend = std::min(svs[i].svStart + halfSize, (int32_t) hdr[file_c]->target_len[refIndex]);
	    }
	    for(uint32_t k = mstart; ((k < (uint32_t) mend) && (k < hdr[file_c]->target_len[refIndex])); ++k) covbase += covBases[k];
	    covMap[file_c][svs[i].id].rc = covbase;
	  
	    // Right region
	    covbase = 0;
	    int32_t rstart = svs[i].svEnd;
	    int32_t rend = std::min(svs[i].svEnd + halfSize, (int32_t) hdr[file_c]->target_len[refIndex]);
	    if ((_translocation(svs[i].svt)) || (svs[i].svt == 4)) {
	      rstart = svs[i].svStart;
	      rend = std::min(svs[i].svStart + halfSize, (int32_t) hdr[file_c]->target_len[refIndex]);
	    }
	    for(uint32_t k = rstart; ((k < (uint32_t) rend) && (k < hdr[file_c]->target_len[refIndex])); ++k) covbase += covBases[k];
	    covMap[file_c][svs[i].id].rightRC = covbase;
	  }
	}
      }
      // Clean-up chromosome sequence
      if (seq != NULL) free(seq);
    }
    // Finalize methylation fractions from accumulated read counts
    for (unsigned int fc = 0; fc < c.files.size(); ++fc) {
      for (uint32_t i = 0; i < svs.size(); ++i) {
        finalizeMethylInfo(methylAccum[fc][i], methylMap[fc][i], c.minCpgDepth);
      }
    }

    // Clean-up
    fai_destroy(fai);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      bam_hdr_destroy(hdr[file_c]);	  
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }

}

#endif
