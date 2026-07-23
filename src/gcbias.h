#ifndef GCBIAS_H
#define GCBIAS_H

#include <boost/unordered_map.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/dynamic_bitset.hpp>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "scan.h"
#include "util.h"

namespace torali
{
  struct GcBias {
    int64_t sample;
    int32_t reference;
    double fractionSample;
    double fractionReference;
    double percentileSample;
    double percentileReference;
    double obsexp;
    double coverage;

    GcBias() : sample(0), reference(0), fractionSample(0), fractionReference(0), percentileSample(0), percentileReference(0), obsexp(0), coverage(0) {}
  };

  template<typename TConfig>
  inline std::pair<uint32_t, uint32_t>
  gcBound(TConfig const& c, std::vector<GcBias>& gcbias) {
    uint32_t lowerBound = 0;
    uint32_t upperBound = gcbias.size();
    for(uint32_t i = 0; i < gcbias.size(); ++i) {
      if ((gcbias[i].percentileSample < c.exclgc) || (gcbias[i].percentileReference < c.exclgc)) lowerBound = i;
      if ((gcbias[i].percentileSample + c.exclgc > 1) || (gcbias[i].percentileReference + c.exclgc > 1)) {
	if (i < upperBound) upperBound = i;
      }
    }
    if (lowerBound >= upperBound) upperBound = lowerBound + 1;
    return std::make_pair(lowerBound, upperBound);
  }


  // Map window GC fraction to region
  inline double
  regCorrFactor(std::vector<double> const& regcorr, double const gcfrac) {
    if (regcorr.empty()) return 1.0;
    int32_t b = (int32_t) (gcfrac * (double) (regcorr.size() - 1) + 0.5);
    if (b < 0) b = 0;
    if (b >= (int32_t) regcorr.size()) b = (int32_t) regcorr.size() - 1;
    return (regcorr[b] > 0) ? regcorr[b] : 1.0;
  }

  // Smooth GC curve
  inline void
  smoothFillCurve(std::vector<double>& curve, std::vector<double> const& weight) {
    int32_t n = (int32_t) curve.size();
    if (n < 3) return;
    double last = 0;
    bool have = false;
    for(int32_t i = 0; i < n; ++i) {
      if (weight[i] > 0) { last = curve[i]; have = true; }
      else if (have) curve[i] = last;
    }
    last = 0; have = false;
    for(int32_t i = n - 1; i >= 0; --i) {
      if (weight[i] > 0) { last = curve[i]; have = true; }
      else if (have) curve[i] = last;
    }
    // 3-bin smoothing
    std::vector<double> sm(curve);
    for(int32_t i = 1; i + 1 < n; ++i) {
      double w0 = weight[i-1] + 1.0;
      double w1 = 2.0 * (weight[i] + 1.0);
      double w2 = weight[i+1] + 1.0;
      sm[i] = (curve[i-1] * w0 + curve[i] * w1 + curve[i+1] * w2) / (w0 + w1 + w2);
    }
    curve = sm;
  }

  // Regional GC correction
  template<typename TConfig, typename TGCBound>
  inline void
  estimateRegionalGc(TConfig const& c, TGCBound const& gcbound, std::vector<GcBias> const& gcbias, std::vector<std::vector<ScanWindow> > const& scanCounts, uint32_t const regWin, std::vector<double>& regcorr) {
    uint32_t const nbin = 101;
    regcorr.assign(nbin, 1.0);
    std::vector<std::vector<double> > ratios(nbin);

    // Reference GC
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    faidx_t* faiRef = fai_load(c.genome.string().c_str());
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Estimate regional GC correction" << std::endl;

    uint32_t const sw = (c.scanWindow > 0) ? c.scanWindow : 10000;
    uint32_t const grp = std::max((uint32_t) 1, regWin / sw);

    for (int refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      if (scanCounts[refIndex].empty()) continue;
      std::string tname(hdr->target_name[refIndex]);
      int32_t seqlen = faidx_seq_len(faiRef, tname.c_str());
      if (seqlen == -1) continue;
      else seqlen = -1;
      char* ref = faidx_fetch_seq(faiRef, tname.c_str(), 0, faidx_seq_len(faiRef, tname.c_str()), &seqlen);
      if (ref == NULL) continue;

      // GC content
      std::vector<uint16_t> gcContent(hdr->target_len[refIndex], 0);
      {
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet gcref(hdr->target_len[refIndex], false);
	for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	  if ((ref[i] == 'c') || (ref[i] == 'C') || (ref[i] == 'g') || (ref[i] == 'G')) gcref[i] = 1;
	}
	int32_t halfwin = (int32_t) (c.meanisize / 2);
	int32_t gcsum = 0;
	for(int32_t pos = halfwin; pos < (int32_t) hdr->target_len[refIndex] - halfwin; ++pos) {
	  if (pos == halfwin) { for(int32_t i = pos - halfwin; i <= pos + halfwin; ++i) gcsum += gcref[i]; }
	  else { gcsum -= gcref[pos - halfwin - 1]; gcsum += gcref[pos + halfwin]; }
	  gcContent[pos] = gcsum;
	}
      }
      free(ref);

      // Regional windows
      uint32_t nb = (uint32_t) scanCounts[refIndex].size();
      for(uint32_t g0 = 0; g0 < nb; g0 += grp) {
	uint32_t g1 = std::min(nb, g0 + grp);
	double observed = 0;
	for(uint32_t bi = g0; bi < g1; ++bi) observed += (double) scanCounts[refIndex][bi].cov;
	int32_t rstart = scanCounts[refIndex][g0].start;
	int32_t rend = scanCounts[refIndex][g1 - 1].end;
	if (rend > (int32_t) hdr->target_len[refIndex]) rend = (int32_t) hdr->target_len[refIndex];
	if (rend <= rstart) continue;
	double fineExp = 0;
	double gcnum = 0;
	uint32_t winlen = 0;
	for(int32_t pos = rstart; pos < rend; ++pos) {
	  if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second)) {
	    fineExp += gcbias[gcContent[pos]].coverage;
	    gcnum += gcContent[pos];
	    ++winlen;
	  }
	}
	uint32_t totalBases = (uint32_t) (rend - rstart);
	if ((winlen >= totalBases / 2) && (fineExp > 0) && (observed > 0)) {
	  double obsValid = observed * ((double) winlen / (double) totalBases);
	  double gcfrac = (gcnum / (double) winlen) / (double) c.meanisize;
	  int32_t b = (int32_t) (gcfrac * (double) (nbin - 1) + 0.5);
	  if ((b >= 0) && (b < (int32_t) nbin)) ratios[b].push_back(obsValid / fineExp);
	}
      }
    }
    fai_destroy(faiRef);
    sam_close(samfile);
    bam_hdr_destroy(hdr);

    // Median per GC bin
    std::vector<double> weight(nbin, 0);
    double wsum = 0;
    double wtot = 0;
    for(uint32_t b = 0; b < nbin; ++b) {
      if (ratios[b].size() >= 10) {
	std::sort(ratios[b].begin(), ratios[b].end());
	double med = ratios[b][ratios[b].size() / 2];
	regcorr[b] = med;
	weight[b] = (double) ratios[b].size();
	wsum += med * (double) ratios[b].size();
	wtot += (double) ratios[b].size();
      } else regcorr[b] = 0;
    }
    // Preserve the ploidy baseline
    double mean = (wtot > 0) ? (wsum / wtot) : 1.0;
    if (mean > 0) {
      for(uint32_t b = 0; b < nbin; ++b) {
	if (weight[b] > 0) regcorr[b] /= mean;
      }
    }
    smoothFillCurve(regcorr, weight);
    for(uint32_t b = 0; b < nbin; ++b) {
      if (regcorr[b] <= 0) regcorr[b] = 1.0;
    }
  }

  template<typename TConfig, typename TGCBound>
  inline void
  gcBias(TConfig const& c, std::vector< std::vector<ScanWindow> > const& scanCounts, LibraryInfo const& li, std::vector<GcBias>& gcbias, TGCBound& gcbound) {
    // Load bam file
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Parse bam (contig by contig)
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Estimate GC bias" << std::endl;

    faidx_t* faiRef = fai_load(c.genome.string().c_str());
    for (int refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      if (scanCounts[refIndex].empty()) continue;

      // Bin map
      std::vector<uint16_t> binMap;
      if (c.hasScanFile) {
	// Fill bin map
	binMap.resize(hdr->target_len[refIndex], LAST_BIN);
	for(uint32_t bin = 0;((bin < scanCounts[refIndex].size()) && (bin < LAST_BIN)); ++bin) {
	  for(int32_t k = scanCounts[refIndex][bin].start; k < scanCounts[refIndex][bin].end; ++k) binMap[k] = bin;
	}
      }
      
      // Reference sequence
      std::string tname(hdr->target_name[refIndex]);
      int32_t seqlen = faidx_seq_len(faiRef, tname.c_str());
      if (seqlen == - 1) continue;
      else seqlen = -1;
      char* ref = faidx_fetch_seq(faiRef, tname.c_str(), 0, faidx_seq_len(faiRef, tname.c_str()), &seqlen);
      if (ref == NULL) continue;

      // Get GC content
      std::vector<uint16_t> uniqContent(hdr->target_len[refIndex], 0);
      std::vector<uint16_t> gcContent(hdr->target_len[refIndex], 0);
      {
	// GC map
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet gcref(hdr->target_len[refIndex], false);
	for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	  if ((ref[i] == 'c') || (ref[i] == 'C') || (ref[i] == 'g') || (ref[i] == 'G')) gcref[i] = 1;
	}

	// Sum across fragments
	int32_t halfwin = (int32_t) (c.meanisize / 2);
	int32_t gcsum = 0;
	for(int32_t pos = halfwin; pos < (int32_t) hdr->target_len[refIndex] - halfwin; ++pos) {
	  if (pos == halfwin) {
	    for(int32_t i = pos - halfwin; i<=pos+halfwin; ++i) gcsum += gcref[i];
	  } else {
	    gcsum -= gcref[pos - halfwin - 1];
	    gcsum += gcref[pos + halfwin];
	  }
	  gcContent[pos] = gcsum;
	}
      }

      // Coverage track
      typedef uint16_t TCount;
      uint32_t maxCoverage = std::numeric_limits<TCount>::max();
      typedef std::vector<TCount> TCoverage;
      TCoverage cov(hdr->target_len[refIndex], 0);
      TCoverage covUniq(hdr->target_len[refIndex], 0);
      TCoverage covTot;
      if (!c.basecov) covTot.resize(hdr->target_len[refIndex], 0);
      TCoverage& covMap = (!c.basecov) ? covTot : cov;

      // Mate map
      typedef boost::unordered_map<std::size_t, bool> TMateMap;
      TMateMap mateMap;
      
      // Parse BAM
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      int32_t lastAlignedPos = 0;
      std::set<std::size_t> lastAlignedPosReads;
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	if ((rec->core.flag & BAM_FPAIRED) && ((rec->core.flag & BAM_FMUNMAP) || (rec->core.tid != rec->core.mtid))) continue;
	if (rec->core.qual < c.minQual) continue;
	if (c.basecov) {
	  addBaseCoverage(rec, cov, covUniq, c.mapqUniq, hdr->target_len[refIndex], maxCoverage);
	  continue;
	}

	// Fill covTot
	addBaseCoverage(rec, covTot, covUniq, c.mapqUniq, hdr->target_len[refIndex], maxCoverage);

	int32_t midPoint = rec->core.pos + halfAlignmentLength(rec);
	if (rec->core.flag & BAM_FPAIRED) {
	  // Clean-up the read store for identical alignment positions
	  if (rec->core.pos > lastAlignedPos) {
	    lastAlignedPosReads.clear();
	    lastAlignedPos = rec->core.pos;
	  }

	  // Process pair
	  if ((rec->core.pos < rec->core.mpos) || ((rec->core.pos == rec->core.mpos) && (lastAlignedPosReads.find(hash_string(bam_get_qname(rec))) == lastAlignedPosReads.end()))) {
	    // First read
	    lastAlignedPosReads.insert(hash_string(bam_get_qname(rec)));
	    std::size_t hv = hash_pair(rec);
	    mateMap[hv]= true;
	    continue;
	  } else {
	    // Second read
	    std::size_t hv = hash_pair_mate(rec);
	    if ((mateMap.find(hv) == mateMap.end()) || (!mateMap[hv])) continue; // Mate discarded
	    mateMap[hv] = false;
	  }
	
	  // Insert size filter
	  int32_t isize = (rec->core.pos + alignmentLength(rec)) - rec->core.mpos;
	  if ((li.minNormalISize < isize) && (isize < li.maxNormalISize)) {
	    midPoint = rec->core.mpos + (int32_t) (isize/2);
	  } else {
	    if (rec->core.flag & BAM_FREVERSE) midPoint = rec->core.pos + alignmentLength(rec) - (c.meanisize / 2);
	    else midPoint = rec->core.pos + (c.meanisize / 2);
	  }
	}

	// Count fragment
	if ((midPoint >= 0) && (midPoint < (int32_t) hdr->target_len[refIndex]) && (cov[midPoint] < maxCoverage - 1)) ++cov[midPoint];
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);

      // Callable positions
      for(uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
	bool u;
	if (covMap[pos] == 0) u = ((ref[pos] != 'N') && (ref[pos] != 'n'));
	else u = (2 * (uint32_t) covUniq[pos] >= (uint32_t) covMap[pos]);
	uniqContent[pos] = (u ? (uint16_t) c.meanisize : 0);
      }
      // hom-del or unmappable?
      uint32_t maxHomDel = 1000000;
      uint32_t rstart = 0;
      while (rstart < hdr->target_len[refIndex]) {
	if (covMap[rstart] == 0) {
	  uint32_t rend = rstart;
	  while ((rend < hdr->target_len[refIndex]) && (covMap[rend] == 0)) ++rend;
	  bool leftOK = (rstart > 0) && (uniqContent[rstart - 1] > 0);
	  bool rightOK = (rend < hdr->target_len[refIndex]) && (uniqContent[rend] > 0);
	  if ((!leftOK) || (!rightOK) || (rend - rstart > maxHomDel)) {
	    for(uint32_t k = rstart; k < rend; ++k) uniqContent[k] = 0;
	  }
	  rstart = rend;
	} else ++rstart;
      }
      if (ref != NULL) free(ref);

      // Summarize GC coverage
      for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	bool uniqPos = (uniqContent[i] >= c.fragmentUnique * c.meanisize);
	if (uniqPos) {
	  // Valid bin?
	  int32_t bin = _findScanWindow(c, hdr->target_len[refIndex], binMap, i);
	  if ((bin >= 0) && (scanCounts[refIndex][bin].select)) {
	    ++gcbias[gcContent[i]].reference;
	    gcbias[gcContent[i]].sample += cov[i];
	    gcbias[gcContent[i]].coverage += cov[i];
	  }
	}
      }
    }
    
    // Normalize GC coverage
    for(uint32_t i = 0; i < gcbias.size(); ++i) {
      if (gcbias[i].reference) gcbias[i].coverage /= (double) gcbias[i].reference;
      else gcbias[i].coverage = 0;
    }
    // Smooth
    {
      std::vector<double> cvals(gcbias.size(), 0);
      std::vector<double> weight(gcbias.size(), 0);
      for(uint32_t i = 0; i < gcbias.size(); ++i) { cvals[i] = gcbias[i].coverage; weight[i] = (double) gcbias[i].reference; }
      smoothFillCurve(cvals, weight);
      for(uint32_t i = 0; i < gcbias.size(); ++i) gcbias[i].coverage = cvals[i];
    }

    // Determine percentiles
    uint64_t totalSampleCount = 0;
    uint64_t totalReferenceCount = 0;
    for(uint32_t i = 0; i < gcbias.size(); ++i) {
      totalSampleCount += gcbias[i].sample;
      totalReferenceCount += gcbias[i].reference;
    }
    uint64_t cumSample = 0;
    uint64_t cumReference = 0;
    for(uint32_t i = 0; i < gcbias.size(); ++i) {
      cumSample += gcbias[i].sample;
      cumReference += gcbias[i].reference;
      gcbias[i].fractionSample = (double) gcbias[i].sample / (double) totalSampleCount;
      gcbias[i].fractionReference = (double) gcbias[i].reference / (double) totalReferenceCount;
      gcbias[i].percentileSample = (double) cumSample / (double) totalSampleCount;
      gcbias[i].percentileReference = (double) cumReference / (double) totalReferenceCount;
      gcbias[i].obsexp = 1;
      if (gcbias[i].fractionReference > 0) gcbias[i].obsexp = gcbias[i].fractionSample / gcbias[i].fractionReference;
    }

    // Estimate correctable GC range
    gcbound = gcBound(c, gcbias);

    // Adjust correction to the callable range
    totalSampleCount = 0;
    totalReferenceCount = 0;
    for(uint32_t i = gcbound.first + 1; i < gcbound.second; ++i) {
      totalSampleCount += gcbias[i].sample;
      totalReferenceCount += gcbias[i].reference;
    }
    cumSample = 0;
    cumReference = 0;
    // Re-initialize
    for(uint32_t i = 0; i < gcbias.size(); ++i) {
      gcbias[i].fractionSample = 0;
      gcbias[i].fractionReference = 0;
      gcbias[i].percentileSample = 0;
      gcbias[i].percentileReference = 0;
      gcbias[i].obsexp = 1;
    }
    for(uint32_t i = gcbound.first + 1; i < gcbound.second; ++i) {
      cumSample += gcbias[i].sample;
      cumReference += gcbias[i].reference;
      gcbias[i].fractionSample = (double) gcbias[i].sample / (double) totalSampleCount;
      gcbias[i].fractionReference = (double) gcbias[i].reference / (double) totalReferenceCount;
      gcbias[i].percentileSample = (double) cumSample / (double) totalSampleCount;
      gcbias[i].percentileReference = (double) cumReference / (double) totalReferenceCount;
      gcbias[i].obsexp = 1;
      if (gcbias[i].fractionReference > 0) gcbias[i].obsexp = gcbias[i].fractionSample / gcbias[i].fractionReference;
    }
    
    fai_destroy(faiRef);
    hts_idx_destroy(idx);
    sam_close(samfile);
    bam_hdr_destroy(hdr);
  }

}

#endif
