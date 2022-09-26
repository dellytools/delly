#ifndef SCAN_H
#define SCAN_H

#include <limits>

#include <boost/icl/split_interval_map.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>

#include <htslib/sam.h>

#include "version.h"
#include "util.h"


namespace torali
{

  struct ScanWindow {
    bool select;
    int32_t start;
    int32_t end;
    uint32_t cov;
    uint32_t uniqcov;

    ScanWindow() : select(false), start(0), end(0), cov(0), uniqcov(0) {}
    explicit ScanWindow(int32_t const s) : select(false), start(s), end(s+1), cov(0), uniqcov(0) {}
  };

  template<typename TScanWindow>
  struct SortScanWindow : public std::binary_function<TScanWindow, TScanWindow, bool>
  {
    inline bool operator()(TScanWindow const& sw1, TScanWindow const& sw2) {
      return ((sw1.start<sw2.start) || ((sw1.start == sw2.start) && (sw1.end < sw2.end)));
    }
  };


  template<typename TConfig>
  inline int32_t
  _findScanWindow(TConfig const& c, uint32_t const reflen, std::vector<uint16_t> const& binMap, int32_t const midPoint) {
    if (c.hasScanFile) {
      if (binMap[midPoint] == LAST_BIN) return -1;
      else return binMap[midPoint];
    } else {
      uint32_t bin = midPoint / c.scanWindow;
      uint32_t allbins = reflen / c.scanWindow;
      if (bin >= allbins) return -1;
      else return bin;
    }
    return -1;
  }

  
  inline std::pair<uint32_t, uint32_t>
  estCountBounds(std::vector< std::vector<ScanWindow> > const& scanCounts) {
    std::vector<uint32_t> all;
    for(uint32_t refIndex = 0; refIndex < scanCounts.size(); ++refIndex) {
      for(uint32_t i = 0; i<scanCounts[refIndex].size(); ++i) {
	if (scanCounts[refIndex][i].select) all.push_back(scanCounts[refIndex][i].cov);
      }
    }
    std::sort(all.begin(), all.end());
    uint32_t median = all[all.size() / 2];
    std::vector<uint32_t> absdev;
    for(uint32_t i = 0; i<all.size(); ++i) absdev.push_back(std::abs((int32_t) all[i] - (int32_t) median));
    std::sort(absdev.begin(), absdev.end());
    uint32_t mad = absdev[absdev.size() / 2];
    uint32_t lowerBound = 0;
    if (mad < median) lowerBound = median - mad;
    uint32_t upperBound = median + mad;
    return std::make_pair(lowerBound, upperBound);
  }

  template<typename TConfig>
  inline void
  scan(TConfig const& c, LibraryInfo const& li, std::vector< std::vector<ScanWindow> >& scanCounts) {

    // Load bam file
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Pre-defined scanning windows
    if (c.hasScanFile) {
      typedef boost::icl::interval_set<uint32_t> TChrIntervals;
      typedef std::vector<TChrIntervals> TRegionsGenome;
      TRegionsGenome scanRegions;
      if (!_parseBedIntervals(c.scanFile.string(), c.hasScanFile, hdr, scanRegions)) {
	std::cerr << "Warning: Couldn't parse BED intervals. Do the chromosome names match?" << std::endl;
      }
      for (int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
	for(typename TChrIntervals::iterator it = scanRegions[refIndex].begin(); it != scanRegions[refIndex].end(); ++it) {
	  if (it->lower() < it->upper()) {
	    if (it->upper() < hdr->target_len[refIndex]) {
	      ScanWindow sw;
	      sw.start = it->lower();
	      sw.end = it->upper();
	      sw.select = true;
	      scanCounts[refIndex].push_back(sw);
	    }
	  }
	}
	// Sort scan windows
	sort(scanCounts[refIndex].begin(), scanCounts[refIndex].end(), SortScanWindow<ScanWindow>());
      }
    }
    
    // Parse BAM file
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Scanning Windows" << std::endl;

    // Iterate chromosomes
    uint64_t totalCov = 0;
    faidx_t* faiMap = fai_load(c.mapFile.string().c_str());
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      if (chrNoData(c, refIndex, idx)) continue;
      // Exclude small chromosomes
      if ((hdr->target_len[refIndex] < c.minChrLen) && (totalCov > 1000000)) continue;
      // Exclude sex chromosomes
      if ((std::string(hdr->target_name[refIndex]) == "chrX") || (std::string(hdr->target_name[refIndex]) == "chrY") || (std::string(hdr->target_name[refIndex]) == "X") || (std::string(hdr->target_name[refIndex]) == "Y")) continue;

      // Check presence in mappability map
      std::string tname(hdr->target_name[refIndex]);
      int32_t seqlen = faidx_seq_len(faiMap, tname.c_str());
      if (seqlen == -1) continue;
      else seqlen = -1;
      char* seq = faidx_fetch_seq(faiMap, tname.c_str(), 0, faidx_seq_len(faiMap, tname.c_str()), &seqlen);

      // Get Mappability
      std::vector<uint16_t> uniqContent(hdr->target_len[refIndex], 0);
      {
	// Mappability map
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet uniq(hdr->target_len[refIndex], false);
	for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	  if (seq[i] == 'C') uniq[i] = 1;
	}

	// Sum across fragments
	int32_t halfwin = (int32_t) (c.meanisize / 2);
	int32_t usum = 0;
	for(int32_t pos = halfwin; pos < (int32_t) hdr->target_len[refIndex] - halfwin; ++pos) {
	  if (pos == halfwin) {
	    for(int32_t i = pos - halfwin; i<=pos+halfwin; ++i) usum += uniq[i];
	  } else {
	    usum -= uniq[pos - halfwin - 1];
	    usum += uniq[pos + halfwin];
	  }
	  uniqContent[pos] = usum;
	}
      }

      // Bins on this chromosome
      std::vector<uint16_t> binMap;
      if (!c.hasScanFile) {
	uint32_t allbins = hdr->target_len[refIndex] / c.scanWindow;
	scanCounts[refIndex].resize(allbins, ScanWindow());
	for(uint32_t i = 0; i < allbins; ++i) {
	  scanCounts[refIndex][i].start = i * c.scanWindow;
	  scanCounts[refIndex][i].end = (i+1) * c.scanWindow;
	}
      } else {
	// Fill bin map
	binMap.resize(hdr->target_len[refIndex], LAST_BIN);
	if (scanCounts[refIndex].size() >= LAST_BIN) {
	  std::cerr << "Warning: Too many scan windows on " << hdr->target_name[refIndex] << std::endl;
	}
	for(uint32_t bin = 0;((bin < scanCounts[refIndex].size()) && (bin < LAST_BIN)); ++bin) {
	  for(int32_t k = scanCounts[refIndex][bin].start; k < scanCounts[refIndex][bin].end; ++k) binMap[k] = bin;
	}
      }
	
      // Mate map
      typedef boost::unordered_map<std::size_t, bool> TMateMap;
      TMateMap mateMap;

      // Count reads
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      int32_t lastAlignedPos = 0;
      std::set<std::size_t> lastAlignedPosReads;
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	if ((rec->core.flag & BAM_FPAIRED) && ((rec->core.flag & BAM_FMUNMAP) || (rec->core.tid != rec->core.mtid))) continue;
	if (rec->core.qual < c.minQual) continue;
	if (getSVType(rec) != 2) continue;

	int32_t midPoint = rec->core.pos + halfAlignmentLength(rec);
	if (rec->core.flag & BAM_FPAIRED) {
	  // Clean-up the read store for identical alignment positions
	  if (rec->core.pos > lastAlignedPos) {
	    lastAlignedPosReads.clear();
	    lastAlignedPos = rec->core.pos;
	  }
	
	  if ((rec->core.pos < rec->core.mpos) || ((rec->core.pos == rec->core.mpos) && (lastAlignedPosReads.find(hash_string(bam_get_qname(rec))) == lastAlignedPosReads.end()))) {
	    // First read
	    lastAlignedPosReads.insert(hash_string(bam_get_qname(rec)));
	    std::size_t hv = hash_pair(rec);
	    mateMap[hv] = true;
	    continue;
	  } else {
	    // Second read
	    std::size_t hv = hash_pair_mate(rec);
	    if ((mateMap.find(hv) == mateMap.end()) || (!mateMap[hv])) continue; // Mate discarded
	    mateMap[hv] = false;
	  }

	  // Insert size filter
	  int32_t isize = (rec->core.pos + alignmentLength(rec)) - rec->core.mpos;
	  if ((li.minNormalISize < isize) && (isize < li.maxNormalISize)) midPoint = rec->core.mpos + (int32_t) (isize/2);
	  else continue;
	}

	// Count fragment
	if ((midPoint >= 0) && (midPoint < (int32_t) hdr->target_len[refIndex])) {
	  int32_t bin = _findScanWindow(c, hdr->target_len[refIndex], binMap, midPoint);
	  if (bin >= 0) {
	    ++scanCounts[refIndex][bin].cov;
	    
	    if (uniqContent[midPoint] >= c.fragmentUnique * c.meanisize) ++scanCounts[refIndex][bin].uniqcov;
	    ++totalCov;
	  }
	}
      }
      // Clean-up
      bam_destroy1(rec);
      hts_itr_destroy(iter);
      if (seq != NULL) free(seq);
    }
    
    // clean-up
    fai_destroy(faiMap);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }


  template<typename TConfig>
  inline void
  selectWindows(TConfig const& c, std::vector< std::vector<ScanWindow> >& scanCounts) {
    if (c.noScanWindowSelection) {
      // Select all windows
      for(uint32_t refIndex = 0; refIndex < scanCounts.size(); ++refIndex) {
	for(uint32_t i = 0; i<scanCounts[refIndex].size(); ++i) {
	  scanCounts[refIndex][i].select = true;
	}
      }
    } else {
      // Pre-screen using PE layout, uniqueness and percent identity
      for(uint32_t refIndex = 0; refIndex < scanCounts.size(); ++refIndex) {
	for(uint32_t i = 0; i<scanCounts[refIndex].size(); ++i) {
	  // Uniqueness
	  double uniqratio = 0;
	  if (scanCounts[refIndex][i].cov > 0) uniqratio = (double) scanCounts[refIndex][i].uniqcov / scanCounts[refIndex][i].cov;
	  if (uniqratio > c.uniqueToTotalCovRatio) scanCounts[refIndex][i].select = true;
	  else scanCounts[refIndex][i].select = false;
	}
      }
      
      // Normalize user-defined scan windows to same length (10,000bp)
      if (c.hasScanFile) {
	for(uint32_t refIndex = 0; refIndex < scanCounts.size(); ++refIndex) {
	  for(uint32_t i = 0; i<scanCounts[refIndex].size(); ++i) {
	    double scale = (double) 10000 / (double) (scanCounts[refIndex][i].end - scanCounts[refIndex][i].start);
	    scanCounts[refIndex][i].uniqcov *= scale;
	    scanCounts[refIndex][i].cov *= scale;
	  }
	}
      }
      
      // Get "normal" coverage windows of CN2
      typedef std::pair<uint32_t, uint32_t> TCountBounds;
      TCountBounds cb = estCountBounds(scanCounts);
      
      // Select CN2 windows
      for(uint32_t refIndex = 0; refIndex < scanCounts.size(); ++refIndex) {
	for(uint32_t i = 0; i<scanCounts[refIndex].size(); ++i) {
	  if (scanCounts[refIndex][i].select) {
	    if ((scanCounts[refIndex][i].cov > cb.first) && (scanCounts[refIndex][i].cov < cb.second)) scanCounts[refIndex][i].select = true;
	    else scanCounts[refIndex][i].select = false;
	  }
	}
      }
    }
  }
  
}

#endif
