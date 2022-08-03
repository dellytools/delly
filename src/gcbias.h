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
    int32_t sample;
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
    /*
    // Adjust total
    uint64_t totalSampleCount = 0;
    uint64_t totalReferenceCount = 0;
    for(uint32_t i = lowerBound + 1; i < upperBound; ++i) {
      totalSampleCount += gcbias[i].sample;
      totalReferenceCount += gcbias[i].reference;
    }    
    // Re-estimate observed/expected
    for(uint32_t i = lowerBound + 1; i < upperBound; ++i) {
      gcbias[i].fractionSample = (double) gcbias[i].sample / (double) totalSampleCount;
      gcbias[i].fractionReference = (double) gcbias[i].reference / (double) totalReferenceCount;
      gcbias[i].obsexp = 1;
      if (gcbias[i].fractionReference > 0) gcbias[i].obsexp = gcbias[i].fractionSample / gcbias[i].fractionReference;
    }
    */
    
    return std::make_pair(lowerBound, upperBound);
  }


  inline double
  getPercentIdentity(bam1_t const* rec, char const* seq) {
    // Sequence
    std::string sequence;
    sequence.resize(rec->core.l_qseq);
    uint8_t* seqptr = bam_get_seq(rec);
    for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

    // Reference slice
    std::string refslice = boost::to_upper_copy(std::string(seq + rec->core.pos, seq + lastAlignedPosition(rec)));
	      
    // Percent identity
    uint32_t rp = 0; // reference pointer
    uint32_t sp = 0; // sequence pointer
    uint32_t* cigar = bam_get_cigar(rec);
    int32_t matchCount = 0;
    int32_t mismatchCount = 0;
    for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	// match or mismatch
	for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k) {
	  if (sequence[sp] == refslice[rp]) ++matchCount;
	  else ++mismatchCount;
	  ++sp;
	  ++rp;
	}
      } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	mismatchCount += bam_cigar_oplen(cigar[i]);
	rp += bam_cigar_oplen(cigar[i]);
      } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	mismatchCount += bam_cigar_oplen(cigar[i]);
	sp += bam_cigar_oplen(cigar[i]);
      } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	sp += bam_cigar_oplen(cigar[i]);
      } else if(bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
      } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	mismatchCount += bam_cigar_oplen(cigar[i]);
	rp += bam_cigar_oplen(cigar[i]);
      } else {
	std::cerr << "Unknown Cigar options" << std::endl;
	return 1;
      }
    }
    double percid = 0;
    if (matchCount + mismatchCount > 0) percid = (double) matchCount / (double) (matchCount + mismatchCount);
    return percid;
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

    faidx_t* faiMap = fai_load(c.mapFile.string().c_str());
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
      
      // Check presence in mappability map
      std::string tname(hdr->target_name[refIndex]);
      int32_t seqlen = faidx_seq_len(faiMap, tname.c_str());
      if (seqlen == - 1) continue;
      else seqlen = -1;
      char* seq = faidx_fetch_seq(faiMap, tname.c_str(), 0, faidx_seq_len(faiMap, tname.c_str()), &seqlen);

      // Check presence in reference
      seqlen = faidx_seq_len(faiRef, tname.c_str());
      if (seqlen == - 1) continue;
      else seqlen = -1;
      char* ref = faidx_fetch_seq(faiRef, tname.c_str(), 0, faidx_seq_len(faiRef, tname.c_str()), &seqlen);

      // Get GC and Mappability
      std::vector<uint16_t> uniqContent(hdr->target_len[refIndex], 0);
      std::vector<uint16_t> gcContent(hdr->target_len[refIndex], 0);
      {
	// Mappability map
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet uniq(hdr->target_len[refIndex], false);
	for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	  if (seq[i] == 'C') uniq[i] = true;
	}

	// GC map
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet gcref(hdr->target_len[refIndex], false);
	for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	  if ((ref[i] == 'c') || (ref[i] == 'C') || (ref[i] == 'g') || (ref[i] == 'G')) gcref[i] = 1;
	}

	// Sum across fragments
	int32_t halfwin = (int32_t) (c.meanisize / 2);
	int32_t usum = 0;
	int32_t gcsum = 0;
	for(int32_t pos = halfwin; pos < (int32_t) hdr->target_len[refIndex] - halfwin; ++pos) {
	  if (pos == halfwin) {
	    for(int32_t i = pos - halfwin; i<=pos+halfwin; ++i) {
	      usum += uniq[i];
	      gcsum += gcref[i];
	    }
	  } else {
	    usum -= uniq[pos - halfwin - 1];
	    gcsum -= gcref[pos - halfwin - 1];
	    usum += uniq[pos + halfwin];
	    gcsum += gcref[pos + halfwin];
	  }
	  gcContent[pos] = gcsum;
	  uniqContent[pos] = usum;
	}
      }

      // Coverage track
      typedef uint16_t TCount;
      uint32_t maxCoverage = std::numeric_limits<TCount>::max();
      typedef std::vector<TCount> TCoverage;
      TCoverage cov(hdr->target_len[refIndex], 0);
      
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
      if (seq != NULL) free(seq);
      if (ref != NULL) free(ref);

      // Summarize GC coverage for this chromosome
      for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	if (uniqContent[i] >= c.fragmentUnique * c.meanisize) {
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
    fai_destroy(faiMap);
    hts_idx_destroy(idx);
    sam_close(samfile);
    bam_hdr_destroy(hdr);
  }

}

#endif
