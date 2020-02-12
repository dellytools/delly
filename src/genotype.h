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
#include <boost/progress.hpp>

#include <htslib/sam.h>

#include "util.h"

namespace torali
{

  template<typename TConfig>
  inline int32_t
  scoreThreshold(std::string const& s1, std::string const& s2, TConfig const& c) {
    int32_t l1 = 0;
    int32_t l2 = 0;
    for(uint32_t i = 0; i < s1.size(); ++i) {
      if (s1[i] != '-') ++l1;
    }
    for(uint32_t i = 0; i < s2.size(); ++i) {
      if (s2[i] != '-') ++l2;
    }
    int32_t minSize = std::min(l1, l2);
    return (int32_t) ((c.flankQuality * minSize * c.aliscore.match) + ((1.0 - c.flankQuality) * minSize * c.aliscore.mismatch));
  }
  
  template<typename TScore>
  inline int32_t
  scoreAlign(std::string const& s1, std::string const& s2, TScore const& sc) {
    int32_t score = 0;
    bool inGap = false;
    for(uint32_t i = 0; i < s1.size(); ++i) {
      if ((s1[i] == '-') || (s2[i] == '-')) {
	if (inGap) score += sc.ge;
	else {
	  inGap = true;
	  score += (sc.go + sc.ge);
	}
      } else {
	inGap = false;
	if (s1[i] == s2[i]) score += sc.match;
	else score += sc.mismatch;
      }
    }
    return score;
  }
  

  template<typename TConfig, typename TValidRegion, typename TJunctionMap, typename TReadCountMap>
  inline void
  trackRef(TConfig const& c, TValidRegion const& validRegions, std::vector<StructuralVariantRecord>& svs, TJunctionMap& jctMap, TReadCountMap& covMap) {
    typedef typename TValidRegion::value_type TChrIntervals;
    if (svs.empty()) return;

    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile[0]);
    
    // Parse genome chr-by-chr
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "SV annotation" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );

    // Collect ref. supporting reads
    typedef std::map<uint32_t, uint32_t> TQualMap;
    typedef std::vector<TQualMap> TSvSupport;
    TSvSupport svSup(svs.size(), TQualMap());

    // Coverage distribution
    typedef uint16_t TMaxCoverage;
    uint32_t maxCoverage = std::numeric_limits<TMaxCoverage>::max();
    typedef std::vector<uint32_t> TCovDist;
    TCovDist covDist(maxCoverage, 0);

    // Error rates
    uint64_t matchCount = 0;
    uint64_t mismatchCount = 0;
    uint64_t delCount = 0;
    uint64_t insCount = 0;

    // Read length distribution
    typedef uint16_t TMaxReadLength;
    uint32_t maxReadLength = std::numeric_limits<TMaxReadLength>::max();
    uint32_t rlBinSize = 100;
    typedef std::vector<uint32_t> TReadLengthDist;
    TReadLengthDist rlDist(maxReadLength * rlBinSize, 0);
    
    // Track reference reads
    char* seq = NULL;
    faidx_t* fai = fai_load(c.genome.string().c_str());
    
    // Iterate chromosomes
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      ++show_progress;
      if (validRegions[refIndex].empty()) continue;

      // Coverage track
      typedef std::vector<TMaxCoverage> TBpCoverage;
      TBpCoverage covBases(hdr->target_len[refIndex], 0);

      // Load sequence
      std::string tname(hdr->target_name[refIndex]);
      int32_t seqlen = -1;
      seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);

      // N-content
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet nrun(hdr->target_len[refIndex], false);
      for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	if ((seq[i] == 'n') || (seq[i] == 'N')) nrun[i] = 1;
      }

      // Flag breakpoints
      typedef std::set<int32_t> TIdSet;
      typedef std::map<uint32_t, TIdSet> TBpToIdMap;
      TBpToIdMap bpid;
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet bpOccupied(hdr->target_len[refIndex], false);
      for(uint32_t i = 0; i < svs.size(); ++i) {
	if ((svs[i].chr == refIndex) || (svs[i].chr2 == refIndex)) {
	  if (svs[i].chr == refIndex) {
	    svs[i].alleles = _addAlleles(boost::to_upper_copy(std::string(seq + svs[i].svStart - 1, seq + svs[i].svStart)), std::string(hdr->target_name[svs[i].chr2]), svs[i], svs[i].svt);
	    bpOccupied[svs[i].svStart] = 1;
	    if (bpid.find(svs[i].svStart) == bpid.end()) bpid.insert(std::make_pair(svs[i].svStart, TIdSet()));
	    bpid[svs[i].svStart].insert(svs[i].id);
	  }
	  if (svs[i].chr2 == refIndex) {
	    bpOccupied[svs[i].svEnd] = 1;
	    if (bpid.find(svs[i].svEnd) == bpid.end()) bpid.insert(std::make_pair(svs[i].svEnd, TIdSet()));
	    bpid[svs[i].svEnd].insert(svs[i].id);
	  }
	}
      }

      // Collect reads from all samples
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	// Read alignments
	for(typename TChrIntervals::const_iterator vRIt = validRegions[refIndex].begin(); vRIt != validRegions[refIndex].end(); ++vRIt) {
	  hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, vRIt->lower(), vRIt->upper());
	  bam1_t* rec = bam_init1();
	  while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	    // Keep secondary alignments
	    if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
	    if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;

	    // Get sequence
	    std::string sequence;
	    sequence.resize(rec->core.l_qseq);
	    uint8_t* seqptr = bam_get_seq(rec);
	    for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	    // Read length
	    if (!(rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP))) {
	      uint32_t slen = rec->core.l_qseq;
	      if (!slen) slen = sequenceLength(rec);
	      if (slen < maxReadLength * rlBinSize) ++rlDist[(int32_t) (slen / rlBinSize)];
	    }

	    // Hashing
	    std::size_t seed = hash_lr(rec);
	    uint32_t rp = rec->core.pos; // reference pointer
	    uint32_t sp = 0; // sequence pointer
      
	    // Parse the CIGAR
	    uint32_t* cigar = bam_get_cigar(rec);
	    for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
		// Fetch reference alignments
		int32_t bpHit = 0;
		for(uint32_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k) {
		  if ((rp < hdr->target_len[refIndex]) && (covBases[rp] < maxCoverage - 1)) ++covBases[rp];
		  if (bpOccupied[rp]) bpHit = rp;
		  if (sequence[sp] == seq[rp]) ++matchCount;
		  else ++mismatchCount;
		  ++sp;
		  ++rp;
		}
		if (bpHit) {
		  // Only spanning alignments
		  int32_t bg = std::max((int32_t) 0, bpHit - c.minimumFlankSize);
		  int32_t ed = std::min((int32_t) hdr->target_len[refIndex], bpHit + c.minimumFlankSize);
		  if ((bg < rec->core.pos) || ((uint32_t) ed > rec->core.pos + alignmentLength(rec))) continue;
		  
		  // Re-iterate cigar to track alignment		
		  int32_t rit = rec->core.pos; // reference pointer
		  int32_t sit = 0; // sequence pointer
		  std::string refAlign;
		  std::string altAlign;
		  for (std::size_t ij = 0; ij < rec->core.n_cigar; ++ij) {
		    if (rit >= ed) break;
		    if ((bam_cigar_op(cigar[ij]) == BAM_CMATCH) || (bam_cigar_op(cigar[ij]) == BAM_CEQUAL) || (bam_cigar_op(cigar[ij]) == BAM_CDIFF)) {
		      for(uint32_t k = 0; k < bam_cigar_oplen(cigar[ij]); ++k) {
			if ((rit >= bg) && (rit < ed)) {
			  refAlign.push_back(seq[rit]);
			  altAlign.push_back(sequence[sit]);
			}
			++rit;
			++sit;
		      }
		    } else if ((bam_cigar_op(cigar[ij]) == BAM_CDEL) || (bam_cigar_op(cigar[ij]) == BAM_CREF_SKIP)) {
		      for(uint32_t k = 0; k < bam_cigar_oplen(cigar[ij]); ++k) {
			if ((rit >= bg) && (rit < ed)) {
			  refAlign.push_back(seq[rit]);
			  altAlign.push_back('-');
			}
			++rit;
		      }
		    } else if (bam_cigar_op(cigar[ij]) == BAM_CINS) {
		      for(uint32_t k = 0; k < bam_cigar_oplen(cigar[ij]); ++k) {
			if ((rit >= bg) && (rit < ed)) {
			  refAlign.push_back('-');
			  altAlign.push_back(sequence[sit]);
			}
			++sit;
		      }
		    } else if (bam_cigar_op(cigar[ij]) == BAM_CSOFT_CLIP) {
		      sit += bam_cigar_oplen(cigar[ij]);
		    } else if (bam_cigar_op(cigar[ij]) == BAM_CHARD_CLIP) {
		      // Do nothing
		    }
		  }
		  int32_t score = scoreAlign(boost::to_upper_copy<std::string>(refAlign), boost::to_upper_copy<std::string>(altAlign), c.aliscore);
		  int32_t scth = scoreThreshold(refAlign, altAlign, c);
		  if ((score > 0) && (scth > 0)) {
		    int32_t qual = (int32_t) ((((float) score / (float) (scth)) - 1.0) * 100.0);
		    if ((qual > 0) && (qual > c.minMapQual)) {
		      // Debug
		      //std::cerr << "Qual: " << qual << std::endl;
		      //std::cerr << refAlign << std::endl;
		      //std::cerr << altAlign << std::endl;
		      //std::cerr << std::endl;
		  
		      // Assign read to all SVs at this position
		      for(typename TIdSet::const_iterator it = bpid[bpHit].begin(); it != bpid[bpHit].end(); ++it) {
			if (svSup[*it].find(seed) == svSup[*it].end()) svSup[*it].insert(std::make_pair(seed, qual));
			else {
			  if (svSup[*it][seed] < (uint32_t) qual) svSup[*it][seed] = (uint32_t) qual;
			}
		      }
		    }
		  }
		}
	      } else if ((bam_cigar_op(cigar[i]) == BAM_CDEL) || (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP)) {
		++delCount;
		rp += bam_cigar_oplen(cigar[i]);
	      } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
		++insCount;
		sp += bam_cigar_oplen(cigar[i]);
	      } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
		sp += bam_cigar_oplen(cigar[i]);
	      } else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
		// Do nothing
	      } else {
		std::cerr << "Unknown Cigar options" << std::endl;
	      }
	    }
	  }
	  bam_destroy1(rec);
	  hts_itr_destroy(iter);
	}
      }
      if (seq != NULL) free(seq);
      
      // Summarize coverage for this chromosome
      for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	if (!nrun[i]) ++covDist[covBases[i]];
      }
      
      // Assign SV support
      for(uint32_t i = 0; i < svs.size(); ++i) {
	if (svs[i].chr == refIndex) {
	  int32_t halfSize = (svs[i].svEnd - svs[i].svStart)/2;
	  if ((_translocation(svs[i].svt)) || (svs[i].svt == 4)) halfSize = 500;

	  // Left region
	  int32_t lstart = std::max(svs[i].svStart - halfSize, 0);
	  int32_t lend = svs[i].svStart;
	  int32_t covbase = 0;
	  for(uint32_t k = lstart; ((k < (uint32_t) lend) && (k < hdr->target_len[refIndex])); ++k) covbase += covBases[k];
	  covMap[0][svs[i].id].leftRC = covbase;

	  // Actual SV
	  covbase = 0;
	  int32_t mstart = svs[i].svStart;
	  int32_t mend = svs[i].svEnd;
	  if ((_translocation(svs[i].svt)) || (svs[i].svt == 4)) {
	    mstart = std::max(svs[i].svStart - halfSize, 0);
	    mend = std::min(svs[i].svStart + halfSize, (int32_t) hdr->target_len[refIndex]);
	  }
	  for(uint32_t k = mstart; ((k < (uint32_t) mend) && (k < hdr->target_len[refIndex])); ++k) covbase += covBases[k];
	  covMap[0][svs[i].id].rc = covbase;

	  // Right region
	  covbase = 0;
	  int32_t rstart = svs[i].svEnd;
	  int32_t rend = std::min(svs[i].svEnd + halfSize, (int32_t) hdr->target_len[refIndex]);
	  if ((_translocation(svs[i].svt)) || (svs[i].svt == 4)) {
	    rstart = svs[i].svStart;
	    rend = std::min(svs[i].svStart + halfSize, (int32_t) hdr->target_len[refIndex]);
	  }
	  for(uint32_t k = rstart; ((k < (uint32_t) rend) && (k < hdr->target_len[refIndex])); ++k) covbase += covBases[k];
	  covMap[0][svs[i].id].rightRC = covbase;
	}
      }
    }
    fai_destroy(fai);

    // Assign SV support
    for(uint32_t i = 0; i < svs.size(); ++i) {
      for (typename TQualMap::const_iterator it = svSup[svs[i].id].begin(); it != svSup[svs[i].id].end(); ++it) {
	jctMap[0][svs[i].id].ref.push_back(it->second);
      }
    }

    // Output coverage info
    uint64_t totalCovCount = 0;
    for (uint32_t i = 0; i < covDist.size(); ++i) totalCovCount += covDist[i];
    std::vector<uint32_t> covPercentiles(5, 0);  // 5%, 25%, 50%, 75%, 95%
    uint64_t cumCovCount = 0;
    for (uint32_t i = 0; i < covDist.size(); ++i) {
      cumCovCount += covDist[i];
      double frac = (double) cumCovCount / (double) totalCovCount;
      if (frac < 0.05) covPercentiles[0] = i + 1;
      if (frac < 0.25) covPercentiles[1] = i + 1;
      if (frac < 0.5) covPercentiles[2] = i + 1;
      if (frac < 0.75) covPercentiles[3] = i + 1;
      if (frac < 0.95) covPercentiles[4] = i + 1;
    }
    std::cout << "Coverage distribution (^COV)" << std::endl;
    std::cout << "COV\t95% of bases are >= " << covPercentiles[0] << "x" << std::endl;
    std::cout << "COV\t75% of bases are >= " << covPercentiles[1] << "x" << std::endl;
    std::cout << "COV\t50% of bases are >= " << covPercentiles[2] << "x" << std::endl;
    std::cout << "COV\t25% of bases are >= " << covPercentiles[3] << "x" << std::endl;
    std::cout << "COV\t5% of bases are >= " << covPercentiles[4] << "x" << std::endl;

    
    // Output read length info
    uint64_t totalRlCount = 0;
    for (uint32_t i = 0; i < rlDist.size(); ++i) totalRlCount += rlDist[i];
    std::vector<uint32_t> rlPercentiles(5, 0);  // 5%, 25%, 50%, 75%, 95%
    uint64_t cumRlCount = 0;
    for (uint32_t i = 0; i < rlDist.size(); ++i) {
      cumRlCount += rlDist[i];
      double frac = (double) cumRlCount / (double) totalRlCount;
      if (frac < 0.05) rlPercentiles[0] = (i + 1) * rlBinSize;
      if (frac < 0.25) rlPercentiles[1] = (i + 1) * rlBinSize;
      if (frac < 0.5) rlPercentiles[2] = (i + 1) * rlBinSize;
      if (frac < 0.75) rlPercentiles[3] = (i + 1) * rlBinSize;
      if (frac < 0.95) rlPercentiles[4] = (i + 1) * rlBinSize;
    }
    std::cout << "Read-length distribution (^RL)" << std::endl;
    std::cout << "RL\t95% of reads are >= " << rlPercentiles[0] << "bp" << std::endl;
    std::cout << "RL\t75% of reads are >= " << rlPercentiles[1] << "bp" << std::endl;
    std::cout << "RL\t50% of reads are >= " << rlPercentiles[2] << "bp" << std::endl;
    std::cout << "RL\t25% of reads are >= " << rlPercentiles[3] << "bp" << std::endl;
    std::cout << "RL\t5% of reads are >= " << rlPercentiles[4] << "bp" << std::endl;

    // Output sequencing error rates
    uint64_t alignedbases = matchCount + mismatchCount + delCount + insCount;
    std::cout << "Sequencing error rates (^ERR)" << std::endl;
    std::cout << "ERR\tMatchRate\t" << (double) matchCount / (double) alignedbases << std::endl;
    std::cout << "ERR\tMismatchRate\t" << (double) mismatchCount / (double) alignedbases << std::endl;
    std::cout << "ERR\tDeletionRate\t" << (double) delCount / (double) alignedbases << std::endl;
    std::cout << "ERR\tInsertionRate\t" << (double) insCount / (double) alignedbases << std::endl;

    // Clean-up
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }
     
  
}

#endif
