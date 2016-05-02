/*
============================================================================
DELLY: Structural variant discovery by integrated PE mapping and SR analysis
============================================================================
Copyright (C) 2012 Tobias Rausch

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Tobias Rausch (rausch@embl.de)
============================================================================
*/

#ifndef PACBIO_H
#define PACBIO_H

#include <boost/dynamic_bitset.hpp>

namespace torali {

struct Peak {
  double score;
  int32_t support;
  int32_t pos;

Peak(double s, int32_t u, int32_t p) : score(s), support(u), pos(p) {}
};

// Sort peaks
template<typename TRecord>
struct SortPeaks : public std::binary_function<TRecord, TRecord, bool>
{
  inline bool operator()(TRecord const& s1, TRecord const& s2) const {
    return (s1.score > s2.score);
  }
};


inline uint32_t
lastAlignedPosition(bam1_t const* rec) {
  uint32_t* cigar = bam_get_cigar(rec);
  uint32_t alen = 0;
  for (std::size_t i = 0; i < rec->core.n_cigar; ++i)
    if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CDEL)) alen += bam_cigar_oplen(cigar[i]);
  return rec->core.pos + alen;
}


template<typename TConfig, typename TStructuralVariantRecord, typename TTag>
inline bool
findGappedReads(TConfig const& c, std::vector<TStructuralVariantRecord>& svs,  SVType<TTag> svType) {
  typedef std::vector<TStructuralVariantRecord> TSVs;
  
  // Open file handles
  typedef std::vector<samFile*> TSamFile;
  typedef std::vector<hts_idx_t*> TIndex;
  TSamFile samfile;
  TIndex idx;
  samfile.resize(c.files.size());
  idx.resize(c.files.size());
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
    idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
  }
  bam_hdr_t* hdr = sam_hdr_read(samfile[0]);

  // Parse genome
  unsigned int totalSplitReadsAligned = 0;
  kseq_t *seq;
  int l;
  gzFile fp = gzopen(c.genome.string().c_str(), "r");
  seq = kseq_init(fp);
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Multiple sequence alignment" << std::endl;
  boost::progress_display show_progress( hdr->n_targets );
  while ((l = kseq_read(seq)) >= 0) {
    // Find reference index
    for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
      if (std::string(seq->name.s) == std::string(hdr->target_name[refIndex])) {
	++show_progress;
	
	// Iterate all structural variants on this chromosome
	typename TSVs::iterator svIt = svs.begin();
	typename TSVs::iterator svItEnd = svs.end();
	for(;svIt!=svItEnd; ++svIt) {
	  if (svIt->chr == refIndex) {
	    // Get the SV reference
	    std::string svRefStr = boost::to_upper_copy(std::string(seq->seq.s + svIt->svStartBeg, seq->seq.s + svIt->svEndEnd));
	    if (svRefStr.find('N') != std::string::npos) continue; // Ignore regions with nearby Ns in the reference for the time being

	    typedef std::set<std::string> TSplitReadSet;
	    TSplitReadSet splitReadSet;
	    std::vector<uint8_t> mapQV;

#pragma omp parallel for default(shared)
	    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	      hts_itr_t* iter = sam_itr_queryi(idx[file_c], svIt->chr, svIt->svStartBeg, svIt->svEndEnd);
	      bam1_t* rec = bam_init1();
	      while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
		if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
		if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	    
		// Get the sequence
		std::string sequence;
		sequence.resize(rec->core.l_qseq);
		uint8_t* seqptr = bam_get_seq(rec);
		for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	      
		// Get the quality vector
		typedef std::vector<uint8_t> TQuality;
		TQuality quality;
		quality.resize(rec->core.l_qseq);
		uint8_t* qualptr = bam_get_qual(rec);
		for (int i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
	      
		// Get the reference slice
		int32_t rp = 0; // reference pointer
		int32_t sp = 0; // sequence pointer
		int32_t spStart = -1;
		int32_t spEnd = -1;
	      
		// Count gaps and aligned bases
		uint32_t totalcount = 0;
		uint32_t gapcount = 0;
		uint32_t* cigar = bam_get_cigar(rec);
		for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
		  if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
		    int32_t countpos = rec->core.pos + rp;
		    // match or mismatch
		    for(uint32_t k = 0; k<bam_cigar_oplen(cigar[i]);++k, ++countpos) 
		      if ((countpos >= svIt->svStart) && (countpos < svIt->svEnd)) {
			if (spStart == -1) {
			  spStart = sp;
			  spEnd = sp;
			}
			++totalcount;
			++spEnd;
		      }
		    rp += bam_cigar_oplen(cigar[i]);
		    sp += bam_cigar_oplen(cigar[i]);
		  } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
		    int32_t countpos = rec->core.pos + rp;
		    for(uint32_t k = 0; k<bam_cigar_oplen(cigar[i]);++k, ++countpos)
		      if ((countpos >= svIt->svStart) && (countpos < svIt->svEnd)) {
			if (spStart == -1) {
			  spStart = sp;
			  spEnd = sp;
			}
			++gapcount;
		      }
		    rp += bam_cigar_oplen(cigar[i]);
		  } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
		    // ToDo
		    sp += bam_cigar_oplen(cigar[i]);
		  } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
		    // Count leading or trailing soft-clips as Cigar D
		    int32_t countpos = rec->core.pos + rp;
		    double meanqual = 0;
		    for(uint32_t spointer = sp; spointer < sp + bam_cigar_oplen(cigar[i]); ++spointer) meanqual += quality[spointer];
		    meanqual /= (double) bam_cigar_oplen(cigar[i]);
		    double clipfraction = (double) bam_cigar_oplen(cigar[i]) / (double) rec->core.l_qseq;
		    if ((meanqual >= c.minMapQual) && (clipfraction < 0.3)) {
		      if ((i == 0) && (std::abs(svIt->svEnd - countpos) < 25)) {
			for(std::size_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k, --countpos) {
			  if ((countpos >= svIt->svStart) && (countpos < svIt->svEnd)) {
			    if (spStart == -1) {
			      spStart = sp;
			      spEnd = sp;
			    }
			    ++gapcount;
			  }
			}
		      } else if ((i == rec->core.n_cigar -1) && (std::abs(svIt->svStart - countpos) < 25)) {
			for(std::size_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k, ++countpos) {
			  if ((countpos >= svIt->svStart) && (countpos < svIt->svEnd)) {
			    if (spStart == -1) {
			      spStart = sp;
			      spEnd = sp;
			    }
			    ++gapcount;
			  }
			}
		      }
		    }
		    sp += bam_cigar_oplen(cigar[i]);
		  }
		}

		spStart = std::max(spStart - 500, 0);
		spEnd = std::min(spEnd + 500, (int32_t) sequence.size());
		sequence = boost::to_upper_copy(sequence.substr(spStart, spEnd - spStart));
		if (sequence.find('N') == std::string::npos) {
		  uint32_t withinSvSize = gapcount + totalcount;
		  if (withinSvSize > 0) {
		    double gapRate = (double) gapcount / (double) (withinSvSize);
		    if (gapRate >= 0.75) {
		      unsigned int qualSum = 0;
		      for(unsigned int pq = spStart; pq < spEnd; ++pq) qualSum += quality[pq];
#pragma omp critical
		      {
			splitReadSet.insert(sequence);
			mapQV.push_back((uint8_t)(qualSum / (spEnd - spStart)));
		      }
		    }
		  }
		}
	      }
	      bam_destroy1(rec);
	      hts_itr_destroy(iter);
	    }
	    totalSplitReadsAligned += splitReadSet.size();

	    // Process reads for this SV
	    if (splitReadSet.size() > 1) {
	      // Set median base quality
	      std::sort(mapQV.begin(), mapQV.end());
	      svIt->peMapQuality = mapQV[mapQV.size()/2];

	      // Calculate MSA
	      svIt->srSupport = msa(c, splitReadSet, svIt->consensus);
	      
	      //std::cerr << hdr->target_name[svIt->chr] << ',' << svIt->svStart << ',' << svIt->svEnd << ',' << (svIt->svEnd - svIt->svStart) << ',' << svIt->srSupport << ',' << svIt->consensus << std::endl;
	      
	      // Search true split in candidates
	      if (!alignConsensus(c, *svIt, svRefStr, svType)) {
		svIt->consensus = ""; 
		svIt->srSupport = 0; 
	      }
	    }
	  }
	}
      }
    }
  }
  kseq_destroy(seq);
  gzclose(fp);

  // Clean-up
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    hts_idx_destroy(idx[file_c]);
    sam_close(samfile[file_c]);
  }

  return (totalSplitReadsAligned>0);
}

template<typename TConfig, typename TSVType>
inline int pacbioRun(TConfig const& c, TSVType svType) {

#ifdef PROFILE
  ProfilerStart("delly.prof");
#endif

  // Collect all promising structural variants
  typedef std::vector<StructuralVariantRecord> TVariants;
  TVariants svs;
  
  // Clique id counter
  unsigned int clique_count = 1;

  // Scanning probe size
  int32_t probe_size = 15;

  // Open file handles
  typedef std::vector<samFile*> TSamFile;
  typedef std::vector<hts_idx_t*> TIndex;
  TSamFile samfile;
  TIndex idx;
  samfile.resize(c.files.size());
  idx.resize(c.files.size());
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
    idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
  }
  bam_hdr_t* hdr = sam_hdr_read(samfile[0]);
  typedef std::vector<uint32_t> TRefLength;
  TRefLength reflen;
  for (int i = 0; i<hdr->n_targets; ++i) reflen.push_back(hdr->target_len[i]);

  // Exclude intervals
  typedef boost::icl::interval_set<int32_t> TChrIntervals;
  typedef typename TChrIntervals::interval_type TIVal;
  typedef std::vector<TChrIntervals> TRegionsGenome;
  TRegionsGenome validRegions;
  validRegions.resize(hdr->n_targets);
  TRegionsGenome exclg;
  exclg.resize(hdr->n_targets);
  std::vector<bool> validChr;
  validChr.resize(hdr->n_targets, true);
  if (c.hasExcludeFile) {
    std::ifstream chrFile(c.exclude.string().c_str(), std::ifstream::in);
    if (chrFile.is_open()) {
      while (chrFile.good()) {
	std::string chrFromFile;
	getline(chrFile, chrFromFile);
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t,;");
	Tokenizer tokens(chrFromFile, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter!=tokens.end()) {
	  std::string chrName = *tokIter++;
	  int32_t tid = bam_name2id(hdr, chrName.c_str());
	  if (tid >= 0) {
	    if (tokIter!=tokens.end()) {
	      int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
	      int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
	      exclg[tid].insert(TIVal::right_open(start, end));
	    } else validChr[tid] = false; // Exclude entire chromosome
	  }
	}
      }
      chrFile.close();
    }
  }

  // Create the valid regions
  for (int i = 0; i<hdr->n_targets; ++i) {
    int32_t istart = 0;
    for(typename TChrIntervals::iterator it = exclg[i].begin(); it != exclg[i].end(); ++it) {
      if (istart + 1 < it->lower()) validRegions[i].insert(TIVal::right_open(istart, it->lower() - 1));
      istart = it->upper();
    }
    if (istart + 1 < (int32_t) hdr->target_len[i]) validRegions[i].insert(TIVal::right_open(istart, hdr->target_len[i]));
  }
  exclg.clear();

  // Debug valid regions
  //for (int i = 0; i<hdr->n_targets; ++i) {
  //for(typename TChrIntervals::const_iterator vRIt = validRegions[i].begin(); vRIt != validRegions[i].end(); ++vRIt) {
  //std::string chrName = hdr->target_name[i];
  //std::cerr << chrName << '\t' << vRIt->lower() << '\t' << vRIt->upper() << std::endl;
  //}
  //}

  // Parse genome, process chromosome by chromosome
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Scanning PacBio Alignments" << std::endl;
  boost::progress_display show_progress( hdr->n_targets );
  if (!c.hasVcfFile) {
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      ++show_progress;
      if (!validChr[refIndex]) continue;
      
      // Collect all gaps
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet gaps(hdr->target_len[refIndex]);

      // Iterate all samples
#pragma omp parallel for default(shared)
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	// Read alignments
	for(typename TChrIntervals::const_iterator vRIt = validRegions[refIndex].begin(); vRIt != validRegions[refIndex].end(); ++vRIt) {
	  int32_t interval_size = vRIt->upper() - vRIt->lower();
	  std::vector<uint8_t> totalcount(interval_size, 0);
	  std::vector<uint8_t> gapcount(interval_size, 0);
	  std::vector<Peak> peaks;
	  
	  hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, vRIt->lower(), vRIt->upper());
	  bam1_t* rec = bam_init1();
	  while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	    if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	    if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	    
	    // Small indel detection
	    if ((c.indels) && (_smallIndelDetection(svType))) {
	      // Get the sequence
	      std::string sequence;
	      sequence.resize(rec->core.l_qseq);
	      uint8_t* seqptr = bam_get_seq(rec);
	      for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	      
	      // Get the quality vector
	      typedef std::vector<uint8_t> TQuality;
	      TQuality quality;
	      quality.resize(rec->core.l_qseq);
	      uint8_t* qualptr = bam_get_qual(rec);
	      for (int i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
	      
	      // Get the reference slice
	      int32_t rp = 0; // reference pointer
	      int32_t sp = 0; // sequence pointer
	      
	      uint32_t* cigar = bam_get_cigar(rec);
	      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
		if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
		  int32_t countpos = rec->core.pos - vRIt->lower() + rp;
		  // match or mismatch
		  for(uint32_t k = 0; k<bam_cigar_oplen(cigar[i]);++k, ++countpos) 
		    if ((countpos >= 0) && (countpos < interval_size) && (totalcount[countpos] < 255)) ++totalcount[countpos];
		  rp += bam_cigar_oplen(cigar[i]);
		  sp += bam_cigar_oplen(cigar[i]);
		} else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
		  int32_t countpos = rec->core.pos - vRIt->lower() + rp;
		  for(uint32_t k = 0; k<bam_cigar_oplen(cigar[i]);++k, ++countpos)
		    if ((countpos >= 0) && (countpos < interval_size) && (gapcount[countpos] < 255)) ++gapcount[countpos];
		  rp += bam_cigar_oplen(cigar[i]);
		} else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
		  // ToDo
		  sp += bam_cigar_oplen(cigar[i]);
		} else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
		  /*
		  // Count leading or trailing soft-clips as Cigar D
		  int32_t countpos = rec->core.pos - vRIt->lower() + rp;
		  double meanqual = 0;
		  for(uint32_t spointer = sp; spointer < sp + bam_cigar_oplen(cigar[i]); ++spointer) meanqual += quality[spointer];
		  meanqual /= (double) bam_cigar_oplen(cigar[i]);
		  double clipfraction = (double) bam_cigar_oplen(cigar[i]) / (double) rec->core.l_qseq;
		  if ((meanqual >= c.minMapQual) && (clipfraction < 0.3)) {
		    if (i == 0) {
		      for(std::size_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k, --countpos) {
			if ((countpos >= 0) && (countpos < interval_size) && (gapcount[countpos] < 255)) ++gapcount[countpos];
		      }
		    } else if (i == rec->core.n_cigar -1) {
		      for(std::size_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k, ++countpos) {
			if ((countpos >= 0) && (countpos < interval_size) && (gapcount[countpos] < 255)) ++gapcount[countpos];
		      }
		    }
		  }
		  */
		  sp += bam_cigar_oplen(cigar[i]);
		}
	      }
	    }
	  }
	  bam_destroy1(rec);
	  hts_itr_destroy(iter);

	  /*
	  // Annotate deletions
	  for(uint32_t i = 0; i < gi[refIndex].size(); ++i) {
	  int32_t s = gi[refIndex][i].first;
	  int32_t e = gi[refIndex][i].second;
	  if ((e < vRIt->lower()) || (s > vRIt->upper())) continue;
	  std::string chrName = hdr->target_name[refIndex];
	  int32_t from = std::max(s - vRIt->lower(), 0);
	  int32_t to = std::min(e - vRIt->lower(), interval_size);
	  double mmsum = 0;
	  double gapsum = 0;
	  for(int32_t i = from; i < to; ++i) {
	    mmsum += totalcount[i];
	    gapsum += gapcount[i];
	  }
	  double gapcov = gapsum / (mmsum + gapsum);
	  std::cerr << chrName << '\t' << s << '\t' << e << "\t" << mmsum << "\t" << gapsum << "\t" << gapcov << std::endl;
	  }
	  */

	  // Find deletions
	  double mmsum = 0;
	  double gapsum = 0;
	  double gaprateUth = 0.3;
	  double gaprateLth = 0.1;
	  double xdrop = 0.1;
	  for(int32_t i = 0; ((i<probe_size) && (i<interval_size)); ++i) {
	    mmsum += totalcount[i];
	    gapsum += gapcount[i];
	  }
	  int32_t bestSupport = 0;
	  double bestGapRate = 0;
	  int32_t bestGapIndex = 0;
	  for(int32_t i = probe_size; i<interval_size; ++i) {
	    int32_t support = (int32_t) (gapsum / probe_size);
	    if (support >= 2) {
	      double gaprate = gapsum / (mmsum + gapsum);
	      if (gaprate >= gaprateUth) {
		if (gaprate > bestGapRate) {
		  bestGapRate = gaprate;
		  bestGapIndex = i;
		  bestSupport = support;
		}
	      } else if (gaprate < gaprateLth) {
		if (bestGapIndex) {
		  peaks.push_back(Peak(bestGapRate , bestSupport, bestGapIndex));
		  bestGapIndex = 0;
		  bestGapRate = 0;
		  bestSupport = 0;
		}
	      }
	    }
	    mmsum -= totalcount[i - probe_size];
	    gapsum -= gapcount[i - probe_size];
	    mmsum += totalcount[i];
	    gapsum += gapcount[i];
	  }
	  if (bestGapIndex) peaks.push_back(Peak(bestGapRate , bestSupport, bestGapIndex));
	  std::sort(peaks.begin(), peaks.end(), SortPeaks<Peak>());
	  
	  // Extend seed gaps
	  for(uint32_t i = 0; i < peaks.size(); ++i) {
	    // Initialization
	    double bestGapRate = peaks[i].score;
	    int32_t svStartIdx = peaks[i].pos - probe_size;
	    int32_t svEndIdx = peaks[i].pos;
	    double mmsum = 0;
	    double gapsum = 0;
	    for(int32_t k = peaks[i].pos - probe_size; k < peaks[i].pos; ++k) {
	      mmsum += totalcount[k];
	      gapsum += gapcount[k];
	    }
	    
	    // Extension
	    while (true) {
	      double rightgaprate = 0;
	      if (svEndIdx < interval_size) {
		double rightmmsum = mmsum + totalcount[svEndIdx];
		double rightgapsum = gapsum + gapcount[svEndIdx];
		rightgaprate = rightgapsum / (rightmmsum + rightgapsum);
	      }
	      double leftgaprate = 0;
	      if (svStartIdx - 1 >= 0) {
		double leftmmsum = mmsum + totalcount[svStartIdx - 1];
		double leftgapsum = gapsum + gapcount[svStartIdx - 1];
		leftgaprate = leftgapsum / (leftmmsum + leftgapsum);
	      }
	      
	      double bestextensionrate = std::max(leftgaprate, rightgaprate);
	      if (bestextensionrate + xdrop < bestGapRate) break;
	      else {
		// Accept extension
		if (leftgaprate > rightgaprate) {
		  mmsum += totalcount[svStartIdx - 1];
		  gapsum += gapcount[svStartIdx - 1];
		  --svStartIdx;
		} else {
		  mmsum += totalcount[svEndIdx];
		  gapsum += gapcount[svEndIdx];
		  ++svEndIdx;
		}
		double gaprate = gapsum / (mmsum + gapsum);
		if (gaprate > bestGapRate) bestGapRate = gaprate;
	      }
	    }
	    double sup = (gapsum / (double) (svEndIdx - svStartIdx));
	    if (sup >= 2) {
#pragma omp critical
	      {
		for(TBitSet::size_type pos = vRIt->lower() + svStartIdx; pos < (TBitSet::size_type) (vRIt->lower() + svEndIdx); ++pos) gaps[pos] = 1;
	      }
	    }
	  }
	}
      }

      // Output deletions for this chromosome
      std::string chrName = hdr->target_name[refIndex];
      int32_t svStart = -1;
      int32_t svEnd = -1;
      for(TBitSet::size_type pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
	if (gaps[pos]) {
	  if (svStart == -1) {
	    svStart = pos;
	    svEnd = pos;
	  } else {
	    ++svEnd;
	  }
	} else {
	  if (svStart != -1) {
	    // Output deletion
	    StructuralVariantRecord svRec;
	    svRec.chr = refIndex;
	    svRec.chr2 = refIndex;
	    svRec.svStartBeg = std::max(svStart - 500, 0);
	    svRec.svStart = svStart;
	    svRec.svStartEnd = std::min((uint32_t) svStart + 500, hdr->target_len[refIndex]);
	    svRec.svEndBeg = std::max((int32_t) svEnd - 500, 0);
	    svRec.svEnd = svEnd;
	    svRec.svEndEnd = std::min((uint32_t) svEnd + 500, hdr->target_len[refIndex]);
	    svRec.peSupport = 0;
	    svRec.wiggle = 2*probe_size;
	    svRec.peMapQuality = 0;
	    svRec.srSupport = 2;
	    svRec.srAlignQuality = 0;
	    svRec.precise=true;
	    svRec.ct=_getCT(svType);
	    svRec.insLen = 0;
	    svRec.id = clique_count++;
	    svs.push_back(svRec);
	    svStart = -1;
	    svEnd = -1;
	  }
	}
      }
    }
  }

  // Clean-up
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    hts_idx_destroy(idx[file_c]);
    sam_close(samfile[file_c]);
  }
  bam_hdr_destroy(hdr);

  // Gapped alignment search
  if (!c.hasVcfFile) {
    if (boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome)) {
      if (!svs.empty()) {
	if (_smallIndelDetection(svType)) {
	  findGappedReads(c, svs, svType);
	  
	  // Sort SVs for look-up and by decreasing PE support
	  sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());
	  
	  // Temporary SV container
	  TVariants svc;
	  
	  // Clean-up SV set
	  for(typename TVariants::iterator svIt = svs.begin(); svIt != svs.end(); ++svIt) {
	    // Unresolved gapped alignments
	    if ((svIt->precise) && (svIt->srAlignQuality == 0)) continue;
	    // Add SV
	    svc.push_back(*svIt);
	  }
	  // Swap back
	  svs = svc;
	}
      }
    }
  } else {
    // Read SV records from input vcffile
    vcfParse(c, reflen, probe_size, svs, svType);
  }

  // Any SVs for genotyping
  if (svs.empty()) {
    std::cout << "No structural variants found!" << std::endl;
    std::cout << "Done." << std::endl;
    return 0;
  }

  // Annotate junction reads
  typedef std::pair<std::string, int> TSampleSVPair;
  typedef std::pair<std::vector<uint8_t>, std::vector<uint8_t> > TReadQual;
  typedef boost::unordered_map<TSampleSVPair, TReadQual> TJunctionCountMap;
  TJunctionCountMap junctionCountMap;
  //if (boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))
  //  _annotateJunctionReads(c, sampleLib, svs, junctionCountMap, svType);

  // Annotate spanning coverage
  typedef boost::unordered_map<TSampleSVPair, TReadQual> TCountMap;
  TCountMap spanCountMap;
  //_annotateSpanningCoverage(c, sampleLib, svs, spanCountMap, svType);

  // Annotate coverage
  typedef boost::unordered_map<TSampleSVPair, ReadCount> TRCMap;
  TRCMap rcMap;
  //_annotateCoverage(c, refnames, sampleLib, svs, rcMap, svType);

  // VCF output
  if (svs.size()) vcfOutput(c, svs, junctionCountMap, rcMap, spanCountMap, svType);

#ifdef PROFILE
  ProfilerStop();
#endif

  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
  return 0;
}


}

#endif
