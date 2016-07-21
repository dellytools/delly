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


 template<typename TConfig, typename TValidRegion, typename TSVs, typename TTag>
 inline bool
 findGappedReads(TConfig const& c, TValidRegion const& validRegions, TSVs& svs,  SVType<TTag> svType) {
   // Open file handles
   typedef std::vector<samFile*> TSamFile;
   typedef std::vector<hts_idx_t*> TIndex;
   TSamFile samfile(c.files.size());
   TIndex idx(c.files.size());
   for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
     samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
     idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
   }
   bam_hdr_t* hdr = sam_hdr_read(samfile[0]);

   // Parse genome
   unsigned int totalSplitReadsAligned = 0;
   boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
   std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Multiple sequence alignment" << std::endl;
   boost::progress_display show_progress( hdr->n_targets );

   faidx_t* fai = fai_load(c.genome.string().c_str());
   // Find reference index
   for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
     ++show_progress;
     if (validRegions[refIndex].empty()) continue;
     char* seq = NULL;
	
     // Iterate all structural variants on this chromosome
     for(typename TSVs::iterator svIt = svs.begin(); svIt!=svs.end(); ++svIt) {
       if ((svIt->chr != refIndex) && (svIt->chr2 != refIndex)) continue;
       
       // Lazy loading of reference sequence
       if (seq == NULL) {
	 int32_t seqlen = -1;
	 std::string tname(hdr->target_name[refIndex]);
	 seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
       }

       // Find gapped reads
       if (svIt->chr == refIndex) {
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
		 bool foundN = false;
		 for(uint32_t spointer = sp; spointer < sp + bam_cigar_oplen(cigar[i]); ++spointer) {
		   if (sequence[spointer] == 'N') { foundN = true; break;  }
		 }
		 if (!foundN) {
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
		     } else if ((i + 1 == rec->core.n_cigar) && (std::abs(svIt->svStart - countpos) < 25)) {
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
		   for(int32_t pq = spStart; pq < spEnd; ++pq) qualSum += quality[pq];
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
	  
	   //std::cerr << hdr->target_name[svIt->chr] << ',' << svIt->svStart << ',' << svIt->svEnd << ',' << (svIt->svEnd - svIt->svStart) << ',' << svIt->srSupport << std::endl;
	  
	   // Get the SV reference
	   svIt->svStartBeg = std::max(0, svIt->svStart - (int32_t) (svIt->consensus.size()));
	   svIt->svEndEnd = std::min((int32_t) (hdr->target_len[svIt->chr]), (int32_t) (svIt->svEnd + svIt->consensus.size()));
	   std::string svRefStr = _getSVRef(c, seq, *svIt, refIndex, svType);

	   // Align consensus to reference
	   if (!alignConsensus(c, *svIt, svRefStr, svType)) { svIt->consensus = ""; svIt->srSupport = 0; }

	   //std::cerr << hdr->target_name[svIt->chr] << ',' << svIt->svStart << ',' << svIt->svEnd << ',' << (svIt->svEnd - svIt->svStart) << ',' << svIt->srSupport << std::endl;
	 }
       }
     }
     if (seq != NULL) free(seq);
   }
   // Clean-up
   fai_destroy(fai);
   for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
     hts_idx_destroy(idx[file_c]);
     sam_close(samfile[file_c]);
   }
   bam_hdr_destroy(hdr);
   
   return (totalSplitReadsAligned>0);
 }

 template<typename TConfig, typename TValidRegion, typename TVariants, typename TSampleLib, typename TTag>
 inline void
 longPacBio(TConfig const& c, TValidRegion const& validRegions, TVariants& svs, TSampleLib&, SVType<TTag> svType)
 {
   typedef typename TValidRegion::value_type TChrIntervals;
   
   int32_t probe_size = 15;

   // Open file handles
   typedef std::vector<samFile*> TSamFile;
   typedef std::vector<hts_idx_t*> TIndex;
   TSamFile samfile(c.files.size());
   TIndex idx(c.files.size());
   for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
     samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
    idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
   }
   bam_hdr_t* hdr = sam_hdr_read(samfile[0]);

   // Parse genome, process chromosome by chromosome
   boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
   std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Scanning PacBio Alignments" << std::endl;
   boost::progress_display show_progress( hdr->n_targets );
   for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
     ++show_progress;
     if (validRegions[refIndex].empty()) continue;
      
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
	   if (c.indels) {
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
		 // Count leading or trailing soft-clips as Cigar D
		 int32_t countpos = rec->core.pos - vRIt->lower() + rp;
		 bool foundN = false;
		 for(uint32_t spointer = sp; spointer < sp + bam_cigar_oplen(cigar[i]); ++spointer) {
		   if (sequence[spointer] == 'N') { foundN = true; break;  }
		 }
		 if (!foundN) {
		   double meanqual = 0;
		   for(uint32_t spointer = sp; spointer < sp + bam_cigar_oplen(cigar[i]); ++spointer) meanqual += quality[spointer];
		   meanqual /= (double) bam_cigar_oplen(cigar[i]);
		   double clipfraction = (double) bam_cigar_oplen(cigar[i]) / (double) rec->core.l_qseq;
		   if ((meanqual >= c.minMapQual) && (clipfraction < 0.3)) {
		     if (i == 0) {
		       for(std::size_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k, --countpos) {
			 if ((countpos >= 0) && (countpos < interval_size) && (gapcount[countpos] < 255)) ++gapcount[countpos];
		       }
		     } else if (i + 1 == rec->core.n_cigar) {
		       for(std::size_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k, ++countpos) {
			 if ((countpos >= 0) && (countpos < interval_size) && (gapcount[countpos] < 255)) ++gapcount[countpos];
		       }
		     }
		   }
		 }
		 sp += bam_cigar_oplen(cigar[i]);
	       }
	     }
	   }
	 }
	 bam_destroy1(rec);
	 hts_itr_destroy(iter);

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
	   double gaprate = gapsum / (mmsum + gapsum);
	   if ((support >= 2) && (gaprate >= gaprateUth)) {
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
	     if ((svEnd - svStart) <= (int32_t) c.indelsize) {
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
	       svRec.homLen = 0;
	       svRec.id = svs.size();
	       svs.push_back(svRec);
	       //std::cerr << hdr->target_name[svRec.chr] << ',' << svRec.svStart << ',' << svRec.svEnd << ',' << (svRec.svEnd - svRec.svStart) << std::endl;
	     }
	     svStart = -1;
	     svEnd = -1;
	   }
	 }
       }
     }
   }

   // Gapped alignment search
   if (!svs.empty()) {
     if (c.indels) {
       findGappedReads(c, validRegions, svs, svType);
	  
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

       // Re-number SVs
       uint32_t cliqueCount = 0;
       for(typename TVariants::iterator svIt = svs.begin(); svIt != svs.end(); ++svIt) svIt->id = cliqueCount++;
     }
   }

   // Clean-up
   bam_hdr_destroy(hdr);
   for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
     hts_idx_destroy(idx[file_c]);
     sam_close(samfile[file_c]);
   }
 }


}

#endif
