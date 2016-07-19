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

#ifndef JUNCTION_H
#define JUNCTION_H


#include <boost/multiprecision/cpp_int.hpp>

#include <htslib/faidx.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <stdio.h>
#include "tags.h"
#include "msa.h"
#include "split.h"


namespace torali {


  template<typename TAlign, typename TQualities>
  inline uint32_t
  _getAlignmentQual(TAlign const& align, TQualities const& qual) {
    typedef typename TAlign::index TAIndex;
    uint32_t baseQualSum = 0;
    uint32_t seqPtr = 0;
    uint32_t alignedBases = 0;
    for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
      if (align[1][j] != '-') {
	if (align[0][j] != '-') {
	  ++alignedBases;
	  baseQualSum += qual[seqPtr];
	}
	++seqPtr;
      }
    }
    return (baseQualSum / alignedBases);
  }

  template<typename TPos, typename TTag>
  inline int32_t
  _cutRefStart(TPos const rStart, TPos const rEnd, TPos const offset, unsigned int bpPoint, uint8_t const, SVType<TTag>) {
    if (bpPoint) return rEnd - offset;
    else return rStart - offset;
  }

  template<typename TPos>
  inline int32_t
  _cutRefStart(TPos const rStart, TPos const rEnd, TPos const offset, unsigned int bpPoint, uint8_t const, SVType<DuplicationTag>) {
    if (!bpPoint) return rEnd - offset;
    else return rStart - offset;
  }

  template<typename TPos>
  inline int32_t
  _cutRefStart(TPos const rStart, TPos const rEnd, TPos const offset, unsigned int bpPoint, uint8_t const ct, SVType<TranslocationTag>) {
    if (ct == 3) {
      if (!bpPoint) return rEnd - offset;
      else return rStart - offset;
    } else {
      if (bpPoint) return rEnd - offset;
      else return rStart - offset;
    }
  }

  template<typename TPos, typename TTag>
  inline int32_t
  _cutRefEnd(TPos const rStart, TPos const rEnd, TPos const offset, unsigned int bpPoint, uint8_t const, SVType<TTag>) {
    if (bpPoint) return rEnd + offset;
    else return rStart + offset;
  }

  template<typename TPos>
  inline int32_t
  _cutRefEnd(TPos const rStart, TPos const rEnd, TPos const offset, unsigned int bpPoint, uint8_t const, SVType<DuplicationTag>) {
    if (!bpPoint) return rEnd + offset;
    else return rStart + offset;
  }

  template<typename TPos>
  inline int32_t
  _cutRefEnd(TPos const rStart, TPos const rEnd, TPos const offset, unsigned int bpPoint, uint8_t const ct, SVType<TranslocationTag>) {
    if (ct == 3) {
      if (!bpPoint) return rEnd + offset;
      else return rStart + offset;
    } else {
      if (bpPoint) return rEnd + offset;
      else return rStart + offset;
    }
  }

  template<typename TConfig, typename TSVs, typename TCountMap, typename TTag>
  inline void
  annotateJunctionReads(TConfig const& c, TSVs& svs, TCountMap& countMap, SVType<TTag> svType)
  {
    typedef typename TCountMap::value_type::value_type TCountPair;
    typedef std::vector<uint8_t> TQuality;

    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    for(std::size_t file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
    }      
    bam_hdr_t* hdr = sam_hdr_read(samfile[0]);

    // Sort Structural Variants
    sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());
    
    // Initialize count map
    countMap.resize(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) countMap[file_c].resize(svs.size(), TCountPair());

    // Process chromosome by chromosome
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Junction read annotation" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );
    
    // Get reference probes
    typedef boost::unordered_map<unsigned int, std::string> TProbes;
    TProbes refProbes;

    // Parse genome
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      ++show_progress;
      std::string tname(hdr->target_name[refIndex]);
      int32_t seqlen = -1;
      char* seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
      // Iterate all structural variants
      for(typename TSVs::iterator itSV = svs.begin(); itSV != svs.end(); ++itSV) {
	if (!itSV->precise) continue;
	int consLen = itSV->consensus.size();

	// Create a pseudo structural variant record
	StructuralVariantRecord svRec;
	svRec.chr = itSV->chr;
	svRec.chr2 = itSV->chr2;
	svRec.svStartBeg = std::max(itSV->svStart - consLen, 0);
	svRec.svStart = itSV->svStart;
	svRec.svStartEnd = std::min((uint32_t) itSV->svStart + consLen, hdr->target_len[itSV->chr]);
	svRec.svEndBeg = std::max(itSV->svEnd - consLen, 0);
	svRec.svEnd = itSV->svEnd;
	svRec.svEndEnd = std::min((uint32_t) itSV->svEnd + consLen, hdr->target_len[itSV->chr2]);
	svRec.ct = itSV->ct;
	if ((itSV->chr != itSV->chr2) && (itSV->chr2 == refIndex)) {
	  refProbes[itSV->id] = _getSVRef(c, seq, svRec, refIndex, svType);
	}
	if (itSV->chr == refIndex) {
	  // Get the reference string
	  if (itSV->chr != itSV->chr2) svRec.consensus=refProbes[itSV->id];
	  std::string svRefStr = _getSVRef(c, seq, svRec, refIndex, svType);
	  
	  // Find breakpoint to reference
	  typedef boost::multi_array<char, 2> TAlign;
	  TAlign align;
	  if (!_consRefAlignment(itSV->consensus, svRefStr, align, svType)) continue;
	  AlignDescriptor ad;
	  if (!_findSplit(c, itSV->consensus, svRefStr, align, ad, svType)) continue;
	  
	  // Debug consensus to reference alignment
	  //for(TAIndex i = 0; i<align.shape()[0]; ++i) {
	  //for(TAIndex j = 0; j<align.shape()[1]; ++j) std::cerr << align[i][j];
	  //std::cerr << std::endl;
	  //}
	  //std::cerr << std::endl;
	  
	  
	  // Iterate all samples
#pragma omp parallel for default(shared)
	  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	    unsigned int refAlignedReadCount = 0;
	    
	    for (unsigned int bpPoint = 0; bpPoint<2; ++bpPoint) {
	      int32_t regionChr, regionStart, regionEnd, cutConsStart, cutConsEnd, cutRefStart, cutRefEnd;
	      if (bpPoint) {
		regionChr = itSV->chr2;
		regionStart = std::max(0, itSV->svEnd - c.minimumFlankSize);
		regionEnd = std::min((uint32_t) (itSV->svEnd + c.minimumFlankSize), hdr->target_len[itSV->chr2]);
		cutConsStart = ad.cEnd - ad.homLeft - c.minimumFlankSize;
		cutConsEnd = ad.cEnd + ad.homRight + c.minimumFlankSize;
		cutRefStart = _cutRefStart(ad.rStart, ad.rEnd, ad.homLeft + c.minimumFlankSize, bpPoint, itSV->ct, svType);
		cutRefEnd = _cutRefEnd(ad.rStart, ad.rEnd, ad.homRight + c.minimumFlankSize, bpPoint, itSV->ct, svType);
	      } else {
		regionChr = itSV->chr;
		regionStart = std::max(0, itSV->svStart - c.minimumFlankSize);
		regionEnd = std::min((uint32_t) (itSV->svStart + c.minimumFlankSize), hdr->target_len[itSV->chr]);
		cutConsStart = ad.cStart - ad.homLeft - c.minimumFlankSize;
		cutConsEnd = ad.cStart + ad.homRight + c.minimumFlankSize;
		cutRefStart = _cutRefStart(ad.rStart, ad.rEnd, ad.homLeft + c.minimumFlankSize, bpPoint, itSV->ct, svType);
		cutRefEnd = _cutRefEnd(ad.rStart, ad.rEnd, ad.homRight + c.minimumFlankSize, bpPoint, itSV->ct, svType);
	      }
	      std::string consProbe = itSV->consensus.substr(cutConsStart, (cutConsEnd - cutConsStart));
	      std::string refProbe = svRefStr.substr(cutRefStart, (cutRefEnd - cutRefStart));
	      
	      //std::cerr << "Junction probes" << std::endl;
	      //std::cerr << bpPoint << "," << consProbe << "," << refProbe << std::endl;
	      
	      
	      hts_itr_t* iter = sam_itr_queryi(idx[file_c], regionChr, regionStart, regionEnd);
	      bam1_t* rec = bam_init1();
	      while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
		if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
		// Check read length & quality
		if ((rec->core.l_qseq < (2 * c.minimumFlankSize)) || (rec->core.qual < c.minGenoQual)) continue;
		
		// Valid soft clip or no soft-clip read?
		bool hasSoftClip = false;
		uint32_t* cigar = bam_get_cigar(rec);
		for (std::size_t i = 0; i < rec->core.n_cigar; ++i)
		  if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) hasSoftClip = true;
		if (!hasSoftClip) {
		  if (bpPoint) {
		    if ((rec->core.pos + c.minimumFlankSize + ad.homLeft > itSV->svEnd) || (rec->core.pos + rec->core.l_qseq < itSV->svEnd + c.minimumFlankSize + ad.homRight)) continue;
		  } else {
		    if ((rec->core.pos + c.minimumFlankSize + ad.homLeft > itSV->svStart) || (rec->core.pos + rec->core.l_qseq < itSV->svStart + c.minimumFlankSize + ad.homRight)) continue;
		  }
		}
		
		// Get sequence
		std::string sequence;
		sequence.resize(rec->core.l_qseq);
		uint8_t* seqptr = bam_get_seq(rec);
		for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
		_adjustOrientation(sequence, bpPoint, itSV->ct, svType);
		
		// Compute alignment to alternative haplotype
		TAlign alignAlt;
		DnaScore<int> simple(5, -4, -4, -4);
		AlignConfig<true, false> semiglobal;
		int32_t scoreA = needle(consProbe, sequence, alignAlt, semiglobal, simple);
		int32_t scoreAltThreshold = (int32_t) (c.flankQuality * consProbe.size() * simple.match + (1.0 - c.flankQuality) * consProbe.size() * simple.mismatch);
		double scoreAlt = (double) scoreA / (double) scoreAltThreshold;

		// Compute alignment to reference haplotype
		TAlign alignRef;
		int32_t scoreR = needle(refProbe, sequence, alignRef, semiglobal, simple);
		int32_t scoreRefThreshold = (int32_t) (c.flankQuality * refProbe.size() * simple.match + (1.0 - c.flankQuality) * refProbe.size() * simple.mismatch);
		double scoreRef = (double) scoreR / (double) scoreRefThreshold;
		
		// Any confident alignment?
		if ((scoreRef > 1) || (scoreAlt > 1)) {
		  // Debug alignment to REF and ALT
		  //std::cerr << "Alt:\t" << scoreAlt << "\tRef:\t" << scoreRef << std::endl;
		  //for(TAIndex i = 0; i< (TAIndex) alignAlt.shape()[0]; ++i) {
		  //for(TAIndex j = 0; j< (TAIndex) alignAlt.shape()[1]; ++j) std::cerr << alignAlt[i][j];
		  //std::cerr << std::endl;
		  //}
		  //for(TAIndex i = 0; i< (TAIndex) alignRef.shape()[0]; ++i) {
		  //for(TAIndex j = 0; j< (TAIndex) alignRef.shape()[1]; ++j) std::cerr << alignRef[i][j];
		  //std::cerr << std::endl;
		  //}
		  
		  if (scoreRef > scoreAlt) {
		    // Take only every second read because we sample both breakpoints
		    if (++refAlignedReadCount % 2) {
		      TQuality quality;
		      quality.resize(rec->core.l_qseq);
		      uint8_t* qualptr = bam_get_qual(rec);
		      for (int i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
		      uint32_t rq = _getAlignmentQual(alignRef, quality);
		      if (rq >= c.minGenoQual) {
#pragma omp critical
			{
			  countMap[file_c][itSV->id].first.push_back((uint8_t) std::min(rq, (uint32_t) rec->core.qual));
			}
		      }
		    }
		  } else {
		    TQuality quality;
		    quality.resize(rec->core.l_qseq);
		    uint8_t* qualptr = bam_get_qual(rec);
		    for (int i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
		    uint32_t aq = _getAlignmentQual(alignAlt, quality);
		    if (aq >= c.minGenoQual) {
#pragma omp critical
		      {
			countMap[file_c][itSV->id].second.push_back((uint8_t) std::min(aq, (uint32_t) rec->core.qual));
		      }
		    }
		  }
		}
	      }
	      bam_destroy1(rec);
	      hts_itr_destroy(iter);
	    }
	  }
	}
      }
      if (seqlen) free(seq);
    }

    // Clean-up
    fai_destroy(fai);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
    bam_hdr_destroy(hdr);
  }



}

#endif


