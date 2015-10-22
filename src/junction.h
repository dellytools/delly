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
#include <htslib/kseq.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <stdio.h>
#include "tags.h"
#include "msa.h"
#include "split.h"

KSEQ_INIT(gzFile, gzread)

namespace torali {

  template<typename TBP, typename TPos, typename TSize, typename TSVBp, typename TTag>
  inline bool
  _overlapsSVBreakpoint(TBP const bpPoint, TPos const pStart, TSize const pEnd, TSVBp const start, TSVBp const end, SVType<TTag>) 
  {
    if (bpPoint) {
      if ((pStart > end) || (pEnd < end)) return false;
    } else {
      if ((pStart > start) || (pEnd < start)) return false;
    }
    return true;
  }

  template<typename TBP, typename TPos, typename TSize, typename TSVBp>
  inline bool
  _overlapsSVBreakpoint(TBP const bpPoint, TPos const pStart, TSize const pEnd, TSVBp const start, TSVBp const end, SVType<DuplicationTag>) 
  {
    if (!bpPoint) {
      if ((pStart > end) || (pEnd < end)) return false;
    } else {
      if ((pStart > start) || (pEnd < start)) return false;
    }
    return true;
  }


  template<typename TAlign, typename TAIndex>
  inline int
  _coreAlignScore(TAlign const& align, TAIndex& consStart, TAIndex& consEnd, TAIndex& alignLength) {
    // Ignore leading and trailing gaps
    bool leadingGap = true;
    int score = 0;
    int savedScore = 0;
    consStart = 0;
    TAIndex runningEnd = 0;
    TAIndex consGaps = 0;
    TAIndex totalConsGaps = 0;
    TAIndex seqGaps = 0;
    for(TAIndex j = 0; j< (TAIndex) align.shape()[1]; ++j) {
      if (leadingGap) {
	if (align[1][j] != '-') {
	  ++consStart;
	  ++runningEnd;
	}
	if ((align[0][j] == '-') || (align[1][j] == '-')) continue;
	else leadingGap = false;
      }
      if ((align[0][j] == '-') && (align[1][j] != '-')) {
	++seqGaps;
	++runningEnd;
      } else if ((align[0][j] != '-') && (align[1][j] == '-')) {
	++consGaps;
      } else {
	if (align[0][j] == align[1][j]) score += 5;
	else score += -4;
	if (consGaps) {
	  score += (-10 + (-1) * consGaps);
	  totalConsGaps += consGaps;
	  consGaps = 0;
	}
	if (seqGaps) {
	  score+= (-10 + (-1) * seqGaps);
	  seqGaps = 0;
	}
	++runningEnd;
	savedScore = score;
	consEnd = runningEnd;
      }
    }
    alignLength = (consEnd - consStart) + totalConsGaps;
    return savedScore;
  }

  template<typename TAlign>
  inline int
  _coreAlignScore(TAlign const& align) {
    typedef typename TAlign::index TAIndex;
    TAIndex consStart, consEnd, alignLength;
    return _coreAlignScore(align, consStart, consEnd, alignLength);
  }

  template<typename TConfig, typename TSampleLibrary, typename TSVs, typename TCountMap, typename TTag>
  inline void
  annotateJunctionReads(TConfig const& c, TSampleLibrary& sampleLib, TSVs& svs, TCountMap& countMap, SVType<TTag> svType)
  {
    typedef typename TCountMap::key_type TSampleSVPair;
    typedef typename TCountMap::mapped_type TCountPair;
    typedef typename TCountPair::first_type TQualVector;
    double qualityThres = (double) (c.flankQuality) / 100.0;

    // Open file handles
    typedef std::vector<std::string> TRefNames;
    typedef std::vector<uint32_t> TRefLength;
    TRefNames refnames;
    TRefLength reflen;
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    TSamFile samfile;
    TIndex idx;
    samfile.resize(c.files.size());
    idx.resize(c.files.size());
    for(std::size_t file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
      if (!file_c) {
	bam_hdr_t* hdr = sam_hdr_read(samfile[file_c]);
	for (int i = 0; i<hdr->n_targets; ++i) {
	  refnames.push_back(hdr->target_name[i]);
	  reflen.push_back(hdr->target_len[i]);
	}
	bam_hdr_destroy(hdr);
      }
    }

    // Sort Structural Variants
    sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());
    
    // Initialize count map
    for(typename TSampleLibrary::iterator sIt = sampleLib.begin(); sIt!=sampleLib.end();++sIt) 
      for(typename TSVs::const_iterator itSV = svs.begin(); itSV!=svs.end(); ++itSV) countMap.insert(std::make_pair(std::make_pair(sIt->first, itSV->id), TCountPair()));
    
    // Process chromosome by chromosome
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Junction read annotation" << std::endl;
    boost::progress_display show_progress( refnames.size() );
    
    // Get reference probes
    typedef boost::unordered_map<unsigned int, std::string> TProbes;
    TProbes refProbes;

    // Parse genome
    kseq_t *seq;
    int l;
    gzFile fp = gzopen(c.genome.string().c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
      // Find reference index
      for(int32_t refIndex=0; refIndex < (int32_t) refnames.size(); ++refIndex) {
	if (seq->name.s == refnames[refIndex]) {
	  ++show_progress;
	  
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
	    svRec.svStartEnd = std::min((uint32_t) itSV->svStart + consLen, reflen[itSV->chr]);
	    svRec.svEndBeg = std::max(itSV->svEnd - consLen, 0);
	    svRec.svEnd = itSV->svEnd;
	    svRec.svEndEnd = std::min((uint32_t) itSV->svEnd + consLen, reflen[itSV->chr2]);
	    svRec.ct = itSV->ct;
	    if ((itSV->chr != itSV->chr2) && (itSV->chr2 == refIndex)) {
	      refProbes[itSV->id] = _getSVRef(seq->seq.s, svRec, refIndex, svType);
	    }
	    if (itSV->chr == refIndex) {
	      // Get the reference string
	      if (itSV->chr != itSV->chr2) svRec.consensus=refProbes[itSV->id];
	      std::string svRefStr = _getSVRef(seq->seq.s, svRec, refIndex, svType);

	      // Find breakpoint to reference
	      typedef boost::multi_array<char, 2> TAlign;
	      typedef typename TAlign::index TAIndex;
	      TAlign alignFwd;
	      AlignConfig<true, false> semiglobal;
	      DnaScore<int> sc(5, -4, -5 * c.minimumFlankSize, 0);
	      gotoh(itSV->consensus, svRefStr, alignFwd, semiglobal, sc);
	      TAIndex cStart, cEnd, rStart, rEnd;
	      _findSplit(alignFwd, cStart, cEnd, rStart, rEnd);
	      
	      // Debug consensus to reference alignment
	      //for(TAIndex i = 0; i<alignFwd.shape()[0]; ++i) {
	      //for(TAIndex j = 0; j<alignFwd.shape()[1]; ++j) std::cerr << alignFwd[i][j];
	      //std::cerr << std::endl;
	      //}
	      //std::cerr << cStart << ',' << cEnd << ':' << rStart << ',' << rEnd << std::endl;

	      // Iterate all samples
#pragma omp parallel for default(shared)
	      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
		std::string sampleName(c.files[file_c].stem().string());
		TQualVector altQual;
		TQualVector refQual;
		unsigned int refAlignedReadCount = 0;

		for (unsigned int bpPoint = 0; bpPoint<2; ++bpPoint) {
		  int32_t regionChr = itSV->chr;
		  int regionStart = itSV->svStartBeg;
		  int regionEnd = std::min((uint32_t) (itSV->svStart + itSV->svStartEnd)/2, reflen[itSV->chr]);
		  int csBp = cStart;
		  if (bpPoint) {
		    csBp = cEnd;
		    regionChr = itSV->chr2;
		    regionStart = (itSV->svEndBeg + itSV->svEnd)/2;
		    regionEnd = std::min((uint32_t) itSV->svEndEnd, reflen[itSV->chr2]);
		  }

		  hts_itr_t* iter = sam_itr_queryi(idx[file_c], regionChr, regionStart, regionEnd);
		  bam1_t* rec = bam_init1();
		  while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
		    if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
		    // Check read length & quality
		    if ((rec->core.l_qseq < 35) || (rec->core.qual < c.minGenoQual)) continue;

		    // Valid soft clip or no soft-clip read?
		    bool hasSoftClip = false;
		    uint32_t* cigar = bam_get_cigar(rec);
		    for (std::size_t i = 0; i < rec->core.n_cigar; ++i)
		      if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) hasSoftClip = true;
		    if (hasSoftClip) {
		      int clipSize = 0;
		      int splitPoint = 0;
		      bool leadingSoftClip = false;
		      if (!_validSoftClip(rec, clipSize, splitPoint, leadingSoftClip, c.minMapQual)) continue;
		      if ((splitPoint < regionStart) || (splitPoint > regionEnd)) continue;
		      if (!_validSCOrientation(bpPoint, leadingSoftClip, itSV->ct, svType)) continue;
		    } else {
		      // Check position
		      if (bpPoint) {
			if ((rec->core.pos + c.minimumFlankSize > itSV->svEnd) || (rec->core.pos + rec->core.l_qseq - c.minimumFlankSize < itSV->svEnd)) continue;
		      } else {
			if ((rec->core.pos + c.minimumFlankSize > itSV->svStart) || (rec->core.pos + rec->core.l_qseq - c.minimumFlankSize < itSV->svStart)) continue;
		      }
		    }

		    // Get sequence
		    std::string sequence;
		    sequence.resize(rec->core.l_qseq);
		    uint8_t* seqptr = bam_get_seq(rec);
		    for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
		    _adjustOrientation(sequence, bpPoint, itSV->ct, svType);

		    // Compute alignment to alternative haplotype
		    AlignConfig<true, true> endFreeAlign;
		    TAlign align;
		    gotoh(sequence, itSV->consensus, align, endFreeAlign);
		    TAIndex consStart = 0;
		    TAIndex consEnd = 0;
		    TAIndex consAlignLength = 0;
		    int altScore = _coreAlignScore(align, consStart, consEnd, consAlignLength);
		    int scoreThresholdAlt = (int) (qualityThres * consAlignLength * 5 + (1.0 - qualityThres) * consAlignLength * (-4));

		    // Debug alignment ALT
		    //std::cerr << "Alt: " << altScore << std::endl;
		    //for(TAIndex i = 0; i< (TAIndex) align.shape()[0]; ++i) {
		    //for(TAIndex j = 0; j< (TAIndex) align.shape()[1]; ++j) std::cerr << align[i][j];
		    //std::cerr << std::endl;
		    //}
		    
		    // Compute alignment to reference haplotype
		    gotoh(sequence, svRefStr, align, endFreeAlign);
		    TAIndex rCoreStart = 0;
		    TAIndex rCoreEnd = 0;
		    TAIndex rAlignLength = 0;
		    int refScore = _coreAlignScore(align, rCoreStart, rCoreEnd, rAlignLength);
		    int scoreThresholdRef = (int) (qualityThres * rAlignLength * 5 + (1.0 - qualityThres) * rAlignLength * (-4));

		    // Debug alignment REF
		    //std::cerr << "Ref: " << refScore << std::endl;
		    //for(TAIndex i = 0; i< (TAIndex) align.shape()[0]; ++i) {
		    //for(TAIndex j = 0; j< (TAIndex) align.shape()[1]; ++j) std::cerr << align[i][j];
		    //std::cerr << std::endl;
		    //}
			  
		    // Any confident alignment?
		    if ((refScore > scoreThresholdRef) || (altScore > scoreThresholdAlt)) {
		      if ( (double) refScore / (double) scoreThresholdRef > (double) altScore / (double) scoreThresholdAlt) {
			if (rCoreEnd - rCoreStart < 35) continue;
			if (!_overlapsSVBreakpoint(bpPoint, rCoreStart + c.minimumFlankSize, rCoreEnd - c.minimumFlankSize, rStart, rEnd, svType)) continue;
			
			//std::cerr << sampleName << ',' << altScore << ',' << refScore << ':' << scoreThresholdAlt << ',' << scoreThresholdRef << std::endl;
			// Take only every second read because we sample both breakpoints
			if (++refAlignedReadCount % 2) refQual.push_back(rec->core.qual);
		      } else {
			if (consEnd - consStart < 35) continue;
			if ( ( (consStart + c.minimumFlankSize) > csBp) || ((consEnd - c.minimumFlankSize) < csBp) ) continue;

			//std::cerr << sampleName << ',' << altScore << ',' << refScore << ':' << scoreThresholdAlt << ',' << scoreThresholdRef << std::endl;
			altQual.push_back(rec->core.qual);
		      }
		    }
		  }
		  bam_destroy1(rec);
		  hts_itr_destroy(iter);
		}

		// Insert counts
		TSampleSVPair svSample = std::make_pair(sampleName, itSV->id);
#pragma omp critical
		{
		  typename TCountMap::iterator countMapIt=countMap.find(svSample);
		  countMapIt->second.first=refQual;
		  countMapIt->second.second=altQual;
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
  }



}

#endif


