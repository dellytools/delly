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

  template<typename TConfig, typename TSampleLibrary, typename TSVs, typename TCountMap, typename TTag>
  inline void
  annotateJunctionReads(TConfig const& c, TSampleLibrary& sampleLib, TSVs& svs, TCountMap& countMap, SVType<TTag> svType)
  {
    typedef typename TCountMap::key_type TSampleSVPair;
    typedef typename TCountMap::mapped_type TCountPair;
    typedef typename TCountPair::first_type TQualVector;
    typedef std::vector<uint8_t> TQuality;

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
	      TAIndex cStart, cEnd, rStart, rEnd, gS, gE;
	      double percId = 0;
	      _findSplit(alignFwd, cStart, cEnd, rStart, rEnd, gS, gE, percId, svType);
	      int32_t homLeft = 0;
	      int32_t homRight = 0;
	      _findHomology(alignFwd, gS, gE, homLeft, homRight);
	      if ((homLeft + c.minimumFlankSize > (int32_t) cStart) || ( (int32_t) (itSV->consensus.size() - cEnd) < homRight + c.minimumFlankSize)) continue;
	      if ((homLeft + c.minimumFlankSize > (int32_t) rStart) || ( (int32_t) (svRefStr.size() - rEnd) < homRight + c.minimumFlankSize)) continue;


	      // Debug consensus to reference alignment
	      //for(TAIndex i = 0; i<alignFwd.shape()[0]; ++i) {
	      //for(TAIndex j = 0; j<alignFwd.shape()[1]; ++j) std::cerr << alignFwd[i][j];
	      //std::cerr << std::endl;
	      //}
	      //std::cerr << "Consensus-to-reference: " << homLeft << ',' << homRight << std::endl;
	      //std::cerr << cStart << ',' << cEnd << ',' << rStart << ',' << rEnd << std::endl;
	      //std::cerr << homLeft << ',' << homRight << std::endl;


	      // Iterate all samples
#pragma omp parallel for default(shared)
	      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
		std::string sampleName(c.files[file_c].stem().string());
		TQualVector altQual;
		TQualVector refQual;
		unsigned int refAlignedReadCount = 0;

		for (unsigned int bpPoint = 0; bpPoint<2; ++bpPoint) {
		  int32_t regionChr, regionStart, regionEnd, cutConsStart, cutConsEnd, cutRefStart, cutRefEnd;
		  if (bpPoint) {
		    regionChr = itSV->chr2;
		    regionStart = std::max(0, itSV->svEnd - c.minimumFlankSize);
		    regionEnd = std::min((uint32_t) (itSV->svEnd + c.minimumFlankSize), reflen[itSV->chr2]);
		    cutConsStart = cEnd - homLeft - c.minimumFlankSize;
		    cutConsEnd = cEnd + homRight + c.minimumFlankSize;
		    cutRefStart = _cutRefStart((int32_t) rStart, (int32_t) rEnd, homLeft + c.minimumFlankSize, bpPoint, itSV->ct, svType);
		    cutRefEnd = _cutRefEnd((int32_t) rStart, (int32_t) rEnd, homRight + c.minimumFlankSize, bpPoint, itSV->ct, svType);
		  } else {
		    regionChr = itSV->chr;
		    regionStart = std::max(0, itSV->svStart - c.minimumFlankSize);
		    regionEnd = std::min((uint32_t) (itSV->svStart + c.minimumFlankSize), reflen[itSV->chr]);
		    cutConsStart = cStart - homLeft - c.minimumFlankSize;
		    cutConsEnd = cStart + homRight + c.minimumFlankSize;
		    cutRefStart = _cutRefStart((int32_t) rStart, (int32_t) rEnd, homLeft + c.minimumFlankSize, bpPoint, itSV->ct, svType);
		    cutRefEnd = _cutRefEnd((int32_t) rStart, (int32_t) rEnd, homRight + c.minimumFlankSize, bpPoint, itSV->ct, svType);
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
		    if ((rec->core.l_qseq < 35) || (rec->core.qual < c.minGenoQual)) continue;

		    // Valid soft clip or no soft-clip read?
		    bool hasSoftClip = false;
		    uint32_t* cigar = bam_get_cigar(rec);
		    for (std::size_t i = 0; i < rec->core.n_cigar; ++i)
		      if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) hasSoftClip = true;
		    if (!hasSoftClip) {
		      if (bpPoint) {
			if ((rec->core.pos + c.minimumFlankSize + homLeft > itSV->svEnd) || (rec->core.pos + rec->core.l_qseq < itSV->svEnd + c.minimumFlankSize + homRight)) continue;
		      } else {
			if ((rec->core.pos + c.minimumFlankSize + homLeft > itSV->svStart) || (rec->core.pos + rec->core.l_qseq < itSV->svStart + c.minimumFlankSize + homRight)) continue;
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
		    int32_t scoreA = needle(consProbe, sequence, alignAlt, semiglobal, simple);
		    int32_t scoreAltThreshold = (int32_t) (c.flankQuality * consProbe.size() * 5 + (1.0 - c.flankQuality) * consProbe.size() * (-4));
		    double scoreAlt = (double) scoreA / (double) scoreAltThreshold;

		    // Compute alignment to reference haplotype
		    TAlign alignRef;
		    int32_t scoreR = needle(refProbe, sequence, alignRef, semiglobal, simple);
		    int32_t scoreRefThreshold = (int32_t) (c.flankQuality * refProbe.size() * 5 + (1.0 - c.flankQuality) * refProbe.size() * (-4));
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
			  if (rq >= c.minGenoQual) refQual.push_back((uint8_t) std::min(rq, (uint32_t) rec->core.qual));
			}
		      } else {
			TQuality quality;
			quality.resize(rec->core.l_qseq);
			uint8_t* qualptr = bam_get_qual(rec);
			for (int i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
			uint32_t aq = _getAlignmentQual(alignAlt, quality);
			if (aq >= c.minGenoQual) altQual.push_back((uint8_t) std::min(aq, (uint32_t) rec->core.qual));
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


