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

#ifndef PACBIOJUNCTION_H
#define PACBIOJUNCTION_H


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
#include "junction.h"

namespace torali {


  template<typename TAlign>
  inline double
  _percentIdentity(TAlign const& align) {
    typedef typename TAlign::index TAIndex;
    TAIndex rI=0;
    TAIndex vI=0;
    uint32_t gapMM = 0;
    uint32_t mm = 0;
    uint32_t ma = 0;
    bool inGap=false;
    for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
      //std::cerr << j << ',' << ma << ',' << mm << std::endl;
      if (align[0][j] != '-') ++vI;
      if (align[1][j] != '-') ++rI;
      // Internal gap?
      if ((align[0][j] == '-') || (align[1][j] == '-')) {
	if ((rI>0) && (vI>0)) {
	  if (!inGap) {
	    inGap = true;
	    gapMM = 0;
	  }
	  gapMM += 1;
	}
      } else {
	if (inGap) {
	  mm += gapMM;
	  inGap=false;
	}
	if (align[0][j] == align[1][j]) ma += 1;
	else mm += 1;
      }
    }
    return (double) ma / (double) (ma + mm);
  }



  template<typename TConfig, typename TSVs, typename TCountMap, typename TTag>
  inline void
  annotatePacbioJunctionReads(TConfig const& c, TSVs& svs, TCountMap& countMap, SVType<TTag> svType)
  {
    typedef typename TCountMap::key_type TSampleSVPair;
    typedef typename TCountMap::mapped_type TCountPair;
    typedef typename TCountPair::first_type TQualVector;
    typedef std::vector<uint8_t> TQuality;

    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    TSamFile samfile;
    TIndex idx;
    samfile.resize(c.files.size());
    idx.resize(c.files.size());
    for(std::size_t file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile[0]);

    // Sort Structural Variants
    sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());
    
    // Initialize count map
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c)
      for(typename TSVs::const_iterator itSV = svs.begin(); itSV!=svs.end(); ++itSV) countMap.insert(std::make_pair(std::make_pair(file_c, itSV->id), TCountPair()));

    // Process chromosome by chromosome
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Junction read annotation" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );
    
    // Parse genome
    kseq_t *seq;
    int l;
    gzFile fp = gzopen(c.genome.string().c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
      // Find reference index
      for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
	if (std::string(seq->name.s) == std::string(hdr->target_name[refIndex])) {
	  ++show_progress;
	  
	  // Iterate all structural variants
	  for(typename TSVs::iterator itSV = svs.begin(); itSV != svs.end(); ++itSV) {
	    if (!itSV->precise) continue;
	    int32_t consLen = itSV->consensus.size();

	    if (itSV->chr == refIndex) {
	      // Get the reference string
	      std::string svRefStr = boost::to_upper_copy(std::string(seq->seq.s + itSV->svStart - consLen, seq->seq.s + itSV->svEnd + consLen));

	      // Find breakpoint to reference
	      typedef boost::multi_array<char, 2> TAlign;
	      typedef typename TAlign::index TAIndex;
	      TAlign alignFwd;
	      _consRefAlignment(itSV->consensus, svRefStr, alignFwd, svType);
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
		TQualVector altQual;
		TQualVector refQual;

		// Region
		int32_t regionChr = itSV->chr;
		int32_t regionStart = std::max(0, itSV->svStart - c.minimumFlankSize);
		int32_t regionEnd = std::min((uint32_t) (itSV->svEnd + c.minimumFlankSize), hdr->target_len[itSV->chr2]);
		int32_t cutConsStart = cStart - homLeft - c.minimumFlankSize;
		int32_t cutConsEnd = cEnd + homRight + c.minimumFlankSize;
		std::string consProbe = itSV->consensus.substr(cutConsStart, (cutConsEnd - cutConsStart));
		int32_t cutRefStart = rStart - homLeft - c.minimumFlankSize;
		int32_t cutRefEnd = rStart + homRight + c.minimumFlankSize;
		std::string refProbe = svRefStr.substr(cutRefStart, (cutRefEnd - cutRefStart));

		hts_itr_t* iter = sam_itr_queryi(idx[file_c], regionChr, regionStart, regionEnd);
		bam1_t* rec = bam_init1();
		while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
		  if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
		  // Check read length & quality
		  if ((rec->core.l_qseq < 35) || (rec->core.qual < c.minGenoQual)) continue;

		  // Get sequence
		  std::string sequence;
		  sequence.resize(rec->core.l_qseq);
		  uint8_t* seqptr = bam_get_seq(rec);
		  for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

		  // Compute alignment to alternative haplotype
		  TAlign alignAlt;
		  AlignConfig<true, false> semiglobal;
		  DnaScore<int> simple(5, -4, -4, -4);
		  int32_t scoreA = needle(consProbe, sequence, alignAlt, semiglobal, simple);
		  int32_t scoreAltThreshold = (int32_t) (c.flankQuality * consProbe.size() * simple.match + (1.0 - c.flankQuality) * consProbe.size() * simple.mismatch);
		  double scoreAlt = (double) scoreA / (double) scoreAltThreshold;

		  // Compute alignment to reference haplotype
		  TAlign alignRef;
		  int32_t scoreR = needle(refProbe, sequence, alignRef, semiglobal, simple);
		  int32_t scoreRefThreshold = (int32_t) (c.flankQuality * refProbe.size() * simple.match + (1.0 - c.flankQuality) * refProbe.size() * simple.mismatch);
		  double scoreRef = (double) scoreR / (double) scoreRefThreshold;
		    
		  // Any confident alignment?
		  //if ((scoreRef > c.flankQuality) || (scoreAlt > c.flankQuality)) {
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
		      TQuality quality;
		      quality.resize(rec->core.l_qseq);
		      uint8_t* qualptr = bam_get_qual(rec);
		      for (int i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
		      uint32_t rq = _getAlignmentQual(alignRef, quality);
		      if (rq >= c.minGenoQual) refQual.push_back((uint8_t) std::min(rq, (uint32_t) rec->core.qual));
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
		
		// Insert counts
		TSampleSVPair svSample = std::make_pair(file_c, itSV->id);
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
    bam_hdr_destroy(hdr);
  }



}

#endif


