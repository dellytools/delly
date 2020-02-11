#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include <boost/progress.hpp>

#include <iostream>
#include "msa.h"
#include "split.h"
#include "gotoh.h"
#include "needle.h"

namespace torali
{
  struct SeqSlice {
    int32_t svid;
    int32_t sstart;
    int32_t inslen;
    int32_t qual;  // Only required for junction count map

    SeqSlice(int32_t const sv, int32_t const sst, int32_t const il, int32_t q) : svid(sv), sstart(sst), inslen(il), qual(q) {}
  };


  // Insertions & deletions only
  template<typename TConfig>
  inline bool
  teguaAlignConsensus(TConfig const& c, bam_hdr_t* hdr, char const* seq, StructuralVariantRecord& sv) {

    // Consensus long enough?
    int32_t seqstart = 0;
    int32_t seqend = 0;
    if (sv.svt == 2) {
      if ((int32_t) sv.consensus.size() < (2 * c.minimumFlankSize)) return false;
      seqstart = std::max(sv.svStart - (int32_t) (sv.consensus.size() - c.minimumFlankSize), 0);
      seqend =  std::min(sv.svEnd + (int32_t) (sv.consensus.size() - c.minimumFlankSize), (int32_t) hdr->target_len[sv.chr]);
    } else if (sv.svt == 4) {
      if ((int32_t) sv.consensus.size() < (2 * c.minimumFlankSize + sv.insLen)) return false;
      int32_t bufferSpace = std::max((int32_t) ((sv.consensus.size() - sv.insLen) / 3), c.minimumFlankSize);
      seqstart = std::max(sv.svStart - bufferSpace, 0);
      seqend =  std::min(sv.svEnd + bufferSpace, (int32_t) hdr->target_len[sv.chr]);
    }
    std::string svRefStr = boost::to_upper_copy(std::string(seq + seqstart, seq + seqend));

    // Consensus to reference alignment
    typedef boost::multi_array<char, 2> TAlign;
    TAlign align;
    AlignConfig<true, false> global;
    DnaScore<int> lnsc(5, -4, -4, -4);
    bool reNeedle = false;
    if (sv.svt == 2) {
      reNeedle = longNeedle(sv.consensus, svRefStr, align, global, lnsc);
    } else if (sv.svt == 4) {
      reNeedle = longNeedle(svRefStr, sv.consensus, align, global, lnsc);
      for(uint32_t j = 0; j < align.shape()[1]; ++j) {
	char tmp = align[0][j];
	align[0][j] = align[1][j];
	align[1][j] = tmp;
      }	
    }
    if (!reNeedle) return false;

    // Debug consensus to reference alignment
    std::cerr << "Consensus-to-Reference alignment (full-length reference)" << std::endl;
    for(uint32_t i = 0; i < align.shape()[0]; ++i) {
      for(uint32_t j = 0; j< align.shape()[1]; ++j) {
	std::cerr << align[i][j];
      }
      std::cerr << std::endl;
    }
    
    // Check breakpoint
    AlignDescriptor ad;
    if (!_findSplit(c, sv.consensus, svRefStr, align, ad, sv.svt)) return false;

    // Get the start and end of the structural variant
    uint32_t finalGapStart = seqstart + ad.rStart;
    uint32_t finalGapEnd = seqstart + ad.rEnd;    
    sv.precise=true;
    sv.svStart=finalGapStart;
    sv.svEnd=finalGapEnd;
    sv.srAlignQuality = ad.percId;
    sv.insLen=ad.cEnd - ad.cStart - 1;
    sv.homLen=std::max(0, ad.homLeft + ad.homRight - 2);
    int32_t ci_wiggle = std::max(ad.homLeft, ad.homRight);
    sv.ciposlow = -ci_wiggle;
    sv.ciposhigh = ci_wiggle;
    sv.ciendlow = -ci_wiggle;
    sv.ciendhigh = ci_wiggle;
    return true;
  }
  
  template<typename TConfig, typename TSRStore>
  inline void
  assemble(TConfig const& c, TSRStore& srStore, std::vector<StructuralVariantRecord>& svs) {
    typedef std::set<std::string> TSequences;
    typedef std::vector<TSequences> TSVSequences;
    TSVSequences seqStore(svs.size(), TSequences());
    std::vector<bool> svcons(svs.size(), false);
    uint32_t maxReadPerSV = 20;
    
    samFile* samfile = sam_open(c.files[0].string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.files[0].string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Parse BAM
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Split-read assembly" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );

    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      ++show_progress;

      // Load sequence
      int32_t seqlen = -1;
      std::string tname(hdr->target_name[refIndex]);
      char* seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
    
      // Parse BAM
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, c.chrlen[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	// Only primary alignments with the full sequence information
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;

	std::size_t seed = hash_lr(rec);
	if (srStore.find(seed) != srStore.end()) {
	  for(uint32_t ri = 0; ri < srStore[seed].size(); ++ri) {
	    int32_t svid = srStore[seed][ri].svid;
	    //std::cerr << seed << '\t' << srStore[seed][ri].svid << '\t' << srStore[seed][ri].sstart << '\t' << srStore[seed][ri].inslen << '\t' << sv[srStore[seed][ri].svid].srSupport << '\t' << sv[srStore[seed][ri].svid].svt << std::endl;

	    if (!svcons[svid]) {
	      // Get sequence
	      std::string sequence;
	      sequence.resize(rec->core.l_qseq);
	      uint8_t* seqptr = bam_get_seq(rec);
	      for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	      if (rec->core.flag & BAM_FREVERSE) reverseComplement(sequence);

	      // Extract subsequence
	      int32_t window = 500; 
	      int32_t sPos = srStore[seed][ri].sstart - window;
	      if (sPos < 0) sPos = 0;
	      int32_t ePos = srStore[seed][ri].sstart + srStore[seed][ri].inslen + window;
	      if (ePos > (int32_t) sequence.size()) ePos = sequence.size();
	      // Min. seq length and max insertion size, 10kbp?
	      if (((ePos - sPos) > window) && ((ePos - sPos) <= c.indelsize)) {
		std::string seqalign = sequence.substr(sPos, (ePos - sPos));
		if (rec->core.flag & BAM_FREVERSE) reverseComplement(seqalign);
		seqStore[svid].insert(seqalign);
	      
		// Enough split-reads?
		if ((!_translocation(svs[svid].svt)) && (svs[svid].chr == refIndex)) {
		  if ((seqStore[svid].size() == maxReadPerSV) || ((int32_t) seqStore[svid].size() == svs[svid].srSupport)) {
		    bool msaSuccess = false;
		    if (seqStore[svid].size() > 1) {
		      //std::cerr << svs[svid].svStart << ',' << svs[svid].svEnd << ',' << svs[svid].svt << ',' << svid << " SV" << std::endl;
		      msa(c, seqStore[svid], svs[svid].consensus);
		      if (((svs[svid].svt == 2) || (svs[svid].svt == 4)) && ((svs[svid].svEnd - svs[svid].svStart) <= c.indelsize)) {
			//if (teguaAlignConsensus(c, hdr, seq, svs[svid])) msaSuccess = true;
			if (alignConsensus(c, hdr, seq, NULL, svs[svid])) msaSuccess = true;
		      } else {
			if (alignConsensus(c, hdr, seq, NULL, svs[svid])) msaSuccess = true;
		      }
		    }
		    if (!msaSuccess) {
		      svs[svid].consensus = "";
		      svs[svid].srSupport = 0;
		      svs[svid].srAlignQuality = 0;
		    }
		    seqStore[svid].clear();
		    svcons[svid] = true;
		  }
		}
	      }
	    }
	  }
	}
      }
      if (seq != NULL) free(seq);
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
    fai_destroy(fai);
    
    // Clean-up unfinished SVs
    for(uint32_t svid = 0; svid < svcons.size(); ++svid) {
      if (!svcons[svid]) {
	//std::cerr << "Missing: " << svid << ',' << svs[svid].svt << std::endl;
	svs[svid].consensus = "";
	svs[svid].srSupport = 0;
	svs[svid].srAlignQuality = 0;
      }
    }
	

    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }


}

#endif
