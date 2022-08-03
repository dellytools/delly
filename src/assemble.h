#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include <iostream>

#include "edlib.h"
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

    SeqSlice() : svid(-1), sstart(-1), inslen(-1), qual(-1) {}
    SeqSlice(int32_t const sv, int32_t const sst, int32_t const il, int32_t q) : svid(sv), sstart(sst), inslen(il), qual(q) {}
  };


  template<typename TAlign>
  inline void
  convertAlignment(std::string const& query, TAlign& align, EdlibAlignMode const modeCode, EdlibAlignResult& cigar) {
    // Input alignment
    TAlign alignIn;
    alignIn.resize(boost::extents[align.shape()[0]][align.shape()[1]]);
    for(uint32_t i = 0; i < align.shape()[0]; ++i) {
      for(uint32_t j = 0; j < align.shape()[1]; ++j) {
	alignIn[i][j] = align[i][j];
      }
    }
	
    // Create new alignment
    uint32_t seqPos = alignIn.shape()[0];
    align.resize(boost::extents[alignIn.shape()[0]+1][cigar.alignmentLength]);
    int32_t tIdx = -1;
    int32_t qIdx = -1;
    if (modeCode == EDLIB_MODE_HW) {
        tIdx = cigar.endLocations[0];
        for (int32_t i = 0; i < cigar.alignmentLength; i++) {
            if (cigar.alignment[i] != EDLIB_EDOP_INSERT) tIdx--;
        }
    }
    // target
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_INSERT) {
	for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j] = '-';
      } else {
	++tIdx;
	for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j] = alignIn[seqIdx][tIdx];
      }
    }

    // query
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_DELETE) align[seqPos][j] = '-';
      else align[seqPos][j] = query[++qIdx];
    }
  }

  
  template<typename TAlign>
  inline void
  consensusEdlib(TAlign const& align, std::string& cons) {
    typedef typename TAlign::index TAIndex;

    cons.resize(align.shape()[1]);
    for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
      std::vector<int32_t> count(5, 0); // ACGT-
      for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
	if ((align[i][j] == 'A') || (align[i][j] == 'a')) ++count[0];
	else if ((align[i][j] == 'C') || (align[i][j] == 'c')) ++count[1];
	else if ((align[i][j] == 'G') || (align[i][j] == 'g')) ++count[2];
	else if ((align[i][j] == 'T') || (align[i][j] == 't')) ++count[3];
	else ++count[4];
      }
      uint32_t maxIdx = 0;
      uint32_t sndIdx = 1;
      if (count[maxIdx] < count[sndIdx]) {
	maxIdx = 1;
	sndIdx = 0;
      }
      for(uint32_t i = 2; i<5; ++i) {
	if (count[i] > count[maxIdx]) {
	  sndIdx = maxIdx;
	  maxIdx = i;
	}
	else if (count[i] > count[sndIdx]) {
	  sndIdx = i;
	}
      }
      if (2 * count[sndIdx] < count[maxIdx]) {
	switch (maxIdx) {
	case 0: cons[j] = 'A'; break;
	case 1: cons[j] = 'C'; break;
	case 2: cons[j] = 'G'; break;
	case 3: cons[j] = 'T'; break;
	default: cons[j] = '-'; break;
	}
      } else {
	uint32_t k1 = maxIdx;
	uint32_t k2 = sndIdx;
	if (k1 > k2) {
	  k1 = sndIdx;
	  k2 = maxIdx;
	}
	// ACGT-
	if ((k1 == 0) && (k2 == 1)) cons[j] = 'M';
	else if ((k1 == 0) && (k2 == 2)) cons[j] = 'R';
	else if ((k1 == 0) && (k2 == 3)) cons[j] = 'W';
	else if ((k1 == 0) && (k2 == 4)) cons[j] = 'B';
	else if ((k1 == 1) && (k2 == 2)) cons[j] = 'S';
	else if ((k1 == 1) && (k2 == 3)) cons[j] = 'Y';
	else if ((k1 == 1) && (k2 == 4)) cons[j] = 'D';
	else if ((k1 == 2) && (k2 == 3)) cons[j] = 'K';
	else if ((k1 == 2) && (k2 == 4)) cons[j] = 'E';
	else if ((k1 == 3) && (k2 == 4)) cons[j] = 'F';
	else cons[j] = '-';
      }
    }
  }
  

  template<typename TConfig, typename TSplitReadSet>
  inline int
  msaEdlib(TConfig const& c, TSplitReadSet& sps, std::string& cs) {
    // Pairwise scores
    std::vector<int32_t> edit(sps.size() * sps.size(), 0);
    for(uint32_t i = 0; i < sps.size(); ++i) {
      for(uint32_t j = i + 1; j < sps.size(); ++j) {
	EdlibAlignResult align = edlibAlign(sps[i].c_str(), sps[i].size(), sps[j].c_str(), sps[j].size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
	edit[i * sps.size() + j] = align.editDistance;
	edit[j * sps.size() + i] = align.editDistance;
	edlibFreeAlignResult(align);
      }
    }

    // Find best sequence to start alignment
    uint32_t bestIdx = 0;
    int32_t bestVal = sps[0].size();
    for(uint32_t i = 0; i < sps.size(); ++i) {
      std::vector<int32_t> dist(sps.size());
      for(uint32_t j = 0; j < sps.size(); ++j) dist[j] = edit[i * sps.size() + j];
      std::sort(dist.begin(), dist.end());
      if (dist[sps.size()/2] < bestVal) {
	bestVal = dist[sps.size()/2];
	bestIdx = i;
      }
    }

    // Align to best sequence
    std::vector<std::pair<int32_t, int32_t> > qscores;
    qscores.push_back(std::make_pair(0, bestIdx));
    std::string revc = sps[bestIdx];
    reverseComplement(revc);
    for(uint32_t j = 0; j < sps.size(); ++j) {
      if (j != bestIdx) {
	EdlibAlignResult align = edlibAlign(revc.c_str(), revc.size(), sps[j].c_str(), sps[j].size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
	if (align.editDistance < edit[bestIdx * sps.size() + j]) {
	  reverseComplement(sps[j]);
	  qscores.push_back(std::make_pair(align.editDistance, j));
	} else qscores.push_back(std::make_pair(edit[bestIdx * sps.size() + j], j));
	edlibFreeAlignResult(align);
      }
    }
    std::sort(qscores.begin(), qscores.end());
    
    // Drop poorest 20% and order by centroid
    std::vector<uint32_t> selectedIdx;
    uint32_t lastIdx = (uint32_t) (0.8 * qscores.size());
    if (lastIdx < 3) lastIdx = 3;
    for(uint32_t i = 0; ((i < qscores.size()) && (i < lastIdx)); ++i) selectedIdx.push_back(qscores[i].second);
    
    // Extended IUPAC code
    EdlibEqualityPair additionalEqualities[20] = {{'M', 'A'}, {'M', 'C'}, {'R', 'A'}, {'R', 'G'}, {'W', 'A'}, {'W', 'T'}, {'B', 'A'}, {'B', '-'}, {'S', 'C'}, {'S', 'G'}, {'Y', 'C'}, {'Y', 'T'}, {'D', 'C'}, {'D', '-'}, {'K', 'G'}, {'K', 'T'}, {'E', 'G'}, {'E', '-'}, {'F', 'T'}, {'F', '-'}};

    // Incrementally align sequences    
    typedef boost::multi_array<char, 2> TAlign;
    TAlign align;
    align.resize(boost::extents[1][sps[selectedIdx[0]].size()]);
    uint32_t ind = 0;
    for(typename std::string::const_iterator str = sps[selectedIdx[0]].begin(); str != sps[selectedIdx[0]].end(); ++str) align[0][ind++] = *str;
    for(uint32_t i = 1; i < selectedIdx.size(); ++i) {
      // Convert to consensus
      std::string alignStr;
      consensusEdlib(align, alignStr);
      // Debug MSA
      //std::cerr << "Input MSA" << std::endl;
      //for(uint32_t i = 0; i<align.shape()[0]; ++i) {
      //for(uint32_t j = 0; j<align.shape()[1]; ++j) std::cerr << align[i][j];
      //std::cerr << std::endl;
      //}
      //std::cerr << alignStr << std::endl;

      // Compute alignment
      EdlibAlignResult cigar = edlibAlign(sps[selectedIdx[i]].c_str(), sps[selectedIdx[i]].size(), alignStr.c_str(), alignStr.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 20));
      convertAlignment(sps[selectedIdx[i]], align, EDLIB_MODE_HW, cigar);
      edlibFreeAlignResult(cigar);
    }
    
    // Debug MSA
    //std::cerr << "Output MSA" << std::endl;
    //for(uint32_t i = 0; i<align.shape()[0]; ++i) {
    //for(uint32_t j = 0; j<align.shape()[1]; ++j) std::cerr << align[i][j];
    //std::cerr << std::endl;
    //}

    // Consensus
    std::string gapped;
    consensus(c, align, gapped, cs);
    //std::cerr << "Consensus:" << std::endl;
    //std::cerr << gapped << std::endl;

    // Trim off 10% from either end
    int32_t trim = (int32_t) (0.1 * cs.size());
    int32_t len = (int32_t) (cs.size()) - 2 * trim;
    if (len > 100) cs = cs.substr(trim, len);
    
    // Return split-read support
    return align.shape()[0];
  }


  
  template<typename TConfig, typename TValidRegion, typename TSRStore>
  inline void
    assemble(TConfig const& c, TValidRegion const& validRegions, std::vector<StructuralVariantRecord>& svs, TSRStore& srStore) {
    typedef typename TSRStore::value_type TPosReadSV;
    
    // Sequence store
    typedef std::vector<std::string> TSequences;
    typedef std::vector<TSequences> TSVSequences;
    TSVSequences seqStore(svs.size(), TSequences());

    // SV consensus done
    std::vector<bool> svcons(svs.size(), false);

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

    // Parse BAM
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Split-read assembly" << std::endl;

    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      if (validRegions[refIndex].empty()) continue;
      if (srStore[refIndex].empty()) continue;

      // Load sequence
      int32_t seqlen = -1;
      std::string tname(hdr->target_name[refIndex]);
      char* seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);

      // Collect all split-read pos
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet hits(hdr->target_len[refIndex]);
      for(typename TPosReadSV::const_iterator it = srStore[refIndex].begin(); it != srStore[refIndex].end(); ++it) hits[it->first.first] = 1;
      
      // Collect reads from all samples
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	// Read alignments (full chromosome because primary alignments might be somewhere else)
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  // Only primary alignments with the full sequence information
	  if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
	  if (!hits[rec->core.pos]) continue;

	  std::size_t seed = hash_lr(rec);
	  if (srStore[refIndex].find(std::make_pair(rec->core.pos, seed)) != srStore[refIndex].end()) {
	    // Get sequence
	    std::string sequence;
	    sequence.resize(rec->core.l_qseq);
	    uint8_t* seqptr = bam_get_seq(rec);
	    for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	    int32_t readlen = sequence.size();

	    // Iterate all spanned SVs
	    for(uint32_t ri = 0; ri < srStore[refIndex][std::make_pair(rec->core.pos, seed)].size(); ++ri) {
	      SeqSlice seqsl = srStore[refIndex][std::make_pair(rec->core.pos, seed)][ri];
	      int32_t svid = seqsl.svid;
	      if ((!svcons[svid]) && (seqStore[svid].size() < c.maxReadPerSV)) {
		// Extract subsequence (otherwise MSA takes forever)
		int32_t window = 1000; // MSA should be larger
		int32_t sPos = seqsl.sstart - window;
		int32_t ePos = seqsl.sstart + seqsl.inslen + window;
		if (rec->core.flag & BAM_FREVERSE) {
		  sPos = (readlen - (seqsl.sstart + seqsl.inslen)) - window;
		  ePos = (readlen - seqsl.sstart) + window;
		}
		if (sPos < 0) sPos = 0;
		if (ePos > (int32_t) readlen) ePos = readlen;
		// Min. seq length and max insertion size, 10kbp?
		if (((ePos - sPos) > window) && ((ePos - sPos) <= (c.maxInsertionSize + window))) {
		  std::string seqalign = sequence.substr(sPos, (ePos - sPos));
		  seqStore[svid].push_back(seqalign);

		  // Enough split-reads?
		  if ((!_translocation(svs[svid].svt)) && (svs[svid].chr == refIndex)) {
		    if ((seqStore[svid].size() == c.maxReadPerSV) || ((int32_t) seqStore[svid].size() == svs[svid].srSupport)) {
		      bool msaSuccess = false;
		      if (seqStore[svid].size() > 1) {
			//std::cerr << svs[svid].svStart << ',' << svs[svid].svEnd << ',' << svs[svid].svt << ',' << svid << " SV" << std::endl;
			msaEdlib(c, seqStore[svid], svs[svid].consensus);
			if (alignConsensus(c, hdr, seq, NULL, svs[svid], true)) msaSuccess = true;
			//std::cerr << msaSuccess << std::endl;
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
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }
      // Handle left-overs and translocations
      for(int32_t refIndex2 = 0; refIndex2 <= refIndex; ++refIndex2) {
	char* sndSeq = NULL;
	for(uint32_t svid = 0; svid < svcons.size(); ++svid) {
	  if (!svcons[svid]) {
	    if (seqStore[svid].size() > 1) {
	      bool computeMSA = false;
	      if (_translocation(svs[svid].svt)) {
		if ((refIndex2 != refIndex) && (svs[svid].chr == refIndex) && (svs[svid].chr2 == refIndex2)) {
		  computeMSA = true;
		  // Lazy loading of references
		  if (sndSeq == NULL) {
		    int32_t seqlen = -1;
		    std::string tname(hdr->target_name[refIndex2]);
		    sndSeq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex2], &seqlen);
		  }
		}
	      } else {
		if ((refIndex2 == refIndex) && (svs[svid].chr == refIndex) && (svs[svid].chr2 == refIndex2)) computeMSA = true;
	      }
	      if (computeMSA) {
		bool msaSuccess = false;
		//std::cerr << svs[svid].svStart << ',' << svs[svid].svEnd << ',' << svs[svid].svt << ',' << svid << " SV" << std::endl;
		msaEdlib(c, seqStore[svid], svs[svid].consensus);
		if (alignConsensus(c, hdr, seq, sndSeq, svs[svid], true)) msaSuccess = true;
		//std::cerr << msaSuccess << std::endl;
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
	if (sndSeq != NULL) free(sndSeq);
      }
      // Clean-up
      if (seq != NULL) free(seq);
    }
    // Clean-up
    fai_destroy(fai);
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
    
    // Clean-up unfinished SVs
    for(uint32_t svid = 0; svid < svcons.size(); ++svid) {
      if (!svcons[svid]) {
	//std::cerr << "Missing: " << svid << ',' << svs[svid].svt << std::endl;
	svs[svid].consensus = "";
	svs[svid].srSupport = 0;
	svs[svid].srAlignQuality = 0;
      }
    }
  }


}

#endif
