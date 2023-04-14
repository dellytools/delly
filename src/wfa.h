#ifndef WFA_H
#define WFA_H

#include <iostream>

#include "bindings/cpp/WFAligner.hpp"
using namespace wfa;

namespace torali
{
  /*
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
  */

  /*
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
  */
  

  template<typename TConfig, typename TSplitReadSet>
  inline int
  msaWfa(TConfig const& c, TSplitReadSet& sps, std::string& cs) {
    // Pairwise scores
    WFAlignerGapAffine aligner(4, 6, 2, WFAligner::Alignment, WFAligner::MemoryHigh);
    std::vector<int32_t> dist(sps.size() * sps.size(), 0);
    for(uint32_t i = 0; i < sps.size(); ++i) {
      for(uint32_t j = i + 1; j < sps.size(); ++j) {
	aligner.alignEnd2End(sps[i], sps[j]);
	std::cerr << aligner.getAlignmentScore() << std::endl;
      }
    }

    /*
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
      //std::cerr << "Progressive MSA: " << i << '(' << align.shape()[0] << ':' << align.shape()[1] << ')' << std::endl;
      //for(uint32_t i = 0; i<align.shape()[0]; ++i) {
      //for(uint32_t j = 0; j<align.shape()[1]; ++j) std::cerr << align[i][j];
      //std::cerr << std::endl;
      //}
      //std::cerr << "Consensus: " << std::endl;
      //std::cerr << alignStr << std::endl;
      //std::cerr << "ToBeAligned: " << sps[selectedIdx[i]] << std::endl;
      // Compute alignment
      EdlibAlignResult cigar = edlibAlign(sps[selectedIdx[i]].c_str(), sps[selectedIdx[i]].size(), alignStr.c_str(), alignStr.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, additionalEqualities, 20));
      convertAlignment(sps[selectedIdx[i]], align, EDLIB_MODE_NW, cigar);
      edlibFreeAlignResult(cigar);
    }
    
    // Debug MSA
    //std::cerr << "Output MSA " << '(' << align.shape()[0] << ':' << align.shape()[1] << ')' << std::endl;
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
    int32_t trim = (int32_t) (0.05 * cs.size());
    if (trim > 50) trim = 50;
    int32_t len = (int32_t) (cs.size()) - 2 * trim;
    if (len > 100) cs = cs.substr(trim, len);
    
    // Return split-read support
    return align.shape()[0];
    */
    return 0;
  }


}

#endif
