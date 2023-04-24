#ifndef WFA_H
#define WFA_H

#include <iostream>

#include "bindings/cpp/WFAligner.hpp"


namespace torali
{

  template<typename TAlign>
  inline void
  convertAlignment(std::string const& query, TAlign& align, std::string const& cigar) {
    // Copy input alignment
    uint32_t seqPos = align.shape()[0];
    TAlign alignIn;
    alignIn.resize(boost::extents[seqPos][align.shape()[1]]);
    for(uint32_t i = 0; i < seqPos; ++i) {
      for(uint32_t j = 0; j < align.shape()[1]; ++j) {
	alignIn[i][j] = align[i][j];
      }
    }
	
    // Create new alignment
    align.resize(boost::extents[seqPos+1][cigar.size()]);

    // Alignment
    uint32_t aidx = 0;
    for (uint32_t j = 0; j < cigar.size(); ++j) {
      if (cigar[j] == 'D') {
	for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j] = '-';
      } else {
	for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j] = alignIn[seqIdx][aidx];
	++aidx;
      }
    }

    // Query
    uint32_t tidx = 0;
    for (uint32_t j = 0; j < cigar.size(); ++j) {
      if (cigar[j] == 'I') align[seqPos][j] = '-';
      else align[seqPos][j] = query[tidx++];
    }
  }




  inline char
  _consLetter(uint32_t const maxIdx) {
    switch (maxIdx) {
     case 0: return 'A';
     case 1: return 'C';
     case 2: return 'G';
     case 3: return 'T';
     default: return '-';
    }
    return 'N';
  }

  template<typename TAlign>
  inline void
  consensusWfa2(TAlign const& align, std::string& cons) {
    typedef typename TAlign::index TAIndex;

    std::vector<uint32_t> seqStart(align.shape()[0], align.shape()[1]);
    std::vector<uint32_t> seqEnd(align.shape()[0], 0);
    for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {    
      for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
	if ((align[i][j] != '-') && (j < seqStart[i])) seqStart[i] = j;
	if ((align[i][j] != '-') && (j > seqEnd[i])) seqEnd[i] = j;
      }
    }

    // Consensus
    cons.resize(align.shape()[1]);
    for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
      std::vector<int32_t> count(6, 0); // ACGT-N
      for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
	if ((j >= seqStart[i]) && (j <= seqEnd[i])) {
	  if ((align[i][j] == 'A') || (align[i][j] == 'a')) ++count[0];
	  else if ((align[i][j] == 'C') || (align[i][j] == 'c')) ++count[1];
	  else if ((align[i][j] == 'G') || (align[i][j] == 'g')) ++count[2];
	  else if ((align[i][j] == 'T') || (align[i][j] == 't')) ++count[3];
	  else if (align[i][j] == '-') ++count[4];
	  else ++count[5];
	}
      }
      uint32_t maxIdx = 0;
      uint32_t sndIdx = 1;
      if (count[maxIdx] < count[sndIdx]) {
	maxIdx = 1;
	sndIdx = 0;
      }
      // Ignore Ns
      for(uint32_t i = 2; i<5; ++i) {
	if (count[i] > count[maxIdx]) {
	  sndIdx = maxIdx;
	  maxIdx = i;
	}
	else if (count[i] > count[sndIdx]) {
	  sndIdx = i;
	}
      }
      if (2 * count[sndIdx] <= count[maxIdx]) cons[j] = _consLetter(maxIdx);
      else {
	// Tie?
	if (count[sndIdx] == count[maxIdx]) {
	  if (sndIdx == 4) cons[j] = _consLetter(maxIdx);
	  else if (maxIdx == 4)  cons[j] = _consLetter(sndIdx);
	  else cons[j] = 'N';
	} else {
	  cons[j] = 'N';
	}
      }
    }
  }

  
}

#endif
