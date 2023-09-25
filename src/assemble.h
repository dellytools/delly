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
    int32_t tIdx = -1;
    int32_t qIdx = -1;
    uint32_t missingEnd = 0;
    uint32_t missingStart = 0;
    if (modeCode == EDLIB_MODE_HW) {
      tIdx = cigar.endLocations[0];
      if (tIdx < (int32_t) alignIn.shape()[1]) missingEnd = alignIn.shape()[1] - tIdx - 1;
      for (int32_t i = 0; i < cigar.alignmentLength; i++) {
	if (cigar.alignment[i] != EDLIB_EDOP_INSERT) tIdx--;
      }
      if (tIdx >= 0) missingStart = tIdx + 1;
    }
    align.resize(boost::extents[alignIn.shape()[0]+1][missingStart + cigar.alignmentLength + missingEnd]);

    // infix alignment, fix start
    if (modeCode == EDLIB_MODE_HW) {
      if (missingStart) {
	for (uint32_t j = 0; j < missingStart; ++j) {
	  for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j] = alignIn[seqIdx][j];
	  align[seqPos][j] = '-';
	}
      }
    }
    
    // target
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_INSERT) {
	for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j + missingStart] = '-';
      } else {
	++tIdx;
	for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j + missingStart] = alignIn[seqIdx][tIdx];
      }
    }

    // query
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_DELETE) align[seqPos][j + missingStart] = '-';
      else align[seqPos][j + missingStart] = query[++qIdx];
    }

    // infix alignment, fix end
    if (modeCode == EDLIB_MODE_HW) {
      if (missingEnd) {
	for (uint32_t j = cigar.alignmentLength + missingStart; j < cigar.alignmentLength + missingStart + missingEnd; ++j) {
	  ++tIdx;
	  for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j] = alignIn[seqIdx][tIdx];
	  align[seqPos][j] = '-';
	}
      }
    }
  }

  inline void
  buildSuperstring(std::string const& seqI, std::string const& seqJ, std::string& outStr, EdlibAlignResult& cigar, uint32_t const preI, uint32_t const postI, uint32_t const preJ, uint32_t const postJ) {
    int32_t iIdx = 0;
    int32_t jIdx = 0;
    // prefix
    bool firstSeq = false;
    if (preI > preJ) {
      firstSeq = true;
      for(uint32_t j = 0; j < preI; ++j) outStr += seqI[iIdx++];
      for(uint32_t j = 0; j < preJ; ++j) ++jIdx;
    } else {
      for(uint32_t j = 0; j < preI; ++j) ++iIdx;
      for(uint32_t j = 0; j < preJ; ++j) outStr += seqJ[jIdx++];  
    }
    //outStr += '#';

    // parse alignment
    int32_t bp = cigar.alignmentLength/2;
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (bp == j) {
	firstSeq = !firstSeq;
	//outStr += "#";
      }
      if (cigar.alignment[j] == EDLIB_EDOP_DELETE) {
	if (!firstSeq) outStr += seqJ[jIdx];
	++jIdx;
      } else if (cigar.alignment[j] == EDLIB_EDOP_INSERT) {
	if (firstSeq) outStr += seqI[iIdx];
	++iIdx;
      } else {
	if (firstSeq) outStr += seqI[iIdx];
	else outStr += seqJ[jIdx];
	++iIdx;
	++jIdx;
      }
    }

    // Post sequence
    if (postI > postJ) {
      for(uint32_t j = 0; j < postI; ++j) outStr += seqI[iIdx++];
    } else {
      for(uint32_t j = 0; j < postJ; ++j) outStr += seqJ[jIdx++];
    }
  }

  
  template<typename TAlign>
  inline void
  convertAlignment(std::string const& query, TAlign& align, EdlibAlignMode const modeCode, EdlibAlignResult& cigar, uint32_t const preI, uint32_t const postI, uint32_t const preJ, uint32_t const postJ) {
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
    int32_t tIdx = 0;
    int32_t qIdx = 0;
    align.resize(boost::extents[alignIn.shape()[0]+1][cigar.alignmentLength + preI + postI + preJ + postJ]);
    if (preI) {
      for(uint32_t j = 0; j < preI; ++j) {
	for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j] = alignIn[seqIdx][tIdx];
	++tIdx;
	align[seqPos][j] = '-';
      }
    }
    if (preJ) {
      for(uint32_t j = 0; j < preJ; ++j) {
	for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j+preI] = '-';
	align[seqPos][j+preI] = query[qIdx++];
      }
    }

    // target
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_DELETE) {
	for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j+preI+preJ] = '-';
      } else {
	for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j+preI+preJ] = alignIn[seqIdx][tIdx];
	++tIdx;
      }
    }

    // query
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_INSERT) align[seqPos][j+preI+preJ] = '-';
      else align[seqPos][j+preI+preJ] = query[qIdx++];
    }

    // Post sequence
    if (postI) {
      for(uint32_t j = 0; j < postI; ++j) {
	for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j+preI+preJ+cigar.alignmentLength] = alignIn[seqIdx][tIdx];
	tIdx++;
	align[seqPos][j+preI+preJ+cigar.alignmentLength] = '-';
      }
    }
    if (postJ) {
      for(uint32_t j = 0; j < postJ; ++j) {
	for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j+preI+preJ+cigar.alignmentLength+postI] = '-';
	align[seqPos][j+preI+preJ+cigar.alignmentLength+postI] = query[qIdx++];
      }
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


  template<typename TAlign>
  inline void
  consensusWfa(TAlign const& align, std::string& cons) {
    typedef typename TAlign::index TAIndex;

    // Read positions
    std::vector<uint32_t> readStart(align.shape()[0], align.shape()[1]);
    std::vector<uint32_t> readEnd(align.shape()[0], 0);
    for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
      for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
	if (align[i][j] != '-') {
	  if (j < readStart[i]) readStart[i] = j;
	  if (j > readEnd[i]) readEnd[i] = j;
	}
      }
    }

    // Consensus
    cons.resize(align.shape()[1]);
    for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
      std::vector<int32_t> count(5, 0); // ACGT-
      for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
	if ((j >= readStart[i]) && (j <= readEnd[i])) {
	  if ((align[i][j] == 'A') || (align[i][j] == 'a')) ++count[0];
	  else if ((align[i][j] == 'C') || (align[i][j] == 'c')) ++count[1];
	  else if ((align[i][j] == 'G') || (align[i][j] == 'g')) ++count[2];
	  else if ((align[i][j] == 'T') || (align[i][j] == 't')) ++count[3];
	  else ++count[4];
	}
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

  inline void
  _trimConsensus(std::string const& prefix, std::string const& suffix, std::string& cs) {
    std::string prefixRev = prefix;
    reverseComplement(prefixRev);
    EdlibAlignResult cigarFwd = edlibAlign(prefix.c_str(), prefix.size(), cs.c_str(), cs.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
    int32_t scoreFwd = cigarFwd.editDistance;
    edlibFreeAlignResult(cigarFwd);
    EdlibAlignResult cigarRev = edlibAlign(prefixRev.c_str(), prefixRev.size(), cs.c_str(), cs.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
    int32_t scoreRev = cigarRev.editDistance;
    edlibFreeAlignResult(cigarRev);
    if (scoreFwd > scoreRev) reverseComplement(cs);
    
    // Anchor reference probes
    EdlibAlignResult cigarPrefix = edlibAlign(prefix.c_str(), prefix.size(), cs.c_str(), cs.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    uint32_t csStart = infixStart(cigarPrefix);
    //std::cerr << "Prefix alignment: " << cigarPrefix.editDistance << std::endl;
    //printAlignment(prefix, cs, EDLIB_MODE_HW, cigarPrefix);
    edlibFreeAlignResult(cigarPrefix);
    EdlibAlignResult cigarSuffix = edlibAlign(suffix.c_str(), suffix.size(), cs.c_str(), cs.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    uint32_t csEnd = infixEnd(cigarSuffix);
    //std::cerr << "Suffix alignment: " << cigarSuffix.editDistance << std::endl;
    //printAlignment(suffix, cs, EDLIB_MODE_HW, cigarSuffix);
    edlibFreeAlignResult(cigarSuffix);
    //std::cerr << "Trimming: " << csStart << ',' << csEnd << '(' << cs.size() << ')' << std::endl;
    
    // Trim consensus
    if ((csStart < csEnd) && (csEnd < cs.size())) cs = cs.substr(csStart, (csEnd - csStart));
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
  }

  inline uint32_t
  charToInt(char c) {
    switch(c)
      {
      case 'A':
	return 0;
      case 'C':
	return 1;
      case 'G':
	return 2;
      case 'T':
	return 3;
      case 'B':
	return 0;
      case 'D':
	return 1;
      case 'E':
	return 2;
      case 'F':
	return 3;
      }
    return 0;
  }


  
  inline void
  fillKmerTable(std::string const& s, std::vector<uint32_t>& kmerpos) {
    uint32_t len = s.size();
    kmerpos.resize(std::pow(4, DELLY_KMER + 1));
    std::fill(kmerpos.begin(), kmerpos.end(), 0);
    uint32_t hash = 0;
    for(uint32_t ki = 0; ((ki < len) && (ki < DELLY_KMER)); ++ki) {
	hash *= 4;
	hash += charToInt(s[ki]);
    }
    for(uint32_t ki = DELLY_KMER; ki < len; ++ki) {
      if (kmerpos[hash]) kmerpos[hash] = DELLY_DUPLICATE; // Duplicate, flag with length
      else kmerpos[hash] = (ki - DELLY_KMER + 1);
      hash -= (charToInt(s[ki - DELLY_KMER]) * 4 * 4 * 4 * 4 * 4 * 4);
      hash *= 4;
      hash += charToInt(s[ki]);
    }
    if (kmerpos[hash]) kmerpos[hash] = DELLY_DUPLICATE; // Duplicate, flag with length
    else kmerpos[hash] = (len - DELLY_KMER + 1);
  }

  inline int32_t
  bestDiagonal(std::vector<uint32_t> const& kmerHitI, std::vector<uint32_t> const& kmerHitJ, uint32_t const lenI, uint32_t const lenJ) {
    // Shared unique kmers
    std::vector<uint32_t> diag(lenI + lenJ, 0);
    for(uint32_t k = 0; k < std::pow(4, DELLY_KMER + 1); ++k) {
      if ((kmerHitI[k]) && (kmerHitJ[k]) && (kmerHitI[k] != DELLY_DUPLICATE) && (kmerHitJ[k] != DELLY_DUPLICATE)) ++diag[lenJ + kmerHitI[k] - kmerHitJ[k]];
    }

    // Find best diagonal
    uint32_t window = 20;
    uint32_t windowVal = 0;
    for(uint32_t d = 0; ((d < diag.size()) && (d < window)); ++d) windowVal += diag[d];
    uint32_t bestDiag = window / 2;
    uint32_t bestWindowVal = windowVal;
    for(uint32_t d = window; d < diag.size(); ++d) {
      windowVal -= diag[d - window];
      windowVal += diag[d];
      if (windowVal > bestWindowVal) {
	bestWindowVal = windowVal;
	bestDiag = d - window / 2;
      }
    }
    return ((int) (bestDiag) - (int) (lenJ));
  }
    
  template<typename TConfig, typename TSplitReadSet>
  inline int
  msaWfa(TConfig const& c, TSplitReadSet& sps, std::string& cs, std::string const& prefix, std::string const& suffix) {
    // Pairwise scores
    std::vector<int32_t> edit(sps.size() * sps.size(), 0);
    std::vector<uint32_t> kmerHitI;
    std::vector<uint32_t> kmerHitJ;
    for(uint32_t i = 0; i < sps.size(); ++i) {
      uint32_t lenI = sps[i].size();
      fillKmerTable(sps[i], kmerHitI);
      for(uint32_t j = i + 1; j < sps.size(); ++j) {
	uint32_t lenJ = sps[j].size();
	fillKmerTable(sps[j], kmerHitJ);
	int32_t bestDiag = bestDiagonal(kmerHitI, kmerHitJ, lenI, lenJ);
	std::string seqI;
	std::string seqJ;
	if (bestDiag >= 0) {
	  uint32_t seqlen = std::min(lenI - bestDiag, lenJ);
	  seqI = sps[i].substr(bestDiag, seqlen);
	  seqJ = sps[j].substr(0, seqlen);
	} else {
	  uint32_t seqlen = std::min(lenJ + bestDiag, lenI);
	  seqI = sps[i].substr(0, seqlen);
	  seqJ = sps[j].substr(-1 * bestDiag, seqlen);
	}
	EdlibAlignResult align = edlibAlign(seqI.c_str(), seqI.size(), seqJ.c_str(), seqJ.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
	int32_t score = (align.editDistance * 1000) / std::max(seqI.size(), seqJ.size());
	edit[i * sps.size() + j] = score;
	edit[j * sps.size() + i] = score;
	//std::cerr << "Diagonal: " << bestDiag << ", Score: " << score << std::endl;
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
    uint32_t lenI = revc.size();
    fillKmerTable(revc, kmerHitI);
    for(uint32_t j = 0; j < sps.size(); ++j) {
      if (j != bestIdx) {
	uint32_t lenJ = sps[j].size();
	fillKmerTable(sps[j], kmerHitJ);
	int32_t bestDiag = bestDiagonal(kmerHitI, kmerHitJ, lenI, lenJ);
	std::string seqI;
	std::string seqJ;
	if (bestDiag >= 0) {
	  uint32_t seqlen = std::min(lenI - bestDiag, lenJ);
	  seqI = revc.substr(bestDiag, seqlen);
	  seqJ = sps[j].substr(0, seqlen);
	} else {
	  uint32_t seqlen = std::min(lenJ + bestDiag, lenI);
	  seqI = revc.substr(0, seqlen);
	  seqJ = sps[j].substr(-1 * bestDiag, seqlen);
	}
	// Compute edit distance
	EdlibAlignResult align = edlibAlign(seqI.c_str(), seqI.size(), seqJ.c_str(), seqJ.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
	int32_t score = (align.editDistance * 1000) / std::max(seqI.size(), seqJ.size());	
	if (score < edit[bestIdx * sps.size() + j]) {
	  reverseComplement(sps[j]);
	  qscores.push_back(std::make_pair(score, j));
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

    // Build superstring
    std::string superStr = sps[selectedIdx[0]];
    for(uint32_t i = 1; i < selectedIdx.size(); ++i) {
      // Find alignment diagonal
      uint32_t lenI = superStr.size();
      fillKmerTable(superStr, kmerHitI);
      uint32_t lenJ = sps[selectedIdx[i]].size();
      fillKmerTable(sps[selectedIdx[i]], kmerHitJ);
      int32_t bestDiag = bestDiagonal(kmerHitI, kmerHitJ, lenI, lenJ);
      std::string seqI;
      std::string seqJ;
      uint32_t preI = 0;
      uint32_t postI = 0;
      uint32_t preJ = 0;
      uint32_t postJ = 0;
      uint32_t seqlen = 0;
      if (bestDiag >= 0) {
	seqlen = std::min(lenI - bestDiag, lenJ);
	preI = bestDiag;
	postI = lenI - (bestDiag + seqlen);
	preJ = 0;
	postJ = lenJ - seqlen;
      } else {
	seqlen = std::min(lenJ + bestDiag, lenI);
	preI = 0;
	postI = lenI - seqlen;
	preJ = -1 * bestDiag;
	postJ = lenJ - (-1 * bestDiag + seqlen);
      }
      
      // Can the superstring be extended?
      //std::cerr << "Superstring clips: (" << preI << ',' << postI << ',' << preJ << ',' << postJ << ')' << std::endl;
      if ((preI > preJ)  && (postI > postJ)) {
	// Nested alignment
      } else if ((preJ > preI) && (postJ > postI)) {
	// Nested alignment, new sequence longer
	superStr = sps[selectedIdx[i]];
      } else {
	// Extend superstring
	if (bestDiag >= 0) {
	  seqI = superStr.substr(bestDiag, seqlen);
	  seqJ = sps[selectedIdx[i]].substr(0, seqlen);
	} else {
	  seqI = superStr.substr(0, seqlen);
	  seqJ = sps[selectedIdx[i]].substr(-1 * bestDiag, seqlen);
	}
      
	// Compute alignment
	EdlibAlignResult cigar = edlibAlign(seqI.c_str(), seqI.size(), seqJ.c_str(), seqJ.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
	//printAlignment(seqI, seqJ, EDLIB_MODE_NW, cigar);
	//std::cerr << std::endl;
	//std::cerr << superStr << std::endl;
	//std::cerr << sps[selectedIdx[i]] << std::endl;
	std::string outStr;
	buildSuperstring(superStr, sps[selectedIdx[i]], outStr, cigar, preI, postI, preJ, postJ);
	edlibFreeAlignResult(cigar);

	// Next iteration
	superStr = outStr;
      }
    }

    // Extended IUPAC code
    EdlibEqualityPair additionalEqualities[20] = {{'M', 'A'}, {'M', 'C'}, {'R', 'A'}, {'R', 'G'}, {'W', 'A'}, {'W', 'T'}, {'B', 'A'}, {'B', '-'}, {'S', 'C'}, {'S', 'G'}, {'Y', 'C'}, {'Y', 'T'}, {'D', 'C'}, {'D', '-'}, {'K', 'G'}, {'K', 'T'}, {'E', 'G'}, {'E', '-'}, {'F', 'T'}, {'F', '-'}};

    // Incrementally align sequences    
    typedef boost::multi_array<char, 2> TAlign;
    TAlign align;
    align.resize(boost::extents[1][superStr.size()]);
    uint32_t ind = 0;
    for(typename std::string::const_iterator str = superStr.begin(); str != superStr.end(); ++str) align[0][ind++] = *str;
    for(uint32_t i = 0; i < selectedIdx.size(); ++i) {
      // Convert to consensus
      std::string alignStr;
      consensusWfa(align, alignStr);
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
      EdlibAlignResult cigar = edlibAlign(sps[selectedIdx[i]].c_str(), sps[selectedIdx[i]].size(), alignStr.c_str(), alignStr.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 20));
      convertAlignment(sps[selectedIdx[i]], align, EDLIB_MODE_HW, cigar);
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
    //std::cerr << cs << std::endl;

    // Fix consensus orientation and trim to prefix & suffix
    if ((!prefix.empty()) && (!suffix.empty())) {
      _trimConsensus(prefix, suffix, cs);
    } else {
      // Trim off 10% from either end
      int32_t trim = (int32_t) (0.05 * cs.size());
      if (trim > 50) trim = 50;
      int32_t len = (int32_t) (cs.size()) - 2 * trim;
      if (len > 100) cs = cs.substr(trim, len);
    }
    
    // Return split-read support
    return selectedIdx.size();
  }


  template<typename TConfig, typename TSplitReadSet>
  inline int
  msaWfa(TConfig const& c, TSplitReadSet& sps, std::string& cs) {
    return msaWfa(c, sps, cs, "", "");
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
		int32_t window = c.minConsWindow; // MSA should be larger
		// Add breakpoint uncertainty
		window += std::max(svs[svid].ciposhigh - svs[svid].ciposlow, svs[svid].ciendhigh - svs[svid].ciendlow);
		window += seqsl.inslen;
		int32_t sPos = seqsl.sstart - window;
		int32_t ePos = seqsl.sstart + window;
		if (rec->core.flag & BAM_FREVERSE) {
		  sPos = (readlen - seqsl.sstart) - window;
		  ePos = (readlen - seqsl.sstart) + window;
		}
		if (sPos < 0) sPos = 0;
		if (ePos > (int32_t) readlen) ePos = readlen;
		//std::cerr << bam_get_qname(rec) << ',' << sPos << ',' << ePos << ":" << window << "\t" << sequence.substr(sPos, (ePos - sPos)) << std::endl;
		// Min. seq length and max insertion size, 10kbp?
		if ((ePos - sPos) > window) {
		  seqStore[svid].push_back(sequence.substr(sPos, (ePos - sPos)));

		  // Enough split-reads?
		  if ((!_translocation(svs[svid].svt)) && (svs[svid].chr == refIndex)) {
		    if ((seqStore[svid].size() == c.maxReadPerSV) || ((int32_t) seqStore[svid].size() == svs[svid].srSupport)) {
		      bool msaSuccess = false;
		      if (seqStore[svid].size() > 1) {
			//std::cerr << svs[svid].svStart << ',' << svs[svid].svEnd << ',' << svs[svid].svt << ',' << svid << " SV" << std::endl;
			if (svs[svid].svt != 4) {
			  msaEdlib(c, seqStore[svid], svs[svid].consensus);
			  // Take care of small inversions
			  std::string tmpCons;
			  int32_t offsetTmpCons = 0;
			  int32_t svSize = svs[svid].svEnd - svs[svid].svStart;
			  if (((svs[svid].svt == 0) || (svs[svid].svt == 1)) && (svSize < (int32_t) svs[svid].consensus.size())) {
			    offsetTmpCons = (svs[svid].consensus.size() - svSize) / 2;
			    tmpCons = svs[svid].consensus;
			    svs[svid].consensus = svs[svid].consensus.substr(offsetTmpCons, svSize);
			  }
			  if (alignConsensus(c, hdr, seq, NULL, svs[svid], true)) msaSuccess = true;
			  if (!tmpCons.empty()) {
			    svs[svid].consensus = tmpCons;
			    svs[svid].consBp += offsetTmpCons;
			  }
			} else {
			  std::string prefix = boost::to_upper_copy(std::string(seq + std::max(svs[svid].svStart - (int32_t) c.minConsWindow, 0), seq + svs[svid].svStart));
			  std::string suffix = boost::to_upper_copy(std::string(seq + svs[svid].svStart, seq + std::min(seqlen, svs[svid].svStart + c.minConsWindow)));
			  msaWfa(c, seqStore[svid], svs[svid].consensus, prefix, suffix);
			  if ((int32_t) svs[svid].consensus.size() < svs[svid].insLen + 4 * c.minConsWindow) {
			    if (alignConsensus(c, hdr, seq, NULL, svs[svid], false)) msaSuccess = true;
			  }
			}
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
		if (svs[svid].svt != 4) {
		  msaEdlib(c, seqStore[svid], svs[svid].consensus);
		  // Take care of small inversions
		  std::string tmpCons;
		  int32_t offsetTmpCons = 0;
		  int32_t svSize = svs[svid].svEnd - svs[svid].svStart;
		  if (((svs[svid].svt == 0) || (svs[svid].svt == 1)) && (svSize < (int32_t) svs[svid].consensus.size())) {
		    offsetTmpCons = (svs[svid].consensus.size() - svSize) / 2;
		    tmpCons = svs[svid].consensus;
		    svs[svid].consensus = svs[svid].consensus.substr(offsetTmpCons, svSize);
		  }
		  if (alignConsensus(c, hdr, seq, sndSeq, svs[svid], true)) msaSuccess = true;
		  if (!tmpCons.empty()) {
		    svs[svid].consensus = tmpCons;
		    svs[svid].consBp += offsetTmpCons;
		  }
		} else {
		  std::string prefix = boost::to_upper_copy(std::string(seq + std::max(svs[svid].svStart - (int32_t) c.minConsWindow, 0), seq + svs[svid].svStart));
		  std::string suffix = boost::to_upper_copy(std::string(seq + svs[svid].svStart, seq + std::min(seqlen, svs[svid].svStart + c.minConsWindow)));
		  msaWfa(c, seqStore[svid], svs[svid].consensus, prefix, suffix);
		  if ((int32_t) svs[svid].consensus.size() < svs[svid].insLen + 4 * c.minConsWindow) {
		    if (alignConsensus(c, hdr, seq, NULL, svs[svid], false)) msaSuccess = true;
		  }
		}
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
