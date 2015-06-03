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

#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/tokenizer.hpp>
#include <boost/functional/hash.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include "api/BamReader.h"
#include "api/BamIndex.h"

#ifdef OPENMP
#include <omp.h>
#endif

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif

#include "version.h"
#include "util.h"
#include "dna_score.h"
#include "align_config.h"
#include "align_gotoh.h"
#include "align_nw.h"
#include "align_nw_mat.h"
#include "index.h"
#include "tags.h"
#include "spanning.h"
#include "coverage.h"
#include "junction.h"
#include "fasta_reader.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <stdio.h>

using namespace torali;

// Config arguments
struct Config {
  unsigned short minMapQual;
  unsigned short minGenoQual;
  unsigned short madCutoff;
  unsigned int minimumFlankSize;
  unsigned int minimumSplitRead;
  unsigned int flankQuality;
  unsigned int graphPruning;
  float epsilon;
  float percentAbnormal;
  std::string svType;
  boost::filesystem::path outfile;
  boost::filesystem::path vcffile;
  boost::filesystem::path genome;
  boost::filesystem::path exclude;
  boost::filesystem::path peDump;
  std::vector<boost::filesystem::path> files;
};


// Reduced bam alignment record data structure
struct BamAlignRecord {
  int32_t     Length;
  int32_t     RefID;         
  int32_t     Position;      
  int32_t     MateRefID; 
  int32_t     MatePosition;
  uint16_t    MapQuality;    
  int32_t     Median;
  int32_t     Mad;
  int32_t     maxNormalISize;
  int         libOrient;
  uint32_t    AlignmentFlag; 

  BamAlignRecord(BamTools::BamAlignment const& al, uint16_t pairQuality, int32_t median, int32_t mad, int32_t maxISize, int lO) : Length(al.Length), RefID(al.RefID), Position(al.Position), MateRefID(al.MateRefID), MatePosition(al.MatePosition), MapQuality(pairQuality), Median(median), Mad(mad), maxNormalISize(maxISize), libOrient(lO), AlignmentFlag(al.AlignmentFlag) {}

};

// Sort reduced bam alignment records
template<typename TRecord>
struct SortBamRecords : public std::binary_function<TRecord, TRecord, bool>
{
  inline bool operator()(TRecord const& s1, TRecord const& s2) const {
    if (s1.RefID==s1.MateRefID) {
      return ((std::min(s1.Position,s1.MatePosition) < std::min(s2.Position,s2.MatePosition)) || 
	      ((std::min(s1.Position,s1.MatePosition) == std::min(s2.Position,s2.MatePosition)) && (std::max(s1.Position,s1.MatePosition) < std::max(s2.Position,s2.MatePosition))) ||
	      ((std::min(s1.Position,s1.MatePosition) == std::min(s2.Position,s2.MatePosition)) && (std::max(s1.Position,s1.MatePosition) == std::max(s2.Position,s2.MatePosition)) && (s1.maxNormalISize < s2.maxNormalISize)));
    } else {
      return ((s1.Position < s2.Position) ||
	      ((s1.Position == s2.Position) && (s1.MatePosition < s2.MatePosition)) ||
	      ((s1.Position == s2.Position) && (s1.MatePosition == s2.MatePosition) && (s1.maxNormalISize < s2.maxNormalISize)));
    }
  }
};

// SplitRead struct
struct SplitReadCoord {
  bool forwardRead;
  unsigned int offset;
  int lastKmer;
  int diag;
  std::vector<char> read;
};


template<typename TSequence>
struct OutputRead {
  bool forwardRead;
  unsigned int readOffset;
  TSequence seq;
};

template<typename TRecord>
struct SortSplitReadRecords : public std::binary_function<TRecord, TRecord, bool>
{
  inline bool operator()(TRecord const& s1, TRecord const& s2) {
    if (s1.offset == s2.offset) {
      if (s1.lastKmer == s2.lastKmer) return (s1.diag < s2.diag);
      else return (s1.lastKmer<s2.lastKmer);
    } else return (s1.offset < s2.offset);
  }
};


// Edge struct
template<typename TWeight, typename TVertex>
struct EdgeRecord {
  TVertex source;
  TVertex target;
  TWeight weight;

  EdgeRecord(TVertex s, TVertex t, TWeight w) : source(s), target(t), weight(w) {}
};

// Sort edge records
template<typename TRecord>
struct SortEdgeRecords : public std::binary_function<TRecord, TRecord, bool>
{
  inline bool operator()(TRecord const& e1, TRecord const& e2) const {
    return ((e1.weight < e2.weight) || ((e1.weight == e2.weight) && (e1.source < e2.source)) || ((e1.weight == e2.weight) && (e1.source == e2.source) && (e1.target < e2.target)));
  }
};

// ExcludeInterval
struct ExcludeInterval {
  int32_t RefID;
  int32_t start;
  int32_t end;
  
  ExcludeInterval() {}
  ExcludeInterval(int32_t r, int32_t s, int32_t e) : RefID(r), start(s), end(e) {}
};

// Sort exclude intervals
template<typename TRecord>
struct SortExcludeIntervals : public std::binary_function<TRecord, TRecord, bool>
{
  inline bool operator()(TRecord const& i1, TRecord const& i2) const {
    return ((i1.RefID < i2.RefID) || ((i1.RefID == i2.RefID) && (i1.start < i2.start)) || ((i1.RefID == i2.RefID) && (i1.start == i2.start) && (i1.end < i2.end)));
  }
};

// Deletions
template<typename TUSize, typename TSize>
inline bool
_betterSplit(TUSize, TUSize, TUSize, TSize, TUSize oldOffset, TUSize skippedLength, TUSize initialLength, double epsilon, TUSize support, TUSize bestSupport, int, int, SVType<DeletionTag>) {
  return ((std::abs((((double) oldOffset + skippedLength) / (double) initialLength) - 1.0) <= epsilon) && (support > bestSupport));
}


// Duplications
template<typename TUSize, typename TSize>
inline bool
_betterSplit(TUSize refSize, TUSize refSizeRight, TUSize, TSize oldDiag, TUSize oldOffset, TUSize skippedLength, TUSize initialLength, double epsilon, TUSize support, TUSize bestSupport, int, int, SVType<DuplicationTag>) {
  return ((oldDiag < (TSize) refSizeRight) && ((oldDiag + (TSize) oldOffset) > (TSize) refSizeRight) &&  (std::abs((( (double) refSize - (double) oldOffset + skippedLength) / (double) initialLength) - 1.0) <= epsilon) && (support > bestSupport));
}


// Inversions
template<typename TUSize, typename TSize>
inline bool
_betterSplit(TUSize refSize, TUSize refSizeRight, TUSize refSizeLeft, TSize oldDiag, TUSize oldOffset, TUSize skippedLength, TUSize initialLength, double epsilon, TUSize support, TUSize bestSupport, int ct, int, SVType<InversionTag>) {
  TSize oldSvSize = 0;
  if (ct==1) oldSvSize = (refSizeLeft - oldDiag) + refSize - (oldDiag + oldOffset);
  else oldSvSize = oldDiag + refSizeRight - (refSize - (oldDiag + oldOffset)); 
  return ((oldDiag < (TSize) refSizeLeft) && ((oldDiag + (TSize) oldOffset) > (TSize) refSizeLeft) &&  (std::abs( ( ((double) oldSvSize + skippedLength) / (double) initialLength ) - 1.0) <= epsilon) && (support > bestSupport));
}

// Translocations
template<typename TUSize, typename TSize>
inline bool
_betterSplit(TUSize, TUSize, TUSize refSizeLeft, TSize oldDiag, TUSize oldOffset, TUSize, TUSize, double, TUSize support, TUSize bestSupport, int, int readLen, SVType<TranslocationTag>) {
  return ((oldDiag < (TSize) refSizeLeft) && ((oldDiag + oldOffset + readLen) > refSizeLeft) && (support > bestSupport));
}

// Insertion
template<typename TUSize, typename TSize>
inline bool
_betterSplit(TUSize, TUSize, TUSize, TSize, TUSize, TUSize, TUSize, double, TUSize, TUSize, int, int, SVType<InsertionTag>) {
  //ToDo
  return false;
}




// Deletions
template<typename TUSize, typename TSize, typename TBamRecord>
inline bool
_translateSVCoordinates(TUSize, TUSize, TUSize, TSize diagBeg, TSize diagEnd, TUSize leftRefWalky, TUSize rightRefWalky, TUSize consLen, TBamRecord const& sv, TUSize& finalGapStart, TUSize& finalGapEnd, TUSize& predictedLength, SVType<DeletionTag>) {
  finalGapStart = diagBeg + leftRefWalky + 1;
  finalGapEnd = diagEnd + consLen - (consLen - rightRefWalky);
  predictedLength = finalGapEnd -  finalGapStart;
  TUSize annealed = (sv.svStartEnd - sv.svStartBeg);
  if ((finalGapStart>=annealed) || (finalGapEnd < annealed)) return true;
  else {
    predictedLength+=(sv.svEndBeg - sv.svStartEnd);
    finalGapStart+=sv.svStartBeg;
    finalGapEnd+=(sv.svEndBeg - annealed);
  } 
  return false;
}

// Duplications
template<typename TUSize, typename TSize, typename TBamRecord>
inline bool
_translateSVCoordinates(TUSize refSize, TUSize refSizeRight, TUSize refSizeLeft, TSize diagBeg, TSize diagEnd, TUSize leftRefWalky, TUSize rightRefWalky, TUSize consLen, TBamRecord const& sv, TUSize& finalGapStart, TUSize& finalGapEnd, TUSize& predictedLength, SVType<DuplicationTag>) {
  finalGapStart = diagBeg + leftRefWalky + 1;
  finalGapEnd = diagEnd + consLen - (consLen - rightRefWalky) + 1;
  if ((finalGapEnd <= refSizeRight) || (finalGapStart >= refSizeRight)) return true;
  unsigned int finalGapStartTmp = finalGapStart;
  finalGapStart = (refSizeLeft - (refSize - finalGapEnd)) + sv.svStartBeg - 1;
  finalGapEnd = refSizeLeft + finalGapStartTmp + sv.svStartBeg + (sv.svEndBeg - sv.svStartEnd);
  predictedLength = finalGapEnd - finalGapStart;
  return false;
}

// Inversions
template<typename TUSize, typename TSize, typename TBamRecord>
inline bool
_translateSVCoordinates(TUSize refSize, TUSize refSizeRight, TUSize refSizeLeft, TSize diagBeg, TSize diagEnd, TUSize leftRefWalky, TUSize rightRefWalky, TUSize consLen, TBamRecord const& sv, TUSize& finalGapStart, TUSize& finalGapEnd, TUSize& predictedLength, SVType<InversionTag>) {
  finalGapStart = diagBeg + leftRefWalky + 1;
  finalGapEnd = diagEnd + consLen - (consLen - rightRefWalky) + 1;
  if ((finalGapEnd <= refSizeRight) || (finalGapStart >= refSizeLeft)) return true;
  if (sv.ct==1) {
    finalGapEnd = refSize - finalGapEnd;
    finalGapStart += sv.svStartBeg;
    finalGapEnd += sv.svStartBeg + (sv.svStartEnd - sv.svStartBeg) + (sv.svEndBeg - sv.svStartEnd) + 1;
  } else {
    finalGapEnd = refSizeRight - (refSize - finalGapEnd);
    finalGapStart = refSizeLeft-finalGapStart + sv.svStartBeg;
    finalGapEnd += sv.svStartBeg + (sv.svStartEnd - sv.svStartBeg) + (sv.svEndBeg - sv.svStartEnd) - 1;	  
  }
  predictedLength = finalGapEnd - finalGapStart;
  return false;
}

// Translocations
template<typename TUSize, typename TSize, typename TBamRecord>
inline bool
_translateSVCoordinates(TUSize refSize, TUSize, TUSize refSizeLeft, TSize diagBeg, TSize diagEnd, TUSize leftRefWalky, TUSize rightRefWalky, TUSize consLen, TBamRecord const& sv, TUSize& fGS, TUSize& fGE, TUSize& predictedLength, SVType<TranslocationTag>) {
  TSize finalGapStart = diagBeg + leftRefWalky + 1;
  TSize finalGapEnd = diagEnd + consLen - (consLen - rightRefWalky) + 1;
  predictedLength = sv.svEnd - sv.svStart + 1;  // Useless for translocations
  if (finalGapStart > (TSize) refSizeLeft) return true;
  if (finalGapEnd < (TSize) refSizeLeft) return true;
  if (sv.ct%2==0) {
    finalGapStart += sv.svStartBeg;
    if (sv.ct==2) finalGapEnd = sv.svEndBeg + (finalGapEnd - refSizeLeft) - 1;
    else finalGapEnd = sv.svEndBeg + (refSize - finalGapEnd) + 1;
  } else {
    finalGapStart = sv.svEndBeg + (refSizeLeft - finalGapStart);
    if (sv.ct==1) finalGapEnd = sv.svStartBeg + (finalGapEnd - refSizeLeft) - 1;
    else finalGapEnd = sv.svStartBeg + (refSize - finalGapEnd) + 1;
  }
  if (sv.ct%2!=0) {
    TUSize tmpStart = finalGapStart;
    finalGapStart = finalGapEnd;
    finalGapEnd = tmpStart;
  }
  if ((finalGapStart < 0) || (finalGapEnd<0)) return true;
  fGS=finalGapStart;
  fGE=finalGapEnd;
  return false;
}

// Insertions
template<typename TUSize, typename TSize, typename TBamRecord>
inline bool
_translateSVCoordinates(TUSize, TUSize, TUSize, TSize, TSize, TUSize, TUSize, TUSize, TBamRecord const&, TUSize&, TUSize&, TUSize&, SVType<InsertionTag>) {
  //ToDo
  return false;
}

// Translocations
template<typename TUSize, typename TStructuralVariantRecord>
inline void
_initRefSize(TUSize const refSize, TStructuralVariantRecord const& sv, TUSize& refSizeLeft, TUSize& refSizeRight, SVType<TranslocationTag>) {
  if (sv.ct%2==0) {
    refSizeRight = sv.svEndEnd - sv.svEndBeg;
    refSizeLeft = refSize - refSizeRight;
  } else {
    refSizeRight = sv.svStartEnd - sv.svStartBeg;
    refSizeLeft = refSize - refSizeRight;
  } 
}

// All other SVs
template<typename TUSize, typename TStructuralVariantRecord, typename TTag>
inline void
_initRefSize(TUSize const, TStructuralVariantRecord const&, TUSize&, TUSize&, SVType<TTag>) {
  //Nop
}

template<typename TConfig, typename TStructuralVariantRecord, typename TReadSet, typename TTag>
inline
void searchSplit(TConfig const& c, TStructuralVariantRecord& sv, std::string const& svRefStr, TReadSet const& splitReadSet, SVType<TTag> svType) {
  // Index SV reference
  typedef Index<int64_t, uint64_t,  char, 11, 4> TIndex;
  TIndex index;
  index.indexSequence(svRefStr);

  // Define kmerHitCutoff
  unsigned int kmerHitCutoff = 2;
  if (TIndex::kmer_size < c.minimumFlankSize) kmerHitCutoff = c.minimumFlankSize - TIndex::kmer_size + 1;

  // Define supposed length of split
  unsigned int initialLength = (sv.svEnd - sv.svStart + 1);
  unsigned int skippedLength = (sv.svEndBeg -sv.svStartEnd);
  unsigned int refSize = index.refSequence.size();
  unsigned int refSizeLeft = sv.svStartEnd - sv.svStartBeg;
  unsigned int refSizeRight = refSize - refSizeLeft;
  _initRefSize(refSize, sv, refSizeLeft, refSizeRight, svType);

  // Collect all potential split sites
  typedef std::vector<SplitReadCoord> TReadVec;
  TReadVec readVec;

  // Iterate over all reads
  typedef std::pair<int, unsigned int> TKmerDiagPos;
  typedef std::vector<TKmerDiagPos> TDiag;
  TDiag forward;
  TDiag reverse;
  TDiag forwardReverse;

  typedef std::vector<char> TSequence;
  typename TReadSet::const_iterator splitIter = splitReadSet.begin();
  for(;splitIter!=splitReadSet.end();++splitIter) {
    TSequence read;
    for(std::string::const_iterator readIt=splitIter->begin(); readIt!=splitIter->end(); ++readIt) read.push_back(dna5_encode[int(*readIt)]);
    unsigned int readLen = read.size();

    // Count kmer's 
    forward.clear();
    reverse.clear();
    index.diagonalForwardKmerCount(read, forward);
    index.diagonalReverseKmerCount(read, reverse);
    
    // No hits at all?
    bool alignDir = true;
    if (forward.empty() && reverse.empty()) continue;
    else {
      // Take better alignment
      if (forward.size()>reverse.size()) forwardReverse=forward;
      else {
	forwardReverse=reverse;
	alignDir=false;
      }
    }

    // Get the best split diagonals
    TDiag bestForwardReverse;
    sort(forwardReverse.begin(), forwardReverse.end());
    bool validForwardReverse = _findBestSupportDiagonal(c, forwardReverse, bestForwardReverse, kmerHitCutoff, readLen - TIndex::kmer_size + 1);

    // Get potential split-reads
    if (validForwardReverse) {
      std::sort(bestForwardReverse.begin(), bestForwardReverse.end()); // Sort diagonals
      typename TDiag::const_iterator itBest = bestForwardReverse.begin();
      typename TDiag::const_iterator itBestEnd = bestForwardReverse.end();
      for(;(itBest + 1)!=itBestEnd; ++itBest) {
	typename TDiag::const_iterator itBestNext = itBest + 1;
	for(;itBestNext !=itBestEnd; ++itBestNext) {
	  if (itBestNext->first - itBest->first >= 5) {   //Previously 20 for dups
	    SplitReadCoord sC;
	    sC.offset = (itBestNext->first - itBest->first);
	    sC.lastKmer = itBest->first + itBest->second;
	    sC.diag = itBest->first;
	    sC.read = read;
	    sC.forwardRead = alignDir;
	    readVec.push_back(sC);
	    //std::cout << sC.offset << ',' << sC.lastKmer << ',' << sC.diag << ',' << sC.forwardRead << ',' << *splitIter << ',' << readLen << std::endl;
	  }
	}
      }
    } 
  }

  // All potential split-reads collected, try to get the split
  if (!readVec.empty()) {
    // Sort by offset
    std::sort(readVec.begin(), readVec.end(),  SortSplitReadRecords<SplitReadCoord>());
    
    // Find best supported offset
    typename TReadVec::const_iterator readVecIt = readVec.begin();
    typename TReadVec::const_iterator readVecItEnd = readVec.end();
    int oldDiag = readVecIt->diag;
    int oldKmer = readVecIt->lastKmer;
    unsigned int oldOffset = readVecIt->offset;
    unsigned int support = 1;
    unsigned int bestSupport = 0;
    unsigned int bestBoundS = 0;
    unsigned int bestBoundE = 0;
    unsigned int bound = 1;
    unsigned int boundS = 0;
    int readLenInit = readVecIt->read.size();
    ++readVecIt;
    for(;readVecIt!=readVecItEnd;++readVecIt, ++bound) {
      // Read pairs should support the same offset and the last kmer should be sufficiently close (allow one mismatch in the last kmer)
      if ((readVecIt->offset == oldOffset) && ((readVecIt->lastKmer - oldKmer) <= TIndex::kmer_size)) {
	if (readVecIt->diag != oldDiag) ++support; // Count only unique reads (unique starting pos);
      } else {
	if (_betterSplit(refSize, refSizeRight, refSizeLeft, oldDiag, oldOffset, skippedLength, initialLength, c.epsilon, support, bestSupport, sv.ct, readVecIt->read.size(), svType)) {
	  bestSupport = support;
	  bestBoundS = boundS;
	  bestBoundE = bound;
	}
	boundS = bound;
	oldOffset = readVecIt->offset;
	oldDiag = readVecIt->diag;
	oldKmer = readVecIt->lastKmer;
	support = 1;
      }
    }
    if (_betterSplit(refSize, refSizeRight, refSizeLeft, oldDiag, oldOffset, skippedLength, initialLength, c.epsilon, support, bestSupport, sv.ct, readLenInit, svType)) {
      bestSupport = support;
      bestBoundS = boundS;
      bestBoundE = bound;
    }

    // At least the desired number of splitted reads
    if (bestSupport >= c.minimumSplitRead) {
      // Collect all reads for output
      typedef std::vector<OutputRead<TSequence> > TReadOut;
      TReadOut readOut;
      
      // Build the consensus of all split-reads
      typedef std::vector<std::pair<unsigned int, char> > TConsensus;
      TConsensus letters;
      readVecIt = readVec.begin() + bestBoundS;
      readVecItEnd = readVec.begin() + bestBoundE;
      int smallestDiag = readVecIt->diag;
      for(++readVecIt;readVecIt!=readVecItEnd;++readVecIt) {
	if (readVecIt->diag < smallestDiag) smallestDiag = readVecIt->diag;
      }
      readVecIt = readVec.begin() + bestBoundS;
      readVecItEnd = readVec.begin() + bestBoundE;
      for(;readVecIt!=readVecItEnd;++readVecIt) {
	unsigned int pos = readVecIt->diag - smallestDiag;
	OutputRead<TSequence> oR;
	oR.forwardRead = readVecIt->forwardRead;
	oR.readOffset = pos;
	oR.seq = readVecIt->read;
	readOut.push_back(oR);
	if (readVecIt->forwardRead) {
	  typename TSequence::const_iterator p = readVecIt->read.begin();
	  typename TSequence::const_iterator pEnd = readVecIt->read.end();
	  for(; p!=pEnd; ++p) letters.push_back(std::make_pair(pos++, *p));
	} else {
	  typename TSequence::const_iterator p = readVecIt->read.end();
	  typename TSequence::const_iterator pEnd = readVecIt->read.begin();
	  do {
	    --p;
	    unsigned int ch = 4;
	    switch((unsigned int) *p) {
	    case 0: ch = 3; break;
	    case 1: ch = 2; break;
	    case 2: ch = 1; break;
	    case 3: ch = 0; break;
	    }
	    letters.push_back(std::make_pair(pos++, ch));
	  } while (p!=pEnd);
	}
      }
      
      // Sort letters by position
      std::sort(letters.begin(), letters.end());
      typename TConsensus::const_iterator letIt = letters.begin();
      typename TConsensus::const_iterator letItEnd = letters.end();
      typedef std::vector<char> TConsSeq;
      TConsSeq consSeq;
      unsigned int oldPos = letIt->first;
      unsigned int consC[5];
      unsigned int belowCovThres = 0;
      for(unsigned int i = 0;i<5;++i) consC[i] = 0;
      ++consC[(int) letIt->second];
      for(++letIt;letIt!=letItEnd; ++letIt) {
	if (oldPos == letIt->first) ++consC[(int) letIt->second];
	else {
	  char maxI = 0;
	  unsigned int maxCount = consC[0];
	  unsigned int totalCount = consC[0];
	  for(int i = 1;i<5;++i) {
	    totalCount += consC[i];
	    if (consC[i] > maxCount) {
	      maxCount = consC[i];
	      maxI = i;
	    }
	  }
	  if (totalCount >= c.minimumSplitRead) consSeq.push_back(maxI);
	  else if (consSeq.empty()) ++belowCovThres;
	  oldPos = letIt->first;
	  for(unsigned int i = 0;i<5;++i) consC[i] = 0;
	  ++consC[(int) letIt->second];
	}
      }
      char maxI = 0;
      unsigned int maxCount = consC[0];
      unsigned int totalCount = consC[0];
      for(int i = 1;i<5;++i) {
	totalCount += consC[i];
	if (consC[i] > maxCount) {
	  maxCount = consC[i];
	  maxI = i;
	}
      }
      if (totalCount >= c.minimumSplitRead) consSeq.push_back(maxI);
      else if (consSeq.empty()) ++belowCovThres;

      // Debug output
      //typename TConsSeq::const_iterator csSeqIt = consSeq.begin();
      //typename TConsSeq::const_iterator csSeqItEnd = consSeq.end();
      //for(;csSeqIt != csSeqItEnd; ++csSeqIt) std::cout << dna5_decode[(unsigned int) *csSeqIt];
      //std::cout << std::endl;

      // Align the consensus
      TDiag final;
      index.diagonalForwardKmerCount(consSeq, final);
      unsigned int consLen = consSeq.size();

      // Collect the two best diagonals, maybe more for multiple splits?
      TDiag best;
      std::set<unsigned int> usedKmer;
      TDiag leftOver;
      while (final.size()) {
	sort(final.begin(), final.end());
	typename TDiag::const_iterator itD = final.begin();
	typename TDiag::const_iterator itBef = final.begin();
	typename TDiag::const_iterator itDEnd = final.end();
	
	// Count the diagonals 
	unsigned int currentCount =1;
	unsigned int bestCount =0;
	int bestDiag =0;
	for(++itD;itD!=itDEnd; ++itD, ++itBef) {
	  if (itBef->first == itD->first) ++currentCount;
	  else {
	    if (bestCount < currentCount) {
	      bestDiag = itBef->first;
	      bestCount = currentCount;
	    }
	    currentCount = 1;
	  }
	}
	if (bestCount < currentCount) {
	  bestDiag = itBef->first;
	  bestCount = currentCount;
	}

	// Add the diagonal or break if it is not good enough
	if (bestCount < kmerHitCutoff) break;

	// Collect used kmer positions
	usedKmer.clear();
	itD = final.begin();
	itDEnd = final.end();
	unsigned int lastSeenKmer = 0;
	for(;itD!=itDEnd; ++itD) {
	  if (bestDiag == itD->first) {
	    if (itD->second > lastSeenKmer) lastSeenKmer = itD->second;
	    usedKmer.insert(itD->second);
	  }
	}
	best.push_back(std::make_pair(bestDiag, lastSeenKmer));
	//std::cout << bestDiag << ',' << lastSeenKmer << '(' << bestCount << ')' << ';';
	if (best.size() == 2) break;

	// Keep the left-over diagonals
	leftOver.clear();
	itD = final.begin();
	itDEnd = final.end();
	for(;itD!=itDEnd; ++itD) if (usedKmer.find(itD->second) == usedKmer.end()) leftOver.push_back(*itD);
	final.swap(leftOver);
      }
      //std::cout << std::endl;

      // Compare to the reference   
      if (best.size() == 2) {
	TConsSeq ref1;
	TConsSeq ref2;

	int diagBeg = best.begin()->first;
	int diagEnd = (best.begin() + 1)->first;
	if (diagBeg > diagEnd) { int tmp = diagBeg; diagBeg = diagEnd; diagEnd = tmp; }
	
	typename TIndex::TSequence::const_iterator refSeqIt = index.refSequence.begin();
	typename TIndex::TSequence::const_iterator refSeqItEnd = index.refSequence.end();
	int seqSize = 0;
	for(;refSeqIt != refSeqItEnd; ++refSeqIt) {
	  if ((seqSize >= diagBeg) && (seqSize < (int) (diagBeg + consSeq.size()))) ref1.push_back(dna5_encode[(int) *refSeqIt]);
	  if ((seqSize >= diagEnd) && (seqSize < (int) (diagEnd + consSeq.size()))) ref2.push_back(dna5_encode[(int) *refSeqIt]);
	  ++seqSize;
	}

	// Debug output
	//if (sv.id==2) {
	//  for(typename TConsSeq::const_iterator csSeqIt = ref1.begin();csSeqIt != ref1.end(); ++csSeqIt) std::cout << dna5_decode[(unsigned int) *csSeqIt];
	//  std::cout << std::endl;
	//  for(typename TConsSeq::const_iterator csSeqIt = ref2.begin();csSeqIt != ref2.end(); ++csSeqIt) std::cout << dna5_decode[(unsigned int) *csSeqIt];
	//  std::cout << std::endl;
	//}

	// Calculate forward and reverse dynamic programming matrix
	int matchScore = 2;
	int penalty = -3;
	DnaScore<int> sc(matchScore,penalty,penalty,penalty);
	AlignConfig<true, false, false, false> alConf;
	typedef std::vector<int> TScoreMatrix;
	TScoreMatrix mat;
	globalNwAlignmentMatrix(ref1, consSeq, sc, alConf, mat);
	typedef std::pair<int, unsigned int> TScoreCol;
	typedef std::vector<TScoreCol>  TScoreColVector;
	TScoreColVector bestRowScore;
	for(unsigned int row = 0; row<=consSeq.size(); ++row) {
	  int maxRow = mat[row];
	  unsigned int colInd = 0;
	  for(unsigned int col = 1; col<=ref1.size(); ++col) {
	    if (mat[col * (consSeq.size() + 1) + row] > maxRow) {
	      maxRow = mat[col * (consSeq.size() + 1) + row];
	      colInd = col;
	    }
	  }
	  bestRowScore.push_back(std::make_pair(maxRow, colInd));
	}

	// Debug code
	//if (sv.id==2) {
	//  typedef FastaRecord<std::string, unsigned int, Dna5GapAlphabet, TSequence, void> TFastaRecord;
	//  std::vector<TFastaRecord> align;
	//  globalNwAlignment(align, ref1, consSeq, sc, alConf);
	//  std::cout << std::endl;
	//  for(TSequence::iterator alItTmp = align[0].seq.begin();alItTmp != align[0].seq.end(); ++alItTmp) std::cout << dna5gap_decode[(int) *alItTmp];
	//  std::cout << std::endl;
	//  for(TSequence::iterator alItTmp = align[1].seq.begin();alItTmp != align[1].seq.end(); ++alItTmp) std::cout << dna5gap_decode[(int) *alItTmp];
	//  std::cout << std::endl;
	//}

	reverseComplement(ref2);
	reverseComplement(consSeq);
	mat.clear();
	globalNwAlignmentMatrix(ref2, consSeq, sc, alConf, mat);
	TScoreColVector bestReverseRowScore;
	bestReverseRowScore.resize(bestRowScore.size());
	for(unsigned int row = 0; row<=consSeq.size(); ++row) {
	  int maxRow = mat[row];
	  unsigned int colInd = 0;
	  for(unsigned int col = 1; col<=ref2.size(); ++col) {
	    if (mat[col * (consSeq.size() + 1) + row] > maxRow) {
	      maxRow = mat[col * (consSeq.size() + 1) + row];
	      colInd = col;
	    }
	  }
	  bestReverseRowScore[bestRowScore.size() - (row + 1)] = std::make_pair(maxRow, (ref2.size() + 1) - (colInd + 1));
	}

	// Debug code
	//if (sv.id==2) {
	//  typedef FastaRecord<std::string, unsigned int, Dna5GapAlphabet, TSequence, void> TFastaRecord;
	//  std::vector<TFastaRecord> align;
	//  globalNwAlignment(align, ref2, consSeq, sc, alConf);
	//  std::cout << std::endl;
	//  for(TSequence::iterator alItTmp = align[0].seq.begin();alItTmp != align[0].seq.end(); ++alItTmp) std::cout << dna5gap_decode[(int) *alItTmp];
	//  std::cout << std::endl;
	//  for(TSequence::iterator alItTmp = align[1].seq.begin();alItTmp != align[1].seq.end(); ++alItTmp) std::cout << dna5gap_decode[(int) *alItTmp];
	//  std::cout << std::endl;
	//}

	// Get back to the true orientation
	reverseComplement(consSeq);
	reverseComplement(ref2);

	// Find the best alignment split (allowing for microinsertions)
	int maxScore = 0;
	unsigned int leftWalky = 0;
	unsigned int rightWalky = 0;
	unsigned int leftRefWalky = 0;
	unsigned int rightRefWalky = 0;
	TScoreColVector::iterator fIt = bestRowScore.begin();
	TScoreColVector::iterator fItEnd = bestRowScore.end();
	TScoreColVector::iterator rIt = bestReverseRowScore.begin();
	TScoreColVector::iterator rItEnd = bestReverseRowScore.end();
	for(;fIt!=fItEnd; ++fIt, ++rIt) {
	  for(TScoreColVector::iterator rItGo = rIt; rItGo != rItEnd; ++rItGo) {
	    if ((fIt->first + rItGo->first) > maxScore) {
	      maxScore = fIt->first + rItGo->first;
	      leftWalky = fIt - bestRowScore.begin();
	      leftRefWalky = fIt->second;
	      rightWalky = rItGo - bestReverseRowScore.begin() + 1;
	      rightRefWalky = rItGo->second + 1;
	    }
	  }
	}

	// Check split point
	bool invalidAlignment = false;
	if (leftWalky == 0) {
	  invalidAlignment = true;
	  leftWalky = consSeq.size() / 2;
	  rightWalky = leftWalky + 1;
	} else if (rightWalky > consSeq.size()) {
	  invalidAlignment = true;
	  leftWalky = consSeq.size() / 2;
	  rightWalky = leftWalky + 1;
	} else {
	  --leftWalky; --rightWalky;
	  --leftRefWalky; --rightRefWalky;
	}
	if (invalidAlignment) return;
	
	// Get the start and end of the structural variant
	unsigned int finalGapStart = 0;
	unsigned int finalGapEnd = 0;
	unsigned int predictedLength = 0;
	invalidAlignment=_translateSVCoordinates(refSize, refSizeRight, refSizeLeft, diagBeg, diagEnd, leftRefWalky, rightRefWalky, consLen, sv, finalGapStart, finalGapEnd, predictedLength, svType);

	// Count the final number of aligned reads
	typename TReadOut::const_iterator readIt = readOut.begin();
	typename TReadOut::const_iterator readItEnd = readOut.end();
	unsigned int readAlignCount = 0;
	for(;readIt!=readItEnd; ++readIt) {
	  if ((readIt->readOffset + readIt->seq.size() < (unsigned int) rightWalky + belowCovThres + TIndex::kmer_size) || (readIt->readOffset + TIndex::kmer_size > (unsigned int) leftWalky + belowCovThres)) continue;
	  else ++readAlignCount;
	}

	// Calculate quality
	double quality = (double) maxScore / (double) (matchScore * leftWalky + matchScore * (consLen - rightWalky + 1));

	// Valid breakpoint?
	if ((!invalidAlignment) && (readAlignCount >= c.minimumSplitRead) && (leftWalky >= c.minimumFlankSize) && ((consLen - rightWalky) >= c.minimumFlankSize) && (std::abs(((double) predictedLength / (double) initialLength) - 1.0) <= c.epsilon) && (quality >= (double) c.flankQuality / 100.0)) {
	  sv.precise=true;
	  sv.svStart=finalGapStart;
	  sv.svEnd=finalGapEnd;
	  sv.srSupport=readAlignCount;
	  sv.srAlignQuality=quality;
	  typename TConsSeq::const_iterator csSeqIt = consSeq.begin();
	  typename TConsSeq::const_iterator csSeqItEnd = consSeq.end();
	  for(;csSeqIt != csSeqItEnd; ++csSeqIt) sv.consensus += dna5_decode[(unsigned int) *csSeqIt];

	  // Output the consensus
	  /*
	  readIt = readOut.begin();
	  readItEnd = readOut.end();
	  for(;readIt!=readItEnd; ++readIt) {
	    if ( (readIt->readOffset + readIt->seq.size() < (unsigned int) rightWalky + belowCovThres + TIndex::kmer_size) || (readIt->readOffset + TIndex::kmer_size > (unsigned int) leftWalky + belowCovThres) ) continue;
	    unsigned int pos = 0;
	    for(unsigned int securityCheck = 0; (pos < readIt->readOffset) && (securityCheck<10000); ++pos, ++securityCheck) std::cout << ' ';
	    if (readIt->forwardRead) {
	      typename TSequence::const_iterator seqIt = readIt->seq.begin();
	      typename TSequence::const_iterator seqItEnd = readIt->seq.end();
	      for(;seqIt!=seqItEnd; ++seqIt, ++pos) {
		if (pos == rightWalky + belowCovThres) std::cout << "\t";
		std::cout << dna5_decode[(int) *seqIt];
		if (pos == leftWalky + belowCovThres) std::cout << "\t";
	      }
	    } else {
	      typename TSequence::const_iterator seqIt = readIt->seq.end();
	      typename TSequence::const_iterator seqItEnd = readIt->seq.begin();
	      do {
		if (pos == rightWalky + belowCovThres) std::cout << "\t";
		unsigned int ch = 4;
		switch ((unsigned int) *(--seqIt)) {
		case 0: ch = 3; break;
		case 1: ch = 2; break;
		case 2: ch = 1; break;
		case 3: ch = 0; break;
		}
		std::cout << dna5_decode[ch];
		//std::cout << (char) (32 + (int) dna5_decode[ch]);
		if (pos == leftWalky + belowCovThres) std::cout << "\t";
		++pos;
	      } while (seqIt!=seqItEnd);
	    }
	    std::cout << std::endl;
	  }
	  std::cout << "---------------------------------------------" << std::endl;
	  std::cout << sv.chr << "\t" << finalGapStart << "\t" << finalGapEnd << "\t" << (finalGapEnd - finalGapStart) << "\t" << readAlignCount << "\t" << quality << "\t>" << sv.id << "<\t" << std::endl;
	  std::cout << std::endl;
	  */

	  /*
	  // Debug code
	  typename TConsSeq::const_iterator refItT = ref1.begin();
	  typename TConsSeq::const_iterator refItTEnd = ref1.end();
	  std::cout << "PO: ";
	  for(unsigned int i = 0; i<consLen; ++i) {
	    std::cout << (i%10);
	  }
	  std::cout << std::endl;
	  std::cout << "RL: ";
	  for (;refItT!=refItTEnd; ++refItT) {
	    std::cout << dna5_decode[(int) *refItT];
	  }
	  std::cout << std::endl;
	  std::cout << "RR: ";
	  refItT = ref2.begin();
	  refItTEnd = ref2.end();
	  for (;refItT!=refItTEnd; ++refItT) {
	    std::cout << dna5_decode[(int) *refItT];
	  }
	  std::cout << std::endl;
	  std::cout << "CO: ";
	  refItT = consSeq.begin();
	  refItTEnd = consSeq.end();
	  for (;refItT!=refItTEnd; ++refItT) {
	    std::cout << dna5_decode[(int) *refItT];
	  }
	  std::cout << std::endl;  
	  std::cout << std::endl;
	  std::cout << "ML: ";
	  typename TIncSum::const_iterator leftIt = leftSc.begin();
	  typename TIncSum::const_iterator leftItEnd = leftSc.end();
	  for (;leftIt!=leftItEnd; ++leftIt) {
	    std::cout << (*leftIt % 10);
	  }
	  std::cout << std::endl;  
	  std::cout << "MR: ";
	  typename TIncSum::const_iterator rightIt = rightSc.begin();
	  typename TIncSum::const_iterator rightItEnd = rightSc.end();
	  for (;rightIt!=rightItEnd; ++rightIt) {
	    std::cout << (*rightIt % 10);
	  }
	  std::cout << std::endl;  
	  std::cout << "SU: ";
	  leftIt = leftSc.begin();
	  rightIt = rightSc.begin();
	  for (;rightIt!=rightItEnd; ++leftIt, ++rightIt) {
	    std::cout << ((*rightIt + *leftIt) % 10);
	  }
	  std::cout << std::endl;  
	  */
	}
      }
    }
  }
}


// Deletions
inline int
_decodeOrientation(std::string const, SVType<DeletionTag>) {
  return 1;
}

// Duplications
inline int
_decodeOrientation(std::string const, SVType<DuplicationTag>) {
  return 1;
}

// Inversions
inline int
_decodeOrientation(std::string const ct, SVType<InversionTag>) {
  if (ct=="3to3") return 1;
  else return 0;
}

// Translocations
inline int
_decodeOrientation(std::string const ct, SVType<TranslocationTag>) {
  if (ct=="3to3") return 0;
  else if (ct=="5to5") return 1;
  else if (ct=="3to5") return 2;
  else return 3;
}

// Insertions
inline int
_decodeOrientation(std::string const, SVType<InsertionTag>) {
  return 4;
}

// Deletions
inline std::string
_addID(SVType<DeletionTag>) {
  return "DEL";
}

// Duplications
inline std::string
_addID(SVType<DuplicationTag>) {
  return "DUP";
}

// Inversions
inline std::string
_addID(SVType<InversionTag>) {
  return "INV";
}

// Translocations
inline std::string
_addID(SVType<TranslocationTag>) {
  return "TRA";
}

// Insertion
inline std::string
_addID(SVType<InsertionTag>) {
  return "INS";
}

// Parse Delly vcf file
template<typename TConfig, typename TReferences, typename TSize, typename TStructuralVariantRecord, typename TTag>
inline void
vcfParse(TConfig const& c, TReferences const references, TSize const overallMaxISize, std::vector<TStructuralVariantRecord>& svs, SVType<TTag> svType)
{
  bool refPresent=false;
  if (boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome)) refPresent=true;
  std::ifstream vcfFile(c.vcffile.string().c_str(), std::ifstream::in);
  if (vcfFile.is_open()) {
    typedef boost::unordered_map<std::string, unsigned int> TMapChr;
    TMapChr mapChr;
    typename TReferences::const_iterator itRef = references.begin();
    for(unsigned int i = 0;itRef!=references.end();++itRef, ++i) mapChr[ itRef->RefName ] = i;
    while (vcfFile.good()) {
      std::string vcfLine;
      getline(vcfFile, vcfLine);
      if (vcfLine[0]=='#') continue;
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep("\t");
      Tokenizer tokens(vcfLine, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      if (tokIter!=tokens.end()) {
	std::string chr=*tokIter++;
	TMapChr::const_iterator mapChrIt = mapChr.find(chr);
	if (mapChrIt != mapChr.end()) {
	  if (tokIter!=tokens.end()) {
	    StructuralVariantRecord svRec;
	    svRec.chr = mapChrIt->second;
	    svRec.chr2 = mapChrIt->second;
	    svRec.svStart = boost::lexical_cast<int32_t>(*tokIter++);
	    std::string id = *tokIter++;
	    if (id.substr(0,3)!=_addID(svType)) continue;
	    svRec.id = parseSVid(id);
	    svRec.peSupport=0;
	    svRec.peMapQuality=0;
	    svRec.srSupport=0;
	    svRec.srAlignQuality=0;
	    svRec.wiggle = 0;
	    svRec.controlID = 0;
	    svRec.precise = false;
	    // Ignore ref, alt, qual and filter
	    tokIter++; tokIter++; tokIter++; tokIter++;
	    // Parse info string
	    std::string infoStr = *tokIter++;
	    boost::char_separator<char> sepInfo(";");
	    Tokenizer infoTokens(infoStr, sepInfo);
	    Tokenizer::iterator infoIter = infoTokens.begin();
	    for (;infoIter!=infoTokens.end();++infoIter) {
	      std::string keyValue = *infoIter;
	      std::size_t found = keyValue.find('=');
	      if (found==std::string::npos) {
		if ((keyValue=="PRECISE") && (refPresent)) svRec.precise=true;
		continue;
	      }
	      std::string key = keyValue.substr(0, found);
	      std::string value = keyValue.substr(found+1);
	      if (key == "PE") svRec.peSupport = boost::lexical_cast<int>(value);
	      else if (key == "MAPQ") svRec.peMapQuality = boost::lexical_cast<uint16_t>(value);
	      else if (key == "SR") svRec.srSupport = boost::lexical_cast<int>(value);
	      else if (key == "SRQ") svRec.srAlignQuality = boost::lexical_cast<double>(value);
	      else if (key == "CHR2") {
		TMapChr::const_iterator mapChr2It = mapChr.find(value);
		if (mapChrIt != mapChr.end()) svRec.chr2 = mapChr2It->second;
	      }
	      else if (key == "END") svRec.svEnd = boost::lexical_cast<int32_t>(value);
	      else if (key == "CONSENSUS") svRec.consensus = value;
	      else if (key == "CIPOS") {
		std::size_t foundComma = value.find(',');
		int off1 = boost::lexical_cast<int>(value.substr(0, foundComma));
		int off2 = boost::lexical_cast<int>(value.substr(foundComma+1));
		if (abs(off1)>svRec.wiggle) svRec.wiggle=abs(off1);
		if (abs(off2)>svRec.wiggle) svRec.wiggle=abs(off2);
	      }
	      else if (key == "CIEND") {
		std::size_t foundComma = value.find(',');
		int endOff1 = boost::lexical_cast<int>(value.substr(0, foundComma));
		int endOff2 = boost::lexical_cast<int>(value.substr(foundComma+1));
		if (abs(endOff1)>svRec.wiggle) svRec.wiggle=abs(endOff1);
		if (abs(endOff2)>svRec.wiggle) svRec.wiggle=abs(endOff2);
	      }
	      else if (key == "CT") svRec.ct = _decodeOrientation(value, svType);
	      else if (key == "CONTROL") svRec.controlID = boost::lexical_cast<unsigned int>(value);
	      else continue;
	    }
	    svRec.svStartBeg = std::max(svRec.svStart - 1 - overallMaxISize, 0);
	    svRec.svStartEnd = std::min(svRec.svStart - 1 + overallMaxISize, references[svRec.chr].RefLength);
	    svRec.svEndBeg = std::max(svRec.svEnd - 1 - overallMaxISize, 0);
	    svRec.svEndEnd = std::min(svRec.svEnd - 1 + overallMaxISize, references[svRec.chr2].RefLength);
	    if ((svRec.chr==svRec.chr2) && (svRec.svStartEnd > svRec.svEndBeg)) {
	      unsigned int midPointDel = ((svRec.svEnd - svRec.svStart) / 2) + svRec.svStart;
	      svRec.svStartEnd = midPointDel -1;
	      svRec.svEndBeg = midPointDel;
	    }
	    svs.push_back(svRec);
	  }
	}
      }
    }
    vcfFile.close();
  }
}


// Deletions
inline std::string
_addOrientation(int, SVType<DeletionTag>) {
  return "3to5";
}

// Duplications
inline std::string
_addOrientation(int, SVType<DuplicationTag>) {
  return "5to3";
}

// Inversions
inline std::string
_addOrientation(int ct, SVType<InversionTag>) {
  if (ct==1) return "3to3";
  else return "5to5";
}

// Translocations
inline std::string
_addOrientation(int ct, SVType<TranslocationTag>) {
  if (ct==0) return "3to3";
  else if (ct==1) return "5to5";
  else if (ct==2) return "3to5";
  else return "5to3";
}

// Insertion
inline std::string
_addOrientation(int, SVType<InsertionTag>) {
  return "NtoN";
}

// Insertion length
template<typename TSize, typename TTag>
inline TSize
_addInsertionLength(TSize, SVType<TTag>) {
  return 0;
}

// Insertion length
template<typename TSize>
inline TSize
_addInsertionLength(TSize l, SVType<InsertionTag>) {
  return l;
}


template<typename TConfig, typename TStructuralVariantRecord, typename TJunctionCountMap, typename TReadCountMap, typename TCountMap, typename TTag>
inline void
vcfOutput(TConfig const& c, std::vector<TStructuralVariantRecord> const& svs, TJunctionCountMap const& jctCountMap, TReadCountMap const& readCountMap, TCountMap const& spanCountMap, SVType<TTag> svType) 
{
  // Typedefs
  typedef typename TCountMap::key_type TSampleSVPair;
  typedef typename TCountMap::mapped_type TCountPair;
  typedef typename TCountPair::first_type TMapqVector;

  // Get the references
  BamTools::BamReader reader;
  reader.Open(c.files[0].string());
  BamTools::RefVector references = reader.GetReferenceData();

  // Output all structural variants
  std::ofstream ofile(c.outfile.string().c_str());

  // Print vcf header
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  boost::gregorian::date today = now.date();
  ofile << "##fileformat=VCFv4.1" << std::endl;
  ofile << "##fileDate=" << boost::gregorian::to_iso_string(today) << std::endl;
  ofile << "##ALT=<ID=DEL,Description=\"Deletion\">" << std::endl;
  ofile << "##ALT=<ID=DUP,Description=\"Duplication\">" << std::endl;
  ofile << "##ALT=<ID=INV,Description=\"Inversion\">" << std::endl;
  ofile << "##ALT=<ID=TRA,Description=\"Translocation\">" << std::endl;
  ofile << "##ALT=<ID=INS,Description=\"Insertion\">" << std::endl;
  ofile << "##FILTER=<ID=LowQual,Description=\"PE support below 3 or mapping quality below 20.\">" << std::endl;
  ofile << "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"PE confidence interval around END\">" << std::endl;
  ofile << "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"PE confidence interval around POS\">" << std::endl;
  ofile << "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">" << std::endl;
  ofile << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">" << std::endl;
  ofile << "##INFO=<ID=PE,Number=1,Type=Integer,Description=\"Paired-end support of the structural variant\">" << std::endl;
  ofile << "##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of paired-ends\">" << std::endl;
  ofile << "##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Split-read support\">" << std::endl;
  ofile << "##INFO=<ID=SRQ,Number=1,Type=Float,Description=\"Split-read consensus alignment quality\">" << std::endl;
  ofile << "##INFO=<ID=CONSENSUS,Number=1,Type=String,Description=\"Split-read consensus sequence\">" << std::endl;
  ofile << "##INFO=<ID=CT,Number=1,Type=String,Description=\"Paired-end signature induced connection type\">" << std::endl;
  ofile << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">" << std::endl;
  ofile << "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">" << std::endl;
  ofile << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << std::endl;
  ofile << "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">" << std::endl;
  ofile << "##INFO=<ID=CONTROL,Number=1,Type=Integer,Description=\"Control variant\">" << std::endl;
  ofile << "##INFO=<ID=INSLEN,Number=1,Type=Integer,Description=\"Predicted length of the insertion\">" << std::endl;
  ofile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
  ofile << "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Log10-scaled genotype likelihoods for RR,RA,AA genotypes\">" << std::endl;
  ofile << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" << std::endl;
  ofile << "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Per-sample genotype filter\">" << std::endl;
  ofile << "##FORMAT=<ID=RC,Number=1,Type=Integer,Description=\"Raw high-quality read counts for the SV\">" << std::endl;
  ofile << "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# high-quality reference pairs\">" << std::endl;
  ofile << "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# high-quality variant pairs\">" << std::endl;
  ofile << "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"# high-quality reference junction reads\">" << std::endl;
  ofile << "##FORMAT=<ID=RV,Number=1,Type=Integer,Description=\"# high-quality variant junction reads\">" << std::endl;
  ofile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    std::string sampleName(c.files[file_c].stem().string());
    ofile << "\t" << sampleName;
  }
  ofile << std::endl;

  // Iterate all structural variants
  typedef std::vector<TStructuralVariantRecord> TSVs;
  typename TSVs::const_iterator svIter = svs.begin();
  typename TSVs::const_iterator svIterEnd = svs.end();
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Genotyping" << std::endl;
  boost::progress_display show_progress( svs.size() );
  for(;svIter!=svIterEnd;++svIter) {
    ++show_progress;

    // Output main vcf fields
    std::string filterField="PASS";
    if ((svIter->peSupport < 3) || (svIter->peMapQuality < 20) || ( (svIter->chr != svIter->chr2) && (svIter->peSupport < 5) ) ) {
      filterField="LowQual";
    }
    
    std::stringstream id;
    id << _addID(svType) << std::setw(8) << std::setfill('0') << svIter->id;
    ofile << references[svIter->chr].RefName << "\t" << svIter->svStart << "\t" << id.str() << "\tN\t<" << _addID(svType) << ">\t.\t" <<  filterField << "\t";

    // Add info fields
    if (svIter->precise) ofile << "PRECISE;";
    else ofile << "IMPRECISE;";
    ofile << "CIEND=" << -svIter->wiggle << "," << svIter->wiggle << ";CIPOS=" << -svIter->wiggle << "," << svIter->wiggle << ";";
    ofile << "SVTYPE=" << _addID(svType) << ";";
    ofile << "SVMETHOD=EMBL.DELLYv" << dellyVersionNumber << ";";
    ofile << "CHR2=" << references[svIter->chr2].RefName << ";";
    ofile << "END=" << svIter->svEnd << ";";
    ofile << "CT=" << _addOrientation(svIter->ct, svType) << ";";
    ofile << "INSLEN=" << _addInsertionLength(svIter->insLen, svType) << ";";
    ofile << "PE=" << svIter->peSupport << ";";
    ofile << "MAPQ=" << svIter->peMapQuality;
    if (svIter->precise)  {
      ofile << ";SR=" << svIter->srSupport;
      ofile << ";SRQ=" << svIter->srAlignQuality;
      ofile << ";CONSENSUS=" << svIter->consensus;
    }
    if (svIter->controlID > 0) {
      ofile << ";CONTROL=" << svIter->controlID;
    }

    // Add genotype columns (right bp only across all samples)
    ofile << "\tGT:GL:GQ:FT:RC:DR:DV:RR:RV";
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      // Get the sample name
      std::string sampleName(c.files[file_c].stem().string());
      TSampleSVPair sampleSVPairLeft = std::make_pair(sampleName, svIter->id);
      TSampleSVPair sampleSVPairRight = std::make_pair(sampleName, -svIter->id);
      typename TJunctionCountMap::const_iterator jctCountMapIt=jctCountMap.find(sampleSVPairLeft);
      typename TCountMap::const_iterator spanLeftIt=spanCountMap.find(sampleSVPairLeft);
      typename TCountMap::const_iterator spanRightIt=spanCountMap.find(sampleSVPairRight);

      // Junction counts
      TMapqVector mapqRef = TMapqVector();
      TMapqVector mapqAlt = TMapqVector();
      if (svIter->precise) {
	// Genotyping for precise events uses junction read qualities
	mapqRef = jctCountMapIt->second.first;
	mapqAlt = jctCountMapIt->second.second;
      } else {
	// Genotyping for imprecise events uses spanning pair qualities, always use the minimum of left and right breakpoint
	if (spanLeftIt->second.first.size() <= spanRightIt->second.first.size()) mapqRef=spanLeftIt->second.first;
	else mapqRef=spanRightIt->second.first;
	if (spanLeftIt->second.second.size() <= spanRightIt->second.second.size()) mapqAlt=spanLeftIt->second.second;
	else mapqAlt=spanRightIt->second.second;
      }

      // Compute genotype likelihoods
      //typedef boost::multiprecision::cpp_dec_float_50 FLP;
      typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<25> > FLP;
      FLP gl[3];
      unsigned int glBest=0;
      unsigned int peDepth=mapqRef.size() + mapqAlt.size();
      FLP glBestVal=-std::numeric_limits<FLP>::infinity();
      FLP scaling = -FLP(peDepth) * boost::multiprecision::log10(FLP(2));
      for(unsigned int geno=0; geno<=2; ++geno) {
	FLP refLike=0;
	FLP altLike=0;
	typename TMapqVector::const_iterator mapqRefIt = mapqRef.begin();
	typename TMapqVector::const_iterator mapqRefItEnd = mapqRef.end();
	for(;mapqRefIt!=mapqRefItEnd;++mapqRefIt) {
	  refLike += boost::multiprecision::log10( (FLP(2.0 - geno) * boost::multiprecision::pow(FLP(10), -FLP(*mapqRefIt)/FLP(10)) + FLP(geno) * (FLP(1) - boost::multiprecision::pow(FLP(10), -FLP(*mapqRefIt)/FLP(10) ) ) ) );
	}
	typename TMapqVector::const_iterator mapqAltIt = mapqAlt.begin();
	typename TMapqVector::const_iterator mapqAltItEnd = mapqAlt.end();
	for(;mapqAltIt!=mapqAltItEnd;++mapqAltIt) {
	  altLike += boost::multiprecision::log10( ((FLP(2.0 - geno) * (FLP(1) - boost::multiprecision::pow(FLP(10), -FLP(*mapqAltIt)/FLP(10) ))) + FLP(geno) * boost::multiprecision::pow(FLP(10), -FLP(*mapqAltIt)/FLP(10) ) ) );
	}
	gl[geno]=scaling+refLike+altLike;
	if (gl[geno] > glBestVal) {
	  glBestVal=gl[geno];
	  glBest = geno;
	}
      }
      // Rescale by best genotype and get second best genotype for GQ
      FLP glSecondBestVal=-std::numeric_limits<FLP>::infinity();
      for(unsigned int geno=0; geno<=2; ++geno) {
	if ((gl[geno]>glSecondBestVal) && (gl[geno]<=glBestVal) && (geno!=glBest)) glSecondBestVal=gl[geno];
	gl[geno] -= glBestVal;
      }
      int gqVal = boost::multiprecision::iround(FLP(10) * boost::multiprecision::log10( boost::multiprecision::pow(FLP(10), glBestVal) / boost::multiprecision::pow(FLP(10), glSecondBestVal) ) );
      // Output genotypes
      if (peDepth) {
	if (glBest==0) ofile << "\t1/1:";
	else if (glBest==1) ofile << "\t0/1:";
	else ofile << "\t0/0:";
	ofile << gl[2] << "," << gl[1] << "," << gl[0] << ":" << gqVal << ":";
	if (gqVal<15) ofile << "LowQual:";
	else ofile << "PASS:";
      } else {
	ofile << "\t./.:.,.,.:0:LowQual:";
      }
      typename TReadCountMap::const_iterator readCountMapIt=readCountMap.find(sampleSVPairLeft);
      int rcCount = 0;
      if (readCountMapIt!=readCountMap.end()) rcCount = readCountMapIt->second.second;
      int drCount = 0;
      int dvCount = 0;
      if (spanLeftIt->second.first.size() <= spanRightIt->second.first.size()) drCount=spanLeftIt->second.first.size();
      else drCount=spanRightIt->second.first.size();
      if (spanLeftIt->second.second.size() <= spanRightIt->second.second.size()) dvCount=spanLeftIt->second.second.size();
      else dvCount=spanRightIt->second.second.size();
      int rrCount = 0;
      int rvCount = 0;
      if (jctCountMapIt!=jctCountMap.end()) {
	rrCount = jctCountMapIt->second.first.size();
	rvCount = jctCountMapIt->second.second.size();
      }
      ofile << rcCount << ":" << drCount << ":" << dvCount << ":" << rrCount << ":" << rvCount;
    }
    ofile << std::endl;
  }

  ofile.close();
}

template<typename TConfig, typename TStructuralVariantRecord>
inline bool
findPutativeSplitReads(TConfig const&, std::vector<TStructuralVariantRecord>&,  SVType<InsertionTag>) 
{
  return false;
}

template<typename TConfig, typename TStructuralVariantRecord, typename TTag>
inline bool
findPutativeSplitReads(TConfig const& c, std::vector<TStructuralVariantRecord>& svs,  SVType<TTag> svType) 
{
  typedef std::vector<TStructuralVariantRecord> TSVs;

  // Parse reference information
  BamTools::BamReader reader;
  if ( ! reader.Open(c.files[0].string())) return -1;
  BamTools::RefVector references = reader.GetReferenceData();

  // Parse genome, no single-anchored reads anymore only soft-clipped reads
  unsigned int totalSplitReadsAligned = 0;
  kseq_t *seq;
  int l;
  gzFile fp = gzopen(c.genome.string().c_str(), "r");
  seq = kseq_init(fp);
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Split-read alignment" << std::endl;
  boost::progress_display show_progress( references.size() );
  while ((l = kseq_read(seq)) >= 0) {
    // Find reference index
    BamTools::RefVector::const_iterator itRef = references.begin();
    for(int32_t refIndex=0;itRef!=references.end();++itRef, ++refIndex) {
      if (seq->name.s == references[refIndex].RefName) {
	++show_progress;

	// Iterate all structural variants on this chromosome
	typename TSVs::iterator svIt = svs.begin();
	typename TSVs::iterator svItEnd = svs.end();
	for(;svIt!=svItEnd; ++svIt) {
	  if ((svIt->chr != svIt->chr2) && (svIt->chr2 == refIndex)) {
	    // For translocations temporarily store the first reference part in the consensus string
	    svIt->consensus = _getSVRef(seq->seq.s, *svIt, refIndex, svType);
	  }
	  if (svIt->chr == refIndex) {
	    // Get the SV reference
	    std::string svRefStr = _getSVRef(seq->seq.s, *svIt, refIndex, svType);
	    svIt->consensus = "";
	    typedef std::set<std::string> TSplitReadSet;
	    TSplitReadSet splitReadSet;

	    // Find putative split reads in all samples
            #pragma omp parallel for default(shared)
	    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	      // Initialize bam file
	      BamTools::BamReader reader;
	      reader.Open(c.files[file_c].string());
	      reader.LocateIndex();

	      BamTools::BamAlignment al;
	      for (unsigned int bpPoint = 0; bpPoint<2; ++bpPoint) {
		int32_t regionChr = svIt->chr;
		int regionStart = (svIt->svStartBeg + svIt->svStart)/2;
		int regionEnd = (svIt->svStart + svIt->svStartEnd)/2;
		if (bpPoint == 1) {
		  regionChr = svIt->chr2;
		  regionStart = (svIt->svEndBeg + svIt->svEnd)/2;
		  regionEnd = (svIt->svEnd + svIt->svEndEnd)/2;
		}
		if (reader.SetRegion(regionChr, regionStart, regionChr, regionEnd)) {
		  while (reader.GetNextAlignmentCore(al)) {
		    if (al.RefID != refIndex) break;
		    if (!(al.AlignmentFlag & 0x0004) && !(al.AlignmentFlag & 0x0100) && !(al.AlignmentFlag & 0x0200) && !(al.AlignmentFlag & 0x0400) && !(al.AlignmentFlag & 0x0800)) {
		      // Clipped read
		      typedef std::vector<int> TIntVector;
		      TIntVector clipSizes;
		      TIntVector readPositions;
		      TIntVector genomePositions;
		      bool hasSoftClips = al.GetSoftClips(clipSizes, readPositions, genomePositions, false);
		      // Exactly one soft clip
		      if ((hasSoftClips) && (clipSizes.size()==1)) {
			// Get the sequence
			al.BuildCharData();
			  
			// Minimum clip size length fulfilled?
			int minClipSize = (int) (log10(al.QueryBases.size()) * 10);
#pragma omp critical
			{
			  if (clipSizes[0]>=minClipSize) splitReadSet.insert(al.QueryBases);
			}
		      }
		    }
		  }
		}
	      }
	    }
	    //std::cout << "Num. split reads: " << splitReadSet.size() << std::endl;

	    // Compare the split-reads against the SV reference
	    if (!splitReadSet.empty()) {
	      // Limit to at most 1000 split reads
	      std::vector<std::string> splitCandidates;
	      typename TSplitReadSet::const_iterator itSplit=splitReadSet.begin();
	      typename TSplitReadSet::const_iterator itSplitEnd=splitReadSet.end();
	      for(unsigned int splitCount=0; ((itSplit!=itSplitEnd) && (splitCount<1000));++itSplit, ++splitCount) {
		// Minimum read-length
		if (itSplit->size()>=35) splitCandidates.push_back(*itSplit);
	      }
	      
	      totalSplitReadsAligned += splitCandidates.size();

	      // Search true split in candidates
	      searchSplit(c, *svIt, svRefStr, splitCandidates, svType);
	    }
	  }
	}
      }
    }
  }
  kseq_destroy(seq);
  gzclose(fp);

  return (totalSplitReadsAligned>0);
}


// Initialize clique, deletions
template<typename TBamRecordIterator, typename TSize>
inline void
_initClique(TBamRecordIterator const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<DeletionTag>) {
  svStart = std::min(el->Position + el->Length, el->MatePosition + el->Length);
  svEnd = std::max(el->Position, el->MatePosition);
  wiggle =  abs(el->Position - el->MatePosition) + el->Length - el->maxNormalISize -(svEnd -svStart);
}

// Initialize clique, insertions
template<typename TBamRecordIterator, typename TSize>
inline void
_initClique(TBamRecordIterator const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<InsertionTag>) {
  svStart = std::min(el->Position + el->Length, el->MatePosition + el->Length);
  svEnd = std::max(el->Position, el->MatePosition);
  wiggle = -(svEnd - svStart);
}

// Initialize clique, duplications
template<typename TBamRecordIterator, typename TSize>
inline void
_initClique(TBamRecordIterator const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<DuplicationTag>) {
  svStart = std::min(el->Position, el->MatePosition);
  svEnd = std::max(el->Position + el->Length, el->MatePosition + el->Length);
  wiggle =  abs(el->Position - el->MatePosition) - el->Length + el->maxNormalISize -(svEnd -svStart);
}

// Initialize clique, inversions
template<typename TBamRecordIterator, typename TSize>
inline void
_initClique(TBamRecordIterator const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<InversionTag>) {
  int ct=_getSpanOrientation(*el, el->libOrient, SVType<InversionTag>());
  if (ct==1) {
    svStart = std::min(el->Position, el->MatePosition) + el->Length;
    svEnd = std::max(el->Position, el->MatePosition) + el->Length;
  } else {
    svStart = std::min(el->Position, el->MatePosition);
    svEnd = std::max(el->Position, el->MatePosition);
  }
  wiggle=el->maxNormalISize - el->Length;
}

// Initialize clique, translocations
template<typename TBamRecordIterator, typename TSize>
inline void
_initClique(TBamRecordIterator const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<TranslocationTag>) {
  int ct=_getSpanOrientation(*el, el->libOrient, SVType<TranslocationTag>());
  if (ct%2==0) {
    svStart = el->Position + el->Length;
    if (ct>=2) svEnd = el->MatePosition;
    else svEnd = el->MatePosition + el->Length;
  } else {
    svStart = el->Position;
    if (ct>=2) svEnd = el->MatePosition + el->Length;
    else svEnd = el->MatePosition;
  }
  wiggle=el->maxNormalISize - el->Length;
}


// Update clique, deletions
template<typename TBamRecordIterator, typename TSize>
inline bool 
_updateClique(TBamRecordIterator const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<DeletionTag>) 
{
  TSize newSvStart = std::max(svStart, std::min(el->Position + el->Length, el->MatePosition + el->Length));
  TSize newSvEnd = std::min(svEnd, std::max(el->Position, el->MatePosition));
  TSize newWiggle = abs(el->Position - el->MatePosition) + el->Length - el->maxNormalISize -(newSvEnd - newSvStart);
  TSize wiggleChange = wiggle + (svEnd-svStart) - (newSvEnd - newSvStart);
  if (wiggleChange > newWiggle) newWiggle=wiggleChange;

  // Does the new deletion size agree with all pairs
  if ((newSvStart < newSvEnd) && (newWiggle<=0)) {
    svStart = newSvStart;
    svEnd = newSvEnd;
    wiggle = newWiggle;
    return true;
  }
  return false;
}

// Update clique, insertions
template<typename TBamRecordIterator, typename TSize>
inline bool 
_updateClique(TBamRecordIterator const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<InsertionTag>) 
{
  TSize newSvStart = std::max(svStart, std::min(el->Position + el->Length, el->MatePosition + el->Length));
  TSize newSvEnd = std::min(svEnd, std::max(el->Position, el->MatePosition));
  TSize newWiggle = -(newSvEnd - newSvStart);

  // Does the new insertion size agree with all pairs
  if ((newSvStart < newSvEnd) && (newWiggle<=0)) {
    svStart = newSvStart;
    svEnd = newSvEnd;
    wiggle = newWiggle;
    return true;
  }
  return false;
}

// Update clique, duplications
template<typename TBamRecordIterator, typename TSize>
inline bool 
_updateClique(TBamRecordIterator const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<DuplicationTag>) 
{
  TSize newSvStart = std::min(svStart, std::min(el->Position, el->MatePosition));
  TSize newSvEnd = std::max(svEnd, std::max(el->Position + el->Length, el->MatePosition + el->Length));
  TSize newWiggle = abs(el->Position - el->MatePosition) - el->Length + el->maxNormalISize -(newSvEnd - newSvStart);
  TSize wiggleChange = wiggle - ((newSvEnd - newSvStart) - (svEnd-svStart));
  if (wiggleChange < newWiggle) newWiggle=wiggleChange;

  // Does the new duplication size agree with all pairs
  if ((newSvStart < newSvEnd) && (newWiggle>=0)) {
    svStart = newSvStart;
    svEnd = newSvEnd;
    wiggle = newWiggle;
    return true;
  }
  return false;
}

// Update clique, inversions
template<typename TBamRecordIterator, typename TSize>
inline bool 
_updateClique(TBamRecordIterator const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<InversionTag>) 
{
  int ct=_getSpanOrientation(*el, el->libOrient, SVType<InversionTag>());
  TSize newSvStart;
  TSize newSvEnd;
  TSize newWiggle;
  TSize wiggleChange;
  if (ct==1) {
    newSvStart = std::max(svStart, std::min(el->Position, el->MatePosition) + el->Length);
    newSvEnd = std::max(svEnd, std::max(el->Position, el->MatePosition) + el->Length);
    newWiggle = std::min(el->maxNormalISize - (newSvStart - std::min(el->Position, el->MatePosition)), el->maxNormalISize - (newSvEnd - std::max(el->Position, el->MatePosition)));
    wiggleChange = wiggle - std::max(newSvStart - svStart, newSvEnd - svEnd);
  } else {
    newSvStart = std::min(svStart, std::min(el->Position, el->MatePosition));
    newSvEnd = std::min(svEnd, std::max(el->Position, el->MatePosition));
    newWiggle = std::min(el->maxNormalISize - (std::min(el->Position, el->MatePosition) + el->Length - newSvStart), el->maxNormalISize - (std::max(el->Position, el->MatePosition) + el->Length - newSvEnd));
    wiggleChange = wiggle - std::max(svStart - newSvStart, svEnd - newSvEnd);
  }
  if (wiggleChange < newWiggle) newWiggle=wiggleChange;

  // Does the new inversion size agree with all pairs
  if ((newSvStart < newSvEnd) && (newWiggle>=0)) {
    svStart = newSvStart;
    svEnd = newSvEnd;
    wiggle = newWiggle;
    return true;
  }
  return false;
}


// Update clique, translocations
template<typename TBamRecordIterator, typename TSize>
inline bool 
_updateClique(TBamRecordIterator const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<TranslocationTag>) 
{
  int ct = _getSpanOrientation(*el, el->libOrient, SVType<TranslocationTag>());
  TSize newSvStart;
  TSize newSvEnd;
  TSize newWiggle = wiggle;
  if (ct%2==0) {
    newSvStart = std::max(svStart, el->Position + el->Length);
    newWiggle -= (newSvStart - svStart);
    if (ct>=2) {
      newSvEnd = std::min(svEnd, el->MatePosition);
      newWiggle -= (svEnd - newSvEnd);
    } else  {
      newSvEnd = std::max(svEnd, el->MatePosition + el->Length);
      newWiggle -= (newSvEnd - svEnd);
    }
  } else {
    newSvStart = std::min(svStart, el->Position);
    newWiggle -= (svStart - newSvStart);
    if (ct>=2) {
      newSvEnd = std::max(svEnd, el->MatePosition + el->Length);
      newWiggle -= (newSvEnd - svEnd);
    } else {
      newSvEnd = std::min(svEnd, el->MatePosition);
      newWiggle -= (svEnd - newSvEnd);
    }
  }
  // Is this still a valid translocation cluster?
  if (newWiggle>0) {
    svStart = newSvStart;
    svEnd = newSvEnd;
    wiggle = newWiggle;
    return true;
  }
  return false;
}

template<typename TSize>
inline bool
_svSizeCheck(TSize const s, TSize const e, SVType<DeletionTag>) {
  return (( e - s ) >= 300);
}

template<typename TSize>
inline bool
_svSizeCheck(TSize const s, TSize const e, SVType<DuplicationTag>) {
  return (( e - s ) >= 100);
}

template<typename TSize>
inline bool
_svSizeCheck(TSize const s, TSize const e, SVType<InversionTag>) {
  return (( e - s ) >= 100);
}

template<typename TSize>
inline bool
_svSizeCheck(TSize const s, TSize const e, SVType<InsertionTag>) {
  return (( e - s ) >= 0);
}

template<typename TSize>
inline bool
_svSizeCheck(TSize const, TSize const, SVType<TranslocationTag>) {
  return true;
}


template<typename TConfig, typename TSampleLibrary, typename TSVs, typename TCountMap, typename TTag>
inline void
_annotateJunctionReads(TConfig const& c, TSampleLibrary& sampleLib, TSVs& svs, TCountMap& junctionCountMap, SVType<TTag> svType) 
{
  annotateJunctionReads(c.files, c.genome, c.minGenoQual, sampleLib, svs, junctionCountMap, svType);
}


template<typename TConfig, typename TSampleLibrary, typename TSVs, typename TCountMap, typename TTag>
inline void
_annotateCoverage(TConfig const& c, TSampleLibrary& sampleLib, TSVs& svs, TCountMap& countMap, SVType<TTag>) 
{
  annotateCoverage(c.files, c.minGenoQual, sampleLib, svs, countMap, BpLevelType<NoBpLevelCount>(), CoverageType<RedundancyFilterTag>());
}

template<typename TConfig, typename TSampleLibrary, typename TSVs, typename TCountMap>
inline void
_annotateCoverage(TConfig const&, TSampleLibrary&, TSVs&, TCountMap&, SVType<TranslocationTag>) 
{
  //Nop
}

template<typename TConfig, typename TSampleLibrary, typename TSVs, typename TCountMap, typename TTag>
inline void
_annotateSpanningCoverage(TConfig const& c, TSampleLibrary& sampleLib, TSVs& svs, TCountMap& spanCountMap, SVType<TTag> svType) 
{
  annotateSpanningCoverage(c.files, c.minGenoQual, sampleLib, svs, spanCountMap, svType);
}


template<typename TSVType>
inline int run(Config const& c, TSVType svType) {

#ifdef PROFILE
  ProfilerStart("delly.prof");
#endif

  // Collect all promising structural variants
  typedef std::vector<StructuralVariantRecord> TVariants;
  TVariants svs;

  // Clique id counter
  unsigned int clique_count = 1;

  // Create library objects
  typedef boost::unordered_map<std::string, LibraryInfo> TLibraryMap;
  typedef boost::unordered_map<std::string, TLibraryMap> TSampleLibrary;
  TSampleLibrary sampleLib;
  int overallMaxISize = 0;

  // Check that all input bam files exist and are indexed
  BamTools::RefVector references;
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    BamTools::BamReader reader;
    if ( ! reader.Open(c.files[file_c].string()) ) {
      std::cerr << "Could not open input bam file: " << c.files[file_c].string() << std::endl;
      reader.Close();
      return -1;
    }
    reader.LocateIndex();
    if ( !reader.HasIndex() ) {
      std::cerr << "Missing bam index file: " << c.files[file_c].string() << std::endl;
      reader.Close();
      return -1;
    }

    // Get references
    if (file_c==0) references = reader.GetReferenceData();
  }

  // Get library parameters
#pragma omp parallel for default(shared)
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    // Get a sample name
    std::string sampleName(c.files[file_c].stem().string());

    // Get library parameters and overall maximum insert size
    TLibraryMap libInfo;
    getLibraryParams(c.files[file_c], libInfo, c.percentAbnormal, c.madCutoff);
#pragma omp critical
    {
      TLibraryMap::const_iterator libIter=libInfo.begin();
      for(;libIter!=libInfo.end();++libIter) 
	if (libIter->second.maxNormalISize > overallMaxISize) overallMaxISize = libIter->second.maxNormalISize;
      sampleLib.insert(std::make_pair(sampleName, libInfo));
    }
  }

  // Dump PE
  bool dumpPe=false;
  std::ofstream dumpPeFile;
  if (c.peDump.string().size()) {
    dumpPe=true;
    dumpPeFile.open(c.peDump.string().c_str());
    dumpPeFile << "#id\tchr\tpos\tmatechr\tmatepos\tmapq" << std::endl;
  }

  // Parse exclude interval list
  BamTools::RefVector::const_iterator itRef = references.begin();
  std::vector<bool> validChr;
  typedef std::vector<ExcludeInterval> TExclInterval;
  TExclInterval exclIntervals;
  validChr.resize(references.size());
  std::fill(validChr.begin(), validChr.end(), true);
  if (boost::filesystem::exists(c.exclude) && boost::filesystem::is_regular_file(c.exclude) && boost::filesystem::file_size(c.exclude)) {
    typedef boost::unordered_map<std::string, unsigned int> TMapChr;
    TMapChr mapChr;
    for(unsigned int i = 0;itRef!=references.end();++itRef, ++i) mapChr[ itRef->RefName ] = i;
    std::ifstream chrFile(c.exclude.string().c_str(), std::ifstream::in);
    if (chrFile.is_open()) {
      while (chrFile.good()) {
	std::string chrFromFile;
	getline(chrFile, chrFromFile);
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t,;");
	Tokenizer tokens(chrFromFile, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter!=tokens.end()) {
	  std::string chr=*tokIter++;
	  TMapChr::const_iterator mapChrIt = mapChr.find(chr);
	  if (mapChrIt != mapChr.end()) {
	    if (tokIter!=tokens.end()) {
	      int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
	      int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
	      exclIntervals.push_back(ExcludeInterval(mapChrIt->second, start, end));
	    } else validChr[mapChrIt->second]=false; // Exclude entire chromosome
	  }
	}
      }
      chrFile.close();
    }
  }
  std::sort(exclIntervals.begin(), exclIntervals.end(), SortExcludeIntervals<ExcludeInterval>());

  // Do we have an input vcffile 
  bool peMapping=true;
  if (boost::filesystem::exists(c.vcffile) && boost::filesystem::is_regular_file(c.vcffile) && boost::filesystem::file_size(c.vcffile)) peMapping=false;

  // Qualities
  typedef boost::unordered_map<unsigned int, uint16_t> TQualities;
  std::vector<TQualities> qualities;
  qualities.resize(c.files.size());

  // Process chromosome by chromosome
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Paired-end clustering" << std::endl;
  boost::progress_display show_progress( (references.end() - references.begin()) );
  itRef = references.begin();
  for(int refIndex=0; (itRef!=references.end()) && (peMapping); ++itRef, ++refIndex) {
    ++show_progress;
    if (!validChr[refIndex]) continue;
      
    // Create bam alignment record vector
    typedef std::vector<BamAlignRecord> TBamRecord;
    TBamRecord bamRecord;

    // Iterate all samples
    #pragma omp parallel for default(shared)
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      // Get a sample name
      std::string sampleName(c.files[file_c].stem().string());
      TSampleLibrary::iterator sampleIt=sampleLib.find(sampleName);

      // Initialize bam file
      BamTools::BamReader reader;
      reader.Open(c.files[file_c].string());
      reader.LocateIndex();

      // Unique pairs for the given sample
      typedef boost::container::flat_set<int32_t> TUniquePairs;
      TUniquePairs unique_pairs;

      // Read alignments
      int32_t oldAlignPos=-1;
      BamTools::BamAlignment al;
      if ( reader.Jump(refIndex, 0) ) {
	while( reader.GetNextAlignmentCore(al) ) {
	  if (al.RefID!=refIndex) break; // Stop when we hit the next chromosome
	  if ((al.AlignmentFlag & 0x0001) && !(al.AlignmentFlag & 0x0004) && !(al.AlignmentFlag & 0x0008) && !(al.AlignmentFlag & 0x0100) && !(al.AlignmentFlag & 0x0200) && !(al.AlignmentFlag & 0x0400) && !(al.AlignmentFlag & 0x0800) && (al.MapQuality >= c.minMapQual) && (al.RefID>=0) && (al.MateRefID>=0)) {
	    // Mapping positions valid?
	    if (_mappingPos(al.RefID, al.MateRefID, al.Position, al.MatePosition, svType)) continue;

	    // Is the read or its mate in a black-masked region
	    if (!exclIntervals.empty()) {
	      typename TExclInterval::const_iterator itBlackMask = std::lower_bound(exclIntervals.begin(), exclIntervals.end(), ExcludeInterval(al.RefID, al.Position, 0), SortExcludeIntervals<ExcludeInterval>());
	      if (itBlackMask!=exclIntervals.begin()) --itBlackMask;
	      if ((itBlackMask->RefID==al.RefID) && (itBlackMask->start <= al.Position) && (al.Position<=itBlackMask->end)) continue;
	      itBlackMask = std::lower_bound(exclIntervals.begin(), exclIntervals.end(), ExcludeInterval(al.MateRefID, al.MatePosition, 0), SortExcludeIntervals<ExcludeInterval>());
	      if (itBlackMask!=exclIntervals.begin()) --itBlackMask;
	      if ((itBlackMask->RefID==al.MateRefID) && (itBlackMask->start <= al.MatePosition) && (al.MatePosition<=itBlackMask->end)) continue;
	    }

	    // Is this a discordantly mapped paired-end?
	    al.BuildCharData();
	    std::string rG = "DefaultLib";
	    al.GetTag("RG", rG);
	    TLibraryMap::iterator libIt=sampleIt->second.find(rG);
	    if (libIt==sampleIt->second.end()) std::cerr << "Missing read group: " << rG << std::endl;
	    if (libIt->second.median == 0) continue; // Single-end library
	    if (_acceptedInsertSize(libIt->second, abs(al.InsertSize), svType)) continue;  // Normal paired-end (for deletions/insertions only)
	    if (_acceptedOrientation(libIt->second.defaultOrient, getStrandIndependentOrientation(al), svType)) continue;  // Orientation disagrees with SV type

	    // Get or store the mapping quality for the partner
	    if (_firstPairObs(al.RefID, al.MateRefID, al.Position, al.MatePosition, svType)) {
	      // Hash the quality
	      unsigned int index=((al.Position % (int)boost::math::pow<14>(2))<<14) + (al.MatePosition % (int)boost::math::pow<14>(2));
	      qualities[file_c][index]=al.MapQuality;
	    } else {
	      // Get the two mapping qualities
	      unsigned int index=((al.MatePosition % (int)boost::math::pow<14>(2))<<14) + (al.Position % (int)boost::math::pow<14>(2));
	      uint16_t pairQuality = std::min(qualities[file_c][index], al.MapQuality);
	      qualities[file_c][index]=0;

	      // Pair quality
	      if (pairQuality < c.minMapQual) continue;

	      // Store the paired-end
	      if (al.Position!=oldAlignPos) {
		oldAlignPos=al.Position;
		unique_pairs.clear();
	      }
	      if (unique_pairs.insert(al.MatePosition).second) {
#pragma omp critical
		{
		  bamRecord.push_back(BamAlignRecord(al, pairQuality, libIt->second.median, libIt->second.mad, libIt->second.maxNormalISize, libIt->second.defaultOrient));
		}
		++libIt->second.abnormal_pairs;
	      } else {
		++libIt->second.non_unique_abnormal_pairs;
	      }
	    }
	  }
	}
      }

      // Clean-up qualities
      _resetQualities(qualities[file_c], svType);
    }
    
    // Sort BAM records according to position
    std::sort(bamRecord.begin(), bamRecord.end(), SortBamRecords<BamAlignRecord>());
    //for(TBamRecord::const_iterator bamIt = bamRecord.begin(); bamIt!=bamRecord.end(); ++bamIt) std::cerr << bamIt->Position << ',' << bamIt->MatePosition << ',' << bamIt->maxNormalISize << ',' << bamIt->Length << ',' << bamIt->Median << ',' << bamIt->Mad << ',' << bamIt->libOrient << ',' << bamIt->AlignmentFlag << std::endl;

    // Define an undirected graph g
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, TBamRecord::const_iterator, boost::property<boost::edge_weight_t, unsigned short> > Graph;
    Graph g;
      
    // Define the reverse map
    typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef boost::unordered_map<unsigned int, Vertex> TNameVertexMap;
    TNameVertexMap nameFrag;
      
    // Define the edge property map
    typedef boost::property_map<Graph, boost::edge_weight_t>::type edge_map_t;
    edge_map_t weightMap = get(boost::edge_weight, g);
      
    // Iterate the chromosome range
    TBamRecord::const_iterator vecBeg = bamRecord.begin();
    TBamRecord::const_iterator vecEnd = bamRecord.end();
    for(;vecBeg!=vecEnd; ++vecBeg) {
      int32_t const minCoord = _minCoord(vecBeg->Position, vecBeg->MatePosition, svType);
      int32_t const maxCoord = _maxCoord(vecBeg->Position, vecBeg->MatePosition, svType);
      TBamRecord::const_iterator vecNext = vecBeg + 1;
      for(; ((vecNext != vecEnd) && (abs(_minCoord(vecNext->Position, vecNext->MatePosition, svType) + vecNext->Length - minCoord) <= overallMaxISize)) ; ++vecNext) {
	// Check that mate chr agree (only for translocations)
	if (vecBeg->MateRefID!=vecNext->MateRefID) continue;

	// Check combinability of pairs
	if (_pairsDisagree(minCoord, maxCoord, vecBeg->Length, vecBeg->maxNormalISize, _minCoord(vecNext->Position, vecNext->MatePosition, svType), _maxCoord(vecNext->Position, vecNext->MatePosition, svType), vecNext->Length, vecNext->maxNormalISize, _getSpanOrientation(*vecBeg, vecBeg->libOrient, svType), _getSpanOrientation(*vecNext, vecNext->libOrient, svType), svType)) continue;
	
	TNameVertexMap::iterator pos;
	bool inserted;
	
	// Add vertex 1
	Vertex u;
	boost::tie(pos, inserted) = nameFrag.insert(TNameVertexMap::value_type(vecBeg-bamRecord.begin(), Vertex()));
	if (inserted) {
	  u = add_vertex(g);
	  pos->second = u;
	  g[u]=vecBeg;
	} else {
	  u = pos->second;
	}
	
	// Add vertex 2
	Vertex v;
	boost::tie(pos, inserted) = nameFrag.insert(TNameVertexMap::value_type(vecNext-bamRecord.begin(), Vertex()));
	if (inserted) {
	  v = add_vertex(g);
	  pos->second = v;
	  g[v]=vecNext;
	} else {
	  v = pos->second;
	}
	
	// Add the edge
	if ((out_degree(u, g) <= c.graphPruning) || (out_degree(v, g) <= c.graphPruning)) {
	  boost::graph_traits<Graph>::edge_descriptor e;
	  tie(e, inserted) = add_edge(u,v, g);
	  if (inserted) {
	    int locEdgeWeight=abs( abs( (_minCoord(vecNext->Position, vecNext->MatePosition, svType) - minCoord) - (_maxCoord(vecNext->Position, vecNext->MatePosition, svType) - maxCoord) ) - abs(vecBeg->Median - vecNext->Median) );
	    weightMap[e] = (locEdgeWeight > 30000) ? 30000 : locEdgeWeight;
	  }
	}
      }
    }
    nameFrag.clear();
  
    // Compute the connected components
    std::vector<int> my_comp(num_vertices(g));
    int numComp = boost::connected_components(g, &my_comp[0]);
    
    // Count the number of vertices for each component
    typedef std::vector<unsigned int> TCompSize;
    TCompSize compSize(numComp);
    std::fill(compSize.begin(), compSize.end(), 0);
    boost::graph_traits<Graph>::vertex_iterator vIt, vItEnd;
    for(boost::tie(vIt, vItEnd) = boost::vertices(g); vIt != vItEnd; ++vIt) ++compSize[my_comp[*vIt]];
    //for(TCompSize::const_iterator compIt = compSize.begin(); compIt!=compSize.end(); ++compIt) std::cerr << *compIt << std::endl;

    // Iterate each component
#pragma omp parallel for default(shared)
    for(int compIt = 0; compIt < numComp; ++compIt) {
      if (compSize[compIt]<2) continue;
      typedef boost::graph_traits<Graph>::vertex_descriptor TVertexDescriptor;
      typedef EdgeRecord<unsigned short, TVertexDescriptor> TEdgeRecord;
      typedef std::vector<TEdgeRecord> TWeightEdge;
      TWeightEdge wEdge;
      boost::graph_traits<Graph>::edge_iterator eIt, eItEnd;
      for(boost::tie(eIt, eItEnd) = boost::edges(g); eIt != eItEnd; ++eIt) {
	if ((my_comp[boost::source(*eIt, g)] == compIt) && (my_comp[boost::source(*eIt, g)] == my_comp[boost::target(*eIt,g)])) {
	  wEdge.push_back(TEdgeRecord(boost::source(*eIt, g), boost::target(*eIt, g), weightMap[*eIt]));
	}
      }

      // Sort edges by weight
      std::sort(wEdge.begin(), wEdge.end(), SortEdgeRecords<TEdgeRecord>());
      //for(TWeightEdge::const_iterator itWEdge = wEdge.begin(); itWEdge!=wEdge.end(); ++itWEdge) std::cerr << refIndex << ',' << compIt << ',' << itWEdge->source << ',' << itWEdge->target << ',' << itWEdge->weight << std::endl;
      
      // Find a large clique
      TWeightEdge::const_iterator itWEdge = wEdge.begin();
      TWeightEdge::const_iterator itWEdgeEnd = wEdge.end();
      typedef std::set<TVertexDescriptor> TCliqueMembers;
      TCliqueMembers clique;
      TCliqueMembers incompatible;
      int svStart, svEnd, wiggle;
      int32_t clusterRefID=g[itWEdge->source]->RefID;
      int32_t clusterMateRefID=g[itWEdge->source]->MateRefID;
      _initClique(g[itWEdge->source], svStart, svEnd, wiggle, svType);
      int connectionType = _getSpanOrientation(*g[itWEdge->source], g[itWEdge->source]->libOrient, svType);
      if ((clusterRefID==clusterMateRefID) && (svStart >= svEnd))  continue;
      clique.insert(itWEdge->source);
      
      // Grow the clique from the seeding edge
      bool cliqueGrow=true;
      while ((cliqueGrow) && (clique.size() < compSize[compIt])) {
	itWEdge = wEdge.begin();
	cliqueGrow = false;
	for(;(!cliqueGrow) && (itWEdge != itWEdgeEnd);++itWEdge) {
	  TVertexDescriptor v;
	  if ((clique.find(itWEdge->source) == clique.end()) && (clique.find(itWEdge->target) != clique.end())) v = itWEdge->source;
	  else if ((clique.find(itWEdge->source) != clique.end()) && (clique.find(itWEdge->target) == clique.end())) v = itWEdge->target;
	  else continue;
	  if (incompatible.find(v) != incompatible.end()) continue;
	  boost::graph_traits<Graph>::adjacency_iterator vi, vi_end;
	  unsigned int cliqSize = 0;
	  for(boost::tie(vi, vi_end) = boost::adjacent_vertices(v, g); vi != vi_end; ++vi)
	    if (clique.find(*vi) != clique.end()) ++cliqSize;
	  if (cliqSize == clique.size()) {
	    //std::cerr << refIndex << ',' << compIt << ',' << v << ',' << svStart << ',' << svEnd << ',' << wiggle << std::endl;
	    cliqueGrow = _updateClique(g[v], svStart, svEnd, wiggle, svType);
	    if (cliqueGrow) clique.insert(v);
	    else incompatible.insert(v);
	  }
	}
      }

      if ((clique.size()>1) && (_svSizeCheck(svStart, svEnd, svType))) {
	StructuralVariantRecord svRec;
	svRec.chr = refIndex;
	svRec.chr2 = clusterMateRefID;
	svRec.svStartBeg = std::max((int) svStart - overallMaxISize, 0);
	svRec.svStart = std::min(svStart + 1, references[refIndex].RefLength);
	svRec.svStartEnd = std::min(svStart + overallMaxISize, references[refIndex].RefLength);
	svRec.svEndBeg = std::max((int) svEnd - overallMaxISize, 0);
	svRec.svEnd = std::min(svEnd+1, references[clusterMateRefID].RefLength);
	svRec.svEndEnd = std::min(svEnd + overallMaxISize, references[clusterMateRefID].RefLength);
	svRec.peSupport = clique.size();
	svRec.wiggle = abs(wiggle);
	svRec.controlID = 0;
	std::vector<uint16_t> mapQV;
	for(typename TCliqueMembers::const_iterator itC = clique.begin(); itC!=clique.end(); ++itC) mapQV.push_back(g[*itC]->MapQuality);
	std::sort(mapQV.begin(), mapQV.end());
	svRec.peMapQuality = mapQV[mapQV.size()/2];
	if ((refIndex==clusterMateRefID) && (svRec.svStartEnd > svRec.svEndBeg)) {
	  unsigned int midPointDel = ((svRec.svEnd - svRec.svStart) / 2) + svRec.svStart;
	  svRec.svStartEnd = midPointDel -1;
	  svRec.svEndBeg = midPointDel;
	}
	svRec.srSupport=0;
	svRec.srAlignQuality=0;
	svRec.precise=false;
	svRec.ct=connectionType;
	std::vector<int32_t> inslenV;
	for(typename TCliqueMembers::const_iterator itC = clique.begin(); itC!=clique.end(); ++itC) inslenV.push_back(g[*itC]->Median - (abs(g[*itC]->Position - g[*itC]->MatePosition) + g[*itC]->Length));
	std::sort(inslenV.begin(), inslenV.end());
	svRec.insLen = inslenV[inslenV.size()/2];
#pragma omp critical
	{
	  svRec.id = clique_count++;
	  svs.push_back(svRec);
	}
	
	// Dump PEs
	if (dumpPe) {
#pragma omp critical
	  {
	    for(typename TCliqueMembers::const_iterator itC=clique.begin(); itC!=clique.end(); ++itC) {
	      std::stringstream id;
	      id << _addID(svType) << std::setw(8) << std::setfill('0') << svRec.id;
	      dumpPeFile << id.str() << "\t" << references[g[*itC]->RefID].RefName << "\t" << g[*itC]->Position << "\t" <<  references[g[*itC]->MateRefID].RefName << "\t" << g[*itC]->MatePosition << "\t" << g[*itC]->MapQuality << std::endl;
	    }
	  }
	}
      }
    }
  }

  // Close dump PE file
  if (dumpPe) dumpPeFile.close();

  // Split-read search
  if (peMapping) {
    if (boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome)) 
      if (!svs.empty()) 
	findPutativeSplitReads(c, svs, svType);
  } else {
    // Read SV records from input vcffile
    vcfParse(c, references, overallMaxISize, svs, svType);
  }

  // Debug output
  //for (TVariants::const_iterator s = svs.begin();s!=svs.end();++s) std::cerr << s->svStart << ',' << s->svEnd << ',' <<  s->svStartBeg << ',' << s->svStartEnd << ',' << s->svEndBeg << ',' << s->svEndEnd << ',' << s->peSupport << ',' << s->srSupport << ',' << s->wiggle << ',' << s->srAlignQuality << ',' << s->precise << ',' << s->ct << ',' << s->peMapQuality << ',' << s->chr << ',' << s->chr2 << ',' << s->consensus << std::endl;

  // Any SVs for genotyping
  if (svs.empty()) {
    std::cout << "No structural variants found!" << std::endl;
    std::cout << "Done." << std::endl;
    return 0;
  }

  // Annotate junction reads
  typedef std::pair<std::string, int> TSampleSVPair;
  typedef std::pair<std::vector<uint16_t>, std::vector<uint16_t> > TReadQual;
  typedef boost::unordered_map<TSampleSVPair, TReadQual> TJunctionCountMap;
  TJunctionCountMap junctionCountMap;
  if (boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))
    _annotateJunctionReads(c, sampleLib, svs, junctionCountMap, svType);

  // Annotate spanning coverage
  typedef boost::unordered_map<TSampleSVPair, TReadQual> TCountMap;
  TCountMap spanCountMap;
  _annotateSpanningCoverage(c, sampleLib, svs, spanCountMap, svType);

  // Annotate coverage
  typedef std::pair<int, int> TBpRead;
  typedef boost::unordered_map<TSampleSVPair, TBpRead> TReadCountMap;
  TReadCountMap readCountMap;
  _annotateCoverage(c, sampleLib, svs, readCountMap, svType);

  // VCF output
  if (svs.size()) {
    vcfOutput(c, svs, junctionCountMap, readCountMap, spanCountMap, svType);
  }


  // Output library statistics
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Library statistics" << std::endl;
  TSampleLibrary::const_iterator sampleIt=sampleLib.begin();
  for(;sampleIt!=sampleLib.end();++sampleIt) {
    std::cout << "Sample: " << sampleIt->first << std::endl;
    TLibraryMap::const_iterator libIt=sampleIt->second.begin();
    for(;libIt!=sampleIt->second.end();++libIt) {
      std::cout << "RG: ID=" << libIt->first << ",Median=" << libIt->second.median << ",MAD=" << libIt->second.mad << ",Orientation=" << (int) libIt->second.defaultOrient << ",InsertSizeCutoff=" << libIt->second.maxNormalISize << ",DuplicateDiscordantPairs=" << libIt->second.non_unique_abnormal_pairs << ",UniqueDiscordantPairs=" << libIt->second.abnormal_pairs << std::endl;
    }
  }

#ifdef PROFILE
  ProfilerStop();
#endif

  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
  return 0;
}


int main(int argc, char **argv) {
  Config c;

  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("type,t", boost::program_options::value<std::string>(&c.svType)->default_value("DEL"), "SV analysis type (DEL, DUP, INV, TRA)")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("sv.vcf"), "SV output file")
    ("exclude,x", boost::program_options::value<boost::filesystem::path>(&c.exclude)->default_value(""), "file with chr to exclude")
    ;

  boost::program_options::options_description pem("PE options");
  pem.add_options()
    ("map-qual,q", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(1), "min. paired-end mapping quality")
    ("mad-cutoff,s", boost::program_options::value<unsigned short>(&c.madCutoff)->default_value(9), "insert size cutoff, median+s*MAD (deletions only)")
    ;

  boost::program_options::options_description breaks("SR options");
  breaks.add_options()
    ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
    ("min-flank,m", boost::program_options::value<unsigned int>(&c.minimumFlankSize)->default_value(13), "minimum flanking sequence size")
    ;

  boost::program_options::options_description geno("Genotyping options");
  geno.add_options()
    ("vcfgeno,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile)->default_value("site.vcf"), "input vcf file for genotyping only")
    ("geno-qual,u", boost::program_options::value<unsigned short>(&c.minGenoQual)->default_value(1), "min. mapping quality for genotyping")
    ;

  // Define hidden options
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
    ("pe-dump,p", boost::program_options::value<boost::filesystem::path>(&c.peDump)->default_value(""), "outfile to dump PE info")
    ("epsilon,e", boost::program_options::value<float>(&c.epsilon)->default_value(0.1), "allowed epsilon deviation of PE vs. SR deletion")
    ("pe-fraction,c", boost::program_options::value<float>(&c.percentAbnormal)->default_value(0.0), "fixed fraction c of discordant PEs, for c=0 MAD cutoff is used")
    ("num-split,n", boost::program_options::value<unsigned int>(&c.minimumSplitRead)->default_value(2), "minimum number of splitted reads")
    ("flanking,f", boost::program_options::value<unsigned int>(&c.flankQuality)->default_value(80), "quality of the aligned flanking region")
    ("pruning,j", boost::program_options::value<unsigned int>(&c.graphPruning)->default_value(100), "PE graph pruning cutoff")
    ("warranty,w", "show warranty")
    ("license,l", "show license")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  // Set the visibility
  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(pem).add(breaks).add(geno).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic).add(pem).add(breaks).add(geno);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);


  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file"))) { 
    printTitle("DELLY");
    if (vm.count("warranty")) {
      displayWarranty();
    } else if (vm.count("license")) {
      gplV3();
    } else {
      std::cout << "Usage: " << argv[0] << " [OPTIONS] <sample1.sort.bam> <sample2.sort.bam> ..." << std::endl;
      std::cout << visible_options << "\n"; 
    }
    return 1; 
  }

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Always ignore reads of mapping quality 0 for genotyping
  if (c.minGenoQual<1) c.minGenoQual=1;

  // Run main program
  if (c.svType == "DEL") return run(c, SVType<DeletionTag>());
  else if (c.svType == "DUP") return run(c, SVType<DuplicationTag>());
  else if (c.svType == "INV") return run(c, SVType<InversionTag>());
  else if (c.svType == "TRA") return run(c, SVType<TranslocationTag>());
  else if (c.svType == "INS") return run(c, SVType<InsertionTag>());
  else {
    std::cerr << "SV analysis type not supported by Delly: " << c.svType << std::endl;
    return -1;
  }
  
}
