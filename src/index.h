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

#ifndef INDEX_H
#define INDEX_H

#include <iostream>

namespace torali {

  // Metaprogramming for power function
  template<long B, unsigned long E>
    struct pow_helper {
      static long const value = B * pow_helper<B, E - 1>::value;
    };
  
  template<long B>
    struct pow_helper<B, 0> {
      static long const value = 1;
    };
  
  template<long B, long E>
    struct power {
      static double const value;
    };

  template<long B, long E>
    double const power<B, E>::value= E < 0 ? 1.0 / pow_helper<B, E < 0 ? -E : 0>::value
    : pow_helper<B, E < 0 ? 0 : E>::value;



  template<typename TGenomicPos, typename TUSize64, typename TAlphabet, uint64_t KMER, uint64_t ALPHSIZE>
    struct Index {
      enum { kmer_size = KMER, alph_size = ALPHSIZE };


      typedef std::vector<TAlphabet> TSequence;
      TGenomicPos* kt;
      TGenomicPos* hit;
      TSequence refSequence;
      
      Index() {
	kt = new TGenomicPos[(TUSize64) power<ALPHSIZE, KMER>::value + 1];
	hit = NULL;
      }

      Index(Index const& ind) {
	kt = new TGenomicPos[(TUSize64) power<ALPHSIZE, KMER>::value + 1];
	for(TUSize64 i = 0; i < (TUSize64) power<ALPHSIZE, KMER>::value + 1; ++i) {
	  kt[i] = ind.kt[i];
	}
	hit = NULL;
	refSequence = ind.refSequence;
      }

      ~Index() {
	delete [] kt;
	delete [] hit;
      }

      inline void indexSequence(std::string const& seq) {
	refSequence = TSequence(seq.begin(), seq.end());

	// Index the reference
	// Kmer table
	for(TUSize64 i = 0; i <= (TUSize64) power<ALPHSIZE, KMER>::value; ++i) kt[i] = 0;

	// Current kmer
	char currentk[KMER];
	for(unsigned int i = 0; i<KMER; ++i) currentk[i] = 0;

	// Count kmer's
	typename TSequence::const_iterator refSeqIt = refSequence.begin();
	typename TSequence::const_iterator refSeqItEnd = refSequence.end();
	TUSize64 bucket = 0;
	TUSize64 kmerlen = 0;
	for(TUSize64 seqIndex=0;refSeqIt != refSeqItEnd; ++refSeqIt, ++seqIndex) {
	  if (dna5_encode[(int) *refSeqIt] != 4) {
	    if (kmerlen == KMER) {
	      ++kt[bucket];
	      bucket -= ((TUSize64) currentk[seqIndex % KMER] * (TUSize64) power<ALPHSIZE, KMER - 1>::value);
	      bucket *= (TUSize64) ALPHSIZE;
	      bucket += (TUSize64) dna5_encode[(int) *refSeqIt];
	    } else bucket += ((TUSize64) dna5_encode[(int) *refSeqIt] * (TUSize64) std::pow((TUSize64) ALPHSIZE, (TUSize64) KMER - (++kmerlen)));
	    currentk[seqIndex % KMER] = dna5_encode[(int) *refSeqIt];
	  } else {
	    // Handle N's properly
	    if (kmerlen == KMER) ++kt[bucket];
	    kmerlen = bucket = 0;
	  }
	}
	if (kmerlen == KMER) ++kt[bucket];

	// Print count table
	//for (TUSize64 i = 0; i < power<ALPHSIZE, KMER>::value; ++i) std::cout << i << ": " << kt[i] << "; ";
	//std::cout << std::endl;

	// Allocate hit list and set pointers
	TUSize64 lenHitList = 0;
	for (TUSize64 i = 0; i < power<ALPHSIZE, KMER>::value; ++i) lenHitList += kt[i];
	hit = new TGenomicPos[lenHitList];
	TGenomicPos endLocal = kt[0];
	kt[0] = 0;
	for(TUSize64 i = 1; i < (TUSize64) power<ALPHSIZE, KMER>::value; ++i) {
	  TGenomicPos tmp = kt[i];
	  kt[i] = endLocal;
	  endLocal += tmp;
	}
	
	// Populate hit list
	bucket = 0;
	kmerlen = 0;
	refSeqIt = refSequence.begin();
	refSeqItEnd = refSequence.end();
	TUSize64 seqIndex=0;
	for(;refSeqIt != refSeqItEnd; ++refSeqIt, ++seqIndex) {
	  if (dna5_encode[(int) *refSeqIt] != 4) {
	    if (kmerlen == KMER) {
	      hit[kt[bucket]++] = seqIndex - KMER;
	      bucket -= ((TUSize64) currentk[seqIndex % KMER] * (TUSize64) power<ALPHSIZE, KMER - 1>::value);
	      bucket *= (TUSize64) ALPHSIZE;
	      bucket += (TUSize64) dna5_encode[(int) *refSeqIt];
	    } else bucket += ((TUSize64) dna5_encode[(int) *refSeqIt] * (TUSize64) std::pow((TUSize64) ALPHSIZE, (TUSize64) KMER - (++kmerlen)));
	    currentk[seqIndex % KMER] = dna5_encode[(int) *refSeqIt];
	  } else {
	    // Handle N's properly
	    if (kmerlen == KMER) hit[kt[bucket]++] = seqIndex - KMER;
	    kmerlen = bucket = 0;
	  }
	}
	if (kmerlen == KMER) hit[kt[bucket]++] = seqIndex - KMER;

	// Reset all the pointers
	TGenomicPos newPos = kt[0];
	kt[0] = 0;
	for(TUSize64 i = 1; i < (TUSize64) power<ALPHSIZE, KMER>::value; ++i) {
	  TGenomicPos oldPos = kt[i];
	  kt[i] = newPos;
	  newPos = oldPos;
	}
	kt[(TUSize64) power<ALPHSIZE, KMER>::value] = lenHitList;
      }

      template<typename TSeq, typename TDiag>
      inline void
      diagonalForwardKmerCount(TSeq& seq, TDiag& forward) {
	// Current kmer
	char currentk[KMER];
	for(unsigned int i = 0; i<KMER; ++i) currentk[i] = 0;

	// Count kmer's for forward read
	TUSize64 bucket = 0;
	TUSize64 kmerlen = 0;
	typename TSeq::const_iterator p = seq.begin();
	typename TSeq::const_iterator pEnd = seq.end();
	TUSize64 seqIndex=0;
	for(; p!=pEnd; ++p, ++seqIndex) {
	  if ((TUSize64) *p != 4) {
	    if (kmerlen == KMER) {
	      TGenomicPos itBeg = kt[bucket];
	      TGenomicPos itEnd = kt[bucket+1];
	      for(;itBeg != itEnd; ++itBeg) forward.push_back(std::make_pair(hit[itBeg] - (seqIndex - KMER), (seqIndex - KMER)));
	      bucket -= ((TUSize64) currentk[seqIndex % KMER] * (TUSize64) power<ALPHSIZE, KMER - 1>::value);
	      bucket *= (TUSize64) ALPHSIZE;
	      bucket += (TUSize64) *p;
	    } else bucket += (((TUSize64) *p) * (TUSize64) std::pow((TUSize64) ALPHSIZE, (TUSize64) KMER - (++kmerlen)));
	    currentk[seqIndex % KMER] = *p;
	  } else {
	    // Handle N's properly
	    if (kmerlen == KMER) {
	      TGenomicPos itBeg = kt[bucket];
	      TGenomicPos itEnd = kt[bucket+1];
	      for(;itBeg != itEnd; ++itBeg) forward.push_back(std::make_pair(hit[itBeg] - (seqIndex - KMER), (seqIndex - KMER)));
	    }
	    kmerlen = bucket = 0;
	  }
	}
	if (kmerlen == KMER) {
	  TGenomicPos itBeg = kt[bucket];
	  TGenomicPos itEnd = kt[bucket+1];
	  for(;itBeg != itEnd; ++itBeg) forward.push_back(std::make_pair(hit[itBeg] - (seqIndex - KMER), (seqIndex - KMER)));
	}
      }

      template<typename TSeq, typename TDiag>
      inline void
      diagonalReverseKmerCount(TSeq& seq, TDiag& reverse) {
	// Current kmer
	char currentk[KMER];
	for(unsigned int i = 0; i<KMER; ++i) currentk[i] = 0;

	// Count kmer's for forward read
	TUSize64 bucket = 0;
	TUSize64 seqIndex = 0;
	TUSize64 kmerlen = 0;
	typename TSeq::const_iterator pEnd = seq.begin();
	typename TSeq::const_iterator p = seq.end();
	do {
	  --p;
	  char ch = 0;
	  switch((TUSize64) *p) {
	  case 0: ch = 3; break;
	  case 1: ch = 2; break;
	  case 2: ch = 1; break;
	  case 3: ch = 0; break;
	  }
	  if ((TUSize64) ch != 4) {
	    if (kmerlen == KMER) {
	      TGenomicPos itBeg = kt[bucket];
	      TGenomicPos itEnd = kt[bucket+1];
	      for(;itBeg != itEnd; ++itBeg) reverse.push_back(std::make_pair(hit[itBeg] - (seqIndex - KMER), (seqIndex - KMER)));
	      bucket -= ((TUSize64) currentk[seqIndex % KMER] * (TUSize64) power<ALPHSIZE, KMER - 1>::value);
	      bucket *= (TUSize64) ALPHSIZE;
	      bucket += (TUSize64) ch;
	    } else bucket += ((TUSize64) ch * (TUSize64) std::pow((TUSize64) ALPHSIZE, (TUSize64) KMER - (++kmerlen)));
	    currentk[seqIndex % KMER] = ch;
	  } else {
	    // Handle N's properly
	    if (kmerlen == KMER) {
	      TGenomicPos itBeg = kt[bucket];
	      TGenomicPos itEnd = kt[bucket+1];
	      for(;itBeg != itEnd; ++itBeg) reverse.push_back(std::make_pair(hit[itBeg] - (seqIndex - KMER), (seqIndex - KMER)));
	    }
	    kmerlen = bucket = 0;
	  }
	  ++seqIndex;
	} while (p!=pEnd);
	if (kmerlen == KMER) {
	  TGenomicPos itBeg = kt[bucket];
	  TGenomicPos itEnd = kt[bucket+1];
	  for(;itBeg != itEnd; ++itBeg) reverse.push_back(std::make_pair(hit[itBeg] - (seqIndex - KMER), (seqIndex - KMER)));
	}
      }


    };




  // Global helper function


  template<typename TConfig, typename TDiag, typename TSize>
    inline bool
    _findBestSupportDiagonal(TConfig const&, TDiag& final, TDiag& best, TSize kmerHitCutoff, TSize totalKmer) {
    // Pick diagonals one by one, input must be sorted!!!
    std::vector<unsigned int> score;
    std::set<unsigned int> usedKmer;
    TDiag leftOver;
    while (final.size()) {
      typename TDiag::const_iterator itD = final.begin();
      typename TDiag::const_iterator itBef = final.begin();
      typename TDiag::const_iterator itDEnd = final.end();
      
      // Count the diagonals 
      unsigned int currentCount = 1;
      unsigned int bestCount = 0;
      int bestDiag =0;
      for(++itD;itD!=itDEnd; ++itD, ++itBef) {
	if (itBef->first == itD->first) ++currentCount;
	else {
	  if (bestCount < currentCount) { bestDiag = itBef->first; bestCount = currentCount; }
	  currentCount = 1;
	}
      }
      if (bestCount < currentCount) { bestDiag = itBef->first; bestCount = currentCount; }
    
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
	  usedKmer.insert(itD->second);  // Exclude all kmers of the read matched to the first diagonal
	}
      }
      best.push_back(std::make_pair(bestDiag, lastSeenKmer));
      score.push_back(bestCount);
      //std::cout << bestDiag << ',' << lastSeenKmer << '(' << bestCount << ')' << ';';
      
      // Process the left-over diagonals
      leftOver.clear();
      itD = final.begin();
      itDEnd = final.end();
      for(;itD!=itDEnd; ++itD) {
	if (usedKmer.find(itD->second) == usedKmer.end()) leftOver.push_back(*itD);
      }
      final.swap(leftOver);
    }
    //std::cout << std::endl;
    return ((score.size() > 1) && (score[0] + score[1] > totalKmer / 2)) ? true : false;
  }


}

#endif

