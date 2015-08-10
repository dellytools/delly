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

#ifndef MSA_H
#define MSA_H

#include <boost/multi_array.hpp>
#include "needle.h"
#include "gotoh.h"

namespace torali {

  inline int
  lcs(std::string s1, std::string s2) {
    int m = s1.size();
    int n = s2.size();
    typedef boost::multi_array<int, 2> T2DArray;
    typedef T2DArray::index T2DIndex;
    T2DArray L(boost::extents[m+1][n+1]);
    for(T2DIndex i=0; i<=m; ++i) {
      for(T2DIndex j=0; j<=n; ++j) {
	if ((i==0) || (j==0)) L[i][j]=0;
	else if (s1[i-1] == s2[j-1]) L[i][j] = L[i-1][j-1] + 1;
	else L[i][j] = (L[i-1][j] > L[i][j-1]) ? L[i-1][j] : L[i][j-1];
      }
    }
    return L[m][n];
  }

  template<typename TSplitReadSet, typename TDistArray>
  inline void
  distanceMatrix(TSplitReadSet const& sps, TDistArray& d) {
    typedef typename TDistArray::index TDIndex;
    typename TSplitReadSet::const_iterator sIt1 = sps.begin();
    for (TDIndex i = 0; sIt1 != sps.end(); ++sIt1, ++i) {
      typename TSplitReadSet::const_iterator sIt2 = sIt1;
      ++sIt2;
      for (TDIndex j = i+1; sIt2 != sps.end(); ++sIt2, ++j) {
	d[i][j] = (lcs(*sIt1, *sIt2) * 100) / std::min(sIt1->size(), sIt2->size());
      }
    }
  }

  template<typename TDistArray, typename TDIndex>
  inline int
  closestPair(TDistArray const& d, TDIndex num, TDIndex& dI, TDIndex& dJ) {
    int dMax = -1;
    for (TDIndex i = 0; i<num; ++i) {
      for (TDIndex j = i+1; j<num; ++j) {
	if (d[i][j]>dMax) {
	  dMax = d[i][j];
	  dI = i;
	  dJ = j;
	}
      }
    }
    return dMax;
  }

  template<typename TDistArray, typename TPhylogeny, typename TDIndex>
  inline void
  updateDistanceMatrix(TDistArray& d, TPhylogeny const& p, TDIndex num, TDIndex& dI, TDIndex& dJ) {
    for (TDIndex i = 0; i < num; ++i) 
      if (p[i][0] == -1) 
	d[i][num] = (((dI < i) ? d[dI][i] : d[i][dI]) + ((dJ < i) ? d[dJ][i] : d[i][dJ])) / 2;
    for (TDIndex i = 0; i<dI; ++i) d[i][dI] = -1;
    for (TDIndex i = dI+1; i<num+1; ++i) d[dI][i] = -1;
    for (TDIndex i = 0; i<dJ; ++i) d[i][dJ] = -1;
    for (TDIndex i = dJ+1; i<num+1; ++i) d[dJ][i] = -1;
  }

  template<typename TDistArray, typename TPhylogeny, typename TDIndex>
  inline TDIndex
  upgma(TDistArray& d, TPhylogeny& p, TDIndex num) {
    TDIndex nn = num;
    for(;nn<2*num+1; ++nn) {
      TDIndex dI = 0;
      TDIndex dJ = 0;
      if (closestPair(d, nn, dI, dJ) == -1) break;
      p[dI][0] = nn;
      p[dJ][0] = nn;
      p[nn][1] = dI;
      p[nn][2] = dJ;
      updateDistanceMatrix(d, p, nn, dI, dJ);
    }
    return (nn > 0) ? (nn - 1) : 0;
  }

  template<typename TSplitReadSet, typename TPhylogeny, typename TDIndex, typename TAlign>
  inline void
  palign(TSplitReadSet const& sps, TPhylogeny const& p, TDIndex root, TAlign& align) {
    typedef typename TAlign::index TAIndex;
    if ((p[root][1] == -1) && (p[root][2] == -1)) {
      typename TSplitReadSet::const_iterator sIt = sps.begin();
      if (root) std::advance(sIt, root);
      align.resize(boost::extents[1][sIt->size()]);
      TAIndex ind = 0;
      for(typename std::string::const_iterator str = sIt->begin(); str != sIt->end(); ++str) align[0][ind++] = *str;
    } else {
      TAlign align1;
      palign(sps, p, p[root][1], align1);
      TAlign align2;
      palign(sps, p, p[root][2], align2);
      gotoh(align1, align2, align);
    }
  }

  template<typename TAlign>
  inline void
  sprealign(TAlign& align) {
    typedef typename TAlign::index TAIndex;
    for(TAIndex i = 0; i<align.shape()[0]; ++i) {
      // Precompute sub-alignment shapes
      int seqLength = 0;
      std::set<TAIndex> gaps;
      for(TAIndex j = 0; j<align.shape()[1]; ++j) {	
	if (align[i][j] != '-') {
	  ++seqLength;
	  bool gapOnly = true;
	  for(TAIndex k = 0; k<align.shape()[0]; ++k) {
	    if ((k!=i) && (align[k][j] != '-')) {
	      gapOnly = false;
	      break;
	    }
	  }
	  if (gapOnly) gaps.insert(j);
	}
      }

      // Create sub-alignments
      TAlign align1;
      TAIndex aind1 = 0;
      align1.resize(boost::extents[1][seqLength]);
      TAlign align2;
      TAIndex aind2 = 0;
      align2.resize(boost::extents[align.shape()[0] - 1][align.shape()[1] - gaps.size()]);
      for(TAIndex j = 0; j<align.shape()[1]; ++j) {
	if (align[i][j] != '-') align1[0][aind1++] = align[i][j];
	if (gaps.find(j) == gaps.end()) {
	  TAIndex kr = 0;
	  for(TAIndex k = 0; k<align.shape()[0]; ++k)
	    if (k!=i) align2[kr++][aind2] = align[k][j];
	  ++aind2;
	}
      }
      
      // Re-align sequence to profile
      gotoh(align1, align2, align);
    }
  }


  template<typename TAlign>
  inline void
  consensus(TAlign const& align, std::string& cs) {
    typedef typename TAlign::index TAIndex;

    // Calculate coverage
    typedef boost::multi_array<bool, 2> TFlag;
    TFlag fl;
    fl.resize(boost::extents[align.shape()[0]][align.shape()[1]]);
    typedef std::vector<int> TCoverage;
    TCoverage cov;
    cov.resize(align.shape()[1], 0);
    for(TAIndex i = 0; i<align.shape()[0]; ++i) {
      int start = 0;
      int end = -1;
      for(TAIndex j = 0; j<align.shape()[1]; ++j) {
	fl[i][j] = false;
	if (align[i][j] != '-') end = j;
	else if (end == -1) start = j + 1;
      }
      for(TAIndex j = start; j<=end; ++j) {
	++cov[j];
	fl[i][j] = true;
      }
    }
    
    int covThreshold = 3;
    TAIndex j = 0;
    std::vector<char> cons;
    for(typename TCoverage::const_iterator itCov = cov.begin(); itCov != cov.end(); ++itCov, ++j) {
      if (*itCov >= covThreshold) {
	// Get consensus letter
	int countA = 0;
	int countC = 0;
	int countG = 0;
	int countT = 0;
	for(TAIndex i = 0; i<align.shape()[0]; ++i) {
	  if (fl[i][j]) {
	    if ((align[i][j] == 'A') || (align[i][j] == 'a')) ++countA;
	    else if ((align[i][j] == 'C') || (align[i][j] == 'c')) ++countC;
	    else if ((align[i][j] == 'G') || (align[i][j] == 'g')) ++countG;
	    else if ((align[i][j] == 'T') || (align[i][j] == 't')) ++countT;
	  }
	}
	if (countA > (*itCov / 2)) cons.push_back('A');
	else if (countC > (*itCov / 2)) cons.push_back('C');
	else if (countG > (*itCov / 2)) cons.push_back('G');
	else if (countT > (*itCov / 2)) cons.push_back('T');
      }
    }
    cs = std::string(cons.begin(), cons.end());
  }


  template<typename TSplitReadSet>
  inline void
  msa(TSplitReadSet const& sps, std::string& cs) {
    // Compute distance matrix
    typedef boost::multi_array<int, 2> TDistArray;
    typedef typename TDistArray::index TDIndex;
    TDIndex num = sps.size();
    TDistArray d(boost::extents[2*num+1][2*num+1]);
    for (TDIndex i = 0; i<(2*num+1); ++i) 
      for (TDIndex j = i+1; j<(2*num+1); ++j) 
	d[i][j]=-1;
    distanceMatrix(sps, d);

    // UPGMA
    typedef boost::multi_array<int, 2> TPhylogeny;
    TPhylogeny p(boost::extents[2*num+1][3]);
    for(TDIndex i = 0; i<(2*num+1); ++i) 
      for (TDIndex j = 0; j<3; ++j) p[i][j] = -1;
    TDIndex root = upgma(d, p, num);

    // Debug guide tree
    //std::cerr << "Phylogeny" << std::endl;
    //std::cerr << "#Sequences: " << sps.size() << std::endl;
    //std::cerr << "Root: " << root << std::endl;
    //std::cerr << "Node:Parent\tLeftChild\tRightChild" << std::endl;
    //for(TDIndex i = 0; i<(2*num+1); ++i) {
    //std::cerr << i << ':' << '\t';
    //for (TDIndex j = 0; j<3; ++j) {
    //std::cerr << p[i][j] << '\t';
    //}
    //std::cerr << std::endl;
    //}
    
    // Progressive Alignment
    typedef boost::multi_array<char, 2> TAlign;
    TAlign align;
    palign(sps, p, root, align);

    // Debug MSA
    typedef typename TAlign::index TAIndex;
    for(TAIndex i = 0; i<align.shape()[0]; ++i) {
      for(TAIndex j = 0; j<align.shape()[1]; ++j) {
	std::cerr << align[i][j];
      }
      std::cerr << std::endl;
    }

    // Sequence to profile re-alignment
    //sprealign(align);

    // Consensus calling
    consensus(align, cs);
    std::cerr << cs << std::endl;
    std::cerr << std::endl;
  }

}

#endif
