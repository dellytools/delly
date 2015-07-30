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
  inline void
  upgma(TDistArray& d, TPhylogeny& p, TDIndex num) {
    for(TDIndex nn = num; nn<2*num+1; ++nn) {
      TDIndex dI = 0;
      TDIndex dJ = 0;
      if (closestPair(d, nn, dI, dJ) == -1) break;
      p[dI][0] = nn;
      p[dJ][0] = nn;
      p[nn][1] = dI;
      p[nn][2] = dJ;
      updateDistanceMatrix(d, p, nn, dI, dJ);
    }
  }

  template<typename TSplitReadSet>
  inline void
  msa(TSplitReadSet const& sps) {
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
    upgma(d, p, num);
    
    // Debug
    std::cerr << "Phylogeny" << std::endl;
    for(TDIndex i = 0; i<(2*num+1); ++i) {
      std::cerr << i << ':' << '\t';
      for (TDIndex j = 0; j<3; ++j) {
	std::cerr << p[i][j] << '\t';
      }
      std::cerr << std::endl;
    }

  }

}

#endif
