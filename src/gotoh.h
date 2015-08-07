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

#ifndef GOTOH_H
#define GOTOH_H

#include <iostream>
#include <list>


namespace torali
{

  template<typename TAlign1, typename TAlign2, typename TProfile, typename TAIndex>
  inline int
  _score(TAlign1 const& a1, TAlign2 const& a2, TProfile const& p1, TProfile const& p2, TAIndex row, TAIndex col)
  {
    typedef int TScoreValue;
    TScoreValue match = 5;
    TScoreValue mismatch = -4;
    if ((a1.shape()[0] == 1) && (a2.shape()[0] == 1)) {
      if (a1[0][row] == a2[0][col]) return match;
      else return mismatch;
    } else {
      typedef typename TProfile::index TPIndex;
      double score = 0;
      for(TPIndex k1 = 0; k1<5; ++k1) 
	for(TPIndex k2 = 0; k2<5; ++k2) 
	  score += p1[k1][row] * p2[k2][col] * ((k1 == k2) ? match : mismatch);
      return ((int) score);
    }
  }

  template<typename TAlign1, typename TAlign2, typename TAlign>
  inline int
  gotoh(TAlign1 const& a1, TAlign2 const& a2, TAlign& align)
  {
    typedef typename TAlign::index TAIndex;
    typedef int TScoreValue;
    TScoreValue inf = 1000000;
    TScoreValue go = -10;
    TScoreValue ge = -1;

    // DP Matrix
    typedef boost::multi_array<TScoreValue, 2> TMatrix;
    std::size_t m = a1.shape()[1];
    std::size_t n = a2.shape()[1];
    TMatrix s(boost::extents[m+1][n+1]);
    TMatrix h(boost::extents[m+1][n+1]);
    TMatrix v(boost::extents[m+1][n+1]);

    // Create profile
    typedef boost::multi_array<double, 2> TProfile;
    typedef typename TProfile::index TPIndex;
    TProfile p1;
    TProfile p2;
    if ((a1.shape()[0] != 1) || (a2.shape()[0] != 1)) {
      p1.resize(boost::extents[6][a1.shape()[1]]);   // 'A', 'C', 'G', 'T', 'N', '-'
      for (TAIndex j = 0; j<a1.shape()[1]; ++j) {
	for(TPIndex k = 0; k<6; ++k) p1[k][j] = 0;
	int sum = 0;
	for(TAIndex i = 0; i<a1.shape()[0]; ++i) {
	  ++sum;
	  if ((a1[i][j] == 'A') || (a1[i][j] == 'a')) p1[0][j] += 1;
	  else if ((a1[i][j] == 'C') || (a1[i][j] == 'c')) p1[1][j] += 1;
	  else if ((a1[i][j] == 'G') || (a1[i][j] == 'g')) p1[2][j] += 1;
	  else if ((a1[i][j] == 'T') || (a1[i][j] == 't')) p1[3][j] += 1;
	  else if ((a1[i][j] == 'N') || (a1[i][j] == 'n')) p1[4][j] += 1;
	  else if (a1[i][j] == '-') p1[5][j] += 1;
	  else --sum;
	}
	for(TPIndex k = 0; k<6; ++k) p1[k][j] /= sum;
      }
      p2.resize(boost::extents[6][a2.shape()[1]]);  
      for (TAIndex j = 0; j<a2.shape()[1]; ++j) {
	for(TPIndex k = 0; k<6; ++k) p2[k][j] = 0;
	int sum = 0;
	for(TAIndex i = 0; i<a2.shape()[0]; ++i) {
	  ++sum;
	  if ((a2[i][j] == 'A') || (a2[i][j] == 'a')) p2[0][j] += 1;
	  else if ((a2[i][j] == 'C') || (a2[i][j] == 'c')) p2[1][j] += 1;
	  else if ((a2[i][j] == 'G') || (a2[i][j] == 'g')) p2[2][j] += 1;
	  else if ((a2[i][j] == 'T') || (a2[i][j] == 't')) p2[3][j] += 1;
	  else if ((a2[i][j] == 'N') || (a2[i][j] == 'n')) p2[4][j] += 1;
	  else if (a2[i][j] == '-') p2[5][j] += 1;
	  else --sum;
	}
	for(TPIndex k = 0; k<6; ++k) p2[k][j] /= sum;
      }
    }

    // Initialization
    for(std::size_t col = 1; col <= n; ++col) {
      v[0][col] = -inf;
      s[0][col] = go/2 + col * ge;
      h[0][col] = go/2 + col * ge;
    }
    for(std::size_t row = 1; row <= m; ++row) {
      h[row][0] = -inf;
      s[row][0] = go/2 + row * ge;
      v[row][0] = go/2 + row * ge;
    }
    s[0][0] = 0;
    v[0][0] = -inf;
    h[0][0] = -inf;

    // Recursion
    for(std::size_t col = 1; col <= n; ++col) {
      for(std::size_t row = 1; row <= m; ++row) {
	h[row][col] = std::max(s[row][col-1] + go + ge, h[row][col-1] + ge);
	v[row][col] = std::max(s[row-1][col] + go + ge, v[row-1][col] + ge);
	s[row][col] = std::max(std::max(s[row-1][col-1] + _score(a1, a2, p1, p2, row-1, col-1), h[row][col]), v[row][col]);
      }
    }

    // Trace-back
    std::size_t row = m;
    std::size_t col = n;
    char lastMatrix = 's';
    typedef std::vector<char> TTrace;
    TTrace trace;
    while ((row>0) || (col>0)) {
      if (lastMatrix == 's') {
	if (s[row][col] == h[row][col]) lastMatrix = 'h';
	else if (s[row][col] == v[row][col]) lastMatrix = 'v';
	else {
	  --row;
	  --col;
	  trace.push_back('s');
	  //std::cerr << a1[0][row] << a2[0][col] << std::endl;
	}
      } else if (lastMatrix == 'h') {
	if (h[row][col] != h[row][col-1] + ge) lastMatrix = 's';
	--col;
	trace.push_back('h');
	//std::cerr << '-' << a2[0][col] << std::endl;
      } else if (lastMatrix == 'v') {
	if (v[row][col] != v[row-1][col] + ge) lastMatrix = 's';
	--row;
	trace.push_back('v');
	//std::cerr << a1[0][row] << '-' << std::endl;
      }
    }

    // Create alignment
    align.resize(boost::extents[a1.shape()[0] + a2.shape()[0]][trace.size()]);
    row=0;
    col=0;
    TAIndex numN = a1.shape()[0]; 
    TAIndex numM = a2.shape()[0]; 
    TAIndex ai = 0;
    for(typename TTrace::const_reverse_iterator itT = trace.rbegin(); itT != trace.rend(); ++itT, ++ai) {
      if (*itT == 's') {
	for(TAIndex i = 0; i<numN; ++i) align[i][ai] = a1[i][row];
	for(TAIndex i = 0; i<numM; ++i) align[numN + i][ai] = a2[i][col];
	++row;
	++col;
      } else if (*itT =='h') {
	for(TAIndex i = 0; i<numN; ++i) align[i][ai] = '-';
	for(TAIndex i = 0; i<numM; ++i) align[numN + i][ai] = a2[i][col];
	++col;
      } else {
	for(TAIndex i = 0; i<numN; ++i) align[i][ai] = a1[i][row];
	for(TAIndex i = 0; i<numM; ++i) align[numN + i][ai] = '-';
	++row;
      }
    }
    return s[m][n];
  }

}

#endif
