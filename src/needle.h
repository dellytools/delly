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

#ifndef NEEDLE_H
#define NEEDLE_H

#include <iostream>
#include "align.h"

namespace torali
{

  template<typename TAlign1, typename TAlign2, typename TAlign, typename TMatrix>
  inline int
  _needle(TAlign1 const& a1, TAlign2 const& a2, TAlign& align, TMatrix& mat)
  {
    typedef int TScoreValue;
    TScoreValue gap = -1;

    // DP Matrix
    std::size_t m = a1.shape()[1];
    std::size_t n = a2.shape()[1];
    mat.resize(boost::extents[m+1][n+1]);

    // Create profile
    typedef boost::multi_array<double, 2> TProfile;
    TProfile p1;
    TProfile p2;
    if ((a1.shape()[0] != 1) || (a2.shape()[0] != 1)) {
      _createProfile(a1, p1);
      _createProfile(a2, p2);
    }

    // Initialization
    mat[0][0] = 0;
    for(std::size_t col = 1; col <= n; ++col) mat[0][col] = mat[0][col-1] + gap;
    for(std::size_t row = 1; row <= m; ++row) mat[row][0] = mat[row-1][0] + gap;

    // Recursion
    for(std::size_t col = 1; col <= n; ++col) 
      for(std::size_t row = 1; row <= m; ++row) 
	mat[row][col] = std::max(std::max(mat[row-1][col-1] + _score(a1, a2, p1, p2, row-1, col-1), mat[row-1][col] + gap), mat[row][col-1] + gap);

    // Trace-back
    std::size_t row = m;
    std::size_t col = n;
    typedef std::vector<char> TTrace;
    TTrace trace;
    while ((row>0) || (col>0)) {
      if ((row>0) && (mat[row][col] == mat[row-1][col] + gap)) {
	--row;
	trace.push_back('v');
      } else if ((col>0) && (mat[row][col] == mat[row][col-1] + gap)) {
	--col;
	trace.push_back('h');
      } else {
	--row;
	--col;
	trace.push_back('s');
      }
    }

    // Create alignment
    _createAlignment(trace, a1, a2, align);

    // Score
    return mat[m][n];
  }


  template<typename TAlign1, typename TAlign2, typename TAlign>
  inline int
  needle(TAlign1 const& a1, TAlign2 const& a2, TAlign& align)
  {
    typedef boost::multi_array<int, 2> TMatrix;
    TMatrix mat;
    return _needle(a1, a2, align, mat);
  }

}

#endif
