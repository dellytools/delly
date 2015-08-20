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
#include "align.h"

namespace torali
{
  
  template<typename TAlign1, typename TAlign2, typename TAlign, typename TAlignConfig, typename TScoreObject>
  inline int
    gotoh(TAlign1 const& a1, TAlign2 const& a2, TAlign& align, TAlignConfig const& ac, TScoreObject const& sc)
  {
    typedef typename TScoreObject::TValue TScoreValue;

    // DP Matrix
    typedef boost::multi_array<TScoreValue, 2> TMatrix;
    std::size_t m = _size(a1, 1);
    std::size_t n = _size(a2, 1);
    TMatrix s(boost::extents[m+1][n+1]);
    TMatrix h(boost::extents[m+1][n+1]);
    TMatrix v(boost::extents[m+1][n+1]);

    // Create profile
    typedef boost::multi_array<double, 2> TProfile;
    TProfile p1;
    TProfile p2;
    if ((_size(a1, 0) != 1) || (_size(a2, 0) != 1)) {
      _createProfile(a1, p1);
      _createProfile(a2, p2);
    }

    // Initialization
    for(std::size_t col = 1; col <= n; ++col) {
      v[0][col] = -sc.inf;
      s[0][col] = _horizontalGap(ac, 0, m, sc.go + col * sc.ge);
      h[0][col] = _horizontalGap(ac, 0, m, sc.go + col * sc.ge);
    }
    for(std::size_t row = 1; row <= m; ++row) {
      h[row][0] = -sc.inf;
      s[row][0] = _verticalGap(ac, 0, n, sc.go + row * sc.ge);
      v[row][0] = _verticalGap(ac, 0, n, sc.go + row * sc.ge);
    }
    s[0][0] = 0;
    v[0][0] = -sc.inf;
    h[0][0] = -sc.inf;

    // Recursion
    for(std::size_t col = 1; col <= n; ++col) {
      for(std::size_t row = 1; row <= m; ++row) {
	h[row][col] = std::max(s[row][col-1] + _horizontalGap(ac, row, m, sc.go + sc.ge), h[row][col-1] + _horizontalGap(ac, row, m, sc.ge));
	v[row][col] = std::max( s[row-1][col] + _verticalGap(ac, col, n, sc.go + sc.ge), v[row-1][col] + _verticalGap(ac, col, n, sc.ge));
	s[row][col] = std::max(std::max(s[row-1][col-1] + _score(a1, a2, p1, p2, row-1, col-1, sc), h[row][col]), v[row][col]);
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
	if (h[row][col] != h[row][col-1] + _horizontalGap(ac, row, m, sc.ge)) lastMatrix = 's';
	--col;
	trace.push_back('h');
	//std::cerr << '-' << a2[0][col] << std::endl;
      } else if (lastMatrix == 'v') {
	if (v[row][col] != v[row-1][col] + _verticalGap(ac, col, n, sc.ge)) lastMatrix = 's';
	--row;
	trace.push_back('v');
	//std::cerr << a1[0][row] << '-' << std::endl;
      }
    }

    // Create alignment
    _createAlignment(trace, a1, a2, align);

    // Score
    return s[m][n];
  }

  template<typename TAlign1, typename TAlign2, typename TAlign, typename TAlignConfig>
  inline int
  gotoh(TAlign1 const& a1, TAlign2 const& a2, TAlign& align, TAlignConfig const& ac) 
  {
    DnaScore<int> dnasc;    
    return gotoh(a1, a2, align, ac, dnasc);
  }

  template<typename TAlign1, typename TAlign2, typename TAlign>
  inline int
  gotoh(TAlign1 const& a1, TAlign2 const& a2, TAlign& align)
  {
    AlignConfig<false, false> ac;
    return gotoh(a1, a2, align, ac);
  }

}

#endif
