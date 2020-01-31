#ifndef GOTOH_H
#define GOTOH_H

#include <boost/dynamic_bitset.hpp>

#include <iostream>
#include "align.h"

namespace torali
{

  template<typename TAlign1, typename TAlign2, typename TAlignConfig, typename TScoreObject>
  inline int
  gotohScore(TAlign1 const& a1, TAlign2 const& a2, TAlignConfig const& ac, TScoreObject const& sc)
  {
    typedef typename TScoreObject::TValue TScoreValue;

    // DP variables
    std::size_t m = _size(a1, 1);
    std::size_t n = _size(a2, 1);
    std::vector<TScoreValue> s(n+1, 0);
    std::vector<TScoreValue> v(n+1, 0);
    TScoreValue newhoz = 0;
    TScoreValue prevsub = 0;
    
    // Create profile
    typedef boost::multi_array<float, 2> TProfile;
    TProfile p1;
    TProfile p2;
    if ((_size(a1, 0) != 1) || (_size(a2, 0) != 1)) {
      _createProfile(a1, p1);
      _createProfile(a2, p2);
    }

    // DP
    for(std::size_t row = 0; row <= m; ++row) {
      for(std::size_t col = 0; col <= n; ++col) {
	// Initialization
	if ((row == 0) && (col == 0)) {
	  s[0] = 0;
	  v[0] = -sc.inf;
	  newhoz = -sc.inf;
	} else if (row == 0) {
	  v[col] = -sc.inf;
	  s[col] = _horizontalGap(ac, 0, m, sc.go + col * sc.ge);
	  newhoz = _horizontalGap(ac, 0, m, sc.go + col * sc.ge);
	} else if (col == 0) {
	  newhoz = -sc.inf;
	  s[0] = _verticalGap(ac, 0, n, sc.go + row * sc.ge);
	  if (row - 1 == 0) prevsub = 0;
	  else prevsub = _verticalGap(ac, 0, n, sc.go + (row - 1) * sc.ge);
	  v[0] = _verticalGap(ac, 0, n, sc.go + row * sc.ge);
	} else {
	  // Recursion
	  TScoreValue prevhoz = newhoz;
	  TScoreValue prevver = v[col];
	  TScoreValue prevprevsub = prevsub;
	  prevsub = s[col];
	  newhoz = std::max(s[col-1] + _horizontalGap(ac, row, m, sc.go + sc.ge), prevhoz + _horizontalGap(ac, row, m, sc.ge));
	  v[col] = std::max(prevsub + _verticalGap(ac, col, n, sc.go + sc.ge), prevver + _verticalGap(ac, col, n, sc.ge));
	  s[col] = std::max(std::max(prevprevsub + _score(a1, a2, p1, p2, row-1, col-1, sc), newhoz), v[col]);
	}
      }
    }

    // Score
    return s[n];
  }

  
  template<typename TAlign1, typename TAlign2, typename TAlign, typename TAlignConfig, typename TScoreObject>
  inline int
    gotoh(TAlign1 const& a1, TAlign2 const& a2, TAlign& align, TAlignConfig const& ac, TScoreObject const& sc)
  {
    typedef typename TScoreObject::TValue TScoreValue;

    // DP variables
    std::size_t m = _size(a1, 1);
    std::size_t n = _size(a2, 1);
    std::vector<TScoreValue> s(n+1, 0);
    std::vector<TScoreValue> v(n+1, 0);
    TScoreValue newhoz = 0;
    TScoreValue prevsub = 0;
    
    // Trace Matrix
    std::size_t mf = n+1;
    typedef boost::dynamic_bitset<> TBitSet;
    TBitSet bit1( (m+1) * (n+1), false);
    TBitSet bit2( (m+1) * (n+1), false);
    TBitSet bit3( (m+1) * (n+1), false);
    TBitSet bit4( (m+1) * (n+1), false);

    // Create profile
    typedef boost::multi_array<float, 2> TProfile;
    TProfile p1;
    TProfile p2;
    if ((_size(a1, 0) != 1) || (_size(a2, 0) != 1)) {
      _createProfile(a1, p1);
      _createProfile(a2, p2);
    }

    // DP
    for(std::size_t row = 0; row <= m; ++row) {
      for(std::size_t col = 0; col <= n; ++col) {
	// Initialization
	if ((row == 0) && (col == 0)) {
	  s[0] = 0;
	  v[0] = -sc.inf;
	  newhoz = -sc.inf;
	  bit1[0] = true;
	  bit2[0] = true;
	} else if (row == 0) {
	  v[col] = -sc.inf;
	  s[col] = _horizontalGap(ac, 0, m, sc.go + col * sc.ge);
	  newhoz = _horizontalGap(ac, 0, m, sc.go + col * sc.ge);
	  bit3[col] = true;
	} else if (col == 0) {
	  newhoz = -sc.inf;
	  s[0] = _verticalGap(ac, 0, n, sc.go + row * sc.ge);
	  if (row - 1 == 0) prevsub = 0;
	  else prevsub = _verticalGap(ac, 0, n, sc.go + (row - 1) * sc.ge);
	  v[0] = _verticalGap(ac, 0, n, sc.go + row * sc.ge);
	  bit4[row * mf] = true;
	} else {
	  // Recursion
	  TScoreValue prevhoz = newhoz;
	  TScoreValue prevver = v[col];
	  TScoreValue prevprevsub = prevsub;
	  prevsub = s[col];
	  newhoz = std::max(s[col-1] + _horizontalGap(ac, row, m, sc.go + sc.ge), prevhoz + _horizontalGap(ac, row, m, sc.ge));
	  v[col] = std::max(prevsub + _verticalGap(ac, col, n, sc.go + sc.ge), prevver + _verticalGap(ac, col, n, sc.ge));
	  s[col] = std::max(std::max(prevprevsub + _score(a1, a2, p1, p2, row-1, col-1, sc), newhoz), v[col]);

	  // Trace
	  if (s[col] == newhoz) bit3[row * mf + col] = true;
	  else if (s[col] == v[col]) bit4[row * mf + col] = true;
	  if (newhoz != prevhoz + _horizontalGap(ac, row, m, sc.ge)) bit1[row * mf + col] = true;
	  if (v[col] != prevver + _verticalGap(ac, col, n, sc.ge)) bit2[row * mf + col] = true;
	}
      }
    }

    // Trace-back using pointers
    std::size_t row = m;
    std::size_t col = n;
    char lastMatrix = 's';
    typedef std::vector<char> TTrace;
    TTrace btr;
    while ((row>0) || (col>0)) {
      if (lastMatrix == 's') {
	if (bit3[row * mf + col]) lastMatrix = 'h';
	else if (bit4[row * mf + col]) lastMatrix = 'v';
	else {
	  --row;
	  --col;
	  btr.push_back('s');
	}
      } else if (lastMatrix == 'h') {
	if (bit1[row * mf + col]) lastMatrix = 's';
	--col;
	btr.push_back('h');
      } else if (lastMatrix == 'v') {
	if (bit2[row * mf + col]) lastMatrix = 's';
	--row;
	btr.push_back('v');
      }
    }

    // Create alignment
    _createAlignment(btr, a1, a2, align);

    // Score
    return s[n];
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
