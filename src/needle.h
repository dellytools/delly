#ifndef NEEDLE_H
#define NEEDLE_H

#define BOOST_DISABLE_ASSERTS
#include <boost/dynamic_bitset.hpp>
#include <boost/multi_array.hpp>
#include <iostream>
#include "align.h"

namespace torali
{

  inline int32_t
  longestHomology(std::string const& s1, std::string const& s2, int32_t scoreThreshold)  {
    // DP Matrix
    typedef boost::multi_array<int32_t, 2> TMatrix;
    int32_t m = s1.size();
    int32_t n = s2.size();
    TMatrix mat(boost::extents[m+1][n+1]);

    // Initialization
    int32_t k = std::abs(scoreThreshold);
    mat[0][0] = 0;
    for(int32_t col = 1; col <= k; ++col) mat[0][col] = mat[0][col-1] - 1;
    for(int32_t row = 1; row <= k; ++row) mat[row][0] = mat[row-1][0] - 1;

    // Edit distance
    for(int32_t row = 1; row <= m; ++row) {
      int32_t bestCol = scoreThreshold - 1;
      for(int32_t h = -k; h <= k; ++h) {
	int32_t col = row + h;
	if ((col >= 1) && (col <= n)) {
	  mat[row][col] = mat[row-1][col-1] + (s1[row-1] == s2[col-1] ? 0 : -1);
	  if ((row - 1 - col >= -k) && (row - 1 - col <= k)) mat[row][col] = std::max(mat[row][col], mat[row-1][col] - 1);
	  if ((row - col + 1 >= -k) && (row - col + 1 <= k)) mat[row][col] = std::max(mat[row][col], mat[row][col-1] - 1);
	  if (mat[row][col] > bestCol) bestCol = mat[row][col];
	}
      }
      if (bestCol < scoreThreshold) return row - 1;
    }
    return 0;
  }


  template<typename TAlign, typename TAlignConfig, typename TScoreObject>
  inline bool
  longNeedle(std::string const& s1, std::string const& s2, TAlign& align, TAlignConfig const& ac, TScoreObject const& sc)
  {
    typedef typename TScoreObject::TValue TScoreValue;
    typedef typename TAlign::index TAIndex;

    // DP Matrix
    typedef boost::multi_array<int16_t, 2> TMatrix;
    std::size_t m = s1.size();
    std::size_t n = s2.size();
    TMatrix mat(boost::extents[m+1][n+1]);

    // Initialization
    mat[0][0] = 0;
    for(std::size_t col = 1; col <= n; ++col) mat[0][col] = mat[0][col-1] + _horizontalGap(ac, 0, m, sc.ge);
    for(std::size_t row = 1; row <= m; ++row) mat[row][0] = mat[row-1][0] + _verticalGap(ac, 0, n, sc.ge);

    // Forward alignment
    for(std::size_t row = 1; row <= m; ++row)
      for(std::size_t col = 1; col <= n; ++col)
	mat[row][col] = std::max(std::max(mat[row-1][col-1] + (s1[row-1] == s2[col-1] ? sc.match : sc.mismatch), mat[row-1][col] + _verticalGap(ac, col, n, sc.ge)), mat[row][col-1] + _horizontalGap(ac, row, m, sc.ge));

    // Reverse input sequences
    std::string sRev1 = s1;
    reverseComplement(sRev1);
    std::string sRev2 = s2;
    reverseComplement(sRev2);

    // Reverse alignment
    TMatrix rev(boost::extents[m+1][n+1]);
    rev[0][0] = 0;
    for(std::size_t col = 1; col <= n; ++col) rev[0][col] = rev[0][col-1] + _horizontalGap(ac, 0, m, sc.ge);
    for(std::size_t row = 1; row <= m; ++row) rev[row][0] = rev[row-1][0] + _verticalGap(ac, 0, n, sc.ge);
    for(std::size_t row = 1; row <= m; ++row)
      for(std::size_t col = 1; col <= n; ++col)
	rev[row][col] = std::max(std::max(rev[row-1][col-1] + (sRev1[row-1] == sRev2[col-1] ? sc.match : sc.mismatch), rev[row-1][col] + _verticalGap(ac, col, n, sc.ge)), rev[row][col-1] + _horizontalGap(ac, row, m, sc.ge));

    if (mat[m][n] != rev[m][n]) {
      //std::cerr << "Warning: Alignment scores disagree!" << std::endl;
      return false;
    } else {
      // Find best join
      TMatrix bestMat(boost::extents[m+1][n+1]);
      for(std::size_t row = 0; row <= m; ++row) {
	bestMat[row][0] = mat[row][0];
	for(std::size_t col = 1; col <= n; ++col) {
	  if (mat[row][col] > bestMat[row][col-1]) bestMat[row][col] = mat[row][col];
	  else bestMat[row][col] = bestMat[row][col-1];
	}
      }
      TMatrix bestRev(boost::extents[m+1][n+1]);
      for(std::size_t row = 0; row <= m; ++row) {
	bestRev[row][0] = rev[row][0];
	for(std::size_t col = 1; col <= n; ++col) {
	  if (rev[row][col] > bestRev[row][col-1]) bestRev[row][col] = rev[row][col];
	  else bestRev[row][col] = bestRev[row][col-1];
	}
      }
      TScoreValue bestScore = mat[m][n];
      std::size_t consLeft = 0;
      std::size_t refLeft = 0;
      for(std::size_t row = 0; row<=m; ++row) {
	for(std::size_t col = 0; col<=n; ++col) {
	  if (bestMat[row][col]+bestRev[m-row][n-col] > bestScore) {
	    bestScore=bestMat[row][col]+bestRev[m-row][n-col];
	    consLeft = row;
	    refLeft = col;
	  }
	}
      }
      std::size_t consRight = m - consLeft;
      std::size_t refRight = 0;
      // Find right bound
      for(std::size_t right = 0; right<=(n-refLeft); ++right) {
	if (mat[consLeft][refLeft] + rev[consRight][right] == bestScore) {
	  refRight = right;
	}
      }

      // Debug best join
      /*
      TScoreValue bScore = mat[m][n];
      std::size_t cLeft = 0;
      std::size_t cRight = 0;
      std::size_t rLeft = 0;
      std::size_t rRight = 0;
      for(std::size_t fwdcut = 0; fwdcut<=m; ++fwdcut) {
	std::size_t revcut = m - fwdcut;
	// Iterate all valid collinear alignments on the reference
	for(std::size_t left = 0; left<=n; ++left) {
	  for(std::size_t right = 0; right<=(n-left); ++right) {
	    if (mat[fwdcut][left] + rev[revcut][right] > bScore) {
	      bScore = mat[fwdcut][left] + rev[revcut][right];
	      cLeft = fwdcut;
	      cRight = revcut;
	      rLeft = left;
	      rRight = right;
	    }
	  }
	}
      }
      std::cerr << mat[m][n] << ',' << bScore << ';' << s1.size() << ',' << s2.size() << ';' << cLeft << ',' << (m-cRight) << ',' << rLeft << ',' << (n-rRight) << std::endl;
      std::cerr << mat[m][n] << ',' << bestScore << ';' << s1.size() << ',' << s2.size() << ';' << consLeft << ',' << (m-consRight) << ',' << refLeft << ',' << (n-refRight) << std::endl;
      */

      // Better split found?
      if (bestScore == mat[m][n]) return false; // No split found

      // Trace-back fwd
      std::size_t rr = consLeft;
      std::size_t cc = refLeft;
      typedef std::vector<char> TTrace;
      TTrace trace;
      while ((rr>0) || (cc>0)) {
	if ((rr>0) && (mat[rr][cc] == mat[rr-1][cc] + _verticalGap(ac, cc, n, sc.ge))) {
	  --rr;
	  trace.push_back('v');
	} else if ((cc>0) && (mat[rr][cc] == mat[rr][cc-1] + _horizontalGap(ac, rr, m, sc.ge))) {
	  --cc;
	  trace.push_back('h');
	} else {
	  --rr;
	  --cc;
	  trace.push_back('s');
	}
      }
      TAlign fwd;
      _createAlignment(trace, s1.substr(0, consLeft), s2.substr(0, refLeft), fwd);

      // Trace-back rev
      rr = consRight;
      cc = refRight;
      typedef std::vector<char> TTrace;
      TTrace rtrace;
      while ((rr>0) || (cc>0)) {
	if ((rr>0) && (rev[rr][cc] == rev[rr-1][cc] + _verticalGap(ac, cc, n, sc.ge))) {
	  --rr;
	  rtrace.push_back('v');
	} else if ((cc>0) && (rev[rr][cc] == rev[rr][cc-1] + _horizontalGap(ac, rr, m, sc.ge))) {
	  --cc;
	  rtrace.push_back('h');
	} else {
	  --rr;
	  --cc;
	  rtrace.push_back('s');
	}
      }
      TAlign rvs;
      _createAlignment(rtrace, sRev1.substr(0, consRight), sRev2.substr(0, refRight), rvs);

      // Concat alignments
      std::size_t gapref = (n-refRight) - refLeft;
      std::size_t alilen = fwd.shape()[1] + rvs.shape()[1] + gapref;
      align.resize(boost::extents[2][alilen]);
      TAIndex jEnd = rvs.shape()[1];
      for(TAIndex i = 0; i < (TAIndex) fwd.shape()[0]; ++i) {
	TAIndex alicol = 0;
	for(;alicol < (TAIndex) fwd.shape()[1]; ++alicol) align[i][alicol]=fwd[i][alicol];
	for(TAIndex j = refLeft; j < (TAIndex) (n-refRight); ++j, ++alicol) {
	  if (i==0) align[i][alicol] = '-';
	  else align[i][alicol] = s2[j];
	}
	for(TAIndex j = 0; j < (TAIndex) rvs.shape()[1]; ++j, ++alicol) {
	  switch (rvs[i][jEnd-j-1]) {
	  case 'A': align[i][alicol] = 'T'; break;
	  case 'C': align[i][alicol] = 'G'; break;
	  case 'G': align[i][alicol] = 'C'; break;
	  case 'T': align[i][alicol] = 'A'; break;
	  case 'N': align[i][alicol] = 'N'; break;
	  case '-': align[i][alicol] = '-'; break;
	  default: break;
	  }
	}
      } 
    }
    return true;
  }
  
  
  template<typename TAlign1, typename TAlign2, typename TAlignConfig, typename TScoreObject>
  inline int
  needleScore(TAlign1 const& a1, TAlign2 const& a2, TAlignConfig const& ac, TScoreObject const& sc)
  {
    typedef typename TScoreObject::TValue TScoreValue;

    // DP Matrix
    std::size_t m = _size(a1, 1);
    std::size_t n = _size(a2, 1);
    std::vector<TScoreValue> s(n+1, 0);
    TScoreValue prevsub = 0;

    // Create profile
    typedef boost::multi_array<double, 2> TProfile;
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
	  prevsub = 0;
	} else if (row == 0) {
	  s[col] = _horizontalGap(ac, 0, m, col * sc.ge);
	} else if (col == 0) {
	  s[0] = _verticalGap(ac, 0, n, row * sc.ge);
	  if (row - 1 == 0) prevsub = 0;
	  else prevsub = _verticalGap(ac, 0, n, (row - 1) * sc.ge);
	} else {
	  // Recursion
	  TScoreValue prevprevsub = prevsub;
	  prevsub = s[col];
	  s[col] = std::max(std::max(prevprevsub + _score(a1, a2, p1, p2, row-1, col-1, sc), prevsub + _verticalGap(ac, col, n, sc.ge)), s[col-1] + _horizontalGap(ac, row, m, sc.ge));
	}
      }
    }

    // Score
    return s[n];
  }


  template<typename TAlignConfig, typename TScoreObject>
  inline int32_t
  needleBanded(std::string const& s1, std::string const& s2, TAlignConfig const& ac, TScoreObject const& sc)
  {
    typedef typename TScoreObject::TValue TScoreValue;

    // DP Matrix
    int32_t m = s1.size();
    int32_t n = s2.size();
    int32_t band = 100;
    int32_t lowBand = band;
    int32_t highBand = band;
    if (m < n) highBand += n - m;
    else lowBand += m - n;
    std::vector<TScoreValue> s(n+1, 0);
    TScoreValue prevsub = 0;
    TScoreValue prevprevsub = 0;

    // DP
    for(int32_t row = 0; row <= m; ++row) {
      for(int32_t col = std::max(0, row - lowBand); col <= std::min(n, row + highBand); ++col) {
	// Initialization
	if ((row == 0) && (col == 0)) {
	  s[0] = 0;
	  prevsub = 0;
	} else if (row == 0) {
	  s[col] = _horizontalGap(ac, 0, m, col * sc.ge);
	} else if (col == 0) {
	  s[0] = _verticalGap(ac, 0, n, row * sc.ge);
	  if (row - 1 == 0) prevsub = 0;
	  else prevsub = _verticalGap(ac, 0, n, (row - 1) * sc.ge);
	} else {
	  // Recursion
	  prevprevsub = prevsub;
	  prevsub = s[col];
	  if (col == row - lowBand) {
	    prevprevsub = s[col-1];
	    s[col - 1] = DELLY_OUTOFBAND;
	  } else if (col == row + highBand) prevsub = DELLY_OUTOFBAND;
	  s[col] = std::max(std::max(prevprevsub + (s1[row-1] == s2[col-1] ? sc.match : sc.mismatch), prevsub + _verticalGap(ac, col, n, sc.ge)), s[col-1] + _horizontalGap(ac, row, m, sc.ge));
	}
      }
    }
	
    // Score
    return s[n];
  }



  
  template<typename TAlign1, typename TAlign2, typename TAlign, typename TAlignConfig, typename TScoreObject>
  inline int
  needle(TAlign1 const& a1, TAlign2 const& a2, TAlign& align, TAlignConfig const& ac, TScoreObject const& sc)
  {
    typedef typename TScoreObject::TValue TScoreValue;

    // DP Matrix
    std::size_t m = _size(a1, 1);
    std::size_t n = _size(a2, 1);
    std::vector<TScoreValue> s(n+1, 0);
    TScoreValue prevsub = 0;

    // Trace Matrix
    std::size_t mf = n+1;
    typedef boost::dynamic_bitset<> TBitSet;
    TBitSet bit3( (m+1) * (n+1), false);
    TBitSet bit4( (m+1) * (n+1), false);
    
    // Create profile
    typedef boost::multi_array<double, 2> TProfile;
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
	  prevsub = 0;
	} else if (row == 0) {
	  s[col] = _horizontalGap(ac, 0, m, col * sc.ge);
	  bit3[col] = true;
	} else if (col == 0) {
	  s[0] = _verticalGap(ac, 0, n, row * sc.ge);
	  if (row - 1 == 0) prevsub = 0;
	  else prevsub = _verticalGap(ac, 0, n, (row - 1) * sc.ge);
	  bit4[row * mf] = true;
	} else {
	  // Recursion
	  TScoreValue prevprevsub = prevsub;
	  prevsub = s[col];
	  s[col] = std::max(std::max(prevprevsub + _score(a1, a2, p1, p2, row-1, col-1, sc), prevsub + _verticalGap(ac, col, n, sc.ge)), s[col-1] + _horizontalGap(ac, row, m, sc.ge));

	  // Trace
	  if (s[col] ==  s[col-1] + _horizontalGap(ac, row, m, sc.ge)) bit3[row * mf + col] = true;
	  else if (s[col] == prevsub + _verticalGap(ac, col, n, sc.ge)) bit4[row * mf + col] = true;
	}
      }
    }
	
    // Trace-back using pointers
    std::size_t row = m;
    std::size_t col = n;
    typedef std::vector<char> TTrace;
    TTrace trace;
    while ((row>0) || (col>0)) {
      if (bit3[row * mf + col]) {
	--col;
	trace.push_back('h');
      } else if (bit4[row * mf + col]) {
	--row;
	trace.push_back('v');
      } else {
	--row;
	--col;
	trace.push_back('s');
      }
    }

    // Create alignment
    _createAlignment(trace, a1, a2, align);

    // Score
    return s[n];
  }

  template<typename TAlign1, typename TAlign2, typename TAlign, typename TAlignConfig>
  inline int
  needle(TAlign1 const& a1, TAlign2 const& a2, TAlign& align, TAlignConfig const& ac)
  {
    DnaScore<int> dnasc;
    return needle(a1, a2, align, ac, dnasc);
  }

  template<typename TAlign1, typename TAlign2, typename TAlign>
  inline int
  needle(TAlign1 const& a1, TAlign2 const& a2, TAlign& align)
  {
    AlignConfig<false, false> ac;
    return needle(a1, a2, align, ac);
  }

}

#endif
