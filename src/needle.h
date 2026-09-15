#ifndef NEEDLE_H
#define NEEDLE_H

#define BOOST_DISABLE_ASSERTS
#include <boost/dynamic_bitset.hpp>
#include <boost/multi_array.hpp>
#include <iostream>
#include <cctype>
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


  struct SplitDPMatrix {
    static constexpr int32_t NEG = -(1 << 28);

    std::vector<int32_t> lo;
    std::vector<int32_t> hi;
    std::vector<std::size_t> off;
    std::vector<int32_t> val;

    inline void
    reset(int32_t const m) {
      lo.assign(m + 1, 1);
      hi.assign(m + 1, 0);
      off.assign(m + 1, 0);
      val.clear();
    }

    inline bool
    empty(int32_t const r) const {
      return (lo[r] > hi[r]);
    }

    inline int32_t
    get(int32_t const r, int32_t const c) const {
      if ((c < lo[r]) || (c > hi[r])) return NEG;
      return val[off[r] + (c - lo[r])];
    }

    inline void
    trim(int32_t const r, int32_t const first, int32_t const last, int32_t const thr) {
      std::size_t start = off[r];
      int32_t lokeep = first;
      while ((lokeep <= last) && (val[start + (lokeep - first)] < thr)) ++lokeep;
      int32_t hikeep = last;
      while ((hikeep >= lokeep) && (val[start + (hikeep - first)] < thr)) --hikeep;
      if (lokeep > hikeep) {
	lo[r] = 1;
	hi[r] = 0;
	val.resize(start);
      } else {
	lo[r] = lokeep;
	hi[r] = hikeep;
	off[r] = start + (lokeep - first);
	val.resize(start + (hikeep - first) + 1);
      }
    }
  };


  template<typename TAlignConfig, typename TScoreObject>
  inline void
  _splitFill(std::string const& s1, std::string const& s2, TAlignConfig const& ac, TScoreObject const& sc, int32_t const lowerBound, int32_t const gainMax, bool const prune, SplitDPMatrix& dp) {
    int32_t const NEG = SplitDPMatrix::NEG;
    int32_t m = s1.size();
    int32_t n = s2.size();
    dp.reset(m);
    dp.val.reserve(4 * (n + 1));

    // Initialization
    int32_t thr = (prune) ? (lowerBound - m * gainMax) : NEG;
    dp.off[0] = dp.val.size();
    int32_t v = 0;
    dp.val.push_back(v);
    for(int32_t col = 1; col <= n; ++col) {
      v += _horizontalGap(ac, 0, m, sc.ge);
      dp.val.push_back(v);
    }
    dp.trim(0, 0, n, thr);

    // Alignment
    for(int32_t row = 1; row <= m; ++row) {
      if (dp.empty(row - 1)) break;
      thr = (prune) ? (lowerBound - (m - row) * gainMax) : NEG;
      int32_t plo = dp.lo[row - 1];
      int32_t phi = dp.hi[row - 1];
      std::size_t poff = dp.off[row - 1];
      int32_t hg = _horizontalGap(ac, row, m, sc.ge);
      char const c1 = s1[row - 1];
      dp.off[row] = dp.val.size();
      // Vertical only
      int32_t left = dp.val[poff] + _verticalGap(ac, plo, n, sc.ge);
      dp.val.push_back(left);
      // Diagonal, vertical and horizontal
      for(int32_t col = plo + 1; col <= phi; ++col) {
	int32_t diag = dp.val[poff + (col - 1 - plo)] + ((c1 == s2[col - 1]) ? sc.match : sc.mismatch);
	int32_t up = dp.val[poff + (col - plo)] + _verticalGap(ac, col, n, sc.ge);
	left = std::max(std::max(diag, up), left + hg);
	dp.val.push_back(left);
      }
      int32_t col = phi + 1;
      if (col <= n) {
	// Diagonal and horizontal
	int32_t diag = dp.val[poff + (phi - plo)] + ((c1 == s2[col - 1]) ? sc.match : sc.mismatch);
	left = std::max(diag, left + hg);
	dp.val.push_back(left);
	++col;
	// Horizontal only
	while ((col <= n) && (left + hg >= thr)) {
	  left += hg;
	  dp.val.push_back(left);
	  ++col;
	}
      }
      dp.trim(row, plo, col - 1, thr);
    }
  }


  template<typename TAlignConfig, typename TScoreObject, typename TTrace>
  inline bool
  _splitTraceback(SplitDPMatrix const& dp, int32_t rr, int32_t cc, int32_t const m, int32_t const n, TAlignConfig const& ac, TScoreObject const& sc, TTrace& trace) {
    while ((rr>0) || (cc>0)) {
      int32_t v = dp.get(rr, cc);
      if (v == SplitDPMatrix::NEG) return false;
      if ((rr>0) && (v == dp.get(rr-1, cc) + _verticalGap(ac, cc, n, sc.ge))) {
	--rr;
	trace.push_back('v');
      } else if ((cc>0) && (v == dp.get(rr, cc-1) + _horizontalGap(ac, rr, m, sc.ge))) {
	--cc;
	trace.push_back('h');
      } else {
	if ((rr == 0) || (cc == 0)) return false;
	--rr;
	--cc;
	trace.push_back('s');
      }
    }
    return true;
  }


  template<typename TAlign, typename TAlignConfig, typename TScoreObject>
  inline bool
  longNeedle(std::string const& s1, std::string const& s2, TAlign& align, TAlignConfig const& ac, TScoreObject const& sc)
  {
    typedef typename TAlign::index TAIndex;
    int32_t const NEG = SplitDPMatrix::NEG;
    int32_t m = s1.size();
    int32_t n = s2.size();
    if ((m == 0) || (n == 0)) return false;

    // Reverse input sequences
    std::string sRev1 = s1;
    reverseComplement(sRev1);
    std::string sRev2 = s2;
    reverseComplement(sRev2);

    // Prune
    int32_t gainMax = std::max(std::max(sc.match, sc.mismatch), std::max(sc.ge, 0));
    int32_t penMax = std::max(std::max(-sc.match, -sc.mismatch), std::max(-sc.ge, 0));
    bool canPrune = (sc.ge <= 0);
    for(std::size_t i = 0; ((i < s1.size()) && (canPrune)); ++i) if (std::islower((unsigned char) s1[i])) canPrune = false;
    for(std::size_t i = 0; ((i < s2.size()) && (canPrune)); ++i) if (std::islower((unsigned char) s2[i])) canPrune = false;
    int32_t fullSlack = m * gainMax + (m + n) * penMax;

    SplitDPMatrix fmat;
    SplitDPMatrix rmat;
    std::vector<int32_t> pm;
    int32_t slack = 16;
    while (true) {
      bool prune = ((canPrune) && (slack < fullSlack));
      int32_t lowerBound = m * gainMax - slack;

      // Forward and reverse alignment
      _splitFill(s1, s2, ac, sc, lowerBound, gainMax, prune, fmat);
      _splitFill(sRev1, sRev2, ac, sc, lowerBound, gainMax, prune, rmat);
      if ((!prune) && (fmat.get(m, n) != rmat.get(m, n))) return false;

      // Find best join
      int32_t bestScore = fmat.get(m, n);
      int32_t consLeft = 0;
      int32_t refLeft = 0;
      for(int32_t row = 0; row <= m; ++row) {
	if ((fmat.empty(row)) || (rmat.empty(m - row))) continue;
	int32_t rlo = rmat.lo[m - row];
	int32_t rhi = rmat.hi[m - row];
	std::size_t roff = rmat.off[m - row];
	pm.resize(rhi - rlo + 1);
	pm[0] = rmat.val[roff];
	for(int32_t k = 1; k <= rhi - rlo; ++k) pm[k] = std::max(pm[k - 1], rmat.val[roff + k]);
	std::size_t foff = fmat.off[row];
	for(int32_t col = fmat.lo[row]; col <= fmat.hi[row]; ++col) {
	  int32_t right = n - col;
	  if (right < rlo) break;
	  if (right > rhi) right = rhi;
	  int32_t score = fmat.val[foff + (col - fmat.lo[row])] + pm[right - rlo];
	  if (score > bestScore) {
	    bestScore = score;
	    consLeft = row;
	    refLeft = col;
	  }
	}
      }
      if ((prune) && (bestScore < lowerBound)) {
	// Any split found below the bound is still a valid alignment, so its score is a valid bound
	slack *= 4;
	if ((bestScore > NEG / 2) && (m * gainMax - bestScore < slack)) slack = m * gainMax - bestScore;
	continue;
      }

      // Better split found?
      if (bestScore == fmat.get(m, n)) return false; // No split found

      // Find right bound
      int32_t consRight = m - consLeft;
      int32_t refRight = 0;
      int32_t base = fmat.get(consLeft, refLeft);
      int32_t rEnd = std::min(rmat.hi[consRight], n - refLeft);
      for(int32_t right = rmat.lo[consRight]; right <= rEnd; ++right) {
	if (base + rmat.get(consRight, right) == bestScore) refRight = right;
      }

      // Trace-back fwd and rev
      typedef std::vector<char> TTrace;
      TTrace trace;
      TTrace rtrace;
      if ((!_splitTraceback(fmat, consLeft, refLeft, m, n, ac, sc, trace)) || (!_splitTraceback(rmat, consRight, refRight, m, n, ac, sc, rtrace))) {
	if (!prune) return false;
	slack *= 4;
	continue;
      }
      TAlign fwd;
      _createAlignment(trace, s1.substr(0, consLeft), s2.substr(0, refLeft), fwd);
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
      return true;
    }
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
