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

#ifndef ALIGN_H
#define ALIGN_H

#include <iostream>

namespace torali
{

  template<typename TProfile, typename TAIndex>
  inline int
  _score(std::string const& s1, std::string const& s2, TProfile const&, TProfile const&, TAIndex row, TAIndex col)
  {
    typedef int TScoreValue;
    TScoreValue match = 5;
    TScoreValue mismatch = -4;
    return (s1[row] == s2[col] ? match : mismatch);
  }

  template<typename TChar, typename TProfile, typename TAIndex>
  inline int
  _score(boost::multi_array<TChar, 2> const& a1, boost::multi_array<TChar, 2> const& a2, TProfile const& p1, TProfile const& p2, TAIndex row, TAIndex col)
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


  template<typename TProfile>
  inline void
  _createProfile(std::string const& s, TProfile& p)
  {
    typedef typename TProfile::index TPIndex;
    p.resize(boost::extents[6][s.size()]);   // 'A', 'C', 'G', 'T', 'N', '-'
    for (std::size_t j = 0; j < s.size(); ++j) {
      for(TPIndex k = 0; k < 6; ++k) p[k][j] = 0;
      if ((s[j] == 'A') || (s[j] == 'a')) p[0][j] += 1;
      else if ((s[j] == 'C') || (s[j] == 'c')) p[1][j] += 1;
      else if ((s[j] == 'G') || (s[j] == 'g')) p[2][j] += 1;
      else if ((s[j] == 'T') || (s[j] == 't')) p[3][j] += 1;
      else if ((s[j] == 'N') || (s[j] == 'n')) p[4][j] += 1;
      else if (s[j] == '-') p[5][j] += 1;
    }
  }

  template<typename TChar, typename TProfile>
  inline void
  _createProfile(boost::multi_array<TChar, 2> const& a, TProfile& p)
  {
    typedef typename boost::multi_array<TChar, 2>::index TAIndex;
    typedef typename TProfile::index TPIndex;
    p.resize(boost::extents[6][a.shape()[1]]);   // 'A', 'C', 'G', 'T', 'N', '-'
    for (TAIndex j = 0; j < (TAIndex) a.shape()[1]; ++j) {
      for(TPIndex k = 0; k < 6; ++k) p[k][j] = 0;
      int sum = 0;
      for(TAIndex i = 0; i < (TAIndex) a.shape()[0]; ++i) {
	++sum;
	if ((a[i][j] == 'A') || (a[i][j] == 'a')) p[0][j] += 1;
	else if ((a[i][j] == 'C') || (a[i][j] == 'c')) p[1][j] += 1;
	else if ((a[i][j] == 'G') || (a[i][j] == 'g')) p[2][j] += 1;
	else if ((a[i][j] == 'T') || (a[i][j] == 't')) p[3][j] += 1;
	else if ((a[i][j] == 'N') || (a[i][j] == 'n')) p[4][j] += 1;
	else if (a[i][j] == '-') p[5][j] += 1;
	else --sum;
      }
      for(TPIndex k = 0; k<6; ++k) p[k][j] /= sum;
    }
  }


  template<typename TTrace, typename TAlign>
  inline void
  _createAlignment(TTrace const& trace, std::string const& s1, std::string const& s2, TAlign& align)
  {
    align.resize(boost::extents[2][trace.size()]);
    std::size_t row = 0;
    std::size_t col = 0;
    std::size_t ai = 0;
    for(typename TTrace::const_reverse_iterator itT = trace.rbegin(); itT != trace.rend(); ++itT, ++ai) {
      if (*itT == 's') {
	align[0][ai] = s1[row];
	align[1][ai] = s2[col];
	//std::cerr << s1[row] << s2[col] << std::endl;
	++row;
	++col;
      } else if (*itT =='h') {
	align[0][ai] = '-';
	align[1][ai] = s2[col];
	//std::cerr << '-' << s2[col] << std::endl;
	++col;
      } else {
	align[0][ai] = s1[row];
	align[1][ai] = '-';
	//std::cerr << s1[row] << '-' << std::endl;
	++row;
      }
    }
  }


  template<typename TTrace, typename TChar, typename TAlign>
  inline void
  _createAlignment(TTrace const& trace, boost::multi_array<TChar, 2> const& a1, boost::multi_array<TChar, 2> const& a2, TAlign& align)
  {
    typedef typename TAlign::index TAIndex;
    TAIndex numN = a1.shape()[0];
    TAIndex numM = a2.shape()[0];
    align.resize(boost::extents[numN + numM][trace.size()]);
    TAIndex row=0;
    TAIndex col=0;
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
  }

}

#endif
