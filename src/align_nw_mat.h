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

#ifndef ALIGN_NW_MAT_H
#define ALIGN_NW_MAT_H

#include <iostream>
#include <list>


#include "alphabet.h"

namespace torali
{

  template<typename TString1, typename TString2, typename TScore, typename TAlignConfig, typename TMatrix>
    inline void
    _align_nw_mat(TString1& str1, TString2& str2, TScore& sc, TAlignConfig, TMatrix& mat) {
    // Traceback values
    typedef unsigned int TSize;
    typedef int TScoreValue;

    // The DP Matrix
    typedef std::vector<TScoreValue> TColumn;
    TColumn column;

    // Initialization
    TSize len1 = str1.size();
    TSize len2 = str2.size();
    column.resize(len2+1);
    mat.resize((len1 + 1) * (len2 + 1));

    // Classical DP
    typename TMatrix::iterator it = mat.begin();
    typename TColumn::iterator coit = column.begin();
    typename TColumn::iterator col_end = column.end();
    *coit = 0; *it = 0;
    for(TSize row = 1; row <= len2; ++row) {
      _initFirstColumn(TAlignConfig(), *(++coit), (TScoreValue) (row) * sc.gapExtendVertical(-1, row - 1, str1, str2));
      *(++it) = *coit;
    }

    // Classical DP                                                                                                                                                                       
    TScoreValue diagVal = 0;
    TScoreValue max_diag = 0;
    TScoreValue max_verti = 0;
    TScoreValue max_hori = 0;
    TSize col2 = 0;
    for(TSize col = 0; col < len1; ++col) {
      coit = column.begin();
      diagVal = *coit;
      _initFirstRow(TAlignConfig(), *coit, (TScoreValue) (col+1) * sc.gapExtendHorizontal(col, -1, str1, str2));
      *(++it) = *coit;
      max_verti = *coit;
      col2 = 0;

      for(;++coit != col_end; ) {
	// Get max for vertical, horizontal and diagonal                                                                                                                      
	max_verti += sc.gapExtendVertical(col, col2, str1, str2);
	max_hori = *coit + sc.gapExtendHorizontal(col, col2, str1, str2);
	max_diag = diagVal + sc.align(col, col2++, str1, str2); //compute the maximum in vertiVal
	diagVal = *coit;

	// Choose the max
	if (max_diag >= ((max_verti > max_hori) ? max_verti : max_hori)) max_verti = *coit = max_diag;
	else if (max_hori >= max_verti) max_verti = *coit = max_hori;
	else *coit = max_verti;
	*(++it) = *coit;
      }
    }
  }


  template<typename TString1, typename TString2, typename TScore, typename TAlignConfig, typename TMatrix>
    inline void
    globalNwAlignmentMatrix(TString1& str1, TString2& str2, TScore& sc, TAlignConfig, TMatrix& mat)
  {
    // Do the alignment
    _align_nw_mat(str1, str2, sc, TAlignConfig(), mat);
  }


}

#endif
