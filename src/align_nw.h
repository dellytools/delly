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

#ifndef ALIGN_NW_H
#define ALIGN_NW_H

#include <iostream>
#include <list>


#include "alphabet.h"
#include "record.h"
#include "align_gotoh.h"

namespace torali
{

  template <typename TAlign, typename TString1, typename TString2, typename TTraceValue, typename TIndexPair>
    inline void
    _align_nw_trace(TAlign& align, TString1& str1, TString2& str2, std::vector<TTraceValue>& trace, TIndexPair& overallMaxIndex)
  {
    typedef unsigned int TSize;
    TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;
    TSize len1 = overallMaxIndex.first;
    TSize len2 = overallMaxIndex.second;
    TSize numCols = str1.size();
    TSize numRows = str2.size();
    if (len1 < numCols) _align_trace_print(align, str1, str2, len1, len2, numCols - len1, Horizontal);
    else if (len2 < numRows) _align_trace_print(align, str1, str2, len1, len2, numRows - len2, Vertical);
    if ((len1 != 0) && (len2 !=0)) {
      // Initialize everything
      TTraceValue tv = trace[(len1-1)*numRows + (len2-1)];
      TTraceValue tvOld = tv;

      TSize segLen = 1;
      if (tv == Diagonal) { --len1; --len2; }
      else if (tv == Horizontal) --len1;
      else if (tv == Vertical) --len2;

      // Now follow the trace
      if ((len1 != 0) && (len2 !=0)) {
	do {
	  tv = trace[(len1-1)*numRows + (len2-1)];
	  if (tv == Diagonal) {
	    if (tv != tvOld) {
	      _align_trace_print(align, str1, str2, len1, len2, segLen, tvOld);
	      tvOld = tv; segLen = 1;
	    } else ++segLen;
	    --len1; --len2;
	  } else if (tv == Horizontal) {
	    if (tv != tvOld) {
	      _align_trace_print(align, str1, str2, len1, len2, segLen, tvOld);
	      tvOld = tv; segLen = 1;
	    } else ++segLen;
	    --len1;
	  } else if (tv == Vertical) {
	    if (tv != tvOld) {
	      _align_trace_print(align, str1, str2, len1, len2, segLen, tvOld);
	      tvOld = tv; segLen = 1;
	    } else ++segLen;
	    --len2;
	  }
	} while ((len1 != 0) && (len2 !=0));
      }
      _align_trace_print(align, str1, str2, len1, len2, segLen, tvOld);
    }
    if (len1 != 0) _align_trace_print(align, str1, str2, (TSize) 0, (TSize) 0, (TSize) len1, Horizontal);
    else if (len2 != 0) _align_trace_print(align, str1, str2, (TSize) 0, (TSize) 0, (TSize) len2, Vertical);
  }



  template<typename TString1, typename TString2, typename TScore, typename TIndexPair, typename TTraceValue, typename TAlignConfig>
    inline int
    _align_nw(TString1& str1, TString2& str2, TScore& sc, TIndexPair& overallMaxIndex, std::vector<TTraceValue>& trace, TAlignConfig)
  {
    // Traceback values
    typedef unsigned int TSize;
    typedef int TScoreValue;
    TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;
    std::pair<int, int> overallMaxValue;


    // The DP Matrix
    typedef std::vector<TScoreValue> TColumn;
    TColumn column;

    // Initialization
    TSize len1 = str1.size();
    TSize len2 = str2.size();
    column.resize(len2+1);
    trace.resize(len1 * len2);

    // Classical DP
    typedef std::vector<TTraceValue> TTrace;
    typename TTrace::iterator it = trace.begin();
    overallMaxValue.first = std::numeric_limits<TScoreValue>::min();
    overallMaxValue.second = std::numeric_limits<TScoreValue>::min();
    overallMaxIndex.first = len1;
    overallMaxIndex.second = len2;
    typename TColumn::iterator coit = column.begin();
    typename TColumn::iterator col_end = column.end();
    *coit = 0;
    for(TSize row = 1; row <= len2; ++row) _initFirstColumn(TAlignConfig(), *(++coit), (TScoreValue) (row) * sc.gapExtendVertical(-1, row - 1, str1, str2));
    _lastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, column[len2], 0);

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
      max_verti = *coit;
      col2 = 0;

      for(;++coit != col_end; ++it) {
	// Get max for vertical, horizontal and diagonal                                                                                                                      
	max_verti += sc.gapExtendVertical(col, col2, str1, str2);
	max_hori = *coit + sc.gapExtendHorizontal(col, col2, str1, str2);
	max_diag = diagVal + sc.align(col, col2++, str1, str2); //compute the maximum in vertiVal
	diagVal = *coit;

	// Choose the max
	if (max_diag >= ((max_verti > max_hori) ? max_verti : max_hori)) {
	  *it = Diagonal;
	  max_verti = *coit = max_diag;
	} else if (max_hori >= max_verti) {
	  *it = Horizontal;
	  max_verti = *coit = max_hori;
	} else {
	  *it = Vertical;
	  *coit = max_verti;
	}
      }
      _lastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, max_verti, col+1);
    }
    _lastColumn(TAlignConfig(), overallMaxValue, overallMaxIndex, column);
    return _maxOfAlignment<TScoreValue>(TAlignConfig(), overallMaxValue, overallMaxIndex, len1, len2);
  }


  template<typename TAlign, typename TString1, typename TString2, typename TScore, typename TAlignConfig>
    inline int
    globalNwAlignment(TAlign& align, TString1& str1, TString2& str2, TScore& sc, TAlignConfig)
  {
    // Do the alignment
    std::pair<unsigned int, unsigned int> maxInd;
    std::vector<unsigned char> trace;
    int score = _align_nw(str1, str2, sc, maxInd, trace, TAlignConfig());

    // Do the traceback
    std::pair<std::list<char>, std::list<char> > alignCon;
    _align_nw_trace(alignCon, str1, str2, trace, maxInd);

    // Do the conversion
    _convert_align(align, alignCon);
    return score;
  }


}

#endif
