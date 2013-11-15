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

#ifndef ALIGN_GOTOH_H
#define ALIGN_GOTOH_H

#include <iostream>
#include <list>


#include "alphabet.h"
#include "record.h"

namespace torali
{

  template <typename TFile, typename TAlign>
    inline void
    _convert_align(TFile& file, std::pair<TAlign, TAlign>& alignCon)
  {
    typename TAlign::const_iterator itA = alignCon.first.begin();
    typename TAlign::const_iterator itAEnd = alignCon.first.end();
    for(;itA!=itAEnd; ++itA) file << *itA;
    file << std::endl;
    itA = alignCon.second.begin();
    itAEnd = alignCon.second.end();
    for(;itA!=itAEnd; ++itA) file << *itA;
    file << std::endl;
  }

  template <typename TFastaRecord, typename TAlign>
    inline void
    _convert_align(std::vector<TFastaRecord>& faVec, std::pair<TAlign, TAlign>& alignCon)
  {
    faVec.clear();
    faVec.resize(2);
    typename TAlign::const_iterator itA = alignCon.first.begin();
    typename TAlign::const_iterator itAEnd = alignCon.first.end();
    faVec[0].chrName = "Seq1";
    for(;itA!=itAEnd; ++itA) faVec[0].seq.push_back(*itA);
    faVec[0].chrLen = faVec[0].seq.size();
    faVec[1].chrName = "Seq2";
    itA = alignCon.second.begin();
    itAEnd = alignCon.second.end();
    for(;itA!=itAEnd; ++itA) faVec[1].seq.push_back(*itA);
    faVec[1].chrLen = faVec[1].seq.size();
  }

  template <typename TAlign, typename TString1, typename TString2, typename TPos, typename TTraceValue>
    inline void
    _align_trace_print(std::pair<TAlign, TAlign>& align, TString1& str1, TString2& str2, TPos pos1, TPos pos2, TPos segLen, TTraceValue tv)
  {
    // TraceBack values                                                                                                                                                                   
    TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;

    if (segLen == 0) return;

    if (tv == Horizontal) {
      for (int i = pos1 + segLen - 1; i>= (int) pos1;--i) {
	align.first.push_front(str1[i]);
	align.second.push_front(dna5gap_encode[(int) '-']);
      }
    }
    else if (tv == Vertical) {
      for (int i = pos2 + segLen - 1; i>= (int) pos2;--i) {
	align.first.push_front(dna5gap_encode[(int) '-']);
	align.second.push_front(str2[i]);
      }
    }
    else if (tv == Diagonal) {
      int j = pos2 + segLen - 1;
      for (int i = pos1 + segLen - 1; i>= (int) pos1;--i, --j) {
        align.first.push_front(str1[i]);
        align.second.push_front(str2[j]);
      }
    }
  }


  template <typename TAlign, typename TString1, typename TString2, typename TTraceValue, typename TIndexPair, typename TVal>
    inline void
    _align_gotoh_trace(TAlign& align, TString1& str1, TString2& str2, std::vector<TTraceValue>& trace, TIndexPair& overallMaxIndex, TVal initialDir)
  {
    typedef unsigned int TSize;
    TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;
    TSize len1 = overallMaxIndex.first;
    TSize len2 = overallMaxIndex.second;
    TSize numCols = str1.size();
    TSize numRows = str2.size();
    if (len1 < numCols) _align_trace_print(align, str1, str2, len1, len2, numCols - len1,  Horizontal);
    else if (len2 < numRows) _align_trace_print(align, str1, str2, len1, len2, numRows - len2,  Vertical);
    numRows = (numRows >> 1) + (numRows & 1);

    if ((len1 != 0) && (len2 !=0)) {

      // Initialize everything	
      TTraceValue nextTraceValue = (len2 & 1) ? trace[(len1 - 1)*numRows + ((len2 - 1) >> 1)] >> 4 : trace[(len1 - 1)*numRows + ((len2 - 1) >> 1)];
      TTraceValue tv = Diagonal;
      if (initialDir == Diagonal) tv = (nextTraceValue & 3);
      else if (initialDir == Horizontal) {
	if ((nextTraceValue >> 2) & 1) _align_trace_print(align, str1, str2, --len1, len2, (TSize) 1, Horizontal);
	else tv = Horizontal;
      } else if (initialDir == Vertical) {
	if ((nextTraceValue >> 3) & 1) _align_trace_print(align, str1, str2, len1, --len2, (TSize) 1, Vertical);
	else tv = Vertical;
      }
      TSize segLen = 0;
      TTraceValue tvOld = tv;

      // Now follow the trace
      do {
	nextTraceValue = (len2 & 1) ? trace[(len1 - 1)*numRows + ((len2 - 1) >> 1)] >> 4 : trace[(len1 - 1)*numRows + ((len2 - 1) >> 1)];
	if (tv == Diagonal) tv = (nextTraceValue & 3);
	else if (tv == Horizontal) {
	  if ((nextTraceValue >> 2) & 1) tv = Diagonal; 
	  else tv =  Horizontal;
	} else if (tv == Vertical) {
	  if ((nextTraceValue >> 3) & 1) tv =  Diagonal; 
	  else tv =  Vertical;
	}
	if (tv == Diagonal) {
	  if (tv != tvOld) {
	    if (tvOld == Vertical) --len2; 
	    else --len1;
	    _align_trace_print(align, str1, str2, len1, len2, (segLen + 1), tvOld);
	    tvOld = tv;
	    segLen = 0;
	  } else {
	    ++segLen;
	    --len1; --len2;
	  }
	} else if(tv == Horizontal) {
	  if (tv != tvOld) {
	    _align_trace_print(align, str1, str2, len1, len2, segLen, tvOld);
	    if ((nextTraceValue >> 2) & 1) {
	      _align_trace_print(align, str1, str2, --len1, len2, (TSize) 1,  Horizontal);
	      tv =  Diagonal; segLen = 0;
	    } else {
	      tvOld = tv; segLen = 1;
	      --len1;
	    }
	  } else {
	    ++segLen;
	    --len1;
	  }
	} else if (tv == Vertical) {
	  if (tv != tvOld) {
	    _align_trace_print(align, str1, str2, len1, len2, segLen, tvOld);
	    if ((nextTraceValue >> 3) & 1) {
	      _align_trace_print(align, str1, str2, len1, --len2, (TSize) 1,  Vertical);
	      tv =  Diagonal; segLen = 0;
	    } else {
	      tvOld = tv; segLen = 1;
	      --len2;
	    }
	  } else {
	    ++segLen;
	    --len2;
	  }
	}
      } while ((len1 != 0) && (len2 !=0));
      
      // Process left-overs
      if (segLen) _align_trace_print(align, str1, str2, len1, len2, segLen, tvOld);
    }
    
    // Handle the remaining sequence
    if (len1 != 0) _align_trace_print(align, str1, str2, (TSize) 0, (TSize) 0, (TSize) len1,  Horizontal);
    else if (len2 != 0) _align_trace_print(align, str1, str2, (TSize) 0, (TSize) 0, (TSize) len2,  Vertical);
  }



  template<typename TString1, typename TString2, typename TScore, typename TIndexPair, typename TTraceValue, typename TAlignConfig>
    inline int
    _align_gotoh(TString1& str1, TString2& str2, TScore& sc, TIndexPair& overallMaxIndex, std::vector<TTraceValue>& trace, TAlignConfig, TTraceValue& initialDir)
  {
    // Traceback values
    typedef unsigned int TSize;
    typedef int TScoreValue;
    TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;
    std::pair<int, int> overallMaxValue;


    // The DP Matrix for diagonal walks
    typedef std::vector<TScoreValue> TColumn;
    TColumn mat;
    // The DP Matrix for gaps from the left
    TColumn horizontal;
    // The DP Matrix for gaps from the top
    TScoreValue vert = 0;

    // Initialization
    TSize len1 = str1.size();
    TSize len2 = str2.size();
    mat.resize(len2+1);   // One column for the diagonal matrix
    horizontal.resize(len2+1);   // One column for the horizontal matrix
    trace.resize(len1 * ((len2 >> 1) + (len2 & 1)));
    TTraceValue tvMat = 0;
	
    // Classical DP
    typedef std::vector<TTraceValue> TTrace;
    typename TTrace::iterator it = trace.begin();
    overallMaxValue.first = std::numeric_limits<TScoreValue>::min();
    overallMaxValue.second = std::numeric_limits<TScoreValue>::min();
    overallMaxIndex.first = len1;
    overallMaxIndex.second = len2;
	
    typename TColumn::iterator matIt = mat.begin();
    *matIt = 0;
    typename TColumn::iterator horiIt = horizontal.begin();
    *horiIt = 0;
    
    TScoreValue a = 0;
    TScoreValue b = 0;
    TScoreValue max_val = 0;
    for(TSize row = 1; row <= len2; ++row) {
      _initFirstColumn(TAlignConfig(), *(++matIt), sc.gapOpenVertical(-1, 0, str1, str2) + (row - 1) * sc.gapExtendVertical(-1, row - 1, str1, str2));
      *(++horiIt) = *matIt + sc.gapOpenHorizontal(0, row-1, str1, str2) - sc.gapExtendHorizontal(0, row-1, str1, str2);
    }

    _lastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, *matIt, 0);
    if (overallMaxIndex.first == 0) initialDir = Vertical;
    for(TSize col = 1; col <= len1; ++col) {
      matIt = mat.begin();
      horiIt = horizontal.begin();
      TScoreValue diagValMat = *matIt;
      _initFirstRow(TAlignConfig(), *matIt, sc.gapOpenHorizontal(0, -1, str1, str2) + (col - 1) * sc.gapExtendHorizontal(col-1, -1, str1, str2));
      vert = *matIt + sc.gapOpenVertical(col-1, 0, str1, str2) - sc.gapExtendVertical(col-1, 0, str1, str2);
      TSize row = 1; 
      while(row <= len2) {
	// Get the new maximum for vertical
	a = *matIt + sc.gapOpenVertical(col-1, row-1, str1, str2);
	b = vert + sc.gapExtendVertical(col-1, row-1, str1, str2);
	if (a > b) {vert = a; *it |= 1;}
	else vert = b;

	// Get the new maximum for left
	*it <<= 1;
	a = *(++matIt) + sc.gapOpenHorizontal(col-1, row-1, str1, str2);
	b = *(++horiIt) + sc.gapExtendHorizontal(col-1, row-1, str1, str2);
	if (a > b) {*horiIt = a; *it |= 1;}
	else *horiIt = b;
			
	// Get the new maximum for mat
	*it <<= 2;
	max_val = diagValMat + sc.align(col-1, row-1, str1, str2);
	tvMat = Diagonal;
	if (vert > max_val) {
	  max_val = vert;
	  tvMat = Vertical;
	}
	if (*horiIt > max_val) {
	  max_val = *horiIt;
	  tvMat = Horizontal;
	}
	*it |= tvMat;

	// Assign the new diagonal values
	diagValMat = *matIt;
	*matIt = max_val;

	if (row & 1) *it <<= 1; else ++it;
	++row;
      }
      if (!(row & 1)) {*it <<= 3; ++it; }
      _lastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, *matIt, col);
      // If we got a new index, store direction
      if (overallMaxIndex.first == col) initialDir = tvMat;
    }
    _lastColumn(TAlignConfig(), overallMaxValue, overallMaxIndex, mat);
	
    // If we got a new index, store direction
    if ((overallMaxIndex.second != len2)  && (overallMaxValue.second > overallMaxValue.first)) 
      initialDir = (horizontal[overallMaxIndex.second] == mat[overallMaxIndex.second]) ? Horizontal : Diagonal;
    
    // If we end up in the bottom right corner, get direction
    if ((overallMaxIndex.first == len1) && (overallMaxIndex.second == len2)) {
      initialDir =  Diagonal;
      if (horizontal[len2] == mat[len2]) initialDir =  Horizontal;
      else if (vert == mat[len2]) initialDir =  Vertical;
    }
    return _maxOfAlignment<TScoreValue>(TAlignConfig(), overallMaxValue, overallMaxIndex, len1, len2);
  }


  template<typename TAlign, typename TString1, typename TString2, typename TScore, typename TAlignConfig>
    inline int
    globalGotohAlignment(TAlign& align, TString1& str1, TString2& str2, TScore& sc, TAlignConfig)
  {
    // Do the alignment
    std::pair<unsigned int, unsigned int> maxInd;
    std::vector<unsigned char> trace;
    unsigned char initialDir;
    int score = _align_gotoh(str1, str2, sc, maxInd, trace, TAlignConfig(), initialDir);

    // Do the traceback
    std::pair<std::list<char>, std::list<char> > alignCon;
    _align_gotoh_trace(alignCon, str1, str2, trace, maxInd, initialDir);

    // Do the conversion
    _convert_align(align, alignCon);

    return score;
  }

}

#endif
