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

#ifndef ALIGN_CONFIG_H
#define ALIGN_CONFIG_H

namespace torali
{

  // Align Config

  template<bool TTop = false, bool TLeft = false, bool TRight = false, bool TBottom = false>
    class AlignConfig;

  // 1 config
  template<>
    class AlignConfig<false, false, false, false> {};

  // 2 config
  template<>
    class AlignConfig<false, false, false, true> {};

  // 3 config
  template<>
    class AlignConfig<false, false, true, false> {};

  // 4 config
  template<>
    class AlignConfig<false, false, true, true> {};

  // 5 config
  template<>
    class AlignConfig<false, true, false, false> {};

  // 6 config
  template<>
    class AlignConfig<false, true, false, true> {};

  // 7 config
  template<>
    class AlignConfig<false, true, true, false> {};

  // 8 config
  template<>
    class AlignConfig<false, true, true, true> {};

  // 9 config
  template<>
    class AlignConfig<true, false, false, false> {};

  // 10 config
  template<>
    class AlignConfig<true, false, false, true> {};

  // 11 config
  template<>
    class AlignConfig<true, false, true, false> {};

  // 12 config
  template<>
    class AlignConfig<true, false, true, true> {};

  // 13 config
  template<>
    class AlignConfig<true, true, false, false> {};

  // 14 config
  template<>
    class AlignConfig<true, true, false, true> {};

  // 15 config
  template<>
    class AlignConfig<true, true, true, false> {};

  // 16 config
  template<>
    class AlignConfig<true, true, true, true> {};





  // Functions

  template<bool TTop, bool TRight, bool TBottom, typename TElement, typename TCost>
    inline void
    _initFirstColumn(AlignConfig<TTop, false, TRight, TBottom> const, TElement& element, TCost const cost)
  {
    element = cost;
  }


  template<bool TTop, bool TRight, bool TBottom, typename TSpec, typename TElement, typename TCost>
    inline void
    _initFirstColumn(AlignConfig<TTop, true, TRight, TBottom> const, TElement& element, TCost const)
  {
    element = 0;
  }


  template<bool TLeft, bool TRight, bool TBottom, typename TElement, typename TCost>
    inline void
    _initFirstRow(AlignConfig<false, TLeft, TRight, TBottom> const, TElement& element, TCost const cost)
  {
    element = cost;
  }


  template<bool TLeft, bool TRight, bool TBottom, typename TElement, typename TCost>
    inline void
    _initFirstRow(AlignConfig<true, TLeft, TRight, TBottom> const, TElement& element, TCost const)
  {
    element = 0;
  }


  template<bool TTop, bool TLeft, bool TRight, typename TValue1, typename TIndex1, typename TValue2, typename TIndex2>
    inline void
    _lastRow(AlignConfig<TTop, TLeft, TRight, false> const, TValue1&, TIndex1&, TValue2 const, TIndex2 const)
  {
  }

  template<bool TTop, bool TLeft, bool TRight, typename TValue1, typename TIndex1, typename TValue2, typename TIndex2>
    inline void
    _lastRow(AlignConfig<TTop, TLeft, TRight, true> const, TValue1& maxValue, TIndex1& maxIndex, TValue2 const val, TIndex2 const index)
  {
    if (val > maxValue.first) {
      maxValue.first = val;
      maxIndex.first = index;
    }
  }

  template<bool TTop, bool TLeft, bool TBottom, typename TValue1, typename TIndex1, typename TColumn>
    inline void
    _lastColumn(AlignConfig<TTop, TLeft, false, TBottom> const, TValue1& maxValue, TIndex1&, TColumn const& column)
  {
    maxValue.second = column[column.size() - 1];
  }


  template<bool TTop, bool TLeft, bool TBottom, typename TValue1, typename TIndex1, typename TColumn>
    inline void
    _lastColumn(AlignConfig<TTop, TLeft, true, TBottom> const, TValue1& maxValue, TIndex1& maxIndex, TColumn const& column)
  {
    typename TColumn::const_iterator itCol = column.begin();
    typename TColumn::const_iterator itColEnd = column.end();
    maxValue.second = column[column.size() - 1];
    for(unsigned int i = 0;itCol != itColEnd; ++i, ++itCol) {
      if (*itCol > maxValue.second) {
	maxValue.second = *itCol;
	maxIndex.second = i;
      }
    }
  }


  template<typename TScoreValue, bool TTop, bool TLeft, typename TValue, typename TIndex, typename TSize>
    inline TScoreValue
    _maxOfAlignment(AlignConfig<TTop, TLeft, false, false> const, TValue& maxValue, TIndex&, TSize const, TSize const)
  {
    return maxValue.second;
  }



  template<typename TScoreValue, bool TTop, bool TLeft, typename TValue, typename TIndex, typename TSize>
    inline TScoreValue
    _maxOfAlignment(AlignConfig<TTop, TLeft, true, false> const, TValue& maxValue, TIndex& maxIndex, TSize const len1, TSize const)
  {
    maxIndex.first = len1;
    return maxValue.second;
  }


  template<typename TScoreValue, bool TTop, bool TLeft, typename TValue, typename TIndex, typename TSize>
    inline TScoreValue
    _maxOfAlignment(AlignConfig<TTop, TLeft, false, true> const, TValue& maxValue, TIndex& maxIndex, TSize const, TSize const len2)
  {
    maxIndex.second = len2;
    return maxValue.first;
  }


  template<typename TScoreValue, bool TTop, bool TLeft, typename TValue, typename TIndex, typename TSize>
    inline TScoreValue
    _maxOfAlignment(AlignConfig<TTop, TLeft, true, true> const, TValue& maxValue, TIndex& maxIndex, TSize const len1, TSize const len2)
  {
    if (maxValue.second > maxValue.first) maxIndex.first = len1;
    else maxIndex.second = len2;
    return (maxValue.first > maxValue.second) ? maxValue.first : maxValue.second;
  }

}

#endif


