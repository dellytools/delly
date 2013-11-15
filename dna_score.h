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

#ifndef DNA_SCORE_H
#define DNA_SCORE_H

namespace torali
{

  template <typename TValue>
    class DnaScore {
  public:
    TValue match;
    TValue mismatch;
    TValue gapOpen;
    TValue gapExtend;
    
    DnaScore() : match(5), mismatch(-4), gapOpen(-10), gapExtend(-1) {}

      DnaScore(TValue m, TValue mm, TValue gO, TValue gE) : match(m), mismatch(mm), gapOpen(gO), gapExtend(gE) {}
      

      template <typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
	inline TValue
	gapOpenHorizontal(TPos1, TPos2, TSeq1 const &, TSeq2 const &)
      {
	return gapOpen;
      }

      template <typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
	inline TValue
	gapOpenVertical(TPos1, TPos2, TSeq1 const &, TSeq2 const &)
      {
	return gapOpen;
      }

      template <typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
        inline TValue
        gapExtendHorizontal(TPos1, TPos2, TSeq1 const &, TSeq2 const &)
      {
        return gapExtend;
      }

      template <typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
        inline TValue
        gapExtendVertical(TPos1, TPos2, TSeq1 const &, TSeq2 const &)
      {
        return gapExtend;
      }


      template <typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
	inline TValue
	align(TPos1 pos1, TPos2 pos2, TSeq1 const &seq1, TSeq2 const &seq2)
      {
	if (seq1[pos1] == seq2[pos2]) return match;
	else return mismatch;
      }
  };
}

#endif 





