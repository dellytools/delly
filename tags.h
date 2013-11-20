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

#ifndef TAGS_H
#define TAGS_H

namespace torali {

  // Constants
  #define MAX_CHROM_SIZE 250000000


  // Tags
  struct DeletionTag;
  struct DuplicationTag;

  template<typename SvTag>
    struct SVType {
    };


  // Unique paired-end data structure for single chromosome only
  struct Hit {
    int32_t minPos;
    int32_t maxPos;
  
  Hit(BamTools::BamAlignment const& al) : minPos(std::min(al.Position, al.MatePosition)), maxPos(std::max(al.Position, al.MatePosition)) {}

    bool operator <(Hit const& other) const {
      return ((minPos<other.minPos) || ((minPos==other.minPos) && (maxPos<other.maxPos)));
    }
  };

  // Structural variant record
  struct StructuralVariantRecord {
    int svStartBeg;
    int svStartEnd;
    int svEndBeg;
    int svEndEnd;
    int svStart;
    int svEnd;
    int peSupport;
    int srSupport;
    int wiggle;
    double srAlignQuality;
    unsigned int id;
    bool precise;
    uint16_t peMapQuality;
    std::string chr;
    std::string consensus;
  };

  template<typename TSV>
  struct SortSVs : public std::binary_function<TSV, TSV, bool>
  {
    inline bool operator()(TSV const& sv1, TSV const& sv2) {
      return (sv1.chr<sv2.chr) || ((sv1.chr==sv2.chr) && (sv1.svStart<sv2.svStart));
    }
  };



}

#endif
