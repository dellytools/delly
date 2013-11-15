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

#ifndef SAM_H
#define SAM_H

#include <iostream>

#include "record.h"

namespace torali
{


  template<typename TRecord>
    struct SortSamPairs : public std::binary_function<TRecord, TRecord, bool>
  {
    inline bool operator()(TRecord const& s1, TRecord const& s2) {
      int qname = s1.f0.compare(s2.f0);
      if (qname != 0) return (qname<0);
      int rname = s1.f2.compare(s2.f2);
      if (rname != 0) return (rname<0);
      if (abs(s1.f8) == abs(s2.f8)) return (s1.f1 & 0x0040);
      else return (abs(s1.f8) < abs(s2.f8));
    }
  };






template<typename TRecord>
class Sam {
 public:
  typedef std::vector<TRecord> TVecRecord;
  TVecRecord record;

  Sam() {}

  inline void readSam(std::string const& path) {
    // Clear old content
    record.clear();
    record.reserve(1000000);

    // Read sam file
    Memory_mapped_file map_file(path);
    char buffer[Memory_mapped_file::MAX_LINE_LENGTH];
    while (map_file.left_bytes() > 0) {
      map_file.read_line(buffer);
      // Skip header lines
      if (buffer[0] == '@') continue;
      Tokenizer token(buffer, Memory_mapped_file::MAX_LINE_LENGTH);
      TRecord a;
      addF0(token, a);
      addF1(token, a);
      addF2(token, a);
      addF3(token, a);
      addF4(token, a);
      addF5(token, a);
      addF6(token, a);
      addF7(token, a);
      addF8(token, a);
      addF9(token, a);
      record.push_back(a);
    }
  }

  template<typename TSort>
    inline void sortSam() 
    {
      std::sort(record.begin(), record.end(), TSort());
    }

  
  inline void clear() {
    record.clear();
  }
};


}

#endif
