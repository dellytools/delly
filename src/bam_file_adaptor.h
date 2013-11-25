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

#ifndef BAM_FILE_ADAPTOR_H
#define BAM_FILE_ADAPTOR_H

#include "memory_mapped_file.h"

#include <sstream>

#include "api/BamReader.h"
#include "api/BamWriter.h"


namespace torali
{

  template<typename TCigarVec>
    inline std::string
    cigarString(TCigarVec const& cigarOperations) 
    {
      std::stringstream cigar;
      typename TCigarVec::const_iterator coIter = cigarOperations.begin();
      if (coIter == cigarOperations.end()) cigar << "*";
      else
	for(; coIter != cigarOperations.end(); ++coIter) cigar << coIter->Length << coIter->Type;
      return cigar.str();
    }

  template<typename TCigarVec>
    inline void
    createCigarString(TCigarVec const& cigarOperations, std::stringstream& cigar) 
    {
      typename TCigarVec::const_iterator coIter = cigarOperations.begin();
      if (coIter == cigarOperations.end()) cigar << "*";
      else
	for(; coIter != cigarOperations.end(); ++coIter) cigar << coIter->Length << coIter->Type;
    }


class Bam_file_adaptor {
public:
  typedef unsigned long int TLen;

  BamTools::BamReader reader;
  BamTools::RefVector ref;
  BamTools::BamAlignment al;  

  Bam_file_adaptor(const std::string& in_path) {
    if ( !reader.Open(in_path) ) {
      std::string message = reader.GetErrorString();
      std::cerr << message << std::endl;
      exit(-1);
    }
    ref = reader.GetReferenceData();
  }

  inline TLen read_line(char* buffer, TLen num_bytes = Memory_mapped_file::MAX_LINE_LENGTH) 
    {
      std::stringstream cigar;
      createCigarString(al.CigarData, cigar);
      std::stringstream l;
      l << al.Name << "\t" << al.AlignmentFlag << "\t" << ((al.RefID == -1) ? "*" : ref[al.RefID].RefName) << "\t" << (al.Position + 1) << "\t" << al.MapQuality << "\t" << cigar.str() << "\t" << ((al.MateRefID == -1) ? "*" : ((al.RefID == al.MateRefID) ? "=" : ref[al.MateRefID].RefName)) << "\t" << (al.MatePosition + 1) << "\t" << al.InsertSize << "\t" << al.QueryBases << "\t" << al.Qualities;
      l.getline(buffer, num_bytes);
      return l.str().length();
    }

  inline TLen read_line(std::string& line) 
    {
      std::stringstream cigar;
      createCigarString(al.CigarData, cigar);
      std::stringstream l;
      l << al.Name << "\t" << al.AlignmentFlag << "\t" << ((al.RefID == -1) ? "*" : ref[al.RefID].RefName) << "\t" << (al.Position + 1) << "\t" << al.MapQuality << "\t" << cigar.str() << "\t" << ((al.MateRefID == -1) ? "*" : ((al.RefID == al.MateRefID) ? "=" : ref[al.MateRefID].RefName)) << "\t" << (al.MatePosition + 1) << "\t" << al.InsertSize << "\t" << al.QueryBases << "\t" << al.Qualities;
      line = l.str();
      return line.length();
    }

  inline TLen left_bytes()
    {
      return (reader.GetNextAlignment(al));
    }

  inline void close() 
    {
      reader.Close();
    }
};

}

#endif
