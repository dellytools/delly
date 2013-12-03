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

#ifndef EXTRACT_H
#define EXTRACT_H

#include <iostream>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "version.h"
#include "fasta_reader.h"

#define MAX_CHROM_SIZE 250000000

namespace torali {

  struct ExtractConfig {
    unsigned int start;
    unsigned int end;
    unsigned int linesize;
    unsigned int field_identifier;
    bool breaks;
    bool closed;
    bool revComp;
    std::string genome;
    std::string chr;
    std::string intervals;
    std::ostream& file;

    ExtractConfig() : start(0), end(1), linesize(60), field_identifier(0), breaks(false), closed(false), revComp(true), genome(""), chr(""), intervals(""), file(std::cout) {}

    template<typename TFile>
    ExtractConfig(TFile& f) : start(0), end(1), linesize(60), field_identifier(0), breaks(false), closed(false), revComp(true), genome(""), chr(""), intervals(""), file(f) {}

  };


  template<typename TConfig, typename TFastaRecord, typename TRecords>
    inline 
    void printIntervals(TConfig const& c, Fasta_reader<TFastaRecord> const& rFasta, TRecords const& rec, std::string*) 
  {
    typedef Fasta_reader<TFastaRecord> TFastaReader;
    typename TRecords::const_iterator recIt = rec.begin();
    typename TRecords::const_iterator recItEnd = rec.end();
    for(;recIt!=recItEnd; ++recIt) {
      typename TFastaReader::TFastaVec::const_iterator fastaIt = rFasta.faVec.begin();
      typename TFastaReader::TFastaVec::const_iterator fastaItEnd = rFasta.faVec.end();
      for(;fastaIt!=fastaItEnd;++fastaIt) {
	if (!recIt->f0.compare(fastaIt->chrName)) {	
	  c.file << ">" << recIt->f4 << std::endl;
	  typename std::string::const_iterator sIt = recIt->f3.begin();
	  typename std::string::const_iterator sItEnd = recIt->f3.end();
	  unsigned int start = 0;
	  unsigned int pos = 0;
	  typedef std::vector<unsigned int> TPosition;
	  TPosition position;
	  for(;sIt!=sItEnd;++sIt, ++pos) {
	    if (*sIt == ',') {
	      position.push_back(atoi(recIt->f3.substr(start, pos - start).c_str()));
	      start = pos + 1;
	    }
	  }
	  if (start<pos) position.push_back(atoi(recIt->f3.substr(start, pos - start).c_str()));
	  if (position.empty()) {
	    position.push_back(recIt->f1);
	    position.push_back(recIt->f2);
	  }
	  TPosition::const_iterator posIt = position.begin();
	  TPosition::const_iterator posItEnd = position.end();
	  unsigned int count = 1;
	  for(;posIt!=posItEnd; ++posIt) {
	    typename TFastaRecord::TSequence::const_iterator beg = fastaIt->seq.begin() + *posIt;
	    ++posIt;
	    if (posIt==posItEnd) break;
	    typename TFastaRecord::TSequence::const_iterator begEnd = fastaIt->seq.begin() + *posIt;
	    if (*(posIt-1) <= *posIt) {
	      for(;((beg<fastaIt->seq.end()) && (beg!=(begEnd + (int) c.closed)));++beg,++count) {
		c.file << dna5_decode[(int) *beg];
		if (count % c.linesize == 0) c.file << std::endl;
	      }
	    } else {
	      beg += (int) c.closed;
	      if (beg > fastaIt->seq.end()) beg=fastaIt->seq.end();
	      if (beg > fastaIt->seq.begin()) {
		do {
		  --beg;
		  int letter = (int) *beg;
		  if (c.revComp) {
		    switch (*beg) {
		    case 0: letter = 3; break;
		    case 1: letter = 2; break;
		    case 2: letter = 1; break;
		    case 3: letter = 0; break;
		    default: break;
		    }
		  }
		  c.file << dna5_decode[letter];
		  if (count % c.linesize == 0) c.file << std::endl;
		  ++count;
		} while ((beg>fastaIt->seq.begin()) && (beg!=begEnd));
	      }
	    }
	  }
	  if ((count - 1) % c.linesize != 0) c.file << std::endl;
	  break;
	}
      }
    }
  }
  

  template<typename TConfig, typename TFastaRecord, typename TRecords>
    inline 
    void printIntervals(TConfig const& c, Fasta_reader<TFastaRecord> const& rFasta, TRecords const& rec, void*) 
  {
    typedef Fasta_reader<TFastaRecord> TFastaReader;
    typename TRecords::const_iterator recIt = rec.begin();
    typename TRecords::const_iterator recItEnd = rec.end();
    for(;recIt!=recItEnd; ++recIt) {
      typename TFastaReader::TFastaVec::const_iterator fastaIt = rFasta.faVec.begin();
      typename TFastaReader::TFastaVec::const_iterator fastaItEnd = rFasta.faVec.end();
      for(;fastaIt!=fastaItEnd;++fastaIt) {
	if (!recIt->f0.compare(fastaIt->chrName)) {
	  c.file << ">" << recIt->f4 << std::endl;
	  unsigned int count = 1;
	  typename TFastaRecord::TSequence::const_iterator beg = fastaIt->seq.begin() + recIt->f1;
	  typename TFastaRecord::TSequence::const_iterator begEnd = fastaIt->seq.begin() + recIt->f2;
	  if (recIt->f1 <= recIt->f2) {
	    for(;((beg!=(begEnd + (int) c.closed)) && (beg!=fastaIt->seq.end())) ;++beg,++count) {
	      c.file << dna5_decode[(int) *beg];
	      if (count % c.linesize == 0) c.file << std::endl;
	    }
	  } else {
	    beg += (int) c.closed;
	    if (beg>=fastaIt->seq.end()) beg=fastaIt->seq.end();
	    do {
	      --beg;
	      int letter = (int) *beg;
	      if (c.revComp) {
		switch (*beg) {
		case 0: letter = 3; break;
		case 1: letter = 2; break;
		case 2: letter = 1; break;
		case 3: letter = 0; break;
		default: break;
		}
	      }
	      c.file << dna5_decode[letter];
	      if (count % c.linesize == 0) c.file << std::endl;
	      ++count;
	    } while (beg!=begEnd);
	  }
	  if ((count - 1) % c.linesize != 0) c.file << std::endl;
	  break;
	}
      }
    }
  }



  template<typename TToken, typename TRecord>
    inline
    void addPartialF3(TToken& token, TRecord& line, std::string*) {
    addF3(token, line);
    line.f4 = line.f3;
  }
  
  template<typename TToken, typename TRecord>
    inline
    void addPartialF3(TToken& token, TRecord& line, void*) {
    token.getString(line.f4);
  }
  
  template<typename TBreaks>
    inline 
    void runExtract(ExtractConfig const& c) 
    {
      TBreaks* tmp = 0;
      typedef Record<std::string, unsigned int, unsigned int, TBreaks, std::string, void, void, void, void, void, void, void> TRecord;
      typedef std::vector<TRecord> TVecRecord;
      TVecRecord rec;
      if (!c.intervals.empty()) {
	// Read reference interval file
	Memory_mapped_file map_file(c.intervals.c_str());
	char buffer[Memory_mapped_file::MAX_LINE_LENGTH];
	while (map_file.left_bytes() > 0) {
	  map_file.read_line(buffer);
	  Tokenizer token(buffer, Memory_mapped_file::MAX_LINE_LENGTH);
	  TRecord line;
	  addF0(token, line);
	  addF1(token, line);
	  addF2(token, line);
	  addPartialF3(token, line, tmp);
	  if (c.field_identifier != 0) {
	    if (c.field_identifier == 1) line.f4 = line.f0;
	    else if (c.field_identifier == 2) line.f4 = line.f1;
	    else if (c.field_identifier == 3) line.f4 = line.f2;
	    else if (c.field_identifier > 4) {
	      unsigned int k = 5;
	      for (;k < c.field_identifier;++k) token.skipNextChunk();
	      token.getString(line.f4);
	    }
	  } else {
	    std::stringstream s;
	    s << line.f0 << "_" << line.f1 << "_" << line.f2;
	    line.f4 = s.str();
	  }
	  line.f0.reserve(0);
	  line.f4.reserve(0);
	  rec.push_back(line);
	}
	map_file.close();
      } else {
	TRecord line;
	line.f0 = c.chr;
	line.f1 = c.start;
	line.f2 = c.end;
	std::stringstream s;
	s << line.f0 << "_" << line.f1 << "_" << line.f2;
	line.f4 = s.str();
	rec.push_back(line);
      }
      
      // Read fasta file
      typedef std::vector<char> TSequence;
      typedef FastaRecord<std::string, unsigned long, Dna5Alphabet, TSequence, void> TFastaRecord;
      typedef Fasta_reader<TFastaRecord> TFastaReader;
      TFastaReader rFasta;
      rFasta.read_fasta(c.genome.c_str());
      

      // Print intervals
      printIntervals(c, rFasta, rec, tmp);
    }

}

#endif
