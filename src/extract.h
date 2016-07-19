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
#include <fstream>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>


#include <htslib/faidx.h>

#include <zlib.h>

namespace torali {

  struct ExtractConfig {
    unsigned int start;
    unsigned int end;
    unsigned int linesize;
    bool breaks;
    bool closed;
    std::string chr;
    boost::filesystem::path outfile;
    boost::filesystem::path genome;
    boost::filesystem::path intervals;
  };


  struct GenomicInterval {
    std::string chr;
    std::string id;
    std::vector<unsigned int> breaks;
    
  GenomicInterval(std::string c, std::string i, std::vector<unsigned int> b) : chr(c), id(i), breaks(b) {}
  };


  inline 
    int runExtract(ExtractConfig const& c) 
  {
      typedef std::vector<GenomicInterval> TGenInt;
      TGenInt rec;
      if (boost::filesystem::exists(c.intervals) && boost::filesystem::is_regular_file(c.intervals) && boost::filesystem::file_size(c.intervals)) {
	std::ifstream iFile(c.intervals.string().c_str(), std::ifstream::in);
	if (iFile.is_open()) {
	  while (iFile.good()) {
	    std::string intervalLine;
	    getline(iFile, intervalLine);
	    typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	    boost::char_separator<char> sep(" \t;");
	    Tokenizer tokens(intervalLine, sep);
	    Tokenizer::iterator tokIter = tokens.begin();
	    if (tokIter!=tokens.end()) {
	      std::string chr = *tokIter++;
	      unsigned int start =  boost::lexical_cast<unsigned int>(*tokIter++);
	      unsigned int end =  boost::lexical_cast<unsigned int>(*tokIter++);
	      std::vector<unsigned int> br;
	      br.push_back(start - 1);
	      br.push_back(end - 1);
	      std::stringstream s;
	      s << chr << "_" << start << "_" << end;
	      std::string id = s.str();
	      if (tokIter!=tokens.end()) {
		id = *tokIter++;
		if ((c.breaks) && (tokIter!=tokens.end())) {
		  // Parse breaks
		  br.clear();
		  boost::char_separator<char> split(",");
		  Tokenizer bTok(*tokIter++, split);
		  Tokenizer::iterator bTokIter = bTok.begin();
		  for(;bTokIter!=bTok.end(); ++bTokIter) br.push_back(boost::lexical_cast<unsigned int>(*bTokIter));
		}
	      }
	      rec.push_back(GenomicInterval(chr, id, br));
	    }
	  }
	} 
      } else {
	std::vector<unsigned int> br;
	br.push_back(c.start - 1);
	br.push_back(c.end - 1);
	std::stringstream s;
	s << c.chr << "_" << c.start << "_" << c.end;
	rec.push_back(GenomicInterval(c.chr, s.str(), br));
      }

      // Outfile
      std::ofstream ofile(c.outfile.string().c_str());
      
      // Print genomic intervals
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Extracting Regions" << std::endl;
      boost::progress_display show_progress( rec.size() );
      unsigned int countProcessed = 0;
      faidx_t* fai = fai_load(c.genome.string().c_str());
      for(int32_t refIndex=0; ((refIndex < faidx_nseq(fai)) && (countProcessed<rec.size())); ++refIndex) {
	std::string seqname(faidx_iseq(fai, refIndex));
	int32_t seqlen = -1;
	char* seq = faidx_fetch_seq(fai, seqname.c_str(), 0, faidx_seq_len(fai, seqname.c_str()), &seqlen);
	TGenInt::const_iterator genIntIter = rec.begin();
	for(;genIntIter!=rec.end(); ++genIntIter) {
	  if (genIntIter->chr == seqname) {
	    ofile << ">" << genIntIter->id << std::endl;
	    typedef std::vector<unsigned int> TPosition;
	    TPosition::const_iterator posIt = genIntIter->breaks.begin();
	    TPosition::const_iterator posItEnd = genIntIter->breaks.end();
	    unsigned int count = 1;
	    for(;posIt!=posItEnd; ++posIt) {
	      unsigned int start= *posIt;
	      if (++posIt==posItEnd) break;
	      unsigned int end = *posIt; 
	      if (start <= end) {
		std::string fasta = boost::to_upper_copy(std::string(seq + start, seq + std::min((unsigned int) seqlen, end + (int) c.closed)));
		std::string::iterator itF = fasta.begin();
		std::string::iterator itFEnd = fasta.end();
		for (; itF!=itFEnd; ++itF, ++count) {
		  ofile << *itF;
		  if (count % c.linesize == 0) ofile << std::endl;
		}
	      } else {
		std::string fasta = boost::to_upper_copy(std::string(seq + end, seq + std::min((unsigned int) seqlen, start + (int) c.closed)));
		std::string::reverse_iterator itR = fasta.rbegin();
		std::string::reverse_iterator itREnd = fasta.rend();
		for(; itR!=itREnd; ++itR, ++count) {
		  switch (*itR) {
		  case 'A': ofile << 'T'; break;
		  case 'C': ofile << 'G'; break;
		  case 'G': ofile << 'C'; break;
		  case 'T': ofile << 'A'; break;
		  case 'N': ofile << 'N'; break;
		  default: break;
		  }
		  if (count % c.linesize == 0) ofile << std::endl;
		}
	      }
	    }
	    if ((count - 1) % c.linesize != 0) ofile << std::endl;
	    ++show_progress;
	    ++countProcessed;
	  }
	}
	if (seqlen) free(seq);
      }
      ofile.close();
      fai_destroy(fai);

      // End
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;

      return 0;
  }

}

#endif
