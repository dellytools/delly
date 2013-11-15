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

#include <iostream>
#include <fstream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "api/BamMultiReader.h"

#include "memory_mapped_file.h"
#include "bam_file_adaptor.h"
#include "version.h"
#include "sam.h"
#include "fasta_reader.h"
#include "util.h"

using namespace torali;

#define MAX_CHROM_SIZE 250000000


struct Config {
  unsigned int window_size;
  unsigned int window_offset;
  unsigned int qual_cut;
  bool bp_flag;
  bool avg_flag;
  bool inclCigar;
  boost::filesystem::path outfile;
  boost::filesystem::path intervals_file;
  std::vector<std::string> files;
};


template<typename TChr, typename TPos, typename TCigar>
struct Hit {
  TChr chr;
  TPos pos;
  TCigar cigar;
  
  Hit() {}

  Hit(BamTools::BamAlignment const& al) : chr(al.RefID), pos(al.Position+1), cigar(cigarString(al.CigarData)) {}
};

template<typename TChr, typename TPos>
struct Hit<TChr, TPos, void> {
  TChr chr;
  TPos pos;

  Hit() {}

  Hit(BamTools::BamAlignment const& al) : chr(al.RefID), pos(al.Position+1) {}
};



template<typename TChr, typename TPos, typename TCigar, typename TIterator, typename TArrayType>
inline void
_addReadAndBpCounts(std::vector<Hit<TChr, TPos, TCigar> > const& hit_vector, TIterator const chrSEit, TArrayType* read_count, TArrayType* bp_count)
{
  typedef std::vector<Hit<TChr, TPos, TCigar> > THits;
  typename THits::const_iterator vecBeg = hit_vector.begin() + chrSEit->second.first;
  typename THits::const_iterator vecEnd = hit_vector.begin() + chrSEit->second.second;
  
  // Add read and bp counts
  for(;vecBeg!=vecEnd; ++vecBeg) {
    ++read_count[vecBeg->pos];
    int num_start = 0;
    int num_end = 0;
    TArrayType* bpPoint = &bp_count[vecBeg->pos];
    std::string::const_iterator cig = vecBeg->cigar.begin();
    std::string::const_iterator cigEnd = vecBeg->cigar.end();
    for(;cig!=cigEnd; ++cig, ++num_end) {
      if (((int) *cig >=48) && ((int) *cig <= 57)) continue;
      unsigned int len = atoi(vecBeg->cigar.substr(num_start, (num_end-num_start)).c_str());
      if (*cig == 'M') 
	for(unsigned int i = 0; i<len; ++i, ++bpPoint) ++(*bpPoint);
      else if ((*cig == 'N') || (*cig=='D')) bpPoint += len;
      num_start = num_end + 1;
    }
  }
}

template<typename TChr, typename TPos, typename TIterator, typename TArrayType>
inline void
_addReadAndBpCounts(std::vector<Hit<TChr, TPos, void> > const& hit_vector, TIterator const chrSEit, TArrayType* read_count, TArrayType*)
{
  typedef std::vector<Hit<TChr, TPos, void> > THits;
  typename THits::const_iterator vecBeg = hit_vector.begin() + chrSEit->second.first;
  typename THits::const_iterator vecEnd = hit_vector.begin() + chrSEit->second.second;
  
  // Add read and bp counts
  for(;vecBeg!=vecEnd; ++vecBeg) ++read_count[vecBeg->pos];
}

template<typename THit>
struct SortHits : public std::binary_function<THit, THit, bool>
{
  inline bool operator()(THit const& hit1, THit const& hit2) {
    if (hit1.chr != hit2.chr) return (hit1.chr < hit2.chr);
    else return (hit1.pos < hit2.pos);
  }
};




template<typename TChr, typename TPos, typename TCigar>
inline void
run(Config const& conf)
{

  BamTools::BamMultiReader reader;
  if ( ! reader.Open(conf.files) ) {
    std::cerr << "Could not open input bam files!" << std::endl;
    return;
  }

  // Built chromosome map
  const BamTools::RefVector references = reader.GetReferenceData();

  // Iterate all SAM files and collect all hits
  typedef Hit<TChr, TPos, TCigar> THit;
  typedef std::vector<THit> THits;
  THits hit_vector;  
  BamTools::BamAlignment al;
  while( reader.GetNextAlignmentCore(al) ) {
    if ((al.AlignmentFlag & 0x0004) || (al.AlignmentFlag & 0x0200) || (al.AlignmentFlag & 0x0400)) continue;
    if (al.MapQuality < conf.qual_cut) continue;
    if (al.RefID < (TPos) references.size()) {
      if (conf.inclCigar) hit_vector.push_back(THit(al));
      else {
	// For simple read counting only the primary alignment is counted
	if (!(al.AlignmentFlag & 0x0100)) hit_vector.push_back(THit(al));
      }
    }
  }

  // Sort SAM records by chromosome and position
  sort(hit_vector.begin(), hit_vector.end(), SortHits<THit>());

  // Get all chromosome ranges
  typedef std::pair<TPos, TPos> TStartEnd;
  typedef std::map<TChr, TStartEnd> TChrSE;
  TChrSE chrSE;
  typename THits::const_iterator vecBeg = hit_vector.begin();
  typename THits::const_iterator vecEnd = hit_vector.end();
  if (vecBeg != vecEnd) {
    TChr curChr = vecBeg->chr;
    TPos posStart = 0;
    TPos posEnd = 0;
    for(;vecBeg!=vecEnd; ++vecBeg, ++posEnd) {
      if (curChr != vecBeg->chr) {
	chrSE.insert(std::make_pair(curChr, std::make_pair(posStart, posEnd)));
	curChr = vecBeg->chr;
	posStart=posEnd;
      }
    }
    if (posStart != posEnd) chrSE.insert(std::make_pair(curChr, std::make_pair(posStart, posEnd)));
  }

  // Declare the chromosome array
  typedef unsigned short TArrayType;
  TArrayType* read_count = new TArrayType[MAX_CHROM_SIZE];
  TArrayType* bp_count = new TArrayType[MAX_CHROM_SIZE];
  TArrayType* bp_countEnd = bp_count + MAX_CHROM_SIZE;

  // Process each chromosome
  boost::iostreams::filtering_ostream dataOut;
  dataOut.push(boost::iostreams::gzip_compressor());
  dataOut.push(boost::iostreams::file_sink(conf.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
  std::cout << "Processing chromosomes..." << std::endl;
  typename TChrSE::const_iterator chrSEit = chrSE.begin();
  typename TChrSE::const_iterator chrSEitEnd = chrSE.end();
  for(;chrSEit != chrSEitEnd; ++chrSEit) {
    const TPos chrLen = references[(int) chrSEit->first].RefLength;
    const std::string chrName = references[(int) chrSEit->first].RefName;
    std::cout << chrName << " (Length: " << chrLen << "); HitRange: " << chrSEit->second.first << " - " << chrSEit->second.second << std::endl;

    // Set all values to zero
    std::fill(read_count, read_count + MAX_CHROM_SIZE, 0);
    std::fill(bp_count, bp_countEnd, 0);

    // Iterate all reads of that chromosome
    _addReadAndBpCounts(hit_vector, chrSEit, read_count, bp_count);

    // Write coverage windows
    if (!(boost::filesystem::exists(conf.intervals_file) && boost::filesystem::is_regular_file(conf.intervals_file) && boost::filesystem::file_size(conf.intervals_file))) {
	TPos pos = 0;
	while (pos < chrLen) {
	  TPos window_len = pos+conf.window_size;
	  if (window_len > chrLen) { window_len = chrLen; }
	  unsigned int bp_sum = 0;
	  unsigned int read_sum = 0;
	  TArrayType* bpPoint = &bp_count[pos+1];
	  TArrayType* readPoint = &read_count[pos+1];
	  for(TPos i=pos; i<window_len; ++i, ++bpPoint, ++readPoint) { 
	    bp_sum += *bpPoint;
	    read_sum += *readPoint;
	  }
	  dataOut << chrName << "\t" << pos << "\t" << window_len << "\t";
	  if (conf.avg_flag) dataOut << ( (bp_sum) / (double) (window_len - pos)) << "\t";
	  if (conf.bp_flag) dataOut << (bp_sum) << "\t";
	  dataOut << read_sum << std::endl;
	  pos += conf.window_offset;
	}
      } else {
      Memory_mapped_file interval_file(conf.intervals_file.string().c_str());
      char interval_buffer[Memory_mapped_file::MAX_LINE_LENGTH];
      while (interval_file.left_bytes() > 0) {
	interval_file.read_line(interval_buffer);

	// Read single interval line
	Tokenizer token(interval_buffer, Memory_mapped_file::MAX_LINE_LENGTH);
	std::string interval_rname;
	token.getString(interval_rname);
	if (interval_rname.compare(chrName)) continue;
	unsigned int posStart = token.getUInt();
	unsigned int posEnd = token.getUInt();
	unsigned int bp_sum = 0;
	unsigned int read_sum = 0;
	TArrayType* bpPoint = &bp_count[posStart];
	TArrayType* readPoint = &read_count[posStart];
	for(unsigned int i = posStart; i <= posEnd; ++i, ++bpPoint, ++readPoint) { 
	  bp_sum += *bpPoint; 
	  read_sum += *readPoint;
	}
	dataOut << interval_buffer << "\t";
	if (conf.avg_flag) dataOut << ( (bp_sum) / (double) (posEnd - posStart + 1)) << "\t";
	if (conf.bp_flag) dataOut << (bp_sum) << "\t";
	dataOut << read_sum << std::endl;
      }
      interval_file.close();
    }
  }

 
  // End
  std::cout << '[' << boost::posix_time::second_clock::local_time() << "] Done." << std::endl;
}


int main(int argc, char **argv) {
  Config c;
  c.bp_flag = false;
  c.avg_flag = false;
  c.inclCigar = false;

  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("bp-count,b", "show base pair count")
    ("avg-cov,a", "show average coverage")
    ("quality-cut,q", boost::program_options::value<unsigned int>(&c.qual_cut)->default_value(0), "exclude all alignments with quality < q")
    ("outfile,f", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("cov.gz"), "coverage output file")
    ;

  // Define window options
  boost::program_options::options_description window("Window options");
  window.add_options()
    ("window-size,s", boost::program_options::value<unsigned int>(&c.window_size)->default_value(10000), "window size")
    ("window-offset,o", boost::program_options::value<unsigned int>(&c.window_offset)->default_value(10000), "window offset")
    ;

  // Define interval options
  boost::program_options::options_description interval("Interval options");
  interval.add_options()
    ("interval-file,i", boost::program_options::value<boost::filesystem::path>(&c.intervals_file), "interval file")
    ;

  // Define hidden options
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value< std::vector<std::string> >(&c.files), "input file")
    ("license,l", "show license")
    ("warranty,w", "show warranty")
    ;
  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  // Set the visibility
  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(window).add(interval).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic).add(window).add(interval);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);


  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file"))) { 
    printTitle("Coverage calculation");
    if (vm.count("warranty")) {
      displayWarranty();
    } else if (vm.count("license")) {
      gplV3();
    } else {
      std::cout << "Usage: " << argv[0] << " [OPTIONS] <aligned1.bam> <aligned2.bam> ..." << std::endl;
      std::cout << visible_options << "\n"; 
    }
    return 1; 
  }
  if (vm.count("bp-count")) c.bp_flag = true;
  if (vm.count("avg-cov")) c.avg_flag = true;
  if ((c.bp_flag) || (c.avg_flag)) c.inclCigar = true;

  // Show cmd
  std::cout << '[' << boost::posix_time::second_clock::local_time() << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;
 

  if (isBinary(c.files[0])) {
    if (c.inclCigar) run<int32_t, int32_t, std::string>(c);
    else run<int32_t, int32_t, void>(c);
  } else {
    if (c.inclCigar) run<int32_t, int32_t, std::string>(c);
    else run<int32_t, int32_t, void>(c);
  }
  return 0;
}
