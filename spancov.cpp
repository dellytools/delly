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
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/functional/hash.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include "api/BamMultiReader.h"
#include "api/BamReader.h"

#include "memory_mapped_file.h"
#include "bam_file_adaptor.h"
#include "version.h"
#include "util.h"
#include "tags.h"
#include "intervaltree.h"
#include "fasta_reader.h"


using namespace torali;

#define MAX_CHROM_SIZE 250000000


struct Config {
  bool mapq;
  unsigned short minMapQual;
  int bpWindowOffset;
  boost::filesystem::path int_file;
  boost::filesystem::path outfile;
  std::vector<boost::filesystem::path> files;
};

template<typename TPos, typename TQual>
struct HitInterval {
  TPos start;
  TPos end;
  TQual qual;
  
  HitInterval() {}

  HitInterval(TPos const s, TPos const e, TQual const q) : start(s+1), end(e+1), qual(q) {}
};


template<typename TPos>
struct HitInterval<TPos, void> {
  TPos start;
  TPos end;
  
  HitInterval() {}

  HitInterval(TPos const s, TPos const e) : start(s+1), end(e+1) {}

  template<typename TQual>
  HitInterval(TPos const s, TPos const e, TQual) : start(s+1), end(e+1) {}
};


template<typename THitInterval>
struct SortHits : public std::binary_function<THitInterval, THitInterval, bool>
{
  inline bool operator()(THitInterval const& hit1, THitInterval const& hit2) {
    return (hit1.start < hit2.start) || ((hit1.start == hit2.start) && (hit1.end < hit2.end));
  }
};


template<typename TPos, typename TQual, typename TArrayType>
inline void
_addReadAndBpCounts(std::vector<HitInterval<TPos, TQual> > const& hit_vector, TArrayType* bp_count)
{
  typedef std::vector<HitInterval<TPos, TQual> > THits;
  typename THits::const_iterator vecBeg = hit_vector.begin();
  typename THits::const_iterator vecEnd = hit_vector.end();
  
  // Add bp counts
  for(;vecBeg!=vecEnd; ++vecBeg) {
    TArrayType* bpPoint = &bp_count[vecBeg->start - 1];
    TArrayType* bpPointEnd = &bp_count[vecBeg->end];
    for(;bpPoint!=bpPointEnd; ++bpPoint) ++(*bpPoint);
  }
}


template<typename TPos, typename TString>
inline void
  _buildMAPQString(std::vector< HitInterval<TPos, void> >&, TPos const, TPos const, std::vector<TString>&)
{
  // Nothing to do
}

template<typename TQual, typename TPos, typename TString>
inline void
  _buildMAPQString(std::vector< HitInterval<TPos, TQual> >& hit_vector, TPos const posStart, TPos const posEnd, std::vector<TString>& str)
{
  typedef std::vector<TString> TStringVector;
  typedef HitInterval<TPos, TQual> THit;
  typedef std::vector< THit > THits;
  typename THits::const_iterator vecBeg = hit_vector.begin();
  typename THits::const_iterator vecEnd = hit_vector.end();

  // Initialize result vector
  str.resize(posEnd-posStart);
  std::fill(str.begin(), str.end(), "");

  // Add mapq counts
  int searchRange = posStart - 10000;
  if (searchRange < 0) searchRange=0;
  THit hit;
  hit.start=searchRange;
  hit.end=searchRange;
  typename THits::const_iterator vecIt = std::lower_bound(vecBeg, vecEnd, hit, SortHits<THit>());
  for(;vecIt!=vecEnd; ++vecIt) {
    if (vecIt->end < posStart) continue;
    if (vecIt->start > posEnd) break;
    for (int i = (vecIt->start - 1); i<vecIt->end; ++i) {
      if (i >= posStart && i<posEnd) {
	std::stringstream s;
	s << (unsigned short) vecIt->qual << ',';
	str[i-posStart].append(s.str());
      }
    }
  }
  // Remove trailing ,
  typename TStringVector::iterator itStr = str.begin();
  typename TStringVector::iterator itStrEnd = str.end();
  for(;itStr!=itStrEnd;++itStr) {
    if (itStr->size()) itStr->erase(itStr->size() - 1);
  }
}


template<typename TDefaultOrientation>
inline bool
_mateIsUpstream(TDefaultOrientation defOrient, bool firstRead, bool reverse) {
  if (firstRead) {
    if (reverse) {
      if (defOrient % 2 == 0) return false;
      else return true;
    } else {
      if (defOrient % 2 == 0) return true;
      else return false;
    }
  } else {
    if (reverse) {
      if ((defOrient==1) || (defOrient==2)) return false;
      else return true;
    } else {
      if ((defOrient==1) || (defOrient==2)) return true;
      else return false;
    }
  }
}   

template<typename THit>
inline int
run(Config const& c, THit)
{
  // Valid interval file?
  if (!(boost::filesystem::exists(c.int_file) && boost::filesystem::is_regular_file(c.int_file) && boost::filesystem::file_size(c.int_file))) {
    std::cerr << "Error: " << c.int_file.string() << " does not exist or it is empty." << std::endl;
    return -1;
  }

  // Create library objects
  typedef std::map<std::string, LibraryInfo> TLibraryMap;
  typedef std::map<std::string, TLibraryMap> TSampleLibrary;
  TSampleLibrary sampleLib;


  // Scan libraries first
  BamTools::RefVector references;
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    // Get a sample name
    std::string sampleName(c.files[file_c].stem().string());

    // Check that all input bam files exist
    BamTools::BamReader reader;
    if ( ! reader.Open(c.files[file_c].string()) ) {
      std::cerr << "Could not open input bam file: " << c.files[file_c].string() << std::endl;
      reader.Close();
      return -1;
    }
    
    // Check that all input bam files are indexed
    reader.LocateIndex();
    if ( !reader.HasIndex() ) {
      std::cerr << "Missing bam index file: " << c.files[file_c].string() << std::endl;
      reader.Close();
      return -1;
    }

    // Get references
    if (file_c==0) references = reader.GetReferenceData();

    // Get library parameters and overall maximum insert size
    TLibraryMap libInfo;
    getLibraryParams(c.files[file_c], libInfo, 0, 3);
    sampleLib.insert(std::make_pair(sampleName, libInfo));
  }

  // Read all SV intervals
  typedef std::vector<StructuralVariantRecord> TSVs;
  TSVs svs;
  std::map<unsigned int, std::string> idToName;
  unsigned int intervalCount=0;
  if (isValidFile(c.int_file.string())) {
    Memory_mapped_file interval_file(c.int_file.string().c_str());
    char interval_buffer[Memory_mapped_file::MAX_LINE_LENGTH];
    while (interval_file.left_bytes() > 0) {
      interval_file.read_line(interval_buffer);
      // Read single interval line
      StructuralVariantRecord sv;
      Tokenizer token(interval_buffer, Memory_mapped_file::MAX_LINE_LENGTH);
      std::string interval_rname;
      token.getString(sv.chr);
      sv.svStart = token.getUInt();
      sv.svEnd = token.getUInt();
      std::string svName;
      token.getString(svName);
      idToName.insert(std::make_pair(intervalCount, svName));
      sv.id = intervalCount++;
      svs.push_back(sv);
    }
    interval_file.close();
  }

  // Output file
  boost::iostreams::filtering_ostream dataOut;
  dataOut.push(boost::iostreams::gzip_compressor());
  dataOut.push(boost::iostreams::file_sink(c.outfile.c_str(), std::ios_base::out | std::ios_base::binary));

  // Process chromosome by chromosome
  std::cout << "Breakpoint spanning coverage annotation" << std::endl;
  boost::progress_display show_progress( (references.end() - references.begin()) );
  BamTools::RefVector::const_iterator  itRef = references.begin();
  for(int refIndex=0;itRef!=references.end();++itRef, ++refIndex) {
    ++show_progress;

    // Iterate all samples
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      // Create hit vector
      typedef std::vector<THit> THits;
      THits hit_vector;  

      // Store all spanning ranges
      typedef std::vector<THit> THits;
      THits normalSpan;  
      THits missingSpan;  

      // Get a sample name
      std::string sampleName(c.files[file_c].stem().string());
      TSampleLibrary::iterator sampleIt=sampleLib.find(sampleName);

      // Initialize bam file
      BamTools::BamReader reader;
      if ( ! reader.Open(c.files[file_c].string()) ) return -1;
      reader.LocateIndex();
      if ( !reader.HasIndex() ) return -1;

      // Unique pairs for the given sample
      typedef std::set<Hit> TUniquePairs;
      TUniquePairs unique_pairs;

      // Read alignments and hash qualities
      uint16_t* qualities = new uint16_t[(int)boost::math::pow<28>(2)];
      uint16_t* qualitiesEnd = qualities + (int) boost::math::pow<28>(2);
      std::fill(qualities, qualitiesEnd, 0);
      BamTools::BamAlignment al;
      if (reader.Jump(refIndex, 0)) {
	while( reader.GetNextAlignment(al) ) {
	  if (al.RefID != refIndex) break;
	  if (!(al.AlignmentFlag & 0x0001) || (al.AlignmentFlag & 0x0004) || (al.AlignmentFlag & 0x0008) || (al.AlignmentFlag & 0x0100) || (al.AlignmentFlag & 0x0200) || (al.AlignmentFlag & 0x0400) || (al.Position==al.MatePosition) || (al.RefID!=al.MateRefID) || (al.MapQuality < c.minMapQual)) continue;
	  // Get the library information
	  std::string rG = "DefaultLib";
	  al.GetTag("RG", rG);
	  TLibraryMap::iterator libIt=sampleIt->second.find(rG);
	  if (libIt->second.median == 0) continue; // Single-end library

	  // Get or store the mapping quality for the partner
	  if (al.Position<al.MatePosition) {
	    // Hash the quality
	    boost::hash<std::string> hashStr;
	    unsigned int index=((hashStr(al.Name) % (int)boost::math::pow<4>(2))<<24) + ((al.Position % (int)boost::math::pow<12>(2))<<12) + (al.MatePosition % (int)boost::math::pow<12>(2));

	    qualities[index]=al.MapQuality;
	  } else {
	    // Get the two mapping qualities
	    boost::hash<std::string> hashStr;
	    unsigned int index=((hashStr(al.Name) % (int)boost::math::pow<4>(2))<<24) + ((al.MatePosition % (int)boost::math::pow<12>(2))<<12) + (al.Position % (int)boost::math::pow<12>(2));
	    uint16_t pairQuality = std::min(qualities[index], al.MapQuality);

	    // Is it a unique pair
	    Hit hitPos(al);
	    TUniquePairs::const_iterator pos = unique_pairs.begin();
	    bool inserted;
	    boost::tie(pos, inserted) = unique_pairs.insert(hitPos);
	    if (inserted) {
	      // Insert the interval
	      unsigned int outerISize = (al.Position + al.Length) - al.MatePosition;
	      if ((getStrandIndependentOrientation(al) == libIt->second.defaultOrient) && (outerISize >= libIt->second.minNormalISize) && (outerISize <= libIt->second.maxNormalISize)) {
		// Normal spanning coverage
		normalSpan.push_back(THit(al.MatePosition, al.Position+al.Length, pairQuality));
	      } else if ((getStrandIndependentOrientation(al) != libIt->second.defaultOrient) || (outerISize >= libIt->second.median + 5 * libIt->second.mad)) {
		// Missing spanning coverage
		if (_mateIsUpstream(libIt->second.defaultOrient, (al.AlignmentFlag & 0x0040), (al.AlignmentFlag & 0x0010))) 
		  missingSpan.push_back(THit(al.Position, al.Position + libIt->second.median, pairQuality));
		else
		  missingSpan.push_back(THit(std::max(0, al.Position + al.Length - libIt->second.median), al.Position + al.Length, pairQuality));
		if (_mateIsUpstream(libIt->second.defaultOrient, !(al.AlignmentFlag & 0x0040), (al.AlignmentFlag & 0x0020)))
		  missingSpan.push_back(THit(al.MatePosition, al.MatePosition + libIt->second.median, pairQuality));
		else
		  missingSpan.push_back(THit(std::max(0, al.MatePosition + al.Length - libIt->second.median), al.MatePosition + al.Length, pairQuality));
	      }
	      ++libIt->second.unique_pairs;
	    } else {
	      ++libIt->second.non_unique_pairs;
	    }
	  }
	}
      }
      delete [] qualities;

      // Sort SAM records by start position
      sort(normalSpan.begin(), normalSpan.end(), SortHits<THit>());
      sort(missingSpan.begin(), missingSpan.end(), SortHits<THit>());

      // Declare the chromosome array
      typedef unsigned short TArrayType;
      TArrayType* normalCount = new TArrayType[MAX_CHROM_SIZE];
      TArrayType* missingCount = new TArrayType[MAX_CHROM_SIZE];
      std::fill(normalCount, normalCount + MAX_CHROM_SIZE, 0);
      std::fill(missingCount, missingCount + MAX_CHROM_SIZE, 0);
      _addReadAndBpCounts(normalSpan, normalCount);
      _addReadAndBpCounts(missingSpan, missingCount);

      // Write spanning coverage for all input intervals
      typename TSVs::const_iterator itSV = svs.begin();
      typename TSVs::const_iterator itSVEnd = svs.end();
      for(;itSV!=itSVEnd;++itSV) {
	if (itSV->chr == references[refIndex].RefName) {
	  // First breakpoint
	  int posStart = (itSV->svStart - c.bpWindowOffset < 0) ? 0 : (itSV->svStart - c.bpWindowOffset);
	  int posEnd = (c.bpWindowOffset) ? (itSV->svStart + c.bpWindowOffset) : (itSV->svStart + 1);
	  if (!c.mapq) {
	    TArrayType* normalCountPoint = &normalCount[posStart];
	    TArrayType* missingCountPoint = &missingCount[posStart];
	    for(int i=posStart; i<posEnd; ++i, ++normalCountPoint, ++missingCountPoint) dataOut << idToName.find(itSV->id)->second << "\t" << 0 << "\t" << itSV->chr << "\t" << i << "\t" << *normalCountPoint << "\t" << *missingCountPoint << std::endl;
	  } else {
	    std::vector<std::string> normalStr;
	    std::vector<std::string> missingStr;
	    _buildMAPQString(normalSpan, posStart, posEnd, normalStr);
	    _buildMAPQString(missingSpan, posStart, posEnd, missingStr);
	    for(int i=posStart; i<posEnd; ++i) {
	      dataOut << idToName.find(itSV->id)->second << "\t" << 0 << "\t" << itSV->chr << "\t" << i << "\t" << normalStr[i - posStart] << "\t" << missingStr[i-posStart] << std::endl;
	    }
	  }
	  
	  // Second breakpoint
	  posStart = (itSV->svEnd - c.bpWindowOffset < 0) ? 0 : (itSV->svEnd - c.bpWindowOffset);
	  posEnd = (c.bpWindowOffset) ? (itSV->svEnd + c.bpWindowOffset) : (itSV->svEnd + 1);
	  if (!c.mapq) {
	    TArrayType* normalCountPoint = &normalCount[posStart];
	    TArrayType* missingCountPoint = &missingCount[posStart];
	    for(int i=posStart; i<posEnd; ++i, ++normalCountPoint, ++missingCountPoint) dataOut << idToName.find(itSV->id)->second << "\t" << 1 << "\t" << itSV->chr << "\t" << i << "\t" << *normalCountPoint << "\t" << *missingCountPoint << std::endl;
	  } else {
	    std::vector<std::string> normalStr;
	    std::vector<std::string> missingStr;
	    _buildMAPQString(normalSpan, posStart, posEnd, normalStr);
	    _buildMAPQString(missingSpan, posStart, posEnd, missingStr);
	    for(int i=posStart; i<posEnd; ++i) {
	      dataOut << idToName.find(itSV->id)->second << "\t" << 1 << "\t" << itSV->chr << "\t" << i << "\t" << normalStr[i-posStart] << "\t" << missingStr[i-posStart] << std::endl;
	    }
	  }
	}
      }

      // Clean-up
      delete[] normalCount;
      delete[] missingCount;
    }
  }

  // Output library statistics
  std::cout << "Library statistics" << std::endl;
  TSampleLibrary::const_iterator sampleIt=sampleLib.begin();
  for(;sampleIt!=sampleLib.end();++sampleIt) {
    std::cout << "Sample: " << sampleIt->first << std::endl;
    TLibraryMap::const_iterator libIt=sampleIt->second.begin();
    for(;libIt!=sampleIt->second.end();++libIt) {
      std::cout << "RG: ID=" << libIt->first << ",Median=" << libIt->second.median << ",MAD=" << libIt->second.mad << ",Orientation=" << (int) libIt->second.defaultOrient << ",MinInsertSize=" << libIt->second.minNormalISize << ",MaxInsertSize=" << libIt->second.maxNormalISize << ",DuplicatePairs=" << libIt->second.non_unique_pairs << ",UniquePairs=" << libIt->second.unique_pairs << std::endl;
    }
  }
  
  // End
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
  return 0;
}


int main(int argc, char **argv) {
  Config c;

  // Define required options
  boost::program_options::options_description required("Required options");
  required.add_options()
    ("intervals,i", boost::program_options::value<boost::filesystem::path>(&c.int_file), "SV interval file")
    ;

  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("show-mapq,s", "use list of PE MAPQ instead of PE counts")
    ("quality-cut,q", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(0), "min. paired-end mapping quality")
    ("bp-offset,b", boost::program_options::value<int>(&c.bpWindowOffset)->default_value(1000), "breakpoint offset")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("span.gz"), "spanning coverage output file")
    ;

  // Define hidden options
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
    ;
  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  // Set the visibility
  boost::program_options::options_description cmdline_options;
  cmdline_options.add(required).add(generic).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(required).add(generic);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);


  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("intervals")) || (!vm.count("input-file"))) { 
    printTitle("Spanning coverage calculation");
    std::cout << "Usage: " << argv[0] << " [OPTIONS] <sample1.bam> <sample2.bam> ..." << std::endl;
    std::cout << visible_options << "\n"; 
    return -1; 
  }

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Do we need to include the qualities
  if (vm.count("show-mapq")) c.mapq=true;
  else c.mapq=false;

  // Run spanning coverage
  if (c.mapq) return run(c, HitInterval<int32_t, uint16_t>());
  else return run(c, HitInterval<int32_t, void>());
}
