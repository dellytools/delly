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
#include <boost/progress.hpp>
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
#include <boost/filesystem.hpp>
#include "api/BamReader.h"
#include "api/BamIndex.h"

#include "version.h"
#include "util.h"
#include "tags.h"

using namespace torali;

// Config arguments
struct Config {
  bool iout;
  boost::filesystem::path bam;
  boost::filesystem::path outfile;
  boost::filesystem::path insertfile;
};


struct PairedMapping {
  int32_t     RefID;         
  int32_t     Position;      
  int32_t     MateRefID; 
  int32_t     MatePosition;
};


struct SortPairedMapping : public std::binary_function<PairedMapping,PairedMapping,bool>
{
  inline bool operator()(PairedMapping const& s1, PairedMapping const& s2) {
    if (s1.RefID!=s2.RefID) return (s1.RefID<s2.RefID);
    if (s1.MateRefID!=s2.MateRefID) return (s1.MateRefID<s2.MateRefID);
    return ((s1.Position <= s2.Position) && ((s1.Position < s2.Position) || (s1.MatePosition<s2.MatePosition)));
  }
};

inline int run(Config const& c) {
  typedef uint64_t TCount;

#ifdef PROFILE
  ProfilerStart("delly.prof");
#endif

  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();

  // Check that all input bam files exist and are indexed
  BamTools::BamReader reader;
  if ( ! reader.Open(c.bam.string()) ) {
    std::cerr << "Could not open input bam file: " << c.bam.string() << std::endl;
    reader.Close();
    return -1;
  }
  reader.LocateIndex();
  if ( !reader.HasIndex() ) {
    std::cerr << "Missing bam index file: " << c.bam.string() << std::endl;
    reader.Close();
    return -1;
  }

  // Open output file
  std::ofstream ofile(c.outfile.string().c_str());

  // Open insert output file
  boost::iostreams::filtering_ostream dataOut;
  if (c.iout) {
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.insertfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
  }

  // Get references
  BamTools::RefVector references = reader.GetReferenceData();
  uint64_t genomeLen = 0;
  BamTools::RefVector::const_iterator itRef = references.begin();
  ofile << "SAM file: " << c.bam.string() << std::endl;
  ofile << std::endl;
  ofile << "Reference information " << std::endl;
  for(int refIndex=0; itRef!=references.end(); ++itRef, ++refIndex) {
    ofile << itRef->RefName << " (" << itRef->RefLength << "bp); ";
    genomeLen += itRef->RefLength;
  }
  uint64_t genomeDiploidLen = 2*genomeLen;
  ofile << std::endl;
  ofile << "Haploid genome length: " << genomeLen << std::endl;
  ofile << "Diploid genome length: " << genomeDiploidLen << std::endl;
  ofile << std::endl;


  // Vector of all paired-reads
  typedef std::vector<PairedMapping> TVecPaired;
  TVecPaired vecPairM;

  // Vector of all read sizes
  typedef std::vector<unsigned int> TVecRSize;
  TVecRSize vecRSize;

  // Vector of all insert size
  typedef std::vector<unsigned int> TVecISize;
  TVecISize vecISize;

  // Single-end statistics
  TCount totalReadCount = 0;
  TCount mappedReadCount = 0;
  TCount totalRead1 = 0;
  TCount totalRead2 = 0;
  TCount totalMappedRead1 = 0;
  TCount totalMappedRead2 = 0;

  // Paired-end statistics
  TCount pairedReadCount = 0;
  TCount pairedReadMappedCount = 0;
  TCount pairedReadSameChr = 0;

  // Read alignments
  int32_t oldRefID=-1;
  boost::progress_display show_progress( (references.end() - references.begin()) + 1 );
  BamTools::BamAlignment al;
  while( reader.GetNextAlignment(al) ) {
    if (oldRefID != al.RefID) {
      ++show_progress;
      oldRefID = al.RefID;
    }
    if (!(al.AlignmentFlag & 0x0100) && !(al.AlignmentFlag & 0x0200) && !(al.AlignmentFlag & 0x0400) && !(al.AlignmentFlag & 0x0800)) {

      // Single-end statistics
      ++totalReadCount;
      if (al.AlignmentFlag & 0x0040) ++totalRead1;
      else ++totalRead2;
      if (!(al.AlignmentFlag & 0x0004)) {
	++mappedReadCount;
	if (al.AlignmentFlag & 0x0040) ++totalMappedRead1;
	else ++totalMappedRead2;
	vecRSize.push_back(al.QueryBases.size());
      }
      
      // Paired-end statistics
      if (al.AlignmentFlag & 0x0001) {
	++pairedReadCount;
	if (!((al.AlignmentFlag & 0x0004) || (al.AlignmentFlag & 0x0008))) {
	  ++pairedReadMappedCount;
	  // Same chr?
	  if (al.RefID == al.MateRefID) {
	    ++pairedReadSameChr;
	    if (al.AlignmentFlag & 0x0040) {
	      switch(getStrandSpecificOrientation(al)) {
	      case 0:
		if (c.iout) dataOut << al.Name << "\tFF+\tF+\t" << references[al.RefID].RefName << "\t" << al.Position << "\t" <<  references[al.MateRefID].RefName << "\t" <<  al.MatePosition << "\t" << (al.InsertSize + al.QueryBases.size()) << std::endl;
		vecISize.push_back(al.InsertSize + al.QueryBases.size());
		break;
	      case 1:
		if (c.iout) dataOut << al.Name << "\tFF-\tF-\t" << references[al.RefID].RefName << "\t" << al.Position << "\t" <<  references[al.MateRefID].RefName << "\t" <<  al.MatePosition << "\t" << (-al.InsertSize + al.QueryBases.size()) << std::endl;
		vecISize.push_back(-al.InsertSize + al.QueryBases.size());
		break;
	      case 2:
		if (c.iout) dataOut << al.Name << "\tFR+\tR+\t" << references[al.RefID].RefName << "\t" << al.Position << "\t" <<  references[al.MateRefID].RefName << "\t" <<  al.MatePosition << "\t" << al.InsertSize << std::endl;
		vecISize.push_back(al.InsertSize);
		break;
	      case 3:
		if (c.iout) dataOut << al.Name << "\tFR-\tR-\t" << references[al.RefID].RefName << "\t" << al.Position << "\t" <<  references[al.MateRefID].RefName << "\t" <<  al.MatePosition << "\t" << (-al.InsertSize + (al.QueryBases.size() + al.QueryBases.size())) << std::endl;
		vecISize.push_back(-al.InsertSize + (al.QueryBases.size() + al.QueryBases.size()));
		break;
	      case 4:
		if (c.iout) dataOut << al.Name << "\tRF+\tR+\t" << references[al.RefID].RefName << "\t" << al.Position << "\t" <<  references[al.MateRefID].RefName << "\t" <<  al.MatePosition << "\t" << -al.InsertSize << std::endl;
		vecISize.push_back(-al.InsertSize);
		break;
	      case 5:
		if (c.iout) dataOut << al.Name << "\tRF-\tR-\t" << references[al.RefID].RefName << "\t" << al.Position << "\t" <<  references[al.MateRefID].RefName << "\t" <<  al.MatePosition << "\t" << (al.InsertSize + (al.QueryBases.size() + al.QueryBases.size())) << std::endl;
		vecISize.push_back(al.InsertSize + (al.QueryBases.size() + al.QueryBases.size()));
		break;
	      case 6:
	        if (c.iout) dataOut << al.Name << "\tRR+\tF+\t" << references[al.RefID].RefName << "\t" << al.Position << "\t" <<  references[al.MateRefID].RefName << "\t" <<  al.MatePosition << "\t" << (-al.InsertSize + al.QueryBases.size()) << std::endl;
		vecISize.push_back(-al.InsertSize + al.QueryBases.size());
		break;
	      case 7:
		if (c.iout) dataOut << al.Name << "\tRR-\tF-\t" << references[al.RefID].RefName << "\t" << al.Position << "\t" <<  references[al.MateRefID].RefName << "\t" <<  al.MatePosition << "\t" << (al.InsertSize + al.QueryBases.size()) << std::endl;
		vecISize.push_back(al.InsertSize + al.QueryBases.size());
		break;
	      default:
		std::cerr << "False orientation." << std::endl;
		return -1;
	      }
	    }
	  }
	  // Collect all paired-end mapped reads
	  if (al.AlignmentFlag & 0x0040) {
	    PairedMapping pM;
	    if (al.Position < al.MatePosition) {
	      pM.Position = al.Position;
	      pM.MatePosition = al.MatePosition;
	      pM.RefID = al.RefID;
	      pM.MateRefID = al.MateRefID;
	    } else {
	      pM.Position = al.MatePosition;
	      pM.MatePosition = al.Position;
	      pM.RefID = al.MateRefID;
	      pM.MateRefID = al.RefID;
	    }
	    vecPairM.push_back(pM);
	  }
	}
      }
    }
  }

  // Get redundant paired reads
  TVecISize vecISizeBins; // 1kb windows
  vecISizeBins.resize(10);
  std::fill(vecISizeBins.begin(), vecISizeBins.end(), 0);
  sort(vecPairM.begin(), vecPairM.end(), SortPairedMapping());
  TVecPaired::const_iterator pBeg = vecPairM.begin();
  TVecPaired::const_iterator pPrevBeg = pBeg;
  TVecPaired::const_iterator pEnd = vecPairM.end();
  TCount redundantPairs = 0;
  TCount chromRedundantPairs = 0;
  if (pBeg!=pEnd) {
    unsigned int bin = (unsigned int) (abs((int) pBeg->Position - (int) pBeg->MatePosition) / 1000);
    if (bin >= vecISizeBins.size()) bin = vecISizeBins.size() - 1;
    ++vecISizeBins[bin];
    ++pBeg;
  }
  for(;pBeg!=pEnd;++pBeg, ++pPrevBeg) {
    if ((pBeg->Position == pPrevBeg->Position) && (pBeg->MatePosition == pPrevBeg->MatePosition) && (pBeg->RefID == pPrevBeg->RefID) && (pBeg->MateRefID == pPrevBeg->MateRefID)) {
      ++redundantPairs;
      if ( pBeg->RefID == pBeg->MateRefID ) ++chromRedundantPairs;
    } else if ( pBeg->RefID == pBeg->MateRefID ) {
      unsigned int bin = (unsigned int) (abs((int) pBeg->Position - (int) pBeg->MatePosition) / 1000);
      if (bin >= vecISizeBins.size()) bin = vecISizeBins.size() - 1;
      ++vecISizeBins[bin];
    }
  }
  TCount nonredundantCount = ((pairedReadMappedCount / 2) - redundantPairs) * 2; 
  TCount nonChromRedundantCount = ((pairedReadSameChr / 2) - chromRedundantPairs) * 2; 
  // Output statistics
  ofile << "Single-end statistics" << std::endl;
  ofile << "Total number of reads: " << totalReadCount << std::endl;  
  ofile << "Number of mapped reads: " << mappedReadCount << std::endl;  
  ofile << "Number of mapped reads (percentage): " << (mappedReadCount * 100 / (double) totalReadCount) << "%" << std::endl;  
  ofile << "Number of unmapped reads: " << (totalReadCount - mappedReadCount) << std::endl;  
  ofile << "Number of unmapped reads (percentage): " << ((totalReadCount - mappedReadCount) * 100 / (double) totalReadCount) << "%" << std::endl;  
  ofile << "Total number of read 1: " << totalRead1 << std::endl;  
  ofile << "Number of mapped reads 1: " << totalMappedRead1 << std::endl;  
  ofile << "Number of mapped reads 1 (percentage): " << (totalMappedRead1 * 100 / (double) totalRead1) << "%" << std::endl;  
  ofile << "Number of unmapped reads 1: " << (totalRead1 - totalMappedRead1) << std::endl;  
  ofile << "Number of unmapped reads 1 (percentage): " << ((totalRead1 - totalMappedRead1) * 100 / (double) totalRead1) << "%" << std::endl;  
  ofile << "Total number of read 2: " << totalRead2 << std::endl;  
  ofile << "Number of mapped reads 2: " << totalMappedRead2 << std::endl;  
  ofile << "Number of mapped reads 2 (percentage): " << (totalMappedRead2 * 100 / (double) totalRead2) << "%" << std::endl;
  ofile << "Number of unmapped reads 2: " << (totalRead2 - totalMappedRead2) << std::endl;  
  ofile << "Number of unmapped reads 2 (percentage): " << ((totalRead2 - totalMappedRead2) * 100 / (double) totalRead2) << "%" << std::endl;  
  ofile << std::endl;

  std::sort(vecRSize.begin(),vecRSize.end());
  ofile << "Mapped read statistics (" << mappedReadCount << " reads)" << std::endl;
  ofile << "Minimum read length: " << vecRSize.at(0)  << std::endl;
  ofile << "Median read length: " << vecRSize.at(vecRSize.size() / 2)  << std::endl;
  ofile << "Maximum read length: " << vecRSize.at(vecRSize.size() - 1)  << std::endl;
  ofile << "Haploid sequencing coverage: " << (vecRSize.at(vecRSize.size() - 1) * mappedReadCount) / (double) genomeLen << std::endl;
  ofile << "Diploid sequencing coverage: " << (vecRSize.at(vecRSize.size() - 1) * mappedReadCount) / (double) genomeDiploidLen << std::endl;
  ofile << std::endl;

  ofile << "Paired-end statistics" << std::endl;
  ofile << "Total number of reads paired in sequencing: " << pairedReadCount << std::endl;
  ofile << "Total number of reads where both reads of a pair mapped: " << pairedReadMappedCount << std::endl;
  ofile << "Total number of reads where both reads of a pair mapped (percentage): " << (pairedReadMappedCount * 100 / (double) pairedReadCount) << "%" << std::endl;
  ofile << "Total number of reads where both reads of a pair mapped onto the same chromosome: "<< pairedReadSameChr << std::endl;
  ofile << "Total number of reads where both reads of a pair mapped onto the same chromosome (percentage): " << (pairedReadSameChr * 100 / (double) pairedReadCount) << "%" << std::endl;
  ofile << "Total number of non-redundant reads paired in sequencing: " << nonredundantCount << std::endl;
  ofile << "Total number of non-redundant reads paired in sequencing (percentage): " << (nonredundantCount * 100 / (double) pairedReadCount) << "%" << std::endl;
  ofile << "Total number of non-redundant reads paired in sequencing and mapped to same chromosome: " << nonChromRedundantCount << std::endl;
  ofile << "Total number of non-redundant reads paired in sequencing and mapped to same chromosome (percentage): " << (nonChromRedundantCount * 100 / (double) pairedReadCount) << "%" << std::endl;
  ofile << std::endl;

  // Get medium insert size
  if (!vecISize.empty()) {
    std::sort(vecISize.begin(),vecISize.end());
    ofile << "Insert size statistics of all paired reads mapped onto the same chromosome (" << pairedReadSameChr << " reads, " << pairedReadSameChr / 2 << " pairs)" << std::endl;
    ofile << "Insert size (lower quartile): " << vecISize.at(vecISize.size() / 4) << std::endl;
    ofile << "Insert size (median): " << vecISize.at(vecISize.size() / 2) << std::endl;
    ofile << "Insert size (upper quartile): " << vecISize.at(vecISize.size() - vecISize.size() / 4) << std::endl;
    ofile << "Haploid insert coverage: " << (vecISize.at(vecISize.size() / 2) * (pairedReadSameChr / 2)) / (double) genomeLen << std::endl;
    ofile << "Diploid insert coverage: " << (vecISize.at(vecISize.size() / 2) * (pairedReadSameChr / 2)) / (double) genomeDiploidLen << std::endl;
    ofile << std::endl;
  }

  // Print insert size bins
  typename TVecISize::const_iterator binIt = vecISizeBins.begin();
  typename TVecISize::const_iterator binItEnd = vecISizeBins.end();
  unsigned int sumBin = 0;
  for(;binIt!=binItEnd; ++binIt) sumBin += *binIt;
  ofile << "Insert size counts for all non-redundant read pairs mapped to the same chromosome" << std::endl;
  binIt = vecISizeBins.begin();
  unsigned int incSumBin = 0;
  for(unsigned int binPos = 0;binIt!=binItEnd; ++binIt, ++binPos) {
    ofile << ">=" << binPos << "kb: " << (sumBin - incSumBin) << std::endl;
    incSumBin += *binIt;
  }

  // Close statistics file
  ofile.close();

#ifdef PROFILE
  ProfilerStop();
#endif

  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
  return 0;
}



  

int main(int argc, char **argv) {
  Config c;

  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("stat.txt"), "statistics output file")
    ("insert,i", boost::program_options::value<boost::filesystem::path>(&c.insertfile)->default_value(""), "gzip insert size output file")
    ;

  // Define hidden options
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bam), "input file")
    ("warranty,w", "show warranty")
    ("license,l", "show license")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);


  // Set the visibility
  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);

  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file"))) { 
    printTitle("DELLY");
    if (vm.count("warranty")) {
      displayWarranty();
    } else if (vm.count("license")) {
      gplV3();
    } else {
      std::cout << "Usage: " << argv[0] << " [OPTIONS] <sample.sort.bam>" << std::endl;
      std::cout << visible_options << "\n"; 
    }
    return 1; 
  }
  
  // Output insert sizes
  if (c.insertfile.string().size()) c.iout = true;
  else c.iout = false;

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Run stats
  return run(c);
}




