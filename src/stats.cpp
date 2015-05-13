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
  bool text_flag;
  uint32_t maxsize;
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

  // Get references
  BamTools::RefVector references = reader.GetReferenceData();
  TCount genomeLen = 0;
  BamTools::RefVector::const_iterator itRef = references.begin();
  for(int refIndex=0; itRef!=references.end(); ++itRef, ++refIndex) genomeLen += itRef->RefLength;
  TCount genomeDiploidLen = 2*genomeLen;

  // Vector of all paired-reads
  typedef std::vector<uint32_t> TVecSize;
  typedef std::vector<PairedMapping> TVecPaired;
  TVecPaired vecPairM;

  // Read length statistics
  TCount minReadLength=genomeDiploidLen;
  TCount sumReadLength=0;
  TCount maxReadLength=0;

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

  // Redundant paired read statistics
  TCount redundantPairs = 0;
  TCount chromRedundantPairs = 0;

  // Orientations
  TVecSize strandOrient(8);
  TVecSize orient(4);
  std::fill(strandOrient.begin(), strandOrient.end(), 0);
  std::fill(orient.begin(), orient.end(), 0);

  // Insert size counts by orientation
  int maxInsertSize=c.maxsize;
  TVecSize fPlus(maxInsertSize);
  TVecSize fMinus(maxInsertSize);
  TVecSize rPlus(maxInsertSize);
  TVecSize rMinus(maxInsertSize);
  std::fill(fPlus.begin(), fPlus.end(), 0);
  std::fill(fMinus.begin(), fMinus.end(), 0);
  std::fill(rPlus.begin(), rPlus.end(), 0);
  std::fill(rMinus.begin(), rMinus.end(), 0);

  // Read alignments
  int32_t oldRefID=-1;
  boost::progress_display show_progress( (references.end() - references.begin()) + 1 );
  BamTools::BamAlignment al;
  while( reader.GetNextAlignmentCore(al) ) {
    if ((al.AlignmentFlag & 0x0100) || (al.AlignmentFlag & 0x0200) || (al.AlignmentFlag & 0x0400) || (al.AlignmentFlag & 0x0800)) continue;
    if (oldRefID != al.RefID) {
      ++show_progress;
      oldRefID = al.RefID;

      // Get redundant pair counts for the old chromosome
      sort(vecPairM.begin(), vecPairM.end(), SortPairedMapping());
      TVecPaired::const_iterator pBeg = vecPairM.begin();
      TVecPaired::const_iterator pPrevBeg = pBeg;
      TVecPaired::const_iterator pEnd = vecPairM.end();
      if (pBeg!=pEnd) ++pBeg;
      for(;pBeg!=pEnd;++pBeg, ++pPrevBeg) {
	if ((pBeg->Position == pPrevBeg->Position) && (pBeg->MatePosition == pPrevBeg->MatePosition) && (pBeg->RefID == pPrevBeg->RefID) && (pBeg->MateRefID == pPrevBeg->MateRefID)) {
	  if ( pBeg->RefID == pBeg->MateRefID ) {
	    ++chromRedundantPairs;
	    ++redundantPairs;
	  }
	}
      }
      // Keep only the translocated pairs
      TVecPaired vecPairC;
      for(TVecPaired::const_iterator pIt=vecPairM.begin(); pIt!=vecPairM.end(); ++pIt) {
	if (pIt->RefID!=pIt->MateRefID) vecPairC.push_back(*pIt);
      }
      vecPairM = vecPairC;
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
	sumReadLength+=al.Length;
	if ((TCount) al.Length < minReadLength) minReadLength=al.Length;
	if ((TCount) al.Length > maxReadLength) maxReadLength=al.Length;
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
	      int iSize = 0;
	      switch(getStrandSpecificOrientation(al)) {
	      case 0:
		++strandOrient[0];
		++orient[0];
		iSize = al.InsertSize + al.Length;
		if (iSize < 0) iSize=0;
		if (iSize > maxInsertSize) iSize = maxInsertSize - 1;
		++fPlus[iSize];
		break;
	      case 1:
		++strandOrient[1];
		++orient[1];
		iSize = -al.InsertSize + al.Length;
		if (iSize < 0) iSize=0;
		if (iSize > maxInsertSize) iSize = maxInsertSize - 1;
		++fMinus[iSize];
		break;
	      case 2:
		++strandOrient[2];
		++orient[2];
		iSize = al.InsertSize;
		if (iSize < 0) iSize=0;
		if (iSize > maxInsertSize) iSize = maxInsertSize - 1;
		++rPlus[iSize];
		break;
	      case 3:
		++strandOrient[3];
		++orient[3];
		iSize = -al.InsertSize + al.Length + al.Length;
		if (iSize < 0) iSize=0;
		if (iSize > maxInsertSize) iSize = maxInsertSize - 1;
		++rMinus[iSize];
		break;
	      case 4:
		++strandOrient[4];
		++orient[2];
		iSize = -al.InsertSize;
		if (iSize < 0) iSize=0;
		if (iSize > maxInsertSize) iSize = maxInsertSize - 1;
		++rPlus[iSize];
		break;
	      case 5:
		++strandOrient[5];
		++orient[3];
		iSize = al.InsertSize + al.Length + al.Length;
		if (iSize < 0) iSize=0;
		if (iSize > maxInsertSize) iSize = maxInsertSize - 1;
		++rMinus[iSize];
		break;
	      case 6:
		++strandOrient[6];
		++orient[0];
		iSize = -al.InsertSize + al.Length;
		if (iSize < 0) iSize=0;
		if (iSize > maxInsertSize) iSize = maxInsertSize - 1;
		++fPlus[iSize];
		break;
	      case 7:
		++strandOrient[7];
		++orient[1];
		iSize = al.InsertSize + al.Length;
		if (iSize < 0) iSize=0;
		if (iSize > maxInsertSize) iSize = maxInsertSize - 1;
		++fMinus[iSize];
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
  // Get translocated redundant pairs and counts for the last chromosome
  sort(vecPairM.begin(), vecPairM.end(), SortPairedMapping());
  TVecPaired::const_iterator pBeg = vecPairM.begin();
  TVecPaired::const_iterator pPrevBeg = pBeg;
  TVecPaired::const_iterator pEnd = vecPairM.end();
  if (pBeg!=pEnd) ++pBeg;
  for(;pBeg!=pEnd;++pBeg, ++pPrevBeg) {
    if ((pBeg->Position == pPrevBeg->Position) && (pBeg->MatePosition == pPrevBeg->MatePosition) && (pBeg->RefID == pPrevBeg->RefID) && (pBeg->MateRefID == pPrevBeg->MateRefID)) {
      ++redundantPairs;
      if ( pBeg->RefID == pBeg->MateRefID ) ++chromRedundantPairs;
    }
  }
  TCount nonredundantCount = ((pairedReadMappedCount / 2) - redundantPairs) * 2; 
  TCount nonChromRedundantCount = ((pairedReadSameChr / 2) - chromRedundantPairs) * 2; 


  // Get insert size statistics
  TCount totalISizeCount = 0;
  for (TVecSize::const_iterator itI = fPlus.begin(); itI != fPlus.end(); ++itI) totalISizeCount+=*itI;
  for (TVecSize::const_iterator itI = fMinus.begin(); itI != fMinus.end(); ++itI) totalISizeCount+=*itI;
  for (TVecSize::const_iterator itI = rPlus.begin(); itI != rPlus.end(); ++itI) totalISizeCount+=*itI;
  for (TVecSize::const_iterator itI = rMinus.begin(); itI != rMinus.end(); ++itI) totalISizeCount+=*itI;
  TCount low99QIndex = totalISizeCount / 200; // Ignore 0.5% of the insert sizes
  TCount lowQIndex = totalISizeCount / 4;
  TCount medianIndex = totalISizeCount / 2;
  TCount upperQIndex = totalISizeCount - (totalISizeCount / 4);
  TCount upper99QIndex = totalISizeCount - (totalISizeCount / 200);
  totalISizeCount = 0;
  uint32_t status = 0;
  uint32_t iSize = 0;
  uint32_t low99QISize = 0;
  uint32_t lowQISize = 0;
  uint32_t medianISize = 0;
  uint32_t upperQISize = 0;
  uint32_t upper99QISize = 0;
  TVecSize::const_iterator itFP = fPlus.begin();
  TVecSize::const_iterator itFM = fMinus.begin();
  TVecSize::const_iterator itRP = rPlus.begin();
  TVecSize::const_iterator itRM = rMinus.begin();
  for(;itFP!=fPlus.end();++itFP, ++itFM, ++itRP, ++itRM, ++iSize) {
    totalISizeCount += (*itFP + *itFM + *itRP + *itRM);
    if ((status==0) && (totalISizeCount>=low99QIndex)) {
      low99QISize = iSize;
      ++status;
    }
    if ((status==1) && (totalISizeCount>=lowQIndex)) {
      lowQISize = iSize;
      ++status;
    }
    if ((status==2) && (totalISizeCount>=medianIndex)) {
      medianISize = iSize;
      ++status;
    }
    if ((status==3) && (totalISizeCount>=upperQIndex)) {
      upperQISize = iSize;
      ++status;
    }
    if ((status==4) && (totalISizeCount>=upper99QIndex)) {
      upper99QISize = iSize;
      ++status;
    }
    if (status==5) break;
  }


  // Open output file
  std::ofstream ofile(c.outfile.string().c_str());

  // Output statistics
  if (c.text_flag) {
    ofile << "SAM file: " << c.bam.string() << std::endl;
    ofile << "Haploid genome length: " << genomeLen << std::endl;
    ofile << "Diploid genome length: " << genomeDiploidLen << std::endl;
    ofile << std::endl;
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
    
    // Mapped read statistics
    ofile << "Mapped read statistics (" << mappedReadCount << " reads)" << std::endl;
    ofile << "Minimum read length: " << minReadLength  << std::endl;
    ofile << "Mean read length: " << (sumReadLength / mappedReadCount)  << std::endl;
    ofile << "Maximum read length: " << maxReadLength  << std::endl;
    ofile << "Haploid sequencing coverage: " << (sumReadLength) / (double) genomeLen << std::endl;
    ofile << "Diploid sequencing coverage: " << (sumReadLength) / (double) genomeDiploidLen << std::endl;
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
    
    
    ofile << "Insert size statistics of all paired reads mapped onto the same chromosome (" << pairedReadSameChr << " reads, " << pairedReadSameChr / 2 << " pairs)" << std::endl;
    ofile << "Insert size (lower quartile): " << lowQISize << std::endl;
    ofile << "Insert size (median): " << medianISize << std::endl;
    ofile << "Insert size (upper quartile): " << upperQISize << std::endl;
    ofile << "Haploid insert coverage: " << (medianISize * (pairedReadSameChr / 2)) / (double) genomeLen << std::endl;
    ofile << "Diploid insert coverage: " << (medianISize * (pairedReadSameChr / 2)) / (double) genomeDiploidLen << std::endl;
    ofile << std::endl;

    // Print orientations
    ofile << "Strand specific orientations" << std::endl;
    ofile << "FF+ (0): " << strandOrient[0] << std::endl; 
    ofile << "FF- (1): " << strandOrient[1] << std::endl; 
    ofile << "FR+ (2): " << strandOrient[2] << std::endl; 
    ofile << "FR- (3): " << strandOrient[3] << std::endl; 
    ofile << "RF+ (4): " << strandOrient[4] << std::endl; 
    ofile << "RF- (5): " << strandOrient[5] << std::endl; 
    ofile << "RR+ (6): " << strandOrient[6] << std::endl; 
    ofile << "RR- (7): " << strandOrient[7] << std::endl; 
    ofile << std::endl;
    ofile << "Strand independent orientations" << std::endl;
    ofile << "F+ (0): " << orient[0] << std::endl; 
    ofile << "F- (1): " << orient[1] << std::endl; 
    ofile << "R+ (2): " << orient[2] << std::endl; 
    ofile << "R- (3): " << orient[3] << std::endl; 
    ofile << std::endl;
  } else {
    ofile << "[Input]" << std::endl;
    ofile << "BamFile = " << c.bam.string() << std::endl;
    ofile << "HapGenomeLen = " << genomeLen << std::endl;
    ofile << "DipGenomeLen = " << genomeDiploidLen << std::endl;
    ofile << std::endl;

    ofile << "[Reads]" << std::endl;
    ofile << "MinReadLength = " << minReadLength  << std::endl;
    ofile << "MeanReadLength = " << (sumReadLength / mappedReadCount)  << std::endl;
    ofile << "MaxReadLength = " << maxReadLength  << std::endl;
    ofile << std::endl;

    ofile << "[Mapping]" << std::endl;
    ofile << "NumReads = " << totalReadCount << std::endl;  
    ofile << "NumMappedReads = " << mappedReadCount << std::endl;  
    ofile << "NumMappedReadsPerc = " << (mappedReadCount * 100 / (double) totalReadCount) << std::endl;  
    ofile << "NumUnmappedReads = " << (totalReadCount - mappedReadCount) << std::endl;  
    ofile << "NumUnmappedReadsPerc = " << ((totalReadCount - mappedReadCount) * 100 / (double) totalReadCount) << std::endl;  
    ofile << "NumReads1 = " << totalRead1 << std::endl;  
    ofile << "NumMappedReads1 = " << totalMappedRead1 << std::endl;  
    ofile << "NumMappedReads1Perc = " << (totalMappedRead1 * 100 / (double) totalRead1) << std::endl;  
    ofile << "NumUnmappedReads1 = " << (totalRead1 - totalMappedRead1) << std::endl;  
    ofile << "NumUnmappedReads1Perc = " << ((totalRead1 - totalMappedRead1) * 100 / (double) totalRead1) << std::endl;  
    ofile << "NumReads2 = " << totalRead2 << std::endl;  
    ofile << "NumMappedReads2 = " << totalMappedRead2 << std::endl;  
    ofile << "NumMappedReads2Perc = " << (totalMappedRead2 * 100 / (double) totalRead2) << std::endl;
    ofile << "NumUnmappedReads2 = " << (totalRead2 - totalMappedRead2) << std::endl;  
    ofile << "NumUnmappedReads2Perc = " << ((totalRead2 - totalMappedRead2) * 100 / (double) totalRead2) << std::endl;  
    ofile << std::endl;
    
    // Mapped read statistics
    ofile << "[Paired]" << std::endl;
    ofile << "NumMappedInPairs = " << pairedReadMappedCount << std::endl;
    ofile << "NumMappedInPairsPerc = " << (pairedReadMappedCount * 100 / (double) pairedReadCount) << std::endl;
    ofile << "NumMappedInPairsSameChr = " << pairedReadSameChr << std::endl;
    ofile << "NumMappedInPairsSameChrPerc = " << (pairedReadSameChr * 100 / (double) pairedReadCount) << std::endl;
    ofile << std::endl;

    ofile << "[Redundancy]" << std::endl;
    ofile << "NonRedundantReads = " << nonredundantCount << std::endl;
    ofile << "NonRedundantReadsPerc = " << (nonredundantCount * 100 / (double) pairedReadCount) << std::endl;
    ofile << "NonRedundantReadsSameChr = " << nonChromRedundantCount << std::endl;
    ofile << "NonRedundantReadsSameChrPerc = " << (nonChromRedundantCount * 100 / (double) pairedReadCount) << std::endl;
    ofile << std::endl;
    
    ofile << "[Insertsize]" << std::endl;
    ofile << "InsertSizeLowerQuartile = " << lowQISize << std::endl;
    ofile << "InsertSizeMedian = " << medianISize << std::endl;
    ofile << "InsertSizeUpperQuartile = " << upperQISize << std::endl;
    ofile << std::endl;

    ofile << "[Coverage]" << std::endl;
    ofile << "HapSeqCoverage = " << (sumReadLength) / (double) genomeLen << std::endl;
    ofile << "DipSeqCoverage = " << (sumReadLength) / (double) genomeDiploidLen << std::endl;
    ofile << "HapInsertCoverage = " << (medianISize * (pairedReadSameChr / 2)) / (double) genomeLen << std::endl;
    ofile << "DipInsertCoverage = " << (medianISize * (pairedReadSameChr / 2)) / (double) genomeDiploidLen << std::endl;
    ofile << std::endl;

    // Print orientations
    ofile << "[Orientation]" << std::endl;
    ofile << "F+ = " << orient[0] << std::endl; 
    ofile << "F- = " << orient[1] << std::endl; 
    ofile << "R+ = " << orient[2] << std::endl; 
    ofile << "R- = " << orient[3] << std::endl; 
    ofile << std::endl;
    
    ofile << "[Strandorientation]" <<std::endl;
    ofile << "FF+ = " << strandOrient[0] << std::endl; 
    ofile << "FF- = " << strandOrient[1] << std::endl; 
    ofile << "FR+ = " << strandOrient[2] << std::endl; 
    ofile << "FR- = " << strandOrient[3] << std::endl; 
    ofile << "RF+ = " << strandOrient[4] << std::endl; 
    ofile << "RF- = " << strandOrient[5] << std::endl; 
    ofile << "RR+ = " << strandOrient[6] << std::endl; 
    ofile << "RR- = " << strandOrient[7] << std::endl; 
    ofile << std::endl;
  }

  // Close statistics file
  ofile.close();

  // Insert histogram
  uint32_t binwidth=((upper99QISize-low99QISize)/200)+1;
  uint32_t binFP=0;
  uint32_t binFM=0;
  uint32_t binRP=0;
  uint32_t binRM=0;
  std::ofstream insfile(c.insertfile.string().c_str());
  iSize = 0;
  uint32_t lastBound = 0;
  if (low99QISize > 1000) {
    iSize = low99QISize;
    lastBound = low99QISize;
  }
  itFP = fPlus.begin() + iSize;
  itFM = fMinus.begin() + iSize;
  itRP = rPlus.begin() + iSize;
  itRM = rMinus.begin() + iSize;
  insfile << "Size\tCount\tOrientation" << std::endl;
  for(;((itFP!=fPlus.end()) && (iSize<=upper99QISize));++itFP, ++itFM, ++itRP, ++itRM, ++iSize) {
    if ((iSize % binwidth == 0) && (iSize != lastBound)) {
      insfile << lastBound << "\t" << binFP << "\tF+" << std::endl;
      insfile << lastBound << "\t" << binFM << "\tF-" << std::endl;
      insfile << lastBound << "\t" << binRP << "\tR+" << std::endl;
      insfile << lastBound << "\t" << binRM << "\tR-" << std::endl;      
      binFP=0;
      binFM=0;
      binRP=0;
      binRM=0;
      lastBound=iSize;
    }
    binFP+=*itFP;
    binFM+=*itFM;
    binRP+=*itRP;
    binRM+=*itRM;
  }
  insfile << lastBound << "\t" << binFP << "\tF+" << std::endl;
  insfile << lastBound << "\t" << binFM << "\tF-" << std::endl;
  insfile << lastBound << "\t" << binRP << "\tR+" << std::endl;
  insfile << lastBound << "\t" << binRM << "\tR-" << std::endl;      
  insfile.close();

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
    ("key,k", "use key-value output")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("stat.txt"), "statistics output file")
    ("max-insert,m", boost::program_options::value<unsigned int>(&c.maxsize)->default_value(15000), "max. plotting insert size")
    ("ifile,i", boost::program_options::value<boost::filesystem::path>(&c.insertfile)->default_value("ins.txt"), "insert size histogram")
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
  if (vm.count("key")) c.text_flag = false;
  else c.text_flag = true;
  
  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Run stats
  return run(c);
}




