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

#ifndef DELLY_H
#define DELLY_H


#include <iostream>
#include <fstream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/icl/split_interval_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/functional/hash.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/sam.h>

#include "version.h"
#include "util.h"
#include "bolog.h"
#include "tags.h"
#include "spanning.h"
#include "coverage.h"
#include "junction.h"
#include "msa.h"
#include "split.h"
#include "pacbio.h"
#include "json.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <stdio.h>

namespace torali
{

// Config arguments
struct Config {
  unsigned short minMapQual;
  unsigned short minGenoQual;
  unsigned short madCutoff;
  int32_t minimumFlankSize;
  uint32_t graphPruning;
  uint32_t indelsize;
  float flankQuality;
  float percentAbnormal;
  bool indels;
  bool hasExcludeFile;
  bool hasVcfFile;
  DnaScore<int> aliscore;
  std::string svType;
  std::string technology;
  std::string format;
  boost::filesystem::path outfile;
  boost::filesystem::path vcffile;
  boost::filesystem::path genome;
  boost::filesystem::path exclude;
  std::vector<boost::filesystem::path> files;
  std::vector<std::string> sampleName;
};

// Reduced split alignment record
struct SplitAlignRecord {
  int32_t alignbeg;
  int32_t splitbeg;
  int32_t splitend;
  int32_t alignend;
  uint8_t MapQuality;
  
  SplitAlignRecord(int32_t ab, int32_t sb, int32_t se, int32_t ae, uint8_t mq) : alignbeg(ab), splitbeg(sb), splitend(se), alignend(ae), MapQuality(mq) {}
};

// Sort split alignment records
template<typename TRecord>
struct SortSplitRecords : public std::binary_function<TRecord, TRecord, bool>
{
  inline bool operator()(TRecord const& s1, TRecord const& s2) const {
    return ((s1.splitbeg < s2.splitbeg) || ((s1.splitbeg == s2.splitbeg) && (s1.splitend < s2.splitend)) || ((s1.splitbeg == s2.splitbeg) && (s1.splitend == s2.splitend) && (s1.alignbeg < s2.alignbeg)));
  }
};

// Reduced bam alignment record data structure
struct BamAlignRecord {
  int32_t tid;         
  int32_t pos;
  int32_t mtid; 
  int32_t mpos;
  int32_t alen;
  int32_t malen;
  int32_t Median;
  int32_t Mad;
  int32_t maxNormalISize;
  int libOrient;
  uint32_t flag;
  uint8_t MapQuality;
  
  BamAlignRecord(bam1_t* rec, uint8_t pairQuality, uint16_t a, uint16_t ma, int32_t median, int32_t mad, int32_t maxISize, int lO) : tid(rec->core.tid), pos(rec->core.pos), mtid(rec->core.mtid), mpos(rec->core.mpos), alen(a), malen(ma), Median(median), Mad(mad), maxNormalISize(maxISize), libOrient(lO), flag(rec->core.flag), MapQuality(pairQuality) {}
};

// Sort reduced bam alignment records
template<typename TRecord>
struct SortBamRecords : public std::binary_function<TRecord, TRecord, bool>
{
  inline bool operator()(TRecord const& s1, TRecord const& s2) const {
    if (s1.tid==s1.mtid) {
      return ((std::min(s1.pos, s1.mpos) < std::min(s2.pos, s2.mpos)) || 
	      ((std::min(s1.pos, s1.mpos) == std::min(s2.pos, s2.mpos)) && (std::max(s1.pos, s1.mpos) < std::max(s2.pos, s2.mpos))) ||
	      ((std::min(s1.pos, s1.mpos) == std::min(s2.pos, s2.mpos)) && (std::max(s1.pos, s1.mpos) == std::max(s2.pos, s2.mpos)) && (s1.maxNormalISize < s2.maxNormalISize)));
    } else {
      return ((s1.pos < s2.pos) ||
	      ((s1.pos == s2.pos) && (s1.mpos < s2.mpos)) ||
	      ((s1.pos == s2.pos) && (s1.mpos == s2.mpos) && (s1.maxNormalISize < s2.maxNormalISize)));
    }
  }
};

// Edge struct
template<typename TWeight, typename TVertex>
struct EdgeRecord {
  TVertex source;
  TVertex target;
  TWeight weight;

  EdgeRecord(TVertex s, TVertex t, TWeight w) : source(s), target(t), weight(w) {}
};

// Sort edge records
template<typename TRecord>
struct SortEdgeRecords : public std::binary_function<TRecord, TRecord, bool>
{
  inline bool operator()(TRecord const& e1, TRecord const& e2) const {
    return ((e1.weight < e2.weight) || ((e1.weight == e2.weight) && (e1.source < e2.source)) || ((e1.weight == e2.weight) && (e1.source == e2.source) && (e1.target < e2.target)));
  }
};




template<typename TValue, typename TPosition>
inline void
_movingAverage(std::vector<TValue> const& spp, TPosition const windowSize, TValue& movingAverage, TPosition& lowerBound, TPosition& upperBound) {
  movingAverage = 0;
  for(TPosition i = 0; (i<windowSize) && (i < (TPosition) spp.size()); ++i) movingAverage += spp[i];
  TValue bestAverage = movingAverage;
  TValue bestAverageIndex = windowSize - 1;
  for(std::size_t i = windowSize; i < spp.size() ; ++i) {
    movingAverage -= spp[i-windowSize];
    movingAverage += spp[i];
    if (movingAverage>bestAverage) {
      bestAverage = movingAverage;
      bestAverageIndex = i;
    }
  }
  movingAverage = bestAverage;
  upperBound = bestAverageIndex + 1;
  if (upperBound > windowSize) lowerBound = upperBound - windowSize;
  else lowerBound = 0;
}


template<typename TConfig, typename TStructuralVariantRecord, typename TTag>
inline bool
findPutativeSplitReads(TConfig const& c, std::vector<TStructuralVariantRecord>& svs,  SVType<TTag> svType) 
{
  typedef std::vector<TStructuralVariantRecord> TSVs;

  // Open file handles
  typedef std::vector<samFile*> TSamFile;
  typedef std::vector<hts_idx_t*> TIndex;
  TSamFile samfile;
  TIndex idx;
  samfile.resize(c.files.size());
  idx.resize(c.files.size());
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
    idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
  }
  bam_hdr_t* hdr = sam_hdr_read(samfile[0]);


  // Parse genome, no single-anchored reads anymore only soft-clipped reads
  unsigned int totalSplitReadsAligned = 0;
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Split-read alignment" << std::endl;
  boost::progress_display show_progress( hdr->n_targets );

  faidx_t* fai = fai_load(c.genome.string().c_str());
  // Find reference index
  for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
    ++show_progress;
    std::string tname(hdr->target_name[refIndex]);
    int32_t seqlen = -1;
    char* seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
    for(typename TSVs::iterator svIt = svs.begin();svIt!=svs.end(); ++svIt) {
      if ((svIt->chr != svIt->chr2) && (svIt->chr2 == refIndex)) {
	// For translocations temporarily store the first reference part in the consensus string
	svIt->consensus = _getSVRef(seq, *svIt, refIndex, svType);
      }
      if (svIt->chr == refIndex) {
	// Get the SV reference
	std::string svRefStr = _getSVRef(seq, *svIt, refIndex, svType);
	svIt->consensus = "";
	typedef std::vector<std::pair<int, std::string> > TOffsetSplit;
	typedef std::vector<int> TSplitPoints;
	TOffsetSplit osp0;
	TSplitPoints spp0;
	TOffsetSplit osp1;
	TSplitPoints spp1;
	
	// Find putative split reads in all samples
	for (unsigned int bpPoint = 0; bpPoint<2; ++bpPoint) {
	  int32_t regionChr = svIt->chr;
	  int regionStart = svIt->svStartBeg;
	  int regionEnd = (svIt->svStart + svIt->svStartEnd)/2;
	  if (bpPoint) {
	    regionChr = svIt->chr2;
	    regionStart = (svIt->svEndBeg + svIt->svEnd)/2;
	    regionEnd = svIt->svEndEnd;
	    spp1.resize(regionEnd-regionStart, 0);
	  } else {
	    spp0.resize(regionEnd-regionStart, 0);
	  }
#pragma omp parallel for default(shared)
	  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	    hts_itr_t* iter = sam_itr_queryi(idx[file_c], regionChr, regionStart, regionEnd);
	    bam1_t* rec = bam_init1();
	    while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	      if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	      if (rec->core.pos < regionStart) continue;
	      
	      // Valid soft clip?
	      int clipSize = 0;
	      int splitPoint = 0;
	      bool leadingSoftClip = false;
	      if (_validSoftClip(rec, clipSize, splitPoint, leadingSoftClip, c.minMapQual)) {
		if ((splitPoint >= regionStart) && (splitPoint < regionEnd)) {
		  splitPoint -= regionStart;
		  // Leading or trailing softclip?
		  if (_validSCOrientation(bpPoint, leadingSoftClip, svIt->ct, svType)) {
		    // Get the sequence
		    std::string sequence;
		    sequence.resize(rec->core.l_qseq);
		    uint8_t* seqptr = bam_get_seq(rec);
		    for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
		    
		    // Reverse complement iff necesssary
		    _adjustOrientation(sequence, bpPoint, svIt->ct, svType);
		    
		    if (bpPoint) {
#pragma omp critical
		      {
			++spp1[splitPoint];
			osp1.push_back(std::make_pair(splitPoint, sequence));
		      } 
		    } else {
#pragma omp critical
		      {
			++spp0[splitPoint];
			osp0.push_back(std::make_pair(splitPoint, sequence));
		      } 
		    }
		  }
		}
	      }
	    }
	    bam_destroy1(rec);
	    hts_itr_destroy(iter);		
	  }
	}
	// Collect candidate split reads
	typedef std::set<std::string> TSplitReadSet;
	TSplitReadSet splitReadSet;
	int mvAvg, lBound, uBound;
	_movingAverage(spp0, 5, mvAvg, lBound, uBound);
	if (mvAvg > 0) 
	  for(typename TOffsetSplit::const_iterator itOS = osp0.begin(); itOS != osp0.end(); ++itOS) 
	    if ((itOS->first >= lBound) && (itOS->first < uBound)) 
	      if (splitReadSet.size() < 100) splitReadSet.insert(itOS->second); // Limit to at most 100 split reads
	_movingAverage(spp1, 5, mvAvg, lBound, uBound);
	if (mvAvg > 0) {
	  for(typename TOffsetSplit::const_iterator itOS = osp1.begin(); itOS != osp1.end(); ++itOS) 
	    if ((itOS->first >= lBound) && (itOS->first < uBound)) 
	      if (splitReadSet.size() < 100) splitReadSet.insert(itOS->second); // Limit to at most 100 split reads
	}
	totalSplitReadsAligned += splitReadSet.size();
	
	// MSA
	if (splitReadSet.size() > 1) svIt->srSupport = msa(c, splitReadSet, svIt->consensus);
	
	// Search true split in candidates
	if (!alignConsensus(c, *svIt, svRefStr, svType)) { svIt->consensus = ""; svIt->srSupport = 0; }
      }
    }
    if (seqlen) free(seq);
  }
  // Clean-up
  fai_destroy(fai);
  bam_hdr_destroy(hdr);
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    hts_idx_destroy(idx[file_c]);
    sam_close(samfile[file_c]);
  }

  return (totalSplitReadsAligned>0);
}


// Initialize clique, deletions
template<typename TBamRecord, typename TSize>
inline void
_initClique(TBamRecord const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<DeletionTag>) {
  svStart = el.mpos + el.malen;
  svEnd = el.pos;
  wiggle =  -el.maxNormalISize;
}

// Initialize clique, insertions
template<typename TBamRecord, typename TSize>
inline void
_initClique(TBamRecord const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<InsertionTag>) {
  svStart = el.mpos + el.malen;
  svEnd = el.pos;
  wiggle = -(svEnd - svStart);
}

// Initialize clique, duplications
template<typename TBamRecord, typename TSize>
inline void
_initClique(TBamRecord const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<DuplicationTag>) {
  svStart = el.mpos;
  svEnd = el.pos + el.alen;
  wiggle = el.maxNormalISize;
}

// Initialize clique, inversions
template<typename TBamRecord, typename TSize>
inline void
_initClique(TBamRecord const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<InversionTag>) {
  int ct=_getSpanOrientation(el, el.libOrient, SVType<InversionTag>());
  if (!ct) {
    svStart = el.mpos + el.malen;
    svEnd = el.pos + el.alen;
  } else {
    svStart = el.mpos;
    svEnd = el.pos;
  }
  wiggle = el.maxNormalISize - std::max(el.alen, el.malen);
}

// Initialize clique, translocations
template<typename TBamRecord, typename TSize>
inline void
_initClique(TBamRecord const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<TranslocationTag>) {
  int ct=_getSpanOrientation(el, el.libOrient, SVType<TranslocationTag>());
  if (ct%2==0) {
    svStart = el.pos + el.alen;
    if (ct>=2) svEnd = el.mpos;
    else svEnd = el.mpos + el.malen;
  } else {
    svStart = el.pos;
    if (ct>=2) svEnd = el.mpos + el.malen;
    else svEnd = el.mpos;
  }
  wiggle=el.maxNormalISize;
}


// Update clique, deletions
template<typename TBamRecord, typename TSize>
inline bool 
_updateClique(TBamRecord const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<DeletionTag>) 
{
  TSize newSvStart = std::max(svStart, el.mpos + el.malen);
  TSize newSvEnd = std::min(svEnd, el.pos);
  TSize newWiggle = el.pos + el.alen - el.mpos - el.maxNormalISize - (newSvEnd - newSvStart);
  TSize wiggleChange = wiggle + (svEnd-svStart) - (newSvEnd - newSvStart);
  if (wiggleChange > newWiggle) newWiggle=wiggleChange;

  // Does the new deletion size agree with all pairs
  if ((newSvStart < newSvEnd) && (newWiggle<=0)) {
    svStart = newSvStart;
    svEnd = newSvEnd;
    wiggle = newWiggle;
    return true;
  }
  return false;
}

// Update clique, insertions
template<typename TBamRecord, typename TSize>
inline bool 
_updateClique(TBamRecord const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<InsertionTag>) 
{
  TSize newSvStart = std::max(svStart, el.mpos + el.malen);
  TSize newSvEnd = std::min(svEnd, el.pos);
  TSize newWiggle = -(newSvEnd - newSvStart);

  // Does the new insertion size agree with all pairs
  if ((newSvStart < newSvEnd) && (newWiggle<=0)) {
    svStart = newSvStart;
    svEnd = newSvEnd;
    wiggle = newWiggle;
    return true;
  }
  return false;
}

// Update clique, duplications
template<typename TBamRecord, typename TSize>
inline bool 
_updateClique(TBamRecord const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<DuplicationTag>) 
{
  TSize newSvStart = std::min(svStart, el.mpos);
  TSize newSvEnd = std::max(svEnd, el.pos + el.alen);
  TSize newWiggle = el.pos - (el.mpos + el.malen) + el.maxNormalISize - (newSvEnd - newSvStart);
  TSize wiggleChange = wiggle - ((newSvEnd - newSvStart) - (svEnd-svStart));
  if (wiggleChange < newWiggle) newWiggle = wiggleChange;

  // Does the new duplication size agree with all pairs
  if ((newSvStart < newSvEnd) && (newWiggle>=0)) {
    svStart = newSvStart;
    svEnd = newSvEnd;
    wiggle = newWiggle;
    return true;
  }
  return false;
}

// Update clique, inversions
template<typename TBamRecord, typename TSize>
inline bool 
_updateClique(TBamRecord const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<InversionTag>) 
{
  int ct=_getSpanOrientation(el, el.libOrient, SVType<InversionTag>());
  TSize newSvStart;
  TSize newSvEnd;
  TSize newWiggle;
  TSize wiggleChange;
  if (!ct) {
    newSvStart = std::max(svStart, el.mpos + el.malen);
    newSvEnd = std::max(svEnd, el.pos + el.alen);
    newWiggle = std::min(el.maxNormalISize - (newSvStart - el.mpos), el.maxNormalISize - (newSvEnd - el.pos));
    wiggleChange = wiggle - std::max(newSvStart - svStart, newSvEnd - svEnd);
  } else {
    newSvStart = std::min(svStart, el.mpos);
    newSvEnd = std::min(svEnd, el.pos);
    newWiggle = el.pos  + el.alen - (el.mpos + el.malen) + el.maxNormalISize - (newSvEnd - newSvStart);
    newWiggle = std::min(el.maxNormalISize - (el.mpos + el.malen - newSvStart), el.maxNormalISize - (el.pos + el.alen - newSvEnd));
    wiggleChange = wiggle - std::max(svStart - newSvStart, svEnd - newSvEnd);
  }
  if (wiggleChange < newWiggle) newWiggle=wiggleChange;

  // Does the new inversion size agree with all pairs
  if ((newSvStart < newSvEnd) && (newWiggle>=0)) {
    svStart = newSvStart;
    svEnd = newSvEnd;
    wiggle = newWiggle;
    return true;
  }
  return false;
}


// Update clique, translocations
template<typename TBamRecord, typename TSize>
inline bool 
_updateClique(TBamRecord const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, SVType<TranslocationTag>) 
{
  int ct = _getSpanOrientation(el, el.libOrient, SVType<TranslocationTag>());
  TSize newSvStart;
  TSize newSvEnd;
  TSize newWiggle = wiggle;
  if (ct%2==0) {
    newSvStart = std::max(svStart, el.pos + el.alen);
    newWiggle -= (newSvStart - svStart);
    if (ct>=2) {
      newSvEnd = std::min(svEnd, el.mpos);
      newWiggle -= (svEnd - newSvEnd);
    } else  {
      newSvEnd = std::max(svEnd, el.mpos + el.malen);
      newWiggle -= (newSvEnd - svEnd);
    }
  } else {
    newSvStart = std::min(svStart, el.pos);
    newWiggle -= (svStart - newSvStart);
    if (ct>=2) {
      newSvEnd = std::max(svEnd, el.mpos + el.malen);
      newWiggle -= (newSvEnd - svEnd);
    } else {
      newSvEnd = std::min(svEnd, el.mpos);
      newWiggle -= (svEnd - newSvEnd);
    }
  }
  // Is this still a valid translocation cluster?
  if (newWiggle>0) {
    svStart = newSvStart;
    svEnd = newSvEnd;
    wiggle = newWiggle;
    return true;
  }
  return false;
}

template<typename TSize>
inline bool
_svSizeCheck(TSize const s, TSize const e, SVType<DeletionTag>) {
  return (( e - s ) >= 300);
}

template<typename TSize>
inline bool
_svSizeCheck(TSize const s, TSize const e, SVType<DuplicationTag>) {
  return (( e - s ) >= 100);
}

template<typename TSize>
inline bool
_svSizeCheck(TSize const s, TSize const e, SVType<InversionTag>) {
  return (( e - s ) >= 100);
}

template<typename TSize>
inline bool
_svSizeCheck(TSize const s, TSize const e, SVType<InsertionTag>) {
  return (( e - s ) >= 0);
}

template<typename TSize>
inline bool
_svSizeCheck(TSize const, TSize const, SVType<TranslocationTag>) {
  return true;
}


template<typename TConfig, typename TSVs, typename TCountMap, typename TTag>
inline void
_annotateCoverage(TConfig const& c, bam_hdr_t* hdr, TSVs& svs, TCountMap& countMap, SVType<TTag>) 
{
  // Find Ns in the reference genome
  typedef boost::icl::interval_set<int> TNIntervals;
  typedef std::vector<TNIntervals> TNGenome;
  TNGenome ni;
  ni.resize( hdr->n_targets );

  // Parse Ns
  faidx_t* fai = fai_load(c.genome.string().c_str());
  for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
    std::string tname(hdr->target_name[refIndex]);
    int32_t seqlen = -1;
    char* seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
    bool nrun = false;
    int nstart = seqlen;
    for(int i=0; i<seqlen; ++i) {
      if ((seq[i] != 'n') && (seq[i] != 'N')) {
	if (nrun) {
	  ni[refIndex].add(boost::icl::discrete_interval<int>::right_open(nstart,i));
	  nrun = false;
	}
      } else {
	if (!nrun) {
	  nrun = true;
	  nstart = i;
	}
      }
    }
    if (nrun) ni[refIndex].add(boost::icl::discrete_interval<int>::right_open(nstart,seqlen));
    if (seqlen) free(seq);
  }
  fai_destroy(fai);
  
  // Add control regions
  typedef std::vector<CovRecord> TCovRecord;
  TCovRecord svc;
  uint32_t lastId = svs.size();
  for (typename TSVs::const_iterator itSV = svs.begin(); itSV != svs.end(); ++itSV) {
    int halfSize = (itSV->svEnd - itSV->svStart)/2;

    // Left control region
    CovRecord sLeft;
    sLeft.chr = itSV->chr;
    sLeft.id = lastId + itSV->id;
    sLeft.svStart = std::max(itSV->svStart - halfSize, 0);
    sLeft.svEnd = itSV->svStart;
    sLeft.peSupport = 0;
    typename TNIntervals::const_iterator itO = ni[sLeft.chr].find(boost::icl::discrete_interval<int>::right_open(sLeft.svStart, sLeft.svEnd));
    while (itO != ni[sLeft.chr].end()) {
      sLeft.svStart = std::max(itO->lower() - halfSize, 0);
      sLeft.svEnd = itO->lower();
      itO = ni[sLeft.chr].find(boost::icl::discrete_interval<int>::right_open(sLeft.svStart, sLeft.svEnd));
    }
    svc.push_back(sLeft);

    // Actual SV
    CovRecord sMiddle;
    sMiddle.chr = itSV->chr;
    sMiddle.id = itSV->id;
    sMiddle.svStart = itSV->svStart;
    sMiddle.svEnd = itSV->svEnd;
    sMiddle.peSupport = itSV->peSupport;
    svc.push_back(sMiddle);

    // Right control region
    CovRecord sRight;
    sRight.chr = itSV->chr;
    sRight.id = 2 * lastId + itSV->id;
    sRight.svStart = itSV->svEnd;
    sRight.svEnd = itSV->svEnd + halfSize;
    sRight.peSupport = 0;
    itO = ni[sRight.chr].find(boost::icl::discrete_interval<int>::right_open(sRight.svStart, sRight.svEnd));
    while (itO != ni[sRight.chr].end()) {
      sRight.svStart = itO->upper();
      sRight.svEnd = itO->upper() + halfSize;
      itO = ni[sRight.chr].find(boost::icl::discrete_interval<int>::right_open(sRight.svStart, sRight.svEnd));
    }
    svc.push_back(sRight);
    //std::cerr << itSV->id << ':' << sLeft.svStart << '-' << sLeft.svEnd << ',' << itSV->svStart << '-' << itSV->svEnd << ',' << sRight.svStart << '-' << sRight.svEnd << std::endl;
  }
  
  typedef std::pair<int, int> TBpRead;
  typedef std::vector<TBpRead> TSVReadCount;
  typedef std::vector<TSVReadCount> TSampleSVReadCount;
  TSampleSVReadCount readCountMap;  
  annotateCoverage(c.files, c.minGenoQual, svc, readCountMap, BpLevelType<NoBpLevelCount>());
  countMap.resize(c.files.size());
  for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
    countMap[file_c].resize(svs.size());
    for (uint32_t id = 0; id < svs.size(); ++id) {
      countMap[file_c][id].rc = readCountMap[file_c][id].second;
      countMap[file_c][id].leftRC = readCountMap[file_c][id + lastId].second;
      countMap[file_c][id].rightRC = readCountMap[file_c][id + 2*lastId].second;
    }
  }
}

template<typename TConfig, typename TSVs, typename TCountMap>
inline void
_annotateCoverage(TConfig const& c, bam_hdr_t*, TSVs& svs, TCountMap& countMap, SVType<TranslocationTag>) 
{
  countMap.resize(c.files.size());
  for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) countMap[file_c].resize(svs.size());
}

template<typename TConfig, typename TSVs, typename TCountMap>
inline void
_annotateCoverage(TConfig const& c, bam_hdr_t*, TSVs& svs, TCountMap& countMap, SVType<InsertionTag>) 
{
  countMap.resize(c.files.size());
  for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) countMap[file_c].resize(svs.size());
}

template<typename TCompEdgeList, typename TBamRecord, typename TSVs, typename TSVType>
inline void
_searchCliques(bam_hdr_t* hdr, TCompEdgeList& compEdge, TBamRecord const& bamRecord, TSVs& svs, int const overallMaxISize, TSVType svType) {
  typedef typename TCompEdgeList::mapped_type TEdgeList;
  typedef typename TEdgeList::value_type TEdgeRecord;

  // Iterate all components
  for(typename TCompEdgeList::iterator compIt = compEdge.begin(); compIt != compEdge.end(); ++compIt) {
    // Sort edges by weight
    std::sort(compIt->second.begin(), compIt->second.end(), SortEdgeRecords<TEdgeRecord>());

    // Find a large clique
    typename TEdgeList::const_iterator itWEdge = compIt->second.begin();
    typename TEdgeList::const_iterator itWEdgeEnd = compIt->second.end();
    typedef std::set<std::size_t> TCliqueMembers;

    TCliqueMembers clique;
    TCliqueMembers incompatible;
    int svStart, svEnd, wiggle;
    int32_t clusterRefID=bamRecord[itWEdge->source].tid;
    int32_t clusterMateRefID=bamRecord[itWEdge->source].mtid;
    _initClique(bamRecord[itWEdge->source], svStart, svEnd, wiggle, svType);
    uint8_t connectionType = _getSpanOrientation(bamRecord[itWEdge->source], bamRecord[itWEdge->source].libOrient, svType);
    if ((clusterRefID==clusterMateRefID) && (svStart >= svEnd))  continue;
    clique.insert(itWEdge->source);
    
    // Grow the clique from the seeding edge
    bool cliqueGrow=true;
    while (cliqueGrow) {
      itWEdge = compIt->second.begin();
      cliqueGrow = false;
      for(;(!cliqueGrow) && (itWEdge != itWEdgeEnd);++itWEdge) {
	std::size_t v;
	if ((clique.find(itWEdge->source) == clique.end()) && (clique.find(itWEdge->target) != clique.end())) v = itWEdge->source;
	else if ((clique.find(itWEdge->source) != clique.end()) && (clique.find(itWEdge->target) == clique.end())) v = itWEdge->target;
	else continue;
	if (incompatible.find(v) != incompatible.end()) continue;
	cliqueGrow = _updateClique(bamRecord[v], svStart, svEnd, wiggle, svType);
	if (cliqueGrow) clique.insert(v);
	else incompatible.insert(v);
      }
    }

    if ((clique.size()>1) && (_svSizeCheck(svStart, svEnd, svType))) {
      StructuralVariantRecord svRec;
      svRec.chr = clusterRefID;
      svRec.chr2 = clusterMateRefID;
      svRec.svStartBeg = std::max((int) svStart - overallMaxISize, 0);
      svRec.svStart = std::min((uint32_t) svStart + 1, hdr->target_len[clusterRefID]);
      svRec.svStartEnd = std::min((uint32_t) svStart + overallMaxISize, hdr->target_len[clusterRefID]);
      svRec.svEndBeg = std::max((int) svEnd - overallMaxISize, 0);
      svRec.svEnd = std::min((uint32_t) svEnd+1, hdr->target_len[clusterMateRefID]);
      svRec.svEndEnd = std::min((uint32_t) svEnd + overallMaxISize, hdr->target_len[clusterMateRefID]);
      svRec.peSupport = clique.size();
      svRec.wiggle = abs(wiggle);
      std::vector<uint8_t> mapQV;
      for(typename TCliqueMembers::const_iterator itC = clique.begin(); itC!=clique.end(); ++itC) mapQV.push_back(bamRecord[*itC].MapQuality);
      std::sort(mapQV.begin(), mapQV.end());
      svRec.peMapQuality = mapQV[mapQV.size()/2];
      if ((clusterRefID==clusterMateRefID) && (svRec.svStartEnd > svRec.svEndBeg)) {
	unsigned int midPointDel = ((svRec.svEnd - svRec.svStart) / 2) + svRec.svStart;
	svRec.svStartEnd = midPointDel -1;
	svRec.svEndBeg = midPointDel;
      }
      svRec.srSupport=0;
      svRec.srAlignQuality=0;
      svRec.precise=false;
      svRec.ct=connectionType;
      svRec.insLen = 0;
      svRec.homLen = 0;
      svRec.id = svs.size();
      svs.push_back(svRec);
    }
  }
}

template<typename TIterator, typename TSVs, typename TSVType>
inline void
_processSRCluster(bam_hdr_t* hdr, TIterator itInit, TIterator itEnd, int32_t refIndex, int32_t bpWindowLen, TSVs& svs, TSVType svType) 
{
  typedef typename TSVs::value_type TStructuralVariant;

  int32_t bestDistance = bpWindowLen;
  TIterator bestSplit = itEnd;
  for (TIterator itBeg = itInit; itBeg != itEnd; ++itBeg) {
    TIterator itNext = itBeg;
    ++itNext;
    for(; itNext != itEnd; ++itNext) {
      int32_t distance = std::abs(itNext->splitbeg - itBeg->splitbeg) + std::abs(itNext->splitend - itBeg->splitend);
      if (distance < bestDistance) {
	bestDistance = distance;
	bestSplit = itBeg;
      }
    }
  }
  if (bestDistance < bpWindowLen) {
    int32_t svStart = bestSplit->splitbeg;
    int32_t svEnd = bestSplit->splitend;
    int32_t svStartBeg = bestSplit->alignbeg;
    int32_t svEndEnd = bestSplit->alignend;
    std::vector<uint8_t> mapQV;
    for (TIterator itBeg = itInit; itBeg != itEnd; ++itBeg) {
      if ( (std::abs(itBeg->splitbeg - svStart) + std::abs(itBeg->splitend - svEnd)) < bpWindowLen ) {
	mapQV.push_back(itBeg->MapQuality);
	if (itBeg->alignbeg < svStartBeg) svStartBeg = itBeg->alignbeg;
	if (itBeg->alignend > svEndEnd) svEndEnd = itBeg->alignend;
      }
    }

    // Augment existing SV call or create a new record
    int32_t searchWindow = 50;
    bool svExists = false;
    typename TSVs::iterator itSV = std::lower_bound(svs.begin(), svs.end(), TStructuralVariant(refIndex, std::max(0, svStart - searchWindow), svEnd), SortSVs<TStructuralVariant>());
    for(; ((itSV != svs.end()) && (std::abs(itSV->svStart - svStart) < searchWindow)); ++itSV) {
      if ((!itSV->precise) && ((std::abs(itSV->svStart - svStart) + std::abs(itSV->svEnd - svEnd)) < searchWindow) && (itSV->chr == refIndex) && (itSV->chr2 == refIndex)) {
	if ((itSV->svEnd < svStart) || (svEnd < itSV->svStart)) continue;
	// Augment existing SV call
	itSV->svStartBeg = std::max(svStartBeg, 0);
	itSV->svStart = svStart;
	itSV->svStartEnd = std::min((uint32_t) svStart + bpWindowLen, hdr->target_len[itSV->chr]);
	itSV->svEndBeg = std::max((int32_t) svEnd - bpWindowLen, 0);
	itSV->svEnd = svEnd;
	itSV->svEndEnd = std::min((uint32_t) svEndEnd, hdr->target_len[itSV->chr2]);
	if ((itSV->chr == itSV->chr2) && (itSV->svStartEnd > itSV->svEndBeg)) {
	  unsigned int midPointDel = ((itSV->svEnd - itSV->svStart) / 2) + itSV->svStart;
	  itSV->svStartEnd = midPointDel -1;
	  itSV->svEndBeg = midPointDel;
	}

	// Found match
	svExists = true;
	break;
      }
    }
    if (!svExists) {
      // Create SV record
      StructuralVariantRecord svRec;
      svRec.chr = refIndex;
      svRec.chr2 = refIndex;
      svRec.svStartBeg = std::max(svStartBeg, 0);
      svRec.svStart = svStart;
      svRec.svStartEnd = std::min((uint32_t) svStart + bpWindowLen, hdr->target_len[svRec.chr]);
      svRec.svEndBeg = std::max((int32_t) svEnd - bpWindowLen, 0);
      svRec.svEnd = svEnd;
      svRec.svEndEnd = std::min((uint32_t) svEndEnd, hdr->target_len[svRec.chr2]);
      svRec.peSupport = 0;
      svRec.wiggle = bpWindowLen;
      std::sort(mapQV.begin(), mapQV.end());
      svRec.peMapQuality = mapQV[mapQV.size()/2];
      svRec.srSupport=mapQV.size();
      svRec.srAlignQuality=0;
      svRec.precise=true;
      svRec.ct=_getCT(svType);
      svRec.insLen = 0;
      svRec.homLen = 0;
      svRec.id = svs.size();
      if ((svRec.svStartBeg < svRec.svStart) && (svRec.svEnd < svRec.svEndEnd)) svs.push_back(svRec);
    }
  }
}




template<typename TSVType>
inline int dellyRun(Config const& c, TSVType svType) {
#ifdef PROFILE
  ProfilerStart("delly.prof");
#endif

  // Collect all promising structural variants
  typedef std::vector<StructuralVariantRecord> TVariants;
  TVariants svs;

  // Create library objects
  typedef boost::unordered_map<std::string, LibraryInfo> TLibraryMap;
  typedef std::vector<TLibraryMap> TSampleLibrary;
  TSampleLibrary sampleLib(c.files.size());
  int overallMaxISize = 0;
  getLibraryParams(c.files, sampleLib, c.percentAbnormal, c.madCutoff);
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c)
    for(TLibraryMap::const_iterator libIt=sampleLib[file_c].begin();libIt!=sampleLib[file_c].end();++libIt)
      if (libIt->second.maxNormalISize > overallMaxISize) overallMaxISize = libIt->second.maxNormalISize;
  
  // Open file handles
  typedef std::vector<samFile*> TSamFile;
  typedef std::vector<hts_idx_t*> TIndex;
  TSamFile samfile;
  TIndex idx;
  samfile.resize(c.files.size());
  idx.resize(c.files.size());
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
    idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
  }
  bam_hdr_t* hdr = sam_hdr_read(samfile[0]);

  // Exclude intervals
  typedef boost::icl::interval_set<uint32_t> TChrIntervals;
  typedef typename TChrIntervals::interval_type TIVal;
  typedef std::vector<TChrIntervals> TRegionsGenome;
  TRegionsGenome validRegions;
  validRegions.resize(hdr->n_targets);
  TRegionsGenome exclg;
  exclg.resize(hdr->n_targets);
  std::vector<bool> validChr;
  validChr.resize(hdr->n_targets, true);
  if (c.hasExcludeFile) {
    std::ifstream chrFile(c.exclude.string().c_str(), std::ifstream::in);
    if (chrFile.is_open()) {
      while (chrFile.good()) {
	std::string chrFromFile;
	getline(chrFile, chrFromFile);
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t,;");
	Tokenizer tokens(chrFromFile, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter!=tokens.end()) {
	  std::string chrName = *tokIter++;
	  int32_t tid = bam_name2id(hdr, chrName.c_str());
	  if (tid >= 0) {
	    if (tokIter!=tokens.end()) {
	      int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
	      int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
	      exclg[tid].insert(TIVal::right_open(start, end));
	    } else validChr[tid] = false; // Exclude entire chromosome
	  }
	}
      }
      chrFile.close();
    }
  }
  // Create the valid regions
  for (int i = 0; i<hdr->n_targets; ++i) {
    uint32_t istart = 0;
    for(typename TChrIntervals::iterator it = exclg[i].begin(); it != exclg[i].end(); ++it) {
      if (istart + 1 < it->lower()) validRegions[i].insert(TIVal::right_open(istart, it->lower() - 1));
      istart = it->upper();
    }
    if (istart + 1 < hdr->target_len[i]) validRegions[i].insert(TIVal::right_open(istart, hdr->target_len[i]));
  }
  exclg.clear();

  // Qualities
  typedef boost::unordered_map<std::size_t, uint8_t> TQualities;
  std::vector<TQualities> qualities;
  qualities.resize(c.files.size());
  typedef boost::unordered_map<std::size_t, int32_t> TAlignmentLength;
  std::vector<TAlignmentLength> alen;
  alen.resize(c.files.size());

  // Parse genome, process chromosome by chromosome
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Paired-end clustering" << std::endl;
  boost::progress_display show_progress( hdr->n_targets );

  faidx_t* fai = fai_load(c.genome.string().c_str());
  for(int32_t refIndex=0; ((refIndex < (int32_t) hdr->n_targets) && (!c.hasVcfFile)); ++refIndex) {
    ++show_progress;
    if (!validChr[refIndex]) continue;
    std::string tname(hdr->target_name[refIndex]);
    int32_t seqlen = -1;
    char* seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
      
    // Create bam alignment record vector
    typedef std::vector<BamAlignRecord> TBamRecord;
    TBamRecord bamRecord;

    // Create split alignment record vector
    typedef std::vector<SplitAlignRecord> TSplitRecord;
    TSplitRecord splitRecord;
    
    // Iterate all samples
#pragma omp parallel for default(shared)
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      // Read alignments
      for(typename TChrIntervals::const_iterator vRIt = validRegions[refIndex].begin(); vRIt != validRegions[refIndex].end(); ++vRIt) {
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, vRIt->lower(), vRIt->upper());
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	  if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	  
	  // Small indel detection using soft clips
	  if (c.indels) {
	    int clipSize = 0;
	    int splitPoint = 0;
	    bool leadingSoftClip = false;
	    if (_validSoftClip(rec, clipSize, splitPoint, leadingSoftClip, c.minMapQual)) {
	      // Iterate both possible breakpoints
	      for (int bpPoint = 0; bpPoint < 2; ++bpPoint) {
		// Leading or trailing softclip?
		if (_validSCOrientation(bpPoint, leadingSoftClip, _getCT(svType), svType)) {
		  // Get the sequence
		  std::string sequence;
		  sequence.resize(rec->core.l_qseq);
		  uint8_t* seqptr = bam_get_seq(rec);
		  for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
		  
		  // Check sequence
		  size_t nCount = std::count(sequence.begin(), sequence.end(), 'N');
		  if ((nCount * 100) / sequence.size() == 0) {
		    std::string cstr = compressStr(sequence);
		    double seqComplexity = (double) cstr.size() / (double) sequence.size();
		    if (seqComplexity >= 0.45) {
		      
		      // Adjust orientation if necessary
		      _adjustOrientation(sequence, bpPoint, _getCT(svType), svType);
		      
		      // Align to local reference
		      int32_t localrefStart = 0;
		      int32_t localrefEnd = 0;
		      int32_t seqLeftOver = 0;
		      if (bpPoint) {
			seqLeftOver = sequence.size() - clipSize;
			localrefStart = std::max(0, (int) rec->core.pos - (int) (c.indelsize + clipSize));
			localrefEnd = std::min(rec->core.pos + seqLeftOver + 25, (int) hdr->target_len[refIndex]);
		      } else {
			seqLeftOver = sequence.size() - (splitPoint - rec->core.pos);
			localrefStart = rec->core.pos;
			localrefEnd = std::min(splitPoint + c.indelsize + seqLeftOver, hdr->target_len[refIndex]);
		      }
		      std::string localref = boost::to_upper_copy(std::string(seq + localrefStart, seq + localrefEnd));
		      typedef boost::multi_array<char, 2> TAlign;
		      TAlign align;
		      AlignConfig<true, false> semiglobal;
		      DnaScore<int> sc(5, -4, -1 * c.aliscore.match * 15, 0);
		      int altScore = gotoh(sequence, localref, align, semiglobal, sc);
		      altScore += c.aliscore.match * 15;
		      
		      // Candidate small indel?
		      AlignDescriptor ad;
		      if (_findSplit(c, sequence, localref, align, ad, svType)) {
			int scoreThresholdAlt = (int) (c.flankQuality * (sequence.size() - (ad.cEnd - ad.cStart - 1)) * sc.match + (1.0 - c.flankQuality) * (sequence.size() - (ad.cEnd - ad.cStart - 1)) * sc.mismatch);
			if (altScore > scoreThresholdAlt) {
			  
			  // Debug consensus to reference alignment
			  //for(TAIndex i = 0; i<align.shape()[0]; ++i) {
			  //for(TAIndex j = 0; j<align.shape()[1]; ++j) std::cerr << align[i][j];
			  //std::cerr << std::endl;
			  //}
			  //std::cerr << bpPoint << ',' << cStart << ',' << cEnd << ',' << rStart << ',' << rEnd << std::endl;
			  
#pragma omp critical
			  {
			    splitRecord.push_back(SplitAlignRecord(localrefStart + ad.rStart - ad.cStart, localrefStart + ad.rStart, localrefStart + ad.rEnd, localrefStart + ad.rEnd + seqLeftOver + 25, rec->core.qual));
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	  
	  // Paired-end clustering
	  if (rec->core.flag & BAM_FPAIRED) {
	    // Mate unmapped
	    if ((rec->core.mtid<0) || (rec->core.flag & BAM_FMUNMAP)) continue;
	    
	    // Mapping positions valid?
	    if (_mappingPos(rec->core.tid, rec->core.mtid, rec->core.pos, rec->core.mpos, svType)) continue;
	      
	    // Is this a discordantly mapped paired-end?
	    std::string rG = "DefaultLib";
	    uint8_t *rgptr = bam_aux_get(rec, "RG");
	    if (rgptr) {
	      char* rg = (char*) (rgptr + 1);
	      rG = std::string(rg);
	    }
	    TLibraryMap::iterator libIt = sampleLib[file_c].find(rG);
	    if (libIt == sampleLib[file_c].end()) std::cerr << "Missing read group: " << rG << std::endl;
	    if (_acceptedInsertSize(libIt->second, abs(rec->core.isize), svType)) continue; 
	    if (_acceptedOrientation(libIt->second.defaultOrient, getStrandIndependentOrientation(rec->core), svType)) continue;
	    
	    // Get or store the mapping quality for the partner
	    if (_firstPairObs(rec->core.tid, rec->core.mtid, rec->core.pos, rec->core.mpos, svType)) {
	      uint8_t r2Qual = rec->core.qual;
	      uint8_t* ptr = bam_aux_get(rec, "AS");
	      if (ptr) {
		int score = std::abs((int) bam_aux2i(ptr));
		r2Qual = std::min(r2Qual, (uint8_t) ( (score<255) ? score : 255 ));
	      }
	      std::size_t hv = hash_pair(rec);
	      qualities[file_c][hv]= r2Qual;
	      alen[file_c][hv]= alignmentLength(rec);
	    } else {
	      // Get the two mapping qualities
	      uint8_t r2Qual = rec->core.qual;
	      uint8_t* ptr = bam_aux_get(rec, "AS");
	      if (ptr) {
		int score = std::abs((int) bam_aux2i(ptr));
		r2Qual = std::min(r2Qual, (uint8_t) ( (score<255) ? score : 255 ));
	      }
	      std::size_t hv=hash_pair_mate(rec);
	      uint8_t pairQuality = std::min(qualities[file_c][hv], r2Qual);
	      qualities[file_c][hv]= (uint8_t) 0;
	      
	      // Pair quality
	      if (pairQuality < c.minMapQual) continue;
	      
#pragma omp critical
	      {
		bamRecord.push_back(BamAlignRecord(rec, pairQuality, alignmentLength(rec), alen[file_c][hv], libIt->second.median, libIt->second.mad, libIt->second.maxNormalISize, libIt->second.defaultOrient));
	      }
	      ++libIt->second.abnormal_pairs;
	    }
	  }
	}
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }
      // Clean-up qualities
      _resetQualities(qualities[file_c], alen[file_c], svType);
    }
    
    // Sort BAM records according to position
    std::sort(bamRecord.begin(), bamRecord.end(), SortBamRecords<BamAlignRecord>());
    
    // Components
    typedef std::vector<uint32_t> TComponent;
    TComponent comp;
    comp.resize(bamRecord.size(), 0);
    uint32_t numComp = 0;
    
    // Edge lists for each component
    typedef uint8_t TWeightType;
    typedef EdgeRecord<TWeightType, std::size_t> TEdgeRecord;
    typedef std::vector<TEdgeRecord> TEdgeList;
    typedef std::map<uint32_t, TEdgeList> TCompEdgeList;
    TCompEdgeList compEdge;
    
    // Iterate the chromosome range
    std::size_t lastConnectedNode = 0;
    std::size_t lastConnectedNodeStart = 0;
    std::size_t bamItIndex = 0;
    for(TBamRecord::const_iterator bamIt = bamRecord.begin(); bamIt != bamRecord.end(); ++bamIt, ++bamItIndex) {
      // Safe to clean the graph?
      if (bamItIndex > lastConnectedNode) {
	// Clean edge lists
	if (!compEdge.empty()) {
	  _searchCliques(hdr, compEdge, bamRecord, svs, overallMaxISize, svType);
	  lastConnectedNodeStart = lastConnectedNode;
	  compEdge.clear();
	}
      }
      int32_t const minCoord = _minCoord(bamIt->pos, bamIt->mpos, svType);
      int32_t const maxCoord = _maxCoord(bamIt->pos, bamIt->mpos, svType);
      TBamRecord::const_iterator bamItNext = bamIt;
      ++bamItNext;
      std::size_t bamItIndexNext = bamItIndex + 1;
      for(; ((bamItNext != bamRecord.end()) && (abs(_minCoord(bamItNext->pos, bamItNext->mpos, svType) + bamItNext->alen - minCoord) <= overallMaxISize)) ; ++bamItNext, ++bamItIndexNext) {
	// Check that mate chr agree (only for translocations)
	if (bamIt->mtid != bamItNext->mtid) continue;
	
	// Check combinability of pairs
	if (_pairsDisagree(minCoord, maxCoord, bamIt->alen, bamIt->maxNormalISize, _minCoord(bamItNext->pos, bamItNext->mpos, svType), _maxCoord(bamItNext->pos, bamItNext->mpos, svType), bamItNext->alen, bamItNext->maxNormalISize, _getSpanOrientation(*bamIt, bamIt->libOrient, svType), _getSpanOrientation(*bamItNext, bamItNext->libOrient, svType), svType)) continue;
	
	// Update last connected node
	if (bamItIndexNext > lastConnectedNode ) lastConnectedNode = bamItIndexNext;
	
	// Assign components
	uint32_t compIndex = 0;
	if (!comp[bamItIndex]) {
	  if (!comp[bamItIndexNext]) {
	    // Both vertices have no component
	    compIndex = ++numComp;
	    comp[bamItIndex] = compIndex;
	    comp[bamItIndexNext] = compIndex;
	    compEdge.insert(std::make_pair(compIndex, TEdgeList()));
	  } else {
	    compIndex = comp[bamItIndexNext];
	    comp[bamItIndex] = compIndex;
	  }
	} else {
	  if (!comp[bamItIndexNext]) {
	    compIndex = comp[bamItIndex];
	    comp[bamItIndexNext] = compIndex;
	  } else {
	    // Both vertices have a component
	    if (comp[bamItIndexNext] == comp[bamItIndex]) {
	      compIndex = comp[bamItIndexNext];
	    } else {
	      // Merge components
	      compIndex = comp[bamItIndex];
	      uint32_t otherIndex = comp[bamItIndexNext];
	      if (otherIndex < compIndex) {
		compIndex = comp[bamItIndexNext];
		otherIndex = comp[bamItIndex];
	      }
	      // Re-label other index
	      for(std::size_t i = lastConnectedNodeStart; i <= lastConnectedNode; ++i) {
		if (otherIndex == comp[i]) comp[i] = compIndex;
	      }
	      // Merge edge lists
	      TCompEdgeList::iterator compEdgeIt = compEdge.find(compIndex);
	      TCompEdgeList::iterator compEdgeOtherIt = compEdge.find(otherIndex);
	      compEdgeIt->second.insert(compEdgeIt->second.end(), compEdgeOtherIt->second.begin(), compEdgeOtherIt->second.end());
	      compEdge.erase(compEdgeOtherIt);
	    }
	  }
	}
	
	// Append new edge
	TCompEdgeList::iterator compEdgeIt = compEdge.find(compIndex);
	if (compEdgeIt->second.size() < c.graphPruning) {
	  TWeightType weight = (TWeightType) ( std::log((float) abs( abs( (_minCoord(bamItNext->pos, bamItNext->mpos, svType) - minCoord) - (_maxCoord(bamItNext->pos, bamItNext->mpos, svType) - maxCoord) ) - abs(bamIt->Median - bamItNext->Median) ) + 1 ) / std::log(2) );
	  compEdgeIt->second.push_back(TEdgeRecord(bamItIndex, bamItIndexNext, weight));
	}
      }
    }
    if (!compEdge.empty()) {
      _searchCliques(hdr, compEdge, bamRecord, svs, overallMaxISize, svType);
      compEdge.clear();
    }
    
    // Sort SVs for look-up
    sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());
      
    // Add the soft clip SV records
    if (c.indels) {
      int32_t bpWindowLen = 10;
      int32_t maxLookAhead = 0;
      TSplitRecord::const_iterator splitClusterIt = splitRecord.end();
      std::sort(splitRecord.begin(), splitRecord.end(), SortSplitRecords<SplitAlignRecord>());
      for(TSplitRecord::const_iterator splitIt = splitRecord.begin(); splitIt!=splitRecord.end(); ++splitIt) {
	if ((maxLookAhead) && (splitIt->splitbeg > maxLookAhead)) {
	  // Process split read cluster
	  _processSRCluster(hdr, splitClusterIt, splitIt, refIndex, bpWindowLen, svs, svType);
	  maxLookAhead = 0;
	  splitClusterIt = splitRecord.end();
	}
	if ((!maxLookAhead) || (splitIt->splitbeg < maxLookAhead)) {
	  if (!maxLookAhead) splitClusterIt = splitIt;
	  maxLookAhead = splitIt->splitbeg + bpWindowLen;
	}
      }
      TSplitRecord::const_iterator splitIt = splitRecord.end();
      _processSRCluster(hdr, splitClusterIt, splitIt, refIndex, bpWindowLen, svs, svType);
    }
    if (seqlen) free(seq);
  }
  fai_destroy(fai);
  
  // Split-read search
  if (!c.hasVcfFile) {
    if (!svs.empty()) {
      findPutativeSplitReads(c, svs, svType);
      
      if (c.indels) {
	// Sort SVs for look-up and by decreasing PE support
	sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());
	
	// Temporary SV container
	TVariants svc;
	
	// Clean-up SV set
	for(typename TVariants::iterator svIt = svs.begin(); svIt != svs.end(); ++svIt) {
	  // Unresolved soft clips
	  if ((svIt->precise) && (svIt->srAlignQuality == 0)) continue;
	  
	  // Precise duplicates
	  int32_t searchWindow = 10;
	  bool svExists = false;
	  typename TVariants::iterator itOther = std::lower_bound(svc.begin(), svc.end(), StructuralVariantRecord(svIt->chr, std::max(0, svIt->svStart - searchWindow), svIt->svEnd), SortSVs<StructuralVariantRecord>());
	  for(; ((itOther != svc.end()) && (std::abs(itOther->svStart - svIt->svStart) < searchWindow)); ++itOther) {
	    if (!svIt->precise) continue;
	    if ((svIt->chr != itOther->chr) || (svIt->chr2 != itOther->chr2)) continue;
	    if ((std::abs(svIt->svStart - itOther->svStart) + std::abs(svIt->svEnd - itOther->svEnd)) > searchWindow) continue;
	    if ((svIt->svEnd < itOther->svStart) || (itOther->svEnd < svIt->svStart)) continue;
	    svExists=true;
	    break;
	  }
	  if (svExists) continue;
	    
	    // Add SV
	  svc.push_back(*svIt);
	}
	
	// Final set of precise and imprecise SVs
	svs = svc;
	// Re-number SVs
	uint32_t cliqueCount = 0;
	for(typename TVariants::iterator svIt = svs.begin(); svIt != svs.end(); ++svIt) svIt->id = cliqueCount++;
      }
    }
  } else {
    // Read SV records from input file
    if (c.format == "json.gz") jsonParse(c, hdr, overallMaxISize, svs, svType);
    else vcfParse(c, hdr, overallMaxISize, svs, svType);
  }


  // Debug output
  //for (TVariants::const_iterator s = svs.begin();s!=svs.end();++s) std::cerr << s->svStart << ',' << s->svEnd << ',' <<  s->svStartBeg << ',' << s->svStartEnd << ',' << s->svEndBeg << ',' << s->svEndEnd << ',' << s->peSupport << ',' << s->srSupport << ',' << s->wiggle << ',' << s->srAlignQuality << ',' << s->precise << ',' << (int32_t) s->peMapQuality << ',' << s->chr << ',' << s->chr2 << ",Consensus:" << s->consensus << ',' << s->id << ',' << s->insLen << ',' << s->homLen << ',' << _addOrientation(s->ct) << std::endl;

  // Any SVs for genotyping
  if (svs.empty()) {
    std::cout << "No structural variants found!" << std::endl;
    std::cout << "Done." << std::endl;
    return 0;
  }

  // Annotate junction reads
  typedef std::pair<std::vector<uint8_t>, std::vector<uint8_t> > TReadQual;
  typedef std::vector<TReadQual> TSVJunctionMap;
  typedef std::vector<TSVJunctionMap> TSampleSVJunctionMap;
  TSampleSVJunctionMap junctionCountMap;
  annotateJunctionReads(c, svs, junctionCountMap, svType);

  // Annotate spanning coverage
  typedef std::vector<TReadQual> TSVSpanningMap;
  typedef std::vector<TSVSpanningMap> TSampleSVSpanningMap;
  TSampleSVSpanningMap spanCountMap;
  annotateSpanningCoverage(c.files, c.minGenoQual, sampleLib, svs, spanCountMap, svType);

  // Annotate coverage
  typedef std::vector<ReadCount> TSVReadCount;
  typedef std::vector<TSVReadCount> TSampleSVReadCount;
  TSampleSVReadCount rcMap;
  _annotateCoverage(c, hdr, svs, rcMap, svType);

  // VCF output
  if (svs.size()) {
    if (c.format == "json.gz") jsonOutput(c, svs, junctionCountMap, rcMap, spanCountMap, svType);
    else vcfOutput(c, svs, junctionCountMap, rcMap, spanCountMap, svType);
  }

  // Clean-up
  bam_hdr_destroy(hdr);
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    hts_idx_destroy(idx[file_c]);
    sam_close(samfile[file_c]);
  }

  // Output library statistics
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Library statistics" << std::endl;
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    std::cout << "Sample: " << c.sampleName[file_c] << std::endl;
    for(TLibraryMap::const_iterator libIt = sampleLib[file_c].begin(); libIt != sampleLib[file_c].end(); ++libIt) {
      std::cout << "RG: ID=" << libIt->first << ",Median=" << libIt->second.median << ",MAD=" << libIt->second.mad << ",LibLayout=" << (int) libIt->second.defaultOrient << ",MaxSizeCut=" << libIt->second.maxISizeCutoff << ",MinSizeCut=" << libIt->second.minISizeCutoff << ",UniqueDiscordantPairs=" << libIt->second.abnormal_pairs << std::endl;
    }
  }
#ifdef PROFILE
  ProfilerStop();
#endif
  
  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
  return 0;
}


int delly(int argc, char **argv) {
  Config c;

  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("type,t", boost::program_options::value<std::string>(&c.svType)->default_value("DEL"), "SV type (DEL, DUP, INV, TRA, INS)")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("sv.bcf"), "SV BCF output file")
    ("exclude,x", boost::program_options::value<boost::filesystem::path>(&c.exclude), "file with regions to exclude")
    ("technology,e", boost::program_options::value<std::string>(&c.technology)->default_value("illumina"), "technology (illumina, pacbio)")
    ;

  boost::program_options::options_description pem("PE options");
  pem.add_options()
    ("map-qual,q", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(1), "min. paired-end mapping quality")
    ("mad-cutoff,s", boost::program_options::value<unsigned short>(&c.madCutoff)->default_value(9), "insert size cutoff, median+s*MAD (deletions only)")
    ;

  boost::program_options::options_description breaks("SR options");
  breaks.add_options()
    ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
    ("noindels,n", "no small InDel calling")
    ("indelsize,i", boost::program_options::value<uint32_t>(&c.indelsize)->default_value(500), "max. small InDel size")
    ;

  boost::program_options::options_description geno("Genotyping options");
  geno.add_options()
    ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input VCF/BCF file for re-genotyping")
    ("geno-qual,u", boost::program_options::value<unsigned short>(&c.minGenoQual)->default_value(5), "min. mapping quality for genotyping")
    ;

  // Define hidden options
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
    ("format,f", boost::program_options::value<std::string>(&c.format)->default_value("bcf"), "output format (bcf, json.gz)")
    ("pe-fraction,c", boost::program_options::value<float>(&c.percentAbnormal)->default_value(0.0), "fixed fraction c of discordant PEs, for c=0 MAD cutoff is used")
    ("pruning,j", boost::program_options::value<uint32_t>(&c.graphPruning)->default_value(1000), "PE graph pruning cutoff")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  // Set the visibility
  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(pem).add(breaks).add(geno).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic).add(pem).add(breaks).add(geno);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);


  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) { 
    std::cout << std::endl;
    std::cout << "Usage: delly " << argv[0] << " [OPTIONS] -g <ref.fa> <sample1.sort.bam> <sample2.sort.bam> ..." << std::endl;
    std::cout << visible_options << "\n";
    return 0;
  }

  // Check reference
  if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
    std::cerr << "Reference file is missing: " << c.genome.string() << std::endl;
    return 1;
  } else {
    faidx_t* fai = fai_load(c.genome.string().c_str());
    if (fai == NULL) {
      if (fai_build(c.genome.string().c_str()) == -1) {
	std::cerr << "Fail to open genome fai index for " << c.genome.string() << std::endl;
	return 1;
      } else fai = fai_load(c.genome.string().c_str());
    }
    fai_destroy(fai);
  }

  // Check input files
  c.sampleName.resize(c.files.size());
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    if (!(boost::filesystem::exists(c.files[file_c]) && boost::filesystem::is_regular_file(c.files[file_c]) && boost::filesystem::file_size(c.files[file_c]))) {
      std::cerr << "Alignment file is missing: " << c.files[file_c].string() << std::endl;
      return 1;
    }
    samFile* samfile = sam_open(c.files[file_c].string().c_str(), "r");
    if (samfile == NULL) {
      std::cerr << "Fail to open file " << c.files[file_c].string() << std::endl;
      return 1;
    }
    hts_idx_t* idx = sam_index_load(samfile, c.files[file_c].string().c_str());
    if (idx == NULL) {
      std::cerr << "Fail to open index for " << c.files[file_c].string() << std::endl;
      return 1;
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    if (hdr == NULL) {
      std::cerr << "Fail to open header for " << c.files[file_c].string() << std::endl;
      return 1;
    }
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
      std::string tname(hdr->target_name[refIndex]);
      if (!faidx_has_seq(fai, tname.c_str())) {
	std::cerr << "BAM file chromosome " << hdr->target_name[refIndex] << " is NOT present in your reference file " << c.genome.string() << std::endl;
	return 1;
      }
    }
    fai_destroy(fai);
    std::string sampleName;
    if (!getSMTag(std::string(hdr->text), c.files[file_c].stem().string(), sampleName)) {
      std::cerr << "Only one sample (@RG:SM) is allowed per input BAM file " << c.files[file_c].string() << std::endl;
      return 1;
    } else c.sampleName[file_c] = sampleName;
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }

  // Check exclude file
  if (vm.count("exclude")) {
    if (!(boost::filesystem::exists(c.exclude) && boost::filesystem::is_regular_file(c.exclude) && boost::filesystem::file_size(c.exclude))) {
      std::cerr << "Exclude file is missing: " << c.exclude.string() << std::endl;
      return 1;
    }
    c.hasExcludeFile = true;
  } else c.hasExcludeFile = false;

  // Check input VCF file
  if (vm.count("vcffile")) {
    if (!(boost::filesystem::exists(c.vcffile) && boost::filesystem::is_regular_file(c.vcffile) && boost::filesystem::file_size(c.vcffile))) {
      std::cerr << "Input VCF/BCF file is missing: " << c.vcffile.string() << std::endl;
      return 1;
    }
    if (c.format != "json.gz") {
      htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
      if (ifile == NULL) {
	std::cerr << "Fail to open file " << c.vcffile.string() << std::endl;
	return 1;
      }
      bcf_hdr_t* hdr = bcf_hdr_read(ifile);
      if (hdr == NULL) {
	std::cerr << "Fail to open index file " << c.vcffile.string() << std::endl;
	return 1;
      }
      bcf_hdr_destroy(hdr);
      bcf_close(ifile);
    }
    c.hasVcfFile = true;
  } else c.hasVcfFile = false;

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  std::cout << "delly ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Always ignore reads of mapping quality <5 for genotyping, otherwise het. is more likely!
  if (c.minGenoQual<5) c.minGenoQual=5;

  // Small InDels?
  c.indels = false;
  if (!vm.count("noindels")) {
    if ((c.svType == "DEL") || (c.svType == "INS")) c.indels = true;
  }

  // Run main program
  if (c.technology == "illumina") {
    c.aliscore = DnaScore<int>(5, -4, -10, -1);
    c.flankQuality = 0.95;
    c.minimumFlankSize = 13;
    if (c.svType == "DEL") return dellyRun(c, SVType<DeletionTag>());
    else if (c.svType == "DUP") return dellyRun(c, SVType<DuplicationTag>());
    else if (c.svType == "INV") return dellyRun(c, SVType<InversionTag>());
    else if (c.svType == "TRA") return dellyRun(c, SVType<TranslocationTag>());
    else if (c.svType == "INS") return dellyRun(c, SVType<InsertionTag>());
    else {
      std::cerr << "SV analysis type not supported by Delly: " << c.svType << std::endl;
      return 1;
    }
  } else if (c.technology == "pacbio") {
    //c.aliscore = DnaScore<int>(2, -5, -2, -1);
    c.aliscore = DnaScore<int>(5, -4, -2, -1);
    c.flankQuality = 0.8;
    c.minimumFlankSize = 13;
    if (c.svType == "DEL") return pacbioRun(c, SVType<DeletionTag>());
    else {
      std::cerr << "PacBio SV analysis not yet supported." << std::endl;
      return 1;
    }
  } else {
    std::cerr << "Technology " << c.technology << " is not supported by Delly." << std::endl;
    return 1;
  }
}

}

#endif
