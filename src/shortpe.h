/*
============================================================================
DELLY: Structural variant discovery by integrated PE mapping and SR analysis
============================================================================
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

#ifndef SHORTPE_H
#define SHORTPE_H


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
#include "coverage.h"
#include "msa.h"
#include "split.h"
#include "junction.h"
#include "cluster.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <stdio.h>

namespace torali
{
  
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
    uint32_t flag;
    uint8_t MapQuality;
  
    BamAlignRecord(bam1_t* rec, uint8_t pairQuality, uint16_t a, uint16_t ma, int32_t median, int32_t mad, int32_t maxISize) : tid(rec->core.tid), pos(rec->core.pos), mtid(rec->core.mtid), mpos(rec->core.mpos), alen(a), malen(ma), Median(median), Mad(mad), maxNormalISize(maxISize), flag(rec->core.flag), MapQuality(pairQuality) {}
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

  template<typename TConfig, typename TValidRegion, typename TStructuralVariantRecord>
  inline bool
  findPutativeSplitReads(TConfig const& c, TValidRegion const& validRegions, std::vector<TStructuralVariantRecord>& svs) 
  {
    typedef std::vector<TStructuralVariantRecord> TSVs;

    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
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
    for(int32_t refIndex2 = 0; refIndex2 < hdr->n_targets; ++refIndex2) {
      ++show_progress;
      if (validRegions[refIndex2].empty()) continue;
      char* sndSeq = NULL;
      for(int32_t refIndex = refIndex2; refIndex < hdr->n_targets; ++refIndex) {
	if (validRegions[refIndex].empty()) continue;
	char* seq = NULL;


	for(typename TSVs::iterator svIt = svs.begin();svIt!=svs.end(); ++svIt) {
	  if ((svIt->chr != refIndex) || (svIt->chr2 != refIndex2)) continue;
	  // Lazy loading of reference sequence
	  if (seq == NULL) {
	    int32_t seqlen = -1;
	    std::string tname(hdr->target_name[refIndex]);
	    seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
	  }
	  if ((sndSeq == NULL) && (refIndex != refIndex2)) {
	    int32_t seqlen = -1;
	    std::string tname(hdr->target_name[refIndex2]);
	    sndSeq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex2], &seqlen);
	  }

	  // Set tag alleles
	  svIt->alleles = _addAlleles(boost::to_upper_copy(std::string(seq + svIt->svStart - 1, seq + svIt->svStart)), std::string(hdr->target_name[svIt->chr2]), *svIt, svIt->svt);
	  
	  typedef std::vector<std::pair<int32_t, std::string> > TClipsizeSplit;
	  typedef std::vector<TClipsizeSplit> TOffsetSplit;
	  typedef std::vector<int32_t> TSplitPoints;
	  TOffsetSplit osp0;
	  TSplitPoints spp0;
	  TOffsetSplit osp1;
	  TSplitPoints spp1;
	  
	  // Find putative split reads in all samples
	  Breakpoint bp(*svIt);
	  _initBreakpoint(hdr, bp, svIt->wiggle, svIt->svt);
	  for (unsigned int bpPoint = 0; bpPoint<2; ++bpPoint) {
	    int32_t regionChr = 0;
	    int32_t regionStart = 0;
	    int32_t regionEnd = 0;
	    if (bpPoint) {
	      regionChr = bp.chr2;
	      regionStart = bp.svEndBeg;
	      regionEnd = bp.svEndEnd;
	      spp1.resize(regionEnd-regionStart, 0);
	      osp1.resize(regionEnd-regionStart, TClipsizeSplit());
	    } else {
	      regionChr = bp.chr;
	      regionStart = bp.svStartBeg;
	      regionEnd = bp.svStartEnd;
	      spp0.resize(regionEnd-regionStart, 0);
	      osp0.resize(regionEnd-regionStart, TClipsizeSplit());
	    }
#pragma omp parallel for default(shared)
	    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	      hts_itr_t* iter = sam_itr_queryi(idx[file_c], regionChr, regionStart, regionEnd);
	      bam1_t* rec = bam_init1();
	      while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
		if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
		
		// Valid soft clip?
		int32_t clipSize = 0;
		int32_t splitPoint = 0;
		bool leadingSoftClip = false;
		if (_validSoftClip(rec, clipSize, splitPoint, leadingSoftClip, c.minMapQual, svIt->svt)) {
		  if ((splitPoint >= regionStart) && (splitPoint < regionEnd)) {
		    splitPoint -= regionStart;
		    // Leading or trailing softclip?
		    if (_validSCOrientation(bpPoint, leadingSoftClip, svIt->svt)) {
		      // Get the sequence
		      std::string sequence;
		      sequence.resize(rec->core.l_qseq);
		      uint8_t* seqptr = bam_get_seq(rec);
		      for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
		      
		      // Reverse complement iff necesssary
		      _adjustOrientation(sequence, bpPoint, svIt->svt);
		      
		      if (bpPoint) {
#pragma omp critical
			{
			  ++spp1[splitPoint];
			  osp1[splitPoint].push_back(std::make_pair(clipSize, sequence));
			} 
		      } else {
#pragma omp critical
			{
			  ++spp0[splitPoint];
			  osp0[splitPoint].push_back(std::make_pair(clipSize, sequence));
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
	  // Sort candidate split reads
	  TClipsizeSplit csSplit;
	  int mvAvg, lBound, uBound;
	  _movingAverage(spp0, 5, mvAvg, lBound, uBound);
	  if (mvAvg > 0)
	    for(int32_t k = lBound; k < uBound; ++k) csSplit.insert(csSplit.end(), osp0[k].begin(), osp0[k].end());
	  _movingAverage(spp1, 5, mvAvg, lBound, uBound);
	  if (mvAvg > 0)
	    for(int32_t k = lBound; k < uBound; ++k) csSplit.insert(csSplit.end(), osp1[k].begin(), osp1[k].end());
	  std::sort(csSplit.begin(), csSplit.end());
	    
	  typedef std::set<std::string> TSplitReadSet;
	  TSplitReadSet splitReadSet;
	  for(typename TClipsizeSplit::reverse_iterator rit = csSplit.rbegin(); rit != csSplit.rend(); ++rit)
	    if (splitReadSet.size() < 10) splitReadSet.insert(rit->second);
	    else break;
	  totalSplitReadsAligned += splitReadSet.size();
	  
	  // MSA
	  if (splitReadSet.size() > 1) {
	    svIt->srSupport = msa(c, splitReadSet, svIt->consensus);
	    if (!alignConsensus(c, hdr, seq, sndSeq, *svIt)) {
	      svIt->consensus = "";
	      svIt->srSupport = 0;
	    } else {
	      // Update REF & ALT alleles because of split-read refinement (except for the small InDels)
	      if ((svIt->peSupport != 0) || (!c.indels)) svIt->alleles = _addAlleles(boost::to_upper_copy(std::string(seq + svIt->svStart - 1, seq + svIt->svStart)), std::string(hdr->target_name[svIt->chr2]), *svIt, svIt->svt);
	    }
	  }
	}
	if (seq != NULL) free(seq);
      }
      if (sndSeq != NULL) free(sndSeq);
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
  _initClique(TBamRecord const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, int32_t const svt) {
    if (_translocation(svt)) {
      uint8_t ct = _getSpanOrientation(svt);
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
    } else {
      if (svt == 0) {
	svStart = el.mpos + el.malen;
	svEnd = el.pos + el.alen;
	wiggle = el.maxNormalISize - std::max(el.alen, el.malen);
      } else if (svt == 1) {
	svStart = el.mpos;
	svEnd = el.pos;
	wiggle = el.maxNormalISize - std::max(el.alen, el.malen);
      } else if (svt == 2) {
	svStart = el.mpos + el.malen;
	svEnd = el.pos;
	wiggle =  -el.maxNormalISize;
      } else if (svt == 3) {
	svStart = el.mpos;
	svEnd = el.pos + el.alen;
	wiggle = el.maxNormalISize;
      }
    } 
  }

  // Update clique, deletions
  template<typename TBamRecord, typename TSize>
  inline bool 
  _updateClique(TBamRecord const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, int32_t const svt) 
  {
    if (_translocation(svt)) {
      int ct = _getSpanOrientation(svt);
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
    } else {
      if ((svt == 0) || (svt == 1)) { 
	int ct = _getSpanOrientation(svt);
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
      } else if (svt == 2) {
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
      } else if (svt == 3) {
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
    }
    return false;
  }


  template<typename TCompEdgeList, typename TBamRecord, typename TSVs>
  inline void
  _searchCliques(bam_hdr_t* hdr, TCompEdgeList& compEdge, TBamRecord const& bamRecord, TSVs& svs, int32_t const svt) {
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
      int32_t svStart = -1;
      int32_t svEnd = -1;
      int32_t wiggle = 0;
      int32_t clusterRefID=bamRecord[itWEdge->source].tid;
      int32_t clusterMateRefID=bamRecord[itWEdge->source].mtid;
      _initClique(bamRecord[itWEdge->source], svStart, svEnd, wiggle, svt);
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
	  cliqueGrow = _updateClique(bamRecord[v], svStart, svEnd, wiggle, svt);
	  if (cliqueGrow) clique.insert(v);
	  else incompatible.insert(v);
	}
      }
      
      if ((clique.size()>1) && (_svSizeCheck(svStart, svEnd, svt))) {
	StructuralVariantRecord svRec;
	svRec.chr = clusterRefID;
	svRec.chr2 = clusterMateRefID;
	svRec.svStart = std::min((uint32_t) svStart + 1, hdr->target_len[clusterRefID]);
	svRec.svEnd = std::min((uint32_t) svEnd+1, hdr->target_len[clusterMateRefID]);
	svRec.peSupport = clique.size();
	svRec.wiggle = std::max(abs(wiggle), 50);
	std::vector<uint8_t> mapQV;
	for(typename TCliqueMembers::const_iterator itC = clique.begin(); itC!=clique.end(); ++itC) mapQV.push_back(bamRecord[*itC].MapQuality);
	std::sort(mapQV.begin(), mapQV.end());
	svRec.peMapQuality = mapQV[mapQV.size()/2];
	svRec.srSupport=0;
	svRec.srAlignQuality=0;
	svRec.precise=false;
	svRec.svt = svt;
	svRec.insLen = 0;
	svRec.homLen = 0;
#pragma omp critical
	{
	  svRec.id = svs.size();
	  svs.push_back(svRec);
	}
      }
    }
  }
  
  template<typename TConfig, typename TValidRegion, typename TVariants, typename TSampleLib>
  inline void
  shortPE(TConfig const& c, TValidRegion const& validRegions, TVariants& svs, TSampleLib& sampleLib)
  {
    typedef typename TValidRegion::value_type TChrIntervals;

    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile[0]);

    // Split-read records
    typedef std::vector<SRBamRecord> TSRBamRecord;
    typedef std::vector<TSRBamRecord> TSvtSRBamRecord;
    TSvtSRBamRecord srBR(2 * DELLY_SVT_TRANS, TSRBamRecord());

    // Create bam alignment record vector
    typedef std::vector<BamAlignRecord> TBamRecord;
    typedef std::vector<TBamRecord> TSvtBamRecord;
    TSvtBamRecord bamRecord(2 * DELLY_SVT_TRANS, TBamRecord());
     
    // Parse genome, process chromosome by chromosome
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Paired-end clustering" << std::endl;
    boost::progress_display show_progress( c.files.size() * hdr->n_targets );
    // Iterate all samples
#pragma omp parallel for default(shared)
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      // Maximum insert size
      int32_t overallMaxISize = std::max(sampleLib[file_c].maxISizeCutoff, sampleLib[file_c].rs);

      // Inter-chromosomal mate map and alignment length
      typedef std::pair<uint8_t, int32_t> TQualLen;
      typedef boost::unordered_map<std::size_t, TQualLen> TMateMap;
      std::vector<TMateMap> matetra(c.files.size());

      // Split-read junctions
      typedef std::vector<Junction> TJunctionVector;
      typedef std::map<unsigned, TJunctionVector> TReadBp;
      TReadBp readBp;
      
      // Iterate all chromosomes for that sample
      for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
	++show_progress;

	// Any data?
	if (validRegions[refIndex].empty()) continue;
	bool nodata = true;
	std::string suffix("cram");
	std::string str(c.files[file_c].string());
	if ((str.size() >= suffix.size()) && (str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0)) nodata = false;
	uint64_t mapped = 0;
	uint64_t unmapped = 0;
	hts_idx_get_stat(idx[file_c], refIndex, &mapped, &unmapped);
	if (mapped) nodata = false;
	if (nodata) continue;

	// Intra-chromosomal mate map and alignment length
	TMateMap mateMap;

	// Read alignments
	for(typename TChrIntervals::const_iterator vRIt = validRegions[refIndex].begin(); vRIt != validRegions[refIndex].end(); ++vRIt) {
	  hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, vRIt->lower(), vRIt->upper());
	  bam1_t* rec = bam_init1();
	  int32_t lastAlignedPos = 0;
	  std::set<std::size_t> lastAlignedPosReads;
	  while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	    if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
	    if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;

	    unsigned seed = hash_string(bam_get_qname(rec));
	    
	    // SV detection using single-end read
	    if (c.indels) {
	      uint32_t rp = rec->core.pos; // reference pointer
	      uint32_t sp = 0; // sequence pointer

	      // Parse the CIGAR
	      uint32_t* cigar = bam_get_cigar(rec);
	      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
		if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
		  // match or mismatch
		  for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k) {
		    ++sp;
		    ++rp;
		  }
		} else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
		  if (bam_cigar_oplen(cigar[i]) > c.minRefSep) _insertJunction(readBp, seed, rec, rp, sp, false);
		  rp += bam_cigar_oplen(cigar[i]);
		  if (bam_cigar_oplen(cigar[i]) > c.minRefSep) _insertJunction(readBp, seed, rec, rp, sp, true);
		} else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
		  sp += bam_cigar_oplen(cigar[i]);
		} else if ((bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) {
		  int32_t finalsp = sp;
		  bool scleft = false;
		  if (sp == 0) {
		    finalsp += bam_cigar_oplen(cigar[i]); // Leading soft-clip / hard-clip
		    scleft = true;
		  }
		  sp += bam_cigar_oplen(cigar[i]);
		  if (bam_cigar_oplen(cigar[i]) > c.minClip) _insertJunction(readBp, seed, rec, rp, finalsp, scleft);
		} else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
		  rp += bam_cigar_oplen(cigar[i]);
		} else {
		  std::cerr << "Warning: Unknown Cigar operation!" << std::endl;
		}
	      }
	    }
	    
	    // Paired-end clustering
	    if (rec->core.flag & BAM_FPAIRED) {
	      // Single-end library
	      if (sampleLib[file_c].median == 0) continue; // Single-end library

	      // Secondary/supplementary alignments, mate unmapped or blacklisted chr
	      if (rec->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
	      if ((rec->core.mtid<0) || (rec->core.flag & BAM_FMUNMAP)) continue;
	      if (validRegions[rec->core.mtid].empty()) continue;
	      if ((_translocation(rec)) && (rec->core.qual < c.minTraQual)) continue;

	      // SV type
	      int32_t svt = _isizeMappingPos(rec, overallMaxISize);
	      if (svt == -1) continue;
	      if ((c.svtcmd) && (c.svtset.find(svt) == c.svtset.end())) continue;

	      // Check library-specific insert size for deletions
	      if ((svt == 2) && (sampleLib[file_c].maxISizeCutoff > std::abs(rec->core.isize))) continue;
	      
	      // Clean-up the read store for identical alignment positions
	      if (rec->core.pos > lastAlignedPos) {
		lastAlignedPosReads.clear();
		lastAlignedPos = rec->core.pos;
	      }
	      
	      // Get or store the mapping quality for the partner
	      if (_firstPairObs(rec, lastAlignedPosReads)) {
		// First read
		lastAlignedPosReads.insert(seed);
		std::size_t hv = hash_pair(rec);
		if (_translocation(svt)) matetra[file_c][hv]= std::make_pair((uint8_t) rec->core.qual, alignmentLength(rec));
		else mateMap[hv]= std::make_pair((uint8_t) rec->core.qual, alignmentLength(rec));
	      } else {
		// Second read
		std::size_t hv = hash_pair_mate(rec);
		int32_t alenmate = 0;
		uint8_t pairQuality = 0;
	        if (_translocation(svt)) {
		  // Inter-chromosomal
		  if ((matetra[file_c].find(hv) == matetra[file_c].end()) || (!matetra[file_c][hv].first)) continue; // Mate discarded
		  TQualLen p = matetra[file_c][hv];
		  pairQuality = std::min((uint8_t) p.first, (uint8_t) rec->core.qual);
		  alenmate = p.second;
		  matetra[file_c][hv].first = 0;
		} else {
		  // Intra-chromosomal
		  if ((mateMap.find(hv) == mateMap.end()) || (!mateMap[hv].first)) continue; // Mate discarded
		  TQualLen p = mateMap[hv];
		  pairQuality = std::min((uint8_t) p.first, (uint8_t) rec->core.qual);
		  alenmate = p.second;
		  mateMap[hv].first = 0;
		}
				
#pragma omp critical
		{
		  bamRecord[svt].push_back(BamAlignRecord(rec, pairQuality, alignmentLength(rec), alenmate, sampleLib[file_c].median, sampleLib[file_c].mad, sampleLib[file_c].maxNormalISize));
		}
		++sampleLib[file_c].abnormal_pairs;
	      }
	    }
	  }
	  bam_destroy1(rec);
	  hts_itr_destroy(iter);
	}
      }

      // Process all junctions for this BAM file
      for(typename TReadBp::iterator it = readBp.begin(); it != readBp.end(); ++it) {
	std::sort(it->second.begin(), it->second.end(), SortJunction<Junction>());
      }
	
      // Collect split-read SVs
#pragma omp critical
      {
	selectDeletions(c, readBp, srBR);
	selectDuplications(c, readBp, srBR);
	selectInversions(c, readBp, srBR);
	selectInsertions(c, readBp, srBR);
	selectTranslocations(c, readBp, srBR);
      }
    }

    // Debug abnormal paired-ends and split-reads
    //outputSRBamRecords(c, srBR);

    // Cluster split-read records
    TVariants srSVs;
    for(uint32_t svt = 0; svt < srBR.size(); ++svt) {
      // Sort
      std::sort(srBR[svt].begin(), srBR[svt].end(), SortSRBamRecord<SRBamRecord>());

      // Cluster
      cluster(c, srBR[svt], srSVs);
      for(uint32_t i = 0; i < srSVs.size(); ++i) {
	if (srSVs[i].svt == -1) srSVs[i].svt = svt;
      }
    }

    // Debug SR SVs
    //outputStructuralVariants(c, srSVs);
      
    // Maximum insert size
    int32_t varisize = getVariability(c, sampleLib);
      
#pragma omp parallel for default(shared)
    for(int32_t svt = 0; svt < (int32_t) bamRecord.size(); ++svt) {
      if (bamRecord[svt].empty()) continue;
	
      // Sort BAM records according to position
      std::sort(bamRecord[svt].begin(), bamRecord[svt].end(), SortBamRecords<BamAlignRecord>());
      
      // Components
      typedef std::vector<uint32_t> TComponent;
      TComponent comp;
      comp.resize(bamRecord[svt].size(), 0);
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
      for(TBamRecord::const_iterator bamIt = bamRecord[svt].begin(); bamIt != bamRecord[svt].end(); ++bamIt, ++bamItIndex) {
	// Safe to clean the graph?
	if (bamItIndex > lastConnectedNode) {
	  // Clean edge lists
	  if (!compEdge.empty()) {
	    _searchCliques(hdr, compEdge, bamRecord[svt], svs, svt);
	    lastConnectedNodeStart = lastConnectedNode;
	    compEdge.clear();
	  }
	}
	int32_t const minCoord = _minCoord(bamIt->pos, bamIt->mpos, svt);
	int32_t const maxCoord = _maxCoord(bamIt->pos, bamIt->mpos, svt);
	TBamRecord::const_iterator bamItNext = bamIt;
	++bamItNext;
	std::size_t bamItIndexNext = bamItIndex + 1;
	for(; ((bamItNext != bamRecord[svt].end()) && (abs(_minCoord(bamItNext->pos, bamItNext->mpos, svt) + bamItNext->alen - minCoord) <= varisize)) ; ++bamItNext, ++bamItIndexNext) {
	  // Check that mate chr agree (only for translocations)
	  if (bamIt->mtid != bamItNext->mtid) continue;
	  
	  // Check combinability of pairs
	  if (_pairsDisagree(minCoord, maxCoord, bamIt->alen, bamIt->maxNormalISize, _minCoord(bamItNext->pos, bamItNext->mpos, svt), _maxCoord(bamItNext->pos, bamItNext->mpos, svt), bamItNext->alen, bamItNext->maxNormalISize, svt)) continue;
	  
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
	    TWeightType weight = (TWeightType) ( std::log((double) abs( abs( (_minCoord(bamItNext->pos, bamItNext->mpos, svt) - minCoord) - (_maxCoord(bamItNext->pos, bamItNext->mpos, svt) - maxCoord) ) - abs(bamIt->Median - bamItNext->Median)) + 1) / std::log(2) );
	    compEdgeIt->second.push_back(TEdgeRecord(bamItIndex, bamItIndexNext, weight));
	  }
	}
      }
      if (!compEdge.empty()) {
	_searchCliques(hdr, compEdge, bamRecord[svt], svs, svt);
	compEdge.clear();
      }
    }

    // Sort SVs for look-up
    sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());
    
    // Add the soft clip SV records
    if (c.indels) {
      /*
	for(int32_t svt = 2; svt < 5; svt += 2) {
	int32_t bpWindowLen = 10;
	int32_t maxLookAhead = 0;
	std::sort(splitRecord[svt/2 - 1].begin(), splitRecord[svt/2 - 1].end(), SortSplitRecords<SplitAlignRecord>());
	TSplitRecord::const_iterator splitClusterIt = splitRecord[svt/2 - 1].end();
	for(TSplitRecord::const_iterator splitIt = splitRecord[svt/2 - 1].begin(); splitIt!=splitRecord[svt/2 - 1].end(); ++splitIt) {
	if ((maxLookAhead) && (splitIt->splitbeg > maxLookAhead)) {
	// Process split read cluster
	_processSRCluster(splitClusterIt, splitIt, refIndex, bpWindowLen, svs, svt);
	maxLookAhead = 0;
	splitClusterIt = splitRecord[svt/2 - 1].end();
	}
	if ((!maxLookAhead) || (splitIt->splitbeg < maxLookAhead)) {
	if (!maxLookAhead) splitClusterIt = splitIt;
	maxLookAhead = splitIt->splitbeg + bpWindowLen;
	}
	}
	TSplitRecord::const_iterator splitIt = splitRecord[svt/2 - 1].end();
	_processSRCluster(splitClusterIt, splitIt, refIndex, bpWindowLen, svs, svt);
	}
      */
    }
  
    // Split-read search
    if (!svs.empty()) {
      findPutativeSplitReads(c, validRegions, svs);
      
      if (c.indels) {
	// Sort SVs for look-up and by decreasing PE support
	sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());
	
	// Temporary SV container
	TVariants svc;
	
	// Clean-up SV set
	for(int32_t svt = 0; svt < 10; ++svt) {
	  for(typename TVariants::iterator svIt = svs.begin(); svIt != svs.end(); ++svIt) {
	    if (svIt->svt != svt) continue;

	    // Unresolved soft clips
	    if ((svIt->precise) && (svIt->srAlignQuality == 0)) continue;
	  
	    // Precise duplicates
	    int32_t searchWindow = 10;
	    bool svExists = false;
	    typename TVariants::iterator itOther = std::lower_bound(svc.begin(), svc.end(), StructuralVariantRecord(svIt->chr, std::max(0, svIt->svStart - searchWindow), svIt->svEnd), SortSVs<StructuralVariantRecord>());
	    for(; ((itOther != svc.end()) && (std::abs(itOther->svStart - svIt->svStart) < searchWindow)); ++itOther) {
	      if (itOther->svt != svt) continue;
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
	}
	
	// Final set of precise and imprecise SVs
	svs = svc;
      }
    }

    // Clean-up
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }
  

}

#endif
