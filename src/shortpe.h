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
  
  template<typename TConfig, typename TValidRegion, typename TSRStore, typename TStructuralVariantRecord>
  inline void
  assembleSplitReads(TConfig const& c, TValidRegion const& validRegions, TSRStore const& srStore, std::vector<TStructuralVariantRecord>& svs) 
  {
    typedef typename TValidRegion::value_type TChrIntervals;
    typedef typename TSRStore::value_type TPosReadSV;

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

    // Reads per SV
    typedef std::set<std::string> TSequences;
    typedef std::vector<TSequences> TSVSequences;
    TSVSequences traStore(svs.size(), TSequences());
    uint32_t maxReadPerSV = 20;
    
    // Parse BAM
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Split-read assembly" << std::endl;
    boost::progress_display show_progress( 2 * hdr->n_targets );

    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      ++show_progress;
      if (validRegions[refIndex].empty()) continue;
      if (srStore[refIndex].empty()) continue;

      // Load sequence
      char* seq = NULL;
      int32_t seqlen = -1;
      std::string tname(hdr->target_name[refIndex]);
      seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
      
      // Collect all split-read pos
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet hits(hdr->target_len[refIndex]);
      for(typename TPosReadSV::const_iterator it = srStore[refIndex].begin(); it != srStore[refIndex].end(); ++it) hits[it->first.first] = 1;

      // Sequences
      TSVSequences seqStore(svs.size(), TSequences());
      
      // Collect reads from all samples
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	// Read alignments
	for(typename TChrIntervals::const_iterator vRIt = validRegions[refIndex].begin(); vRIt != validRegions[refIndex].end(); ++vRIt) {
	  hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, vRIt->lower(), vRIt->upper());
	  bam1_t* rec = bam_init1();
	  while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	    if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
	    if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	    if (!hits[rec->core.pos]) continue;

	    // Valid split-read
	    std::size_t seed = hash_string(bam_get_qname(rec));
	    typename TPosReadSV::const_iterator it = srStore[refIndex].find(std::make_pair(rec->core.pos, seed));
	    if (it != srStore[refIndex].end()) {
	      int32_t svid = it->second;

	      // Get the sequence
	      if (it->second == (int32_t) svs[svid].id) {  // Should be always true
		std::string sequence;
		sequence.resize(rec->core.l_qseq);
		uint8_t* seqptr = bam_get_seq(rec);
		for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

		// Adjust orientation
		bool bpPoint = false;
		if (_translocation(rec)) {
		  if (rec->core.tid == svs[svid].chr2) bpPoint = true;
		} else {
		  if (rec->core.pos > svs[svid].svStart) bpPoint = true;
		}
		_adjustOrientation(sequence, bpPoint, svs[svid].svt);
		
		// At most n split-reads
		if (seqStore[it->second].size() < maxReadPerSV) {
		  if (_translocation(svs[svid].svt)) traStore[it->second].insert(sequence);
		  else seqStore[it->second].insert(sequence);
		}
	      }
	    }
	  }
	}
      }

      // Process all SVs on this chromosome
      for(uint32_t svid = 0; svid < seqStore.size(); ++svid) {
	if (_translocation(svs[svid].svt)) continue;
	if (svs[svid].chr != refIndex) continue;

	// MSA
	bool msaSuccess = false;
	if (seqStore[svid].size() > 1) {
	  msa(c, seqStore[svid], svs[svid].consensus);
	  if (alignConsensus(c, hdr, seq, NULL, svs[svid])) msaSuccess = true;
	}
	if (!msaSuccess) {
	  svs[svid].consensus = "";
	  svs[svid].srSupport = 0;
	  svs[svid].srAlignQuality = 0;
	} else {
	  svs[svid].srSupport = seqStore[svid].size();
	}
      }
      // Clean-up
      if (seq != NULL) free(seq);
    }

    // Process translocations
    for(int32_t refIndex2 = 0; refIndex2 < hdr->n_targets; ++refIndex2) {
      ++show_progress;
      if (validRegions[refIndex2].empty()) continue;
      char* sndSeq = NULL;
      for(int32_t refIndex = refIndex2 + 1; refIndex < hdr->n_targets; ++refIndex) {
	if (validRegions[refIndex].empty()) continue;
	char* seq = NULL;

	// Iterate SVs
	for(uint32_t svid = 0; svid < traStore.size(); ++svid) {
	  if (!_translocation(svs[svid].svt)) continue;
	  if ((svs[svid].chr != refIndex) || (svs[svid].chr2 != refIndex2)) continue;

	  bool msaSuccess = false;
	  if (traStore[svid].size() > 1) {
	    // Lazy loading of references
	    if (seq == NULL) {
	      int32_t seqlen = -1;
	      std::string tname(hdr->target_name[refIndex]);
	      seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
	    }
	    if (sndSeq == NULL) {
	      int32_t seqlen = -1;
	      std::string tname(hdr->target_name[refIndex2]);
	      sndSeq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex2], &seqlen);
	    }
	    msa(c, traStore[svid], svs[svid].consensus);
	    if (alignConsensus(c, hdr, seq, sndSeq, svs[svid])) msaSuccess = true;
	  }
	  if (!msaSuccess) {
	    svs[svid].consensus = "";
	    svs[svid].srSupport = 0;
	    svs[svid].srAlignQuality = 0;
	  } else {
	    svs[svid].srSupport = traStore[svid].size();
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
  }

      
  template<typename TConfig, typename TValidRegion, typename TSRStore, typename TSampleLib>
  inline void
  scanPEandSR(TConfig const& c, TValidRegion const& validRegions, std::vector<StructuralVariantRecord>& svs, std::vector<StructuralVariantRecord>& srSVs, TSRStore& srStore, TSampleLib& sampleLib)
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
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Paired-end and split-read scanning" << std::endl;
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
	if ((!c.svtcmd) || (c.svtset.find(2) != c.svtset.end())) selectDeletions(c, readBp, srBR);
	if ((!c.svtcmd) || (c.svtset.find(3) != c.svtset.end())) selectDuplications(c, readBp, srBR);
	if ((!c.svtcmd) || (c.svtset.find(0) != c.svtset.end()) || (c.svtset.find(1) != c.svtset.end())) selectInversions(c, readBp, srBR);
	if ((!c.svtcmd) || (c.svtset.find(4) != c.svtset.end())) selectInsertions(c, readBp, srBR);
	if ((!c.svtcmd) || (c.svtset.find(DELLY_SVT_TRANS) != c.svtset.end()) || (c.svtset.find(DELLY_SVT_TRANS + 1) != c.svtset.end()) || (c.svtset.find(DELLY_SVT_TRANS + 2) != c.svtset.end()) || (c.svtset.find(DELLY_SVT_TRANS + 3) != c.svtset.end())) selectTranslocations(c, readBp, srBR);
      }
    }

    // Debug abnormal paired-ends and split-reads
    //outputSRBamRecords(c, srBR);

    // Cluster split-read records
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Split-read clustering" << std::endl;
    boost::progress_display spSR( srBR.size() );
    for(uint32_t svt = 0; svt < srBR.size(); ++svt) {
      ++spSR;
      if ((c.svtcmd) && (c.svtset.find(svt) == c.svtset.end())) continue;
      if (srBR[svt].empty()) continue;
      
      // Sort
      std::sort(srBR[svt].begin(), srBR[svt].end(), SortSRBamRecord<SRBamRecord>());

      // Cluster
      cluster(c, srBR[svt], srSVs, c.maxReadSep, svt);
    }

    // Debug SR SVs
    //outputStructuralVariants(c, srSVs);

    // Cluster paired-end records
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Paired-end clustering" << std::endl;
    boost::progress_display spPE( bamRecord.size() );

    // Maximum variability in insert size
    int32_t varisize = getVariability(c, sampleLib);      
    for(int32_t svt = 0; svt < (int32_t) bamRecord.size(); ++svt) {
      ++spPE;
      if ((c.svtcmd) && (c.svtset.find(svt) == c.svtset.end())) continue;
      if (bamRecord[svt].empty()) continue;
	
      // Sort BAM records according to position
      std::sort(bamRecord[svt].begin(), bamRecord[svt].end(), SortBamRecords<BamAlignRecord>());

      // Cluster
      cluster(c, bamRecord[svt], svs, varisize, svt);
    }

    // Track split-reads
    for(uint32_t svt = 0; svt < srBR.size(); ++svt) {
      for(uint32_t i = 0; i < srBR[svt].size(); ++i) {
	// Read assigned?
	if (srBR[svt][i].svid != -1) {
	  srStore[srBR[svt][i].chr][std::make_pair(srBR[svt][i].pos, srBR[svt][i].id)] = srBR[svt][i].svid;
	  srStore[srBR[svt][i].chr2][std::make_pair(srBR[svt][i].pos2, srBR[svt][i].id)] = srBR[svt][i].svid;
	}
      }
    }

    // Clean-up
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }


  inline void
  mergeSort(std::vector<StructuralVariantRecord>& pe, std::vector<StructuralVariantRecord>& sr) {
    typedef typename std::vector<StructuralVariantRecord> TVariants;
    // Sort PE records for look-up
    sort(pe.begin(), pe.end(), SortSVs<StructuralVariantRecord>());

    // Sort SR records for look-up
    sort(sr.begin(), sr.end(), SortSVs<StructuralVariantRecord>());
    
    // Augment PE SVs and append missing SR SVs
    for(int32_t svt = 0; svt < 10; ++svt) {
      for(uint32_t i = 0; i < sr.size(); ++i) {
	if (sr[i].svt != svt) continue;
	if ((sr[i].srSupport == 0) || (sr[i].srAlignQuality == 0)) continue; // SR assembly failed

	// Precise duplicates
	int32_t searchWindow = 500;
	bool svExists = false;
	typename TVariants::iterator itOther = std::lower_bound(pe.begin(), pe.end(), StructuralVariantRecord(sr[i].chr, std::max(0, sr[i].svStart - searchWindow), sr[i].svEnd), SortSVs<StructuralVariantRecord>());
	for(; ((itOther != pe.end()) && (std::abs(itOther->svStart - sr[i].svStart) < searchWindow)); ++itOther) {
	  if ((itOther->svt != svt) || (itOther->precise)) continue; 
	  if ((sr[i].chr != itOther->chr) || (sr[i].chr2 != itOther->chr2)) continue;  // Mismatching chr

	  // Breakpoints within PE confidence interval?
	  if ((itOther->svStart + itOther->ciposlow < sr[i].svStart) && (sr[i].svStart < itOther->svStart + itOther->ciposhigh)) {
	    if ((itOther->svEnd + itOther->ciendlow < sr[i].svEnd) && (sr[i].svEnd < itOther->svEnd + itOther->ciendhigh)) {
	      svExists = true;
	      // Augment PE record
	      itOther->svStart = sr[i].svStart;
	      itOther->svEnd = sr[i].svEnd;
	      itOther->ciposlow = sr[i].ciposlow;
	      itOther->ciposhigh = sr[i].ciposhigh;
	      itOther->ciendlow = sr[i].ciendlow;
	      itOther->ciendhigh = sr[i].ciendhigh;
	      itOther->srSupport = sr[i].srSupport;
	      itOther->insLen = sr[i].insLen;
	      itOther->homLen = sr[i].homLen;
	      itOther->srAlignQuality = sr[i].srAlignQuality;
	      itOther->precise = true;
	      itOther->consensus = sr[i].consensus;
	    }
	  }
	}
	
	// SR only SV
	if (!svExists) pe.push_back(sr[i]);
      }
    }
  }
  

}

#endif
