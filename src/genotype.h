#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include <htslib/sam.h>

#include "util.h"

namespace torali
{

  inline int32_t
  _editDistanceNW(std::string const& query, std::string const& target) {
    EdlibAlignResult align = edlibAlign(query.c_str(), query.size(), target.c_str(), target.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
    // Debug: requires EDLIB_TASK_PATH otherwise EDLIB_TASK_DISTANCE
    //printAlignment(query, target, EDLIB_MODE_NW, align);
    int32_t ed = align.editDistance;
    edlibFreeAlignResult(align);
    return ed;
  }

  inline int32_t
  _readStart(bam1_t* rec) {
    uint32_t rp = rec->core.pos;
    uint32_t* cigar = bam_get_cigar(rec);
    if (rec->core.n_cigar) {
      if ((bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP)) {
	if (rp > bam_cigar_oplen(cigar[0])) rp -= bam_cigar_oplen(cigar[0]);
	else rp = 0;
      }
    }
    return rp;
  }

  inline int32_t
  _readEnd(bam1_t* rec) {
    uint32_t rp = rec->core.pos;
    uint32_t* cigar = bam_get_cigar(rec);
    if (rec->core.n_cigar) {
      for (uint32_t i = 0; i < rec->core.n_cigar; ++i) {
	if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF) || (bam_cigar_op(cigar[i]) == BAM_CDEL) || (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP)) rp += bam_cigar_oplen(cigar[i]);
      }
      if ((bam_cigar_op(cigar[rec->core.n_cigar - 1]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[rec->core.n_cigar - 1]) == BAM_CHARD_CLIP)) {
	rp += bam_cigar_oplen(cigar[rec->core.n_cigar - 1]);
      }
    }
    return rp;
  }

  inline int32_t
  _findSeqBp(bam1_t* rec, uint32_t const pos) {
    uint32_t rp = rec->core.pos; // reference pointer
    uint32_t sp = 0; // sequence pointer

    // Parse the CIGAR
    uint32_t* cigar = bam_get_cigar(rec);
    if (rec->core.n_cigar) {
      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	  for(uint32_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k, ++rp, ++sp) {
	    if (rp >= pos) return sp;
	  }
	} else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	  rp += bam_cigar_oplen(cigar[i]);
	  if (rp >= pos) return sp;
	} else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	  sp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	  rp += bam_cigar_oplen(cigar[i]);
	  if (rp >= pos) return sp;
	} else if ((bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) {
	  sp += bam_cigar_oplen(cigar[i]);
	} else {
	  std::cerr << "Unknown Cigar options" << std::endl;
	}
      }
      if ((bam_cigar_op(cigar[rec->core.n_cigar - 1]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[rec->core.n_cigar - 1]) == BAM_CHARD_CLIP)) {
	return sp - bam_cigar_oplen(cigar[rec->core.n_cigar - 1]);
      }
    }
    return -1;
  }

  template<typename TConfig, typename TJunctionMap, typename TReadCountMap>
  inline void
  genotypeLR(TConfig& c, std::vector<StructuralVariantRecord>& svs, TJunctionMap& jctMap, TReadCountMap& covMap) {
    typedef std::vector<StructuralVariantRecord> TSVs;
    typedef std::vector<uint8_t> TQuality;
    if (svs.empty()) return;

    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    typedef std::vector<bam_hdr_t*> THeader;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    THeader hdr(c.files.size());
    for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
      hdr[file_c] = sam_hdr_read(samfile[file_c]);
    }

    // Count aligned reads per SV
    typedef std::vector<uint32_t> TSVReadCount;
    typedef std::vector<TSVReadCount> TSVFileReadCount;
    TSVFileReadCount readSV(c.files.size());
    
    // Ref aligned reads
    typedef std::vector<uint32_t> TRefAlignCount;
    typedef std::vector<TRefAlignCount> TFileRefAlignCount;
    TFileRefAlignCount refAlignedReadCount(c.files.size(), TRefAlignCount());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      readSV[file_c].resize(svs.size(), 0);
      refAlignedReadCount[file_c].resize(svs.size(), 0);
    }

    // Dump file
    boost::iostreams::filtering_ostream dumpOut;
    if (c.hasDumpFile) {
      dumpOut.push(boost::iostreams::gzip_compressor());
      dumpOut.push(boost::iostreams::file_sink(c.dumpfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
      dumpOut << "#svid\tbam\tqname\tchr\tpos\tmapq\ttype" << std::endl;
    }

    // Genotype SVs
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "SV annotation" << std::endl;
    
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < (int32_t) hdr[0]->n_targets; ++refIndex) {
      // Fetch breakpoints
      typedef std::multimap<int32_t, int32_t> TBreakpointMap;
      TBreakpointMap bpMap;
      for(typename TSVs::iterator itSV = svs.begin(); itSV != svs.end(); ++itSV) {
	if (itSV->chr == refIndex) bpMap.insert(std::make_pair(itSV->svStart, itSV->id));
	if (itSV->chr2 == refIndex) bpMap.insert(std::make_pair(itSV->svEnd, itSV->id));
      }
      if (bpMap.empty()) continue;
      
      // Load sequence
      int32_t seqlen = -1;
      std::string tname(hdr[0]->target_name[refIndex]);
      char* seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr[0]->target_len[refIndex], &seqlen);

      // Take care of symbolic ALTs
      for(typename TSVs::iterator itSV = svs.begin(); itSV != svs.end(); ++itSV) {
	if ((itSV->chr == refIndex) && (itSV->alleles.empty())) itSV->alleles = _addAlleles(boost::to_upper_copy(std::string(seq + itSV->svStart - 1, seq + itSV->svStart)), std::string(hdr[0]->target_name[itSV->chr2]), *itSV, itSV->svt);
      }
      
#pragma omp parallel for default(shared)      
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {

	// Coverage track
	typedef uint16_t TMaxCoverage;
	uint32_t maxCoverage = std::numeric_limits<TMaxCoverage>::max();
	typedef std::vector<TMaxCoverage> TBpCoverage;
	TBpCoverage covBases(hdr[file_c]->target_len[refIndex], 0);

	// Parse BAM
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr[file_c]->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  // Coverage track
	  if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
	  if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;

	  // Annotate coverage
	  {
	    uint32_t rp = rec->core.pos; // reference pointer
	    uint32_t* cigar = bam_get_cigar(rec);
	    for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
		for(uint32_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k, ++rp) {
		  if ((rp < hdr[file_c]->target_len[refIndex]) && (covBases[rp] < maxCoverage - 1)) ++covBases[rp];
		}
	      } else if ((bam_cigar_op(cigar[i]) == BAM_CDEL) || (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP)) {
		rp += bam_cigar_oplen(cigar[i]);
	      }
	    }
	  }
	  
	  // Only primary alignments for genotyping (full sequence)
	  if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) continue;
	  if (rec->core.l_qseq < 2 * c.minimumFlankSize) continue;
	  
	  // Overlaps any SV breakpoints?
	  typedef std::set<int32_t> TSVSet;
	  TSVSet process;
	  int32_t rStart = _readStart(rec) + c.minimumFlankSize;
	  int32_t rEnd = _readEnd(rec);
	  if (rEnd > c.minimumFlankSize) {
	    rEnd -= c.minimumFlankSize;
	    if (rStart < rEnd) {
	      TBreakpointMap::const_iterator itBegin = bpMap.lower_bound(rStart);
	      TBreakpointMap::const_iterator itEnd = bpMap.upper_bound(rEnd);
	      for(;((itBegin != itEnd) && (itBegin != bpMap.end())); ++itBegin) process.insert(itBegin->second);
	    }
	  }

	  // Genotype SVs
	  std::string sequence;
	  for(typename TSVSet::const_iterator it = process.begin(); it != process.end(); ++it) {
	    int32_t svid = *it;
	    if ((jctMap[file_c][svid].ref.size() + jctMap[file_c][svid].alt.size()) >= c.maxGenoReadCount) continue;
	    // Enough candidates?
	    if (readSV[file_c][svid] >= c.maxGenoReadCount) continue;
	    ++readSV[file_c][svid];

	    // Which SV breakpoint does the read overlap
	    std::vector<int32_t> candidates;
	    if ((svs[svid].chr == refIndex) && (svs[svid].svStart >= rStart) && (svs[svid].svStart <= rEnd)) candidates.push_back(svs[svid].svStart);
	    if ((svs[svid].chr2 == refIndex) && (svs[svid].svEnd >= rStart) && (svs[svid].svEnd <= rEnd)) candidates.push_back(svs[svid].svEnd);
	    if (candidates.empty()) continue;
	    
	    // Genotype breakpoints
	    std::vector<double> scoreR(candidates.size());
	    std::vector<double> scoreA(candidates.size());
	    for(uint32_t i = 0; i < candidates.size(); ++i) {
	      int32_t pos = candidates[i];
	      int32_t spBp = _findSeqBp(rec, pos);	      
	      int32_t consBp = svs[svid].consBp;
	      if (pos == svs[svid].svEnd) consBp += svs[svid].insLen;

	      // Find flanking sequence offsets
	      int32_t rStartOffset = pos - std::max(0, pos - spBp);
	      int32_t rEndOffset = std::min(pos + rec->core.l_qseq - spBp, (int32_t) hdr[file_c]->target_len[refIndex]) - pos;
	      int32_t cStartOffset = consBp - std::max(0, consBp - spBp);
	      int32_t cEndOffset = std::min(consBp + rec->core.l_qseq - spBp, (int32_t) svs[svid].consensus.size()) - consBp;
	      // Breakpoint should be in the middle so flanking sequences do not bias edit distance
	      int32_t offset = std::min(std::min(rStartOffset, cStartOffset), std::min(rEndOffset, cEndOffset));
	      if (offset < c.minimumFlankSize) continue;
	      if (2 * offset < c.minConsWindow) continue;

	      // Load sequence
	      if (sequence.empty()) {
		sequence.resize(rec->core.l_qseq);
		uint8_t* seqptr = bam_get_seq(rec);
		for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	      }
	      //std::cerr << hdr[file_c]->target_name[refIndex] << ":" << svs[svid].svStart << "-" << svs[svid].svEnd << "," << svs[svid].svt << "," << (int) (rec->core.flag & BAM_FREVERSE) << "," << ((pos == svs[svid].svStart) ? "Start" : "End") << std::endl;
	      std::string ref = boost::to_upper_copy(std::string(seq + pos - offset, seq + pos + offset));
	      std::string alt = svs[svid].consensus.substr(consBp - offset, 2 * offset);
	      std::string probe = sequence.substr(spBp - offset, 2 * offset);
	      
	      // Edit distances
	      int32_t refScore = _editDistanceNW(ref, probe);
	      if ( ((svs[svid].svt == 0) && (pos == svs[svid].svEnd)) ||
		   ((svs[svid].svt == 1) && (pos == svs[svid].svStart)) ||
		   ((svs[svid].svt == 5) && (pos == svs[svid].svEnd)) ||
		   ((svs[svid].svt == 6) && (pos == svs[svid].svStart))
		   ) {
		reverseComplement(probe);
	      }
	      int32_t altScore = _editDistanceNW(alt, probe);
	      
	      scoreA[i] = (1.0 - c.flankQuality) * alt.size();
	      scoreR[i] = (1.0 - c.flankQuality) * ref.size();
	      scoreA[i] = scoreA[i] / (double) altScore;
	      scoreR[i] = scoreR[i] / (double) refScore;
	    }

	    // Pick better alignment
	    double scoreRef = scoreR[0];
	    double scoreAlt = scoreA[0];
	    if (candidates.size() == 2) {
	      scoreAlt = std::max(scoreA[0], scoreA[1]);
	      scoreRef = std::max(scoreR[0], scoreR[1]);
	    }
	    //std::cerr << bam_get_qname(rec) << "," << (int) (rec->core.flag & BAM_FREVERSE) << "," << hdr[file_c]->target_name[refIndex] << ":" << svs[svid].svStart << "-" << svs[svid].svEnd << "(" << svs[svid].svt << ")" << scoreAlt << ',' << scoreRef << std::endl;

	    // Any confident alignment?
	    if ((scoreRef > 0.8) || (scoreAlt > 0.8)) {
	      if (scoreRef > scoreAlt) {
		// Account for reference bias
		if (++refAlignedReadCount[file_c][svid] % 2) {
		  TQuality quality;
		  quality.resize(rec->core.l_qseq);
		  uint8_t* qualptr = bam_get_qual(rec);
		  for (int i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
		  uint32_t rq = scoreRef * 35;
		  if (rq >= c.minGenoQual) {
#pragma omp critical
		    {
		      jctMap[file_c][svid].ref.push_back((uint8_t) std::min(rq, (uint32_t) rec->core.qual));
		    }
		  }
		}
	      } else {
		TQuality quality;
		quality.resize(rec->core.l_qseq);
		uint8_t* qualptr = bam_get_qual(rec);
		for (int i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
		uint32_t aq = scoreAlt * 35;
		if (aq >= c.minGenoQual) {
#pragma omp critical
		  {
		    if (c.hasDumpFile) {
		      std::string svidStr(_addID(svs[svid].svt));
		      std::string padNumber = boost::lexical_cast<std::string>(svid);
		      padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
		      svidStr += padNumber;
		      dumpOut << svidStr << "\t" << c.files[file_c].string() << "\t" << bam_get_qname(rec) << "\t" << hdr[file_c]->target_name[rec->core.tid] << "\t" << rec->core.pos << "\t" << (int32_t) rec->core.qual << "\tSR" << std::endl;
		    }
		    jctMap[file_c][svid].alt.push_back((uint8_t) std::min(aq, (uint32_t) rec->core.qual));
		  }
		}
	      }
	    }
	  }
	}
	// Clean-up
	bam_destroy1(rec);
	hts_itr_destroy(iter);

	// Coverage annotation
	for(uint32_t i = 0; i < svs.size(); ++i) {
	  if (svs[i].chr == refIndex) {
	    int32_t halfSize = (svs[i].svEnd - svs[i].svStart)/2;
	    if ((_translocation(svs[i].svt)) || (svs[i].svt == 4)) halfSize = 500;

	    // Left region
	    int32_t lstart = std::max(svs[i].svStart - halfSize, 0);
	    int32_t lend = svs[i].svStart;
	    int32_t covbase = 0;
	    for(uint32_t k = lstart; ((k < (uint32_t) lend) && (k < hdr[file_c]->target_len[refIndex])); ++k) covbase += covBases[k];
	    covMap[file_c][svs[i].id].leftRC = covbase;
	  
	    // Actual SV
	    covbase = 0;
	    int32_t mstart = svs[i].svStart;
	    int32_t mend = svs[i].svEnd;
	    if ((_translocation(svs[i].svt)) || (svs[i].svt == 4)) {
	      mstart = std::max(svs[i].svStart - halfSize, 0);
	      mend = std::min(svs[i].svStart + halfSize, (int32_t) hdr[file_c]->target_len[refIndex]);
	    }
	    for(uint32_t k = mstart; ((k < (uint32_t) mend) && (k < hdr[file_c]->target_len[refIndex])); ++k) covbase += covBases[k];
	    covMap[file_c][svs[i].id].rc = covbase;
	  
	    // Right region
	    covbase = 0;
	    int32_t rstart = svs[i].svEnd;
	    int32_t rend = std::min(svs[i].svEnd + halfSize, (int32_t) hdr[file_c]->target_len[refIndex]);
	    if ((_translocation(svs[i].svt)) || (svs[i].svt == 4)) {
	      rstart = svs[i].svStart;
	      rend = std::min(svs[i].svStart + halfSize, (int32_t) hdr[file_c]->target_len[refIndex]);
	    }
	    for(uint32_t k = rstart; ((k < (uint32_t) rend) && (k < hdr[file_c]->target_len[refIndex])); ++k) covbase += covBases[k];
	    covMap[file_c][svs[i].id].rightRC = covbase;
	  }
	}
      }
      // Clean-up chromosome sequence
      if (seq != NULL) free(seq);
    }
    // Clean-up
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      bam_hdr_destroy(hdr[file_c]);	  
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }

}

#endif
