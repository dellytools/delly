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

#ifndef COVERAGE_H
#define COVERAGE_H

#include <boost/container/flat_set.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <htslib/sam.h>

#include "tags.h"
#include "util.h"
#include "msa.h"
#include "split.h"


namespace torali {

  struct SpanPoint {
    int32_t bppos;
    int32_t svt;
    uint32_t id;

    SpanPoint() : bppos(0), svt(0), id(0) {}
    explicit SpanPoint(int32_t bp) : bppos(bp), svt(0), id(0) {}
    SpanPoint(int32_t bp, int32_t s, uint32_t identifier) : bppos(bp), svt(s), id(identifier) {}
  };
  
  struct BpRegion {
    int32_t regionStart;
    int32_t regionEnd;
    int32_t bppos;
    int32_t homLeft;
    int32_t homRight;
    int32_t svt;
    uint32_t id;
    uint8_t bpPoint;
    
    BpRegion() : regionStart(0), regionEnd(0), bppos(0), homLeft(0), homRight(0), svt(0), id(0), bpPoint(0) {}
    explicit BpRegion(int32_t bp) : regionStart(0), regionEnd(0), bppos(bp), homLeft(0), homRight(0), svt(0), id(0), bpPoint(0) {}
  BpRegion(int32_t rs, int32_t re, int32_t bpos, int32_t hl, int32_t hr, int32_t s, uint32_t identifier, uint8_t bpp) : regionStart(rs), regionEnd(re), bppos(bpos), homLeft(hl), homRight(hr), svt(s), id(identifier), bpPoint(bpp) {}
  };

  template<typename TRecord>
  struct SortBp : public std::binary_function<TRecord, TRecord, bool> {
    inline bool operator()(TRecord const& s1, TRecord const& s2) const {
      return (s1.bppos < s2.bppos);
    }
  };
  
  // Reduced structural variant record for cov
  struct CovRecord {
    int32_t svStart;
    int32_t svEnd;
    uint32_t id;

    CovRecord() : svStart(0), svEnd(0), id(0) {}
  };

  template<typename TRecord>
  struct SortCovRecord : public std::binary_function<TRecord, TRecord, bool> {
    inline bool operator()(TRecord const& s1, TRecord const& s2) const {
      return ((s1.svStart<s2.svStart) || ((s1.svStart==s2.svStart) && (s1.svEnd<s2.svEnd)) || ((s1.svStart==s2.svStart) && (s1.svEnd==s2.svEnd) && (s1.id < s2.id)));
    }
  };


  struct SpanningCount {
    int32_t refh1;
    int32_t refh2;
    int32_t alth1;
    int32_t alth2;
    std::vector<uint8_t> ref;
    std::vector<uint8_t> alt;

    SpanningCount() : refh1(0), refh2(0), alth1(0), alth2(0) {}
  };

  
  struct JunctionCount {
    int32_t refh1;
    int32_t refh2;
    int32_t alth1;
    int32_t alth2;
    std::vector<uint8_t> ref;
    std::vector<uint8_t> alt;

    JunctionCount() : refh1(0), refh2(0), alth1(0), alth2(0) {}
  };

  template<typename TAlign, typename TQualities>
  inline uint32_t
  _getAlignmentQual(TAlign const& align, TQualities const& qual) {
    typedef typename TAlign::index TAIndex;
    uint32_t baseQualSum = 0;
    uint32_t seqPtr = 0;
    uint32_t alignedBases = 0;
    for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
      if (align[1][j] != '-') {
	if (align[0][j] != '-') {
	  ++alignedBases;
	  baseQualSum += qual[seqPtr];
	}
	++seqPtr;
      }
    }
    return (baseQualSum / alignedBases);
  }

  template<typename TPos>
  inline int32_t
  _cutRefStart(TPos const rStart, TPos const rEnd, TPos const offset, unsigned int bpPoint, int32_t const svt) {
    if (_translocation(svt)) {
      uint8_t ct = _getSpanOrientation(svt);
      if (ct == 3) {
	if (!bpPoint) return rEnd - offset;
	else return rStart - offset;
      } else {
	if (bpPoint) return rEnd - offset;
	else return rStart - offset;
      }
    } else {
      if (svt == 3) {
	if (!bpPoint) return rEnd - offset;
	else return rStart - offset;
      } else {
	if (bpPoint) return rEnd - offset;
	else return rStart - offset;
      }
    }
  }

  template<typename TPos>
  inline int32_t
  _cutRefEnd(TPos const rStart, TPos const rEnd, TPos const offset, unsigned int bpPoint, int32_t const svt) {
    if (_translocation(svt)) {
      uint8_t ct = _getSpanOrientation(svt);
      if (ct == 3) {
	if (!bpPoint) return rEnd + offset;
	else return rStart + offset;
      } else {
	if (bpPoint) return rEnd + offset;
	else return rStart + offset;
      }
    } else {
      if (svt == 3) {
	if (!bpPoint) return rEnd + offset;
	else return rStart + offset;
      } else {
	if (bpPoint) return rEnd + offset;
	else return rStart + offset;
      }
    }
  }

  

template<typename TConfig, typename TSampleLibrary, typename TCovRecord, typename TCoverageCount, typename TSVs, typename TCountMap, typename TSpanMap>
inline void
annotateCoverage(TConfig& c, TSampleLibrary& sampleLib, TCovRecord& ict, TCoverageCount& covCount, TSVs& svs, TCountMap& countMap, TSpanMap& spanMap)
{
  typedef typename TSpanMap::value_type::value_type TSpanPair;
  typedef typename TCountMap::value_type::value_type TCountPair;
  typedef std::vector<uint8_t> TQuality;
  
  // Open file handles
  typedef std::vector<samFile*> TSamFile;
  typedef std::vector<hts_idx_t*> TIndex;
  typedef std::vector<bam_hdr_t*> THeader;
  TSamFile samfile(c.files.size());
  TIndex idx(c.files.size());
  THeader hdr(c.files.size());
  int32_t totalTarget = 0;
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
    hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
    idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
    hdr[file_c] = sam_hdr_read(samfile[file_c]);
    totalTarget += hdr[file_c]->n_targets;
  }

  // Sort Coverage Records
  uint32_t totalSVs = 0;
  for(int32_t refIndex=0; refIndex < hdr[0]->n_targets; ++refIndex) {
    sort(ict[refIndex].begin(), ict[refIndex].end(), SortCovRecord<CovRecord>());
    totalSVs += ict[refIndex].size();
  }

  // Initialize coverage count maps
  covCount.resize(c.files.size());
  countMap.resize(c.files.size());
  spanMap.resize(c.files.size());
  for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
    covCount[file_c].resize(totalSVs);
    countMap[file_c].resize(svs.size(), TCountPair());
    spanMap[file_c].resize(svs.size(), TSpanPair());
  }

  // Preprocess REF and ALT
  boost::posix_time::ptime noww = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(noww) << "] " << "Generate REF and ALT probes" << std::endl;
  boost::progress_display show_progresss( hdr[0]->n_targets );
  
  // Reference and consensus probes
  typedef std::vector<std::string> TProbes;
  typedef std::vector<TProbes> TBreakProbes;
  TBreakProbes refProbeArr(2, TProbes()); // Left and right breakpoint
  TBreakProbes consProbeArr(2, TProbes()); // Left and right breakpoint
  for(uint32_t k = 0; k < 2; ++k) {
    refProbeArr[k].resize(svs.size());
    consProbeArr[k].resize(svs.size());
  }
  typedef std::vector<BpRegion> TBpRegion;
  typedef std::vector<TBpRegion> TGenomicBpRegion;
  TGenomicBpRegion bpRegion(hdr[0]->n_targets, TBpRegion());
  typedef std::vector<uint32_t> TRefAlignCount;
  typedef std::vector<TRefAlignCount> TFileRefAlignCount;
  TFileRefAlignCount refAlignedReadCount(c.files.size(), TRefAlignCount());
  TFileRefAlignCount refAlignedSpanCount(c.files.size(), TRefAlignCount());
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    refAlignedReadCount[file_c].resize(svs.size(), 0);
    refAlignedSpanCount[file_c].resize(svs.size(), 0);
  }
  
  // Iterate all structural variants
  {
    TProbes refProbes(svs.size());
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < (int32_t) hdr[0]->n_targets; ++refIndex) {
      ++show_progresss;
      char* seq = NULL;

      // Iterate all structural variants
      for(typename TSVs::iterator itSV = svs.begin(); itSV != svs.end(); ++itSV) {
	if ((itSV->chr != refIndex) && (itSV->chr2 != refIndex)) continue;

	// Lazy loading of reference sequence
	if (seq == NULL) {
	  int32_t seqlen = -1;
	  std::string tname(hdr[0]->target_name[refIndex]);
	  seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr[0]->target_len[refIndex], &seqlen);
	}

	// Set tag alleles
	if (itSV->chr == refIndex) {
	  itSV->alleles = _addAlleles(boost::to_upper_copy(std::string(seq + itSV->svStart - 1, seq + itSV->svStart)), std::string(hdr[0]->target_name[itSV->chr2]), *itSV, itSV->svt);
	}
	if (!itSV->precise) continue;

	// Get the reference sequence
	if ((itSV->chr != itSV->chr2) && (itSV->chr2 == refIndex)) {
	  Breakpoint bp(*itSV);
	  _initBreakpoint(hdr[0], bp, (int32_t) itSV->consensus.size(), itSV->svt);
	  refProbes[itSV->id] = _getSVRef(seq, bp, refIndex, itSV->svt);
	}
	if (itSV->chr == refIndex) {
	  Breakpoint bp(*itSV);
	  bp.part1 = refProbes[itSV->id];
	  _initBreakpoint(hdr[0], bp, (int32_t) itSV->consensus.size(), itSV->svt);
	  std::string svRefStr = _getSVRef(seq, bp, refIndex, itSV->svt);

	  // Find breakpoint to reference
	  typedef boost::multi_array<char, 2> TAlign;
	  TAlign align;
	  if (!_consRefAlignment(itSV->consensus, svRefStr, align, itSV->svt)) continue;
	  AlignDescriptor ad;
	  if (!_findSplit(c, itSV->consensus, svRefStr, align, ad, itSV->svt)) continue;
	  
	  // Debug consensus to reference alignment
	  //for(TAIndex i = 0; i<align.shape()[0]; ++i) {
	  //for(TAIndex j = 0; j<align.shape()[1]; ++j) std::cerr << align[i][j];
	  //std::cerr << std::endl;
	  //}
	  //std::cerr << std::endl;
	  
	  
	  // Iterate all samples
	  for (unsigned int bpPoint = 0; bpPoint<2; ++bpPoint) {
	    int32_t regionChr, regionStart, regionEnd, cutConsStart, cutConsEnd, cutRefStart, cutRefEnd, bppos;
	    if (bpPoint) {
	      regionChr = itSV->chr2;
	      regionStart = std::max(0, itSV->svEnd - c.minimumFlankSize);
	      regionEnd = std::min((uint32_t) (itSV->svEnd + c.minimumFlankSize), hdr[0]->target_len[itSV->chr2]);
	      cutConsStart = ad.cEnd - ad.homLeft - c.minimumFlankSize;
	      cutConsEnd = ad.cEnd + ad.homRight + c.minimumFlankSize;
	      cutRefStart = _cutRefStart(ad.rStart, ad.rEnd, ad.homLeft + c.minimumFlankSize, bpPoint, itSV->svt);
	      cutRefEnd = _cutRefEnd(ad.rStart, ad.rEnd, ad.homRight + c.minimumFlankSize, bpPoint, itSV->svt);
	      bppos = itSV->svEnd;
	    } else {
	      regionChr = itSV->chr;
	      regionStart = std::max(0, itSV->svStart - c.minimumFlankSize);
	      regionEnd = std::min((uint32_t) (itSV->svStart + c.minimumFlankSize), hdr[0]->target_len[itSV->chr]);
	      cutConsStart = ad.cStart - ad.homLeft - c.minimumFlankSize;
	      cutConsEnd = ad.cStart + ad.homRight + c.minimumFlankSize;
	      cutRefStart = _cutRefStart(ad.rStart, ad.rEnd, ad.homLeft + c.minimumFlankSize, bpPoint, itSV->svt);
	      cutRefEnd = _cutRefEnd(ad.rStart, ad.rEnd, ad.homRight + c.minimumFlankSize, bpPoint, itSV->svt);
	      bppos = itSV->svStart;
	    }
	    consProbeArr[bpPoint][itSV->id] = itSV->consensus.substr(cutConsStart, (cutConsEnd - cutConsStart));
	    refProbeArr[bpPoint][itSV->id] = svRefStr.substr(cutRefStart, (cutRefEnd - cutRefStart));
	    bpRegion[regionChr].push_back(BpRegion(regionStart, regionEnd, bppos, ad.homLeft, ad.homRight, itSV->svt, itSV->id, bpPoint));
	  }
	}
      }
      if (seq != NULL) free(seq);
    }
    // Clean-up
    fai_destroy(fai);
    for(int32_t refIndex=0; refIndex < (int32_t) hdr[0]->n_targets; ++refIndex) {
      // Sort breakpoint regions
      std::sort(bpRegion[refIndex].begin(), bpRegion[refIndex].end(), SortBp<BpRegion>());
    }
  }
  // Debug
  //for(uint32_t k = 0; k < 2; ++k) {
  //for(uint32_t i = 0; i < svs.size(); ++i) {
  //std::cerr << k << ',' << i << ',' << refProbeArr[k][i] << ',' << consProbeArr[k][i] << std::endl;
  //}
  //}
  
  // Iterate all samples
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "SV annotation" << std::endl;
  boost::progress_display show_progress( totalTarget );

  // Dump file
  boost::iostreams::filtering_ostream dumpOut;
  if (c.dumpflag) {
    dumpOut.push(boost::iostreams::gzip_compressor());
    dumpOut.push(boost::iostreams::file_sink(c.srpedump.string().c_str(), std::ios_base::out | std::ios_base::binary));
    dumpOut << "#svid\tbam\tqname\tchr\tpos\tmatechr\tmatepos\tmapq\ttype" << std::endl;
  }

  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) { 
    // Parallelize by chromosome, go go go!
    #pragma omp parallel for default(shared)
    for(int32_t refIndex=0; refIndex < (int32_t) hdr[file_c]->n_targets; ++refIndex) {
      ++show_progress;

      // Pair qualities and features
      typedef boost::unordered_map<std::size_t, uint8_t> TQualities;
      TQualities qualities;
      TQualities qualitiestra;
      typedef boost::unordered_map<std::size_t, bool> TClip;
      TClip clip;
      TClip cliptra;
      
      // Any SVs on this chromosome?
      if (ict[refIndex].empty()) continue;

      // Check we have mapped reads on this chromosome
      bool nodata = true;
      std::string suffix("cram");
      std::string str(c.files[file_c].string());
      if ((str.size() >= suffix.size()) && (str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0)) nodata = false;
      uint64_t mapped = 0;
      uint64_t unmapped = 0;
      #pragma omp critical
      {
        hts_idx_get_stat(idx[file_c], refIndex, &mapped, &unmapped);
      }
      if (mapped) nodata = false;
      if (nodata) continue;

      // Open bam file again per chromosome
      samFile *mysamfile;
      hts_idx_t *myidx;
      #pragma omp critical
      {
        mysamfile = sam_open(c.files[file_c].string().c_str(), "r");
        hts_set_fai_filename(mysamfile, c.genome.string().c_str());
        myidx = sam_index_load(mysamfile, c.files[file_c].string().c_str());
      }
      // Coverage track
      typedef uint16_t TCount;
      uint32_t maxCoverage = std::numeric_limits<TCount>::max();
      typedef std::vector<TCount> TCoverage;
      TCoverage covFragment(hdr[file_c]->target_len[refIndex], 0);
      TCoverage covBases(hdr[file_c]->target_len[refIndex], 0);

      // Flag breakpoint regions
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet bpOccupied(hdr[file_c]->target_len[refIndex]);
      for(uint32_t i = 0; i < bpRegion[refIndex].size(); ++i) {
	for(int32_t k = bpRegion[refIndex][i].regionStart; k < bpRegion[refIndex][i].regionEnd; ++k) {
	  bpOccupied[k] = 1;
	}
      }

      // Flag spanning breakpoints
      typedef std::vector<SpanPoint> TSpanPoint;
      TSpanPoint spanPoint;
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet spanBp(hdr[file_c]->target_len[refIndex]);
      for(typename TSVs::iterator itSV = svs.begin(); itSV != svs.end(); ++itSV) {
	if (itSV->peSupport == 0) continue;
	#pragma omp critical
	{
	  if ((itSV->chr == refIndex) && (itSV->svStart < (int32_t) hdr[file_c]->target_len[refIndex])) {
	    spanBp[itSV->svStart] = 1;
	    spanPoint.push_back(SpanPoint(itSV->svStart, itSV->svt, itSV->id));
	  }
	  if ((itSV->chr2 == refIndex) && (itSV->svEnd < (int32_t) hdr[file_c]->target_len[refIndex])) {
	    spanBp[itSV->svEnd] = 1;
	    spanPoint.push_back(SpanPoint(itSV->svEnd, itSV->svt, itSV->id));
	  }
	}
      }
      std::sort(spanPoint.begin(), spanPoint.end(), SortBp<SpanPoint>());
      
      // Count reads
      int32_t lastAlignedPos = 0;
      std::set<std::size_t> lastAlignedPosReads;
      hts_itr_t* iter = sam_itr_queryi(myidx, refIndex, 0, hdr[file_c]->target_len[refIndex]);
      bam1_t* rec = bam_init1();

      while (sam_itr_next(mysamfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP | BAM_FMUNMAP)) continue;
	if (rec->core.qual < c.minGenoQual) continue;

	// Count aligned basepair (small InDels)
	#pragma omp critical
	{
	  uint32_t rp = 0; // reference pointer
	  uint32_t* cigar = bam_get_cigar(rec);
	  for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	    if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
	      for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k) {
		if ((rec->core.pos + rp < hdr[file_c]->target_len[refIndex]) && (covBases[rec->core.pos + rp] < maxCoverage - 1)) ++covBases[rec->core.pos + rp];
		++rp;
	      }
	    } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	      rp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	      rp += bam_cigar_oplen(cigar[i]);
	    }
	  }
	}

	// Any (leading) soft clip
	bool hasSoftClip = false;
	bool hasClip = false;
	int32_t leadingSC = 0;
	uint32_t* cigar = bam_get_cigar(rec);
	for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	  if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	    hasClip = true;
	    hasSoftClip = true;
	    if (i == 0) leadingSC = bam_cigar_oplen(cigar[i]);
	  } else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) hasClip = true;
	}
	
	// Check read length for junction annotation
	if (rec->core.l_qseq >= (2 * c.minimumFlankSize)) {
	  bool bpvalid = false;
	  int32_t rbegin = std::max(0, rec->core.pos - leadingSC);
	  for(int32_t k = rbegin; ((k < (rec->core.pos + rec->core.l_qseq)) && (k < (int32_t) hdr[file_c]->target_len[refIndex])); ++k) {
	    if (bpOccupied[k]) {
	      bpvalid = true;
	      break;
	    }
	  }
	  if (bpvalid) {
	    // Fetch all relevant SVs
	    typename TBpRegion::iterator itBp = std::lower_bound(bpRegion[refIndex].begin(), bpRegion[refIndex].end(), BpRegion(rbegin), SortBp<BpRegion>());
	    for(; ((itBp != bpRegion[refIndex].end()) && (rec->core.pos + rec->core.l_qseq >= itBp->bppos)); ++itBp) {
	      // Read spans breakpoint?
	      if ((hasSoftClip) || ((!hasClip) && (rec->core.pos + c.minimumFlankSize + itBp->homLeft <= itBp->bppos) &&  (rec->core.pos + rec->core.l_qseq >= itBp->bppos + c.minimumFlankSize + itBp->homRight))) {
		std::string consProbe = consProbeArr[itBp->bpPoint][itBp->id];
		std::string refProbe = refProbeArr[itBp->bpPoint][itBp->id];

		// Get sequence
		std::string sequence;
		sequence.resize(rec->core.l_qseq);
		uint8_t* seqptr = bam_get_seq(rec);
		for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
		_adjustOrientation(sequence, itBp->bpPoint, itBp->svt);
		
		// Compute alignment to alternative haplotype
		typedef boost::multi_array<char, 2> TAlign;
		TAlign alignAlt;
		DnaScore<int> simple(5, -4, -4, -4);
		AlignConfig<true, false> semiglobal;
		int32_t scoreA = needle(consProbe, sequence, alignAlt, semiglobal, simple);
		int32_t scoreAltThreshold = (int32_t) (c.flankQuality * consProbe.size() * simple.match + (1.0 - c.flankQuality) * consProbe.size() * simple.mismatch);
		double scoreAlt = (double) scoreA / (double) scoreAltThreshold;
		
		// Compute alignment to reference haplotype
		TAlign alignRef;
		int32_t scoreR = needle(refProbe, sequence, alignRef, semiglobal, simple);
		int32_t scoreRefThreshold = (int32_t) (c.flankQuality * refProbe.size() * simple.match + (1.0 - c.flankQuality) * refProbe.size() * simple.mismatch);
		double scoreRef = (double) scoreR / (double) scoreRefThreshold;

		// Any confident alignment?
		if ((scoreRef > 1) || (scoreAlt > 1)) {
		  // Debug alignment to REF and ALT
		  //std::cerr << "Alt:\t" << scoreAlt << "\tRef:\t" << scoreRef << std::endl;
		  //for(TAIndex i = 0; i< (TAIndex) alignAlt.shape()[0]; ++i) {
		  //for(TAIndex j = 0; j< (TAIndex) alignAlt.shape()[1]; ++j) std::cerr << alignAlt[i][j];
		  //std::cerr << std::endl;
		  //}
		  //for(TAIndex i = 0; i< (TAIndex) alignRef.shape()[0]; ++i) {
		  //for(TAIndex j = 0; j< (TAIndex) alignRef.shape()[1]; ++j) std::cerr << alignRef[i][j];
		  //std::cerr << std::endl;
		  //}
		  
		  if (scoreRef > scoreAlt) {
		    // Account for reference bias
		    if (++refAlignedReadCount[file_c][itBp->id] % 2) {
		      TQuality quality;
		      quality.resize(rec->core.l_qseq);
		      uint8_t* qualptr = bam_get_qual(rec);
		      for (int i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
		      uint32_t rq = _getAlignmentQual(alignRef, quality);
		      if (rq >= c.minGenoQual) {
			uint8_t* hpptr = bam_aux_get(rec, "HP");
#pragma omp critical
			{
			  countMap[file_c][itBp->id].ref.push_back((uint8_t) std::min(rq, (uint32_t) rec->core.qual));
			  if (hpptr) {
			    c.isHaplotagged = true;
			    int hap = bam_aux2i(hpptr);
			    if (hap == 1) ++countMap[file_c][itBp->id].refh1;
			    else ++countMap[file_c][itBp->id].refh2;
			  }
			}
		      }
		    }
		  } else {
		    TQuality quality;
		    quality.resize(rec->core.l_qseq);
		    uint8_t* qualptr = bam_get_qual(rec);
		    for (int i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
		    uint32_t aq = _getAlignmentQual(alignAlt, quality);
		    if (aq >= c.minGenoQual) {
		      uint8_t* hpptr = bam_aux_get(rec, "HP");
#pragma omp critical
		      {
			if (c.dumpflag) {
			  std::string svid(_addID(itBp->svt));
			  std::string padNumber = boost::lexical_cast<std::string>(itBp->id);
			  padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
			  svid += padNumber;
			  dumpOut << svid << "\t" << c.files[file_c].string() << "\t" << bam_get_qname(rec) << "\t" << hdr[file_c]->target_name[rec->core.tid] << "\t" << rec->core.pos << "\t" << hdr[file_c]->target_name[rec->core.mtid] << "\t" << rec->core.mpos << "\t" << (int32_t) rec->core.qual << "\tSR" << std::endl;
			}
			countMap[file_c][itBp->id].alt.push_back((uint8_t) std::min(aq, (uint32_t) rec->core.qual));
			if (hpptr) {
			  c.isHaplotagged = true;
			  int hap = bam_aux2i(hpptr);
			  if (hap == 1) ++countMap[file_c][itBp->id].alth1;
			  else ++countMap[file_c][itBp->id].alth2;
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	
	// Read-count and spanning annotation
	if ((!(rec->core.flag & BAM_FPAIRED)) || (ict[rec->core.mtid].empty())) continue;

	// Clean-up the read store for identical alignment positions
	if (rec->core.pos > lastAlignedPos) {
	  lastAlignedPosReads.clear();
	  lastAlignedPos = rec->core.pos;
	}

	if (_firstPairObs(rec, lastAlignedPosReads)) {
	  // First read
	  std::size_t hv;
	  #pragma omp critical
	  {
	    lastAlignedPosReads.insert(hash_string(bam_get_qname(rec)));
	    hv = hash_pair(rec);
	    if (rec->core.tid == rec->core.mtid) {
	      qualities[hv] = rec->core.qual;
	      clip[hv] = hasSoftClip;
	    } else {
	      qualitiestra[hv] = rec->core.qual;
	      cliptra[hv] = hasSoftClip;
	    }
	  }
	} else {
	  // Second read
	  std::size_t hv;
	  uint8_t pairQuality = 0;
	  bool pairClip = false;
	  bool skip = false;
	  #pragma omp critical
	  {
	    hv = hash_pair_mate(rec);
	    if (rec->core.tid == rec->core.mtid) {
	      if (qualities.find(hv) == qualities.end()) { skip = true; goto discard; } // Mate discarded
	      pairQuality = std::min((uint8_t) qualities[hv], (uint8_t) rec->core.qual);
	      if ((clip[hv]) || (hasSoftClip)) pairClip = true;
	      qualities[hv] = 0;
	      clip[hv] = false;
	    } else {
	      if (qualitiestra.find(hv) == qualitiestra.end()) { skip = true; goto discard; } // Mate discarded
	      pairQuality = std::min((uint8_t) qualitiestra[hv], (uint8_t) rec->core.qual);
	      if ((cliptra[hv]) || (hasSoftClip)) pairClip = true;
	      qualitiestra[hv] = 0;
	      cliptra[hv] = false;
	    }
discard:
	    ; // nop
	  }
          if (skip) continue;

	  // Pair quality
	  if (pairQuality < c.minGenoQual) continue; // Low quality pair

	  // Read-depth fragment counting
	  if (rec->core.tid == rec->core.mtid) {
	    // Count mid point (fragment counting)
	    int32_t midPoint = rec->core.pos + halfAlignmentLength(rec);
	    if ((midPoint < (int32_t) hdr[file_c]->target_len[refIndex]) && (covFragment[midPoint] < maxCoverage - 1)) ++covFragment[midPoint];
	  }

	  // Spanning counting
	  int32_t outerISize = 0;
	  if (rec->core.pos < rec->core.mpos) outerISize = rec->core.mpos + rec->core.l_qseq - rec->core.pos;
	  else outerISize = rec->core.pos + rec->core.l_qseq - rec->core.mpos;

	  // Get the library information
	  if (sampleLib[file_c].median == 0) continue; // Single-end library or non-valid library

	  // Normal spanning pair
	  if ((!pairClip) && (getSVType(rec->core) == 2) && (outerISize >= sampleLib[file_c].minNormalISize) && (outerISize <= sampleLib[file_c].maxNormalISize) && (rec->core.tid==rec->core.mtid)) {
	    // Take X% of the outerisize as the spanned interval
	    int32_t spanlen = 0.8 * outerISize;
	    int32_t pbegin = std::min(rec->core.pos, rec->core.mpos);
	    int32_t st = pbegin + (outerISize - spanlen) / 2;
	    bool spanvalid = false;
	    for(int32_t i = st; ((i < (st + spanlen)) &&  (i < (int32_t) hdr[file_c]->target_len[refIndex])); ++i) {
	      if (spanBp[i]) {
		spanvalid = true;
		break;
	      }
	    }
	    if (spanvalid) {
	      // Fetch all relevant SVs
	      typename TSpanPoint::iterator itSpan = std::lower_bound(spanPoint.begin(), spanPoint.end(), SpanPoint(st), SortBp<SpanPoint>());
	      for(; ((itSpan != spanPoint.end()) && (st + spanlen >= itSpan->bppos)); ++itSpan) {
		// Account for reference bias
		if (++refAlignedSpanCount[file_c][itSpan->id] % 2) {
		  uint8_t* hpptr = bam_aux_get(rec, "HP");
#pragma omp critical
		  {
		    spanMap[file_c][itSpan->id].ref.push_back(pairQuality);
		    if (hpptr) {
		      c.isHaplotagged = true;
		      int hap = bam_aux2i(hpptr);
		      if (hap == 1) ++spanMap[file_c][itSpan->id].refh1;
		      else ++spanMap[file_c][itSpan->id].refh2;
		    }
		  }
		}
	      }
	    }
	  }

	  // Abnormal spanning coverage
	  if ((getSVType(rec->core) != 2) || (outerISize < sampleLib[file_c].minNormalISize) || (outerISize > sampleLib[file_c].maxNormalISize) || (rec->core.tid!=rec->core.mtid)) {
	    // SV type
	    int32_t svt = _isizeMappingPos(rec, sampleLib[file_c].maxISizeCutoff);
	    if (svt == -1) continue;
	    
	    // Spanning a breakpoint?
	    bool spanvalid = false;
	    int32_t pbegin = rec->core.pos;
	    int32_t pend = std::min(rec->core.pos + sampleLib[file_c].maxNormalISize, (int32_t) hdr[file_c]->target_len[refIndex]);
	    if (rec->core.flag & BAM_FREVERSE) {
	      pbegin = std::max(0, rec->core.pos + rec->core.l_qseq - sampleLib[file_c].maxNormalISize);
	      pend = std::min(rec->core.pos + rec->core.l_qseq, (int32_t) hdr[file_c]->target_len[refIndex]);
	    }
	    for(int32_t i = pbegin; i < pend; ++i) {
	      if (spanBp[i]) {
		spanvalid = true;
		break;
	      }
	    }
	    if (spanvalid) {
	      // Fetch all relevant SVs
	      typename TSpanPoint::iterator itSpan = std::lower_bound(spanPoint.begin(), spanPoint.end(), SpanPoint(pbegin), SortBp<SpanPoint>());
	      for(; ((itSpan != spanPoint.end()) && (pend >= itSpan->bppos)); ++itSpan) {
		if (svt == itSpan->svt) {
		  uint8_t* hpptr = bam_aux_get(rec, "HP");
#pragma omp critical
		  {
		    if (c.dumpflag) {
		      std::string svid(_addID(itSpan->svt));
		      std::string padNumber = boost::lexical_cast<std::string>(itSpan->id);
		      padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
		      svid += padNumber;
		      dumpOut << svid << "\t" << c.files[file_c].string() << "\t" << bam_get_qname(rec) << "\t" << hdr[file_c]->target_name[rec->core.tid] << "\t" << rec->core.pos << "\t" << hdr[file_c]->target_name[rec->core.mtid] << "\t" << rec->core.mpos << "\t" << (int32_t) rec->core.qual << "\tPE" << std::endl;
		    }
		    spanMap[file_c][itSpan->id].alt.push_back(pairQuality);
		    if (hpptr) {
		      c.isHaplotagged = true;
		      int hap = bam_aux2i(hpptr);
		      if (hap == 1) ++spanMap[file_c][itSpan->id].alth1;
		      else ++spanMap[file_c][itSpan->id].alth2;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      // Clean-up
      bam_destroy1(rec);
      hts_itr_destroy(iter);
      qualities.clear();
      clip.clear();

      // Assign fragment and base counts to SVs
      for(uint32_t i = 0; i < ict[refIndex].size(); ++i) {
	int32_t covfrag = 0;
	int32_t covbase = 0;
	for(int32_t k = ict[refIndex][i].svStart; ((k < ict[refIndex][i].svEnd) && (k < (int32_t) hdr[file_c]->target_len[refIndex])); ++k) {
	  covfrag += covFragment[k]; 
	  covbase += covBases[k];
	}
	#pragma omp critical
	{
	  // Store counts
	  covCount[file_c][ict[refIndex][i].id].first = covbase;
	  covCount[file_c][ict[refIndex][i].id].second = covfrag;
	}
      }

      hts_idx_destroy(myidx);
      sam_close(mysamfile);
    }
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
