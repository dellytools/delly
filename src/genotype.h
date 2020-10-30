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
#include <boost/progress.hpp>

#include <htslib/sam.h>

#include "util.h"

namespace torali
{

  struct Geno {
    int32_t svStartPrefix;
    int32_t svStartSuffix;
    int32_t svEndPrefix;
    int32_t svEndSuffix;
    int32_t svStart;
    int32_t svEnd;
    int32_t svt;
    std::string ref;
    std::string alt;

    Geno() : svStartPrefix(-1), svStartSuffix(-1), svEndPrefix(-1), svEndSuffix(-1), svStart(-1), svEnd(-1), svt(-1) {}
  };
  

  template<typename TConfig, typename TSRStore, typename TJunctionMap, typename TReadCountMap>
  inline void
  trackRef(TConfig& c, std::vector<StructuralVariantRecord>& svs, TSRStore& srStore, TJunctionMap& jctMap, TReadCountMap& covMap) {

    typedef uint16_t TMaxCoverage;
    uint32_t maxCoverage = std::numeric_limits<TMaxCoverage>::max();
    
    typedef std::vector<StructuralVariantRecord> TSVs;
    typedef std::vector<uint8_t> TQuality;
    typedef boost::multi_array<char, 2> TAlign;
    if (svs.empty()) return;

    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    typedef std::vector<bam_hdr_t*> THeader;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    THeader hdr(c.files.size());
    int32_t totalTarget = 0;
    for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
      hdr[file_c] = sam_hdr_read(samfile[file_c]);
      totalTarget += hdr[file_c]->n_targets;
    }

    // Parse genome chr-by-chr
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "SV annotation" << std::endl;
    boost::progress_display show_progress( hdr[0]->n_targets );

    // Ref aligned reads
    typedef std::vector<uint32_t> TRefAlignCount;
    typedef std::vector<TRefAlignCount> TFileRefAlignCount;
    TFileRefAlignCount refAlignedReadCount(c.files.size(), TRefAlignCount());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) refAlignedReadCount[file_c].resize(svs.size(), 0);

    // Dump file
    boost::iostreams::filtering_ostream dumpOut;
    if (c.hasDumpFile) {
      dumpOut.push(boost::iostreams::gzip_compressor());
      dumpOut.push(boost::iostreams::file_sink(c.dumpfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
      dumpOut << "#svid\tbam\tqname\tchr\tpos\tmatechr\tmatepos\tmapq\ttype" << std::endl;
    }

    // Iterate chromosomes
    std::vector<std::string> refProbes(svs.size());
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < (int32_t) hdr[0]->n_targets; ++refIndex) {
      ++show_progress;
      char* seq = NULL;

      // Reference and consensus probes for this chromosome
      typedef std::vector<Geno> TGenoRegion;
      TGenoRegion gbp(svs.size(), Geno());
      
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
	  if (_translocation(itSV->svt)) bp.part1 = refProbes[itSV->id];
	  if (itSV->svt ==4) {
	    int32_t bufferSpace = std::max((int32_t) ((itSV->consensus.size() - itSV->insLen) / 3), c.minimumFlankSize);
	    _initBreakpoint(hdr[0], bp, bufferSpace, itSV->svt);
	  } else _initBreakpoint(hdr[0], bp, (int32_t) itSV->consensus.size(), itSV->svt);
	  std::string svRefStr = _getSVRef(seq, bp, refIndex, itSV->svt);
	  
	  // Find breakpoint to reference
	  TAlign align;
	  if (!_consRefAlignment(itSV->consensus, svRefStr, align, itSV->svt)) continue;

	  AlignDescriptor ad;
	  if (!_findSplit(c, itSV->consensus, svRefStr, align, ad, itSV->svt)) continue;

	  // Get exact alleles for INS and DEL
	  if ((itSV->svt == 2) || (itSV->svt == 4)) {
	    std::string refVCF;
	    std::string altVCF;
	    int32_t cpos = 0;
	    bool inSV = false;
	    for(uint32_t j = 0; j<align.shape()[1]; ++j) {
	      if (align[0][j] != '-') {
		++cpos;
		if (cpos == ad.cStart) inSV = true;
	  	else if (cpos == ad.cEnd) inSV = false;
	      }
	      if (inSV) {
		if (align[0][j] != '-') altVCF += align[0][j];
		if (align[1][j] != '-') refVCF += align[1][j];
	      }
	    }
	    itSV->alleles = _addAlleles(refVCF, altVCF);
	  }
	  
	  // Debug consensus to reference alignment
	  //std::cerr << "svid:" << itSV->id << ",consensus-to-reference-alignment" << std::endl;
	  //for(uint32_t i = 0; i<align.shape()[0]; ++i) {
	  //if (i == 0) {
	  //  int32_t cpos = 0;
	  //  for(uint32_t j = 0; j<align.shape()[1]; ++j) {
	  //	if (align[i][j] != '-') ++cpos;
	  //	if (cpos == ad.cStart) std::cerr << '|';
	  //	else if (cpos == ad.cEnd) std::cerr << '|';
	  //	else std::cerr << '#';
	  //  }
	  //  std::cerr << std::endl;
	  //}
	  //for(uint32_t j = 0; j<align.shape()[1]; ++j) std::cerr << align[i][j];
	  //std::cerr << std::endl;
	  //}
	  //std::cerr << std::endl;

	  // Trim aligned sequences
	  std::string altSeq;
	  std::string refSeq;
	  int32_t leadCrop = _trimAlignedSequences(align, altSeq, refSeq);

	  // Allele-tagging probes
	  gbp[itSV->id].svStartPrefix = std::max(ad.cStart - leadCrop, 0);
	  gbp[itSV->id].svStartSuffix = std::max((int32_t) altSeq.size() - gbp[itSV->id].svStartPrefix, 0);
	  gbp[itSV->id].svStart = itSV->svStart;
	  if (itSV->chr2 == refIndex) {
	    gbp[itSV->id].svEndPrefix = std::max(ad.cEnd - leadCrop, 0);
	    gbp[itSV->id].svEndSuffix = std::max((int32_t) altSeq.size() - gbp[itSV->id].svEndPrefix, 0);
	    gbp[itSV->id].svEnd = itSV->svEnd;
	  }
	  gbp[itSV->id].ref = refSeq;
	  gbp[itSV->id].alt = altSeq;
	  gbp[itSV->id].svt = itSV->svt;
	}
      }
      if (seq != NULL) free(seq);

      // Genotype
      // Iterate samples
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	// Check we have mapped reads on this chromosome
	bool nodata = true;
	std::string suffix("cram");
	std::string str(c.files[file_c].string());
	if ((str.size() >= suffix.size()) && (str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0)) nodata = false;
	uint64_t mapped = 0;
	uint64_t unmapped = 0;
	hts_idx_get_stat(idx[file_c], refIndex, &mapped, &unmapped);
	if (mapped) nodata = false;
	if (nodata) continue;

	// Coverage track
	typedef std::vector<TMaxCoverage> TBpCoverage;
	TBpCoverage covBases(hdr[file_c]->target_len[refIndex], 0);

	// Flag breakpoints
	typedef std::set<int32_t> TIdSet;
	typedef std::map<uint32_t, TIdSet> TBpToIdMap;
	TBpToIdMap bpid;
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet bpOccupied(hdr[file_c]->target_len[refIndex], false);
	for(uint32_t i = 0; i < gbp.size(); ++i) {
	  if (gbp[i].svStart != -1) {
	    bpOccupied[gbp[i].svStart] = 1;
	    if (bpid.find(gbp[i].svStart) == bpid.end()) bpid.insert(std::make_pair(gbp[i].svStart, TIdSet()));
	    bpid[gbp[i].svStart].insert(i);
	  }
	  if (gbp[i].svEnd != -1) {
	    bpOccupied[gbp[i].svEnd] = 1;
	    if (bpid.find(gbp[i].svEnd) == bpid.end()) bpid.insert(std::make_pair(gbp[i].svEnd, TIdSet()));
	    bpid[gbp[i].svEnd].insert(i);
	  }
	}

	// Count reads
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr[file_c]->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  // Genotyping only primary alignments
	  if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	  
	  // Read length
	  int32_t readlen = readLength(rec);

	  // Reference and sequence pointer
	  uint32_t rp = rec->core.pos; // reference pointer
	  uint32_t sp = 0; // sequence pointer

	  // All SV hits
	  typedef std::pair<int32_t, int32_t> TRefSeq;
	  typedef std::map<int32_t, TRefSeq> TSVSeqHit;
	  TSVSeqHit genoMap;

	  /*
	  // Any direct SV support
	  std::size_t seed = hash_lr(rec);
	  if (srStore.find(seed) != srStore.end()) {
	    for(uint32_t ri = 0; ri < srStore[seed].size(); ++ri) {
	      int32_t svid = srStore[seed][ri].svid;
	      genoMap.insert(std::make_pair(srStore[seed][ri].svid, std::make_pair(gbp[svid].svStart, srStore[seed][ri].sstart)));
	    }
	  }
	  */
	  
	  // Parse the CIGAR
	  uint32_t* cigar = bam_get_cigar(rec);
	  for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	    if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	      // Fetch reference alignments
	      for(uint32_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k) {
		if ((rp < hdr[file_c]->target_len[refIndex]) && (covBases[rp] < maxCoverage - 1)) ++covBases[rp];
		if (bpOccupied[rp]) {
		  for(typename TIdSet::const_iterator it = bpid[rp].begin(); it != bpid[rp].end(); ++it) {
		    // Ensure fwd alignment and each SV only once
		    if (genoMap.find(*it) == genoMap.end()) {
		      if (rec->core.flag & BAM_FREVERSE) genoMap.insert(std::make_pair(*it, std::make_pair(rp, readlen - sp)));
		      else genoMap.insert(std::make_pair(*it, std::make_pair(rp, sp)));
		    }
		  }
		}
		++sp;
		++rp;
	      }
	    } else if ((bam_cigar_op(cigar[i]) == BAM_CDEL) || (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP)) {
	      for(uint32_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k) {
		if (bpOccupied[rp]) {
		  for(typename TIdSet::const_iterator it = bpid[rp].begin(); it != bpid[rp].end(); ++it) {
		    // Ensure fwd alignment and each SV only once
		    if (genoMap.find(*it) == genoMap.end()) {
		      if (rec->core.flag & BAM_FREVERSE) genoMap.insert(std::make_pair(*it, std::make_pair(rp, readlen - sp)));
		      else genoMap.insert(std::make_pair(*it, std::make_pair(rp, sp)));
		    }
		  }
		}
		++rp;
	      }
	    } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	      sp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	      sp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	      // Do nothing
	    } else {
	      std::cerr << "Unknown Cigar options" << std::endl;
	    }
	  }

	  // Read for genotyping?
	  if (!genoMap.empty()) {
	    // Get sequence
	    std::string sequence;
	    sequence.resize(rec->core.l_qseq);
	    uint8_t* seqptr = bam_get_seq(rec);
	    for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	    // Genotype all SVs covered by this read
	    for(typename TSVSeqHit::iterator git = genoMap.begin(); git != genoMap.end(); ++git) {
	      int32_t svid = git->first;
	      uint32_t maxGenoReadCount = 500;
	      if ((jctMap[file_c][svid].ref.size() + jctMap[file_c][svid].alt.size()) >= maxGenoReadCount) continue;
	      
	      int32_t rpHit = git->second.first;
	      int32_t spHit = git->second.second;

	      // Require spanning reads
	      std::string subseq;
	      if (rpHit == gbp[svid].svStart) {
		if (rec->core.flag & BAM_FREVERSE) {
		  if (spHit < gbp[svid].svStartSuffix) continue;
		  if (readlen < gbp[svid].svStartPrefix + spHit) continue;
		  int32_t st = std::max((readlen - spHit) - gbp[svid].svStartPrefix - c.minimumFlankSize, 0);
		  subseq = sequence.substr(st, gbp[svid].svStartPrefix + gbp[svid].svStartSuffix + 2 * c.minimumFlankSize);
		} else {
		  if (spHit < gbp[svid].svStartPrefix) continue;
		  if (readlen < gbp[svid].svStartSuffix + spHit) continue;
		  int32_t st = std::max(spHit - gbp[svid].svStartPrefix - c.minimumFlankSize, 0);
		  subseq = sequence.substr(st, gbp[svid].svStartPrefix + gbp[svid].svStartSuffix + 2 * c.minimumFlankSize);
		}
	      } else {
		if (rec->core.flag & BAM_FREVERSE) {
		  if (spHit < gbp[svid].svEndSuffix) continue;
		  if (readlen < gbp[svid].svEndPrefix + spHit) continue;
		  int32_t st = std::max((readlen - spHit) - gbp[svid].svEndPrefix - c.minimumFlankSize, 0);
		  subseq = sequence.substr(st, gbp[svid].svEndPrefix + gbp[svid].svEndSuffix + 2 * c.minimumFlankSize);
		} else {
		  if (spHit < gbp[svid].svEndPrefix) continue;
		  if (readlen < gbp[svid].svEndSuffix + spHit) continue;
		  int32_t st = std::max(spHit - gbp[svid].svEndPrefix - c.minimumFlankSize, 0);
		  subseq = sequence.substr(st, gbp[svid].svEndPrefix + gbp[svid].svEndSuffix + 2 * c.minimumFlankSize);
		}
	      }
	    
	      // Compute alignment to alternative haplotype
	      DnaScore<int> simple(c.aliscore.match, c.aliscore.mismatch, c.aliscore.mismatch, c.aliscore.mismatch);
	      AlignConfig<true, false> semiglobal;
	      double scoreAlt = needleBanded(gbp[svid].alt, subseq, semiglobal, simple);
	      scoreAlt /= (double) (c.flankQuality * gbp[svid].alt.size() * simple.match + (1.0 - c.flankQuality) * gbp[svid].alt.size() * simple.mismatch);
	    
	      // Compute alignment to reference haplotype
	      double scoreRef = needleBanded(gbp[svid].ref, subseq, semiglobal, simple);
	      scoreRef /= (double) (c.flankQuality * gbp[svid].ref.size() * simple.match + (1.0 - c.flankQuality) * gbp[svid].ref.size() * simple.mismatch);

	      // Any confident alignment?
	      if ((scoreRef > 1) || (scoreAlt > 1)) {
		if (scoreRef > scoreAlt) {
		  // Account for reference bias
		  if (++refAlignedReadCount[file_c][svid] % 2) {
		    TQuality quality;
		    quality.resize(rec->core.l_qseq);
		    uint8_t* qualptr = bam_get_qual(rec);
		    for (int i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
		    uint32_t rq = scoreRef * 35;
		    if (rq >= c.minGenoQual) {
		      uint8_t* hpptr = bam_aux_get(rec, "HP");
		      jctMap[file_c][svid].ref.push_back((uint8_t) std::min(rq, (uint32_t) rec->core.qual));
		      if (hpptr) {
			c.isHaplotagged = true;
			int hap = bam_aux2i(hpptr);
			if (hap == 1) ++jctMap[file_c][svid].refh1;
			else ++jctMap[file_c][svid].refh2;
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
		    uint8_t* hpptr = bam_aux_get(rec, "HP");
		    if (c.hasDumpFile) {
		      std::string svidStr(_addID(gbp[svid].svt));
		      std::string padNumber = boost::lexical_cast<std::string>(svid);
		      padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
		      svidStr += padNumber;
		      dumpOut << svidStr << "\t" << c.files[file_c].string() << "\t" << bam_get_qname(rec) << "\t" << hdr[file_c]->target_name[rec->core.tid] << "\t" << rec->core.pos << "\t" << hdr[file_c]->target_name[rec->core.mtid] << "\t" << rec->core.mpos << "\t" << (int32_t) rec->core.qual << "\tSR" << std::endl;
		    }
		    jctMap[file_c][svid].alt.push_back((uint8_t) std::min(aq, (uint32_t) rec->core.qual));
		    if (hpptr) {
		      c.isHaplotagged = true;
		      int hap = bam_aux2i(hpptr);
		      if (hap == 1) ++jctMap[file_c][svid].alth1;
		      else ++jctMap[file_c][svid].alth2;
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
      
	// Assign SV support
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
    }
    // Clean-up
    fai_destroy(fai);

    // Clean-up
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      bam_hdr_destroy(hdr[file_c]);	  
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }



  template<typename TConfig, typename TSRStore, typename TJunctionMap, typename TReadCountMap>
  inline void
  genotypeLR(TConfig& c, std::vector<StructuralVariantRecord>& svs, TSRStore& srStore, TJunctionMap& jctMap, TReadCountMap& covMap) {
    typedef std::vector<StructuralVariantRecord> TSVs;
    if (svs.empty()) return;
    
    typedef uint16_t TMaxCoverage;
    uint32_t maxCoverage = std::numeric_limits<TMaxCoverage>::max();
    
    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    typedef std::vector<bam_hdr_t*> THeader;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    THeader hdr(c.files.size());
    int32_t totalTarget = 0;
    for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
      hdr[file_c] = sam_hdr_read(samfile[file_c]);
      totalTarget += hdr[file_c]->n_targets;
    }

    // Parse genome chr-by-chr
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "SV annotation" << std::endl;
    boost::progress_display show_progress( hdr[0]->n_targets );

    // Ref aligned reads
    typedef std::vector<uint32_t> TRefAlignCount;
    typedef std::vector<TRefAlignCount> TFileRefAlignCount;
    TFileRefAlignCount refAlignedReadCount(c.files.size(), TRefAlignCount());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) refAlignedReadCount[file_c].resize(svs.size(), 0);

    // Dump file
    boost::iostreams::filtering_ostream dumpOut;
    if (c.hasDumpFile) {
      dumpOut.push(boost::iostreams::gzip_compressor());
      dumpOut.push(boost::iostreams::file_sink(c.dumpfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
      dumpOut << "#svid\tbam\tqname\tchr\tpos\tmatechr\tmatepos\tmapq\ttype" << std::endl;
    }


    faidx_t* fai = fai_load(c.genome.string().c_str());
    // Iterate chromosomes
    for(int32_t refIndex=0; refIndex < (int32_t) hdr[0]->n_targets; ++refIndex) {
      ++show_progress;
      char* seq = NULL;
      
      // Iterate samples
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	// Check we have mapped reads on this chromosome
	bool nodata = true;
	std::string suffix("cram");
	std::string str(c.files[file_c].string());
	if ((str.size() >= suffix.size()) && (str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0)) nodata = false;
	uint64_t mapped = 0;
	uint64_t unmapped = 0;
	hts_idx_get_stat(idx[file_c], refIndex, &mapped, &unmapped);
	if (mapped) nodata = false;
	if (nodata) continue;

	// Flag breakpoints
	typedef std::set<int32_t> TIdSet;
	typedef std::map<uint32_t, TIdSet> TBpToIdMap;
	TBpToIdMap bpid;
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet bpOccupied(hdr[file_c]->target_len[refIndex], false);
	for(typename TSVs::iterator itSV = svs.begin(); itSV != svs.end(); ++itSV) {
	  if (itSV->chr == refIndex) {
	    bpOccupied[itSV->svStart] = 1;
	    if (bpid.find(itSV->svStart) == bpid.end()) bpid.insert(std::make_pair(itSV->svStart, TIdSet()));
	    bpid[itSV->svStart].insert(itSV->id);
	  }
	  if (itSV->chr2 == refIndex) {
	    bpOccupied[itSV->svEnd] = 1;
	    if (bpid.find(itSV->svEnd) == bpid.end()) bpid.insert(std::make_pair(itSV->svEnd, TIdSet()));
	    bpid[itSV->svEnd].insert(itSV->id);
	  }
	}
	if (bpid.empty()) continue;

	// Lazy loading of reference sequence
	if (seq == NULL) {
	  int32_t seqlen = -1;
	  std::string tname(hdr[0]->target_name[refIndex]);
	  seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr[0]->target_len[refIndex], &seqlen);
	}    
	
	// Coverage track
	typedef std::vector<TMaxCoverage> TBpCoverage;
	TBpCoverage covBases(hdr[file_c]->target_len[refIndex], 0);

	// Count reads
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr[file_c]->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  // Genotyping only primary alignments
	  if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	  
	  // Read length
	  std::size_t seed = hash_lr(rec);

	  // Reference and sequence pointer
	  uint32_t rp = rec->core.pos; // reference pointer
	  uint32_t sp = 0; // sequence pointer

	  // Get sequence
	  std::string sequence;
	  sequence.resize(rec->core.l_qseq);
	  uint8_t* seqptr = bam_get_seq(rec);
	  for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	  
	  // Any REF support
	  std::string refAlign = "";
	  std::string altAlign = "";
	  std::vector<int32_t> hits;
	  uint32_t* cigar = bam_get_cigar(rec);
	  for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	    if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	      // Fetch reference alignments
	      for(uint32_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k) {
		if ((rp < hdr[file_c]->target_len[refIndex]) && (covBases[rp] < maxCoverage - 1)) ++covBases[rp];
		refAlign += seq[rp];
		altAlign += sequence[sp];
		if (bpOccupied[rp]) hits.push_back(rp);
		++sp;
		++rp;
	      }
	    } else if ((bam_cigar_op(cigar[i]) == BAM_CDEL) || (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP)) {
	      for(uint32_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k) {
		refAlign += seq[rp];
		altAlign += "-";
		if (bpOccupied[rp]) hits.push_back(rp);
		++rp;
	      }
	    } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	      for(uint32_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k) {
		refAlign += "-";
		altAlign += sequence[sp];
		++sp;
	      }
	    } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	      sp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	      // Do nothing
	    } else {
	      std::cerr << "Unknown Cigar options" << std::endl;
	    }
	  }

	  // Sufficiently long flank mapping?
	  if (refAlign.size() < (uint32_t) c.minimumFlankSize) continue;

	  // Confident mapping?
	  double score = _score(refAlign, altAlign, c.aliscore);
	  score /= (double) refAlign.size();
	  score /= c.aliscore.match;
	  score *= 10;
	  if (score < c.minGenoQual) continue;

	  // Any ALT support?
	  TIdSet altAssigned;
	  if (srStore.find(seed) != srStore.end()) {
	    for(uint32_t ri = 0; ri < srStore[seed].size(); ++ri) {
	      int32_t svid = srStore[seed][ri].svid;
	      if (svid == -1) continue;
	      altAssigned.insert(svid);
	      uint8_t* hpptr = bam_aux_get(rec, "HP");
	      if (c.hasDumpFile) {
		std::string svidStr(_addID(svs[svid].svt));
		std::string padNumber = boost::lexical_cast<std::string>(svid);
		padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
		svidStr += padNumber;
		dumpOut << svidStr << "\t" << c.files[file_c].string() << "\t" << bam_get_qname(rec) << "\t" << hdr[file_c]->target_name[rec->core.tid] << "\t" << rec->core.pos << "\t" << hdr[file_c]->target_name[rec->core.mtid] << "\t" << rec->core.mpos << "\t" << (int32_t) rec->core.qual << "\tSR" << std::endl;
	      }
	      jctMap[file_c][svid].alt.push_back((uint8_t) std::min((uint32_t) score, (uint32_t) rec->core.qual));
	      if (hpptr) {
		c.isHaplotagged = true;
		int hap = bam_aux2i(hpptr);
		if (hap == 1) ++jctMap[file_c][svid].alth1;
		else ++jctMap[file_c][svid].alth2;
	      }
	    }
	  }

	  // Any REF support
	  if (hits.empty()) continue;

	  /*
	  std::cerr << ">" << hdr[file_c]->target_len[refIndex] << ":" << rec->core.pos << "\t" << score << "\t" << seed << std::endl;
	  uint32_t refpos = rec->core.pos;
	  for(uint32_t i = 0; i < refAlign.size(); ++i) {
	    std::cerr << refpos << "\t" << refAlign[i] << "\t" << altAlign[i];
	    if (refAlign != "-") {
	      if (bpOccupied[refpos]) std::cerr << "\t<---";
	      ++refpos;
	    }
	    std::cerr << std::endl;
	  }
	  */
	  
	  // Iterate all spanned SVs
	  for(uint32_t idx = 0; idx < hits.size(); ++idx) {
	    for(typename TIdSet::const_iterator its = bpid[hits[idx]].begin(); its != bpid[hits[idx]].end(); ++its) {
	      int32_t svid = *its;
	      if (altAssigned.find(svid) != altAssigned.end()) continue; 
	      std::cerr << svs[svid].chr << ',' << svs[svid].svStart << ',' << svs[svid].chr2 << ',' << svs[svid].svEnd << std::endl;

	      // Account for reference bias
	      if (++refAlignedReadCount[file_c][svid] % 2) {
		uint8_t* hpptr = bam_aux_get(rec, "HP");
		jctMap[file_c][svid].ref.push_back((uint8_t) std::min((uint32_t) score, (uint32_t) rec->core.qual));
		if (hpptr) {
		  c.isHaplotagged = true;
		  int hap = bam_aux2i(hpptr);
		  if (hap == 1) ++jctMap[file_c][svid].refh1;
		  else ++jctMap[file_c][svid].refh2;
		}
	      }
	    }
	  }
	}
	// Clean-up
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      
	// Assign SV support
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
      if (seq != NULL) free(seq);
    }
    // Clean-up
    fai_destroy(fai);

    // Clean-up
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      bam_hdr_destroy(hdr[file_c]);	  
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }
     
  
}

#endif
