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
    int32_t svStart;
    int32_t svEnd;
    int32_t svt;
    std::string refStart;
    std::string altStart;
    std::string refEnd;
    std::string altEnd;

    Geno() : svStart(-1), svEnd(-1), svt(-1) {}
  };


  inline void
  printAlignment(std::string const& query, std::string const& target, EdlibAlignMode const modeCode, EdlibAlignResult& align) {
    int32_t tIdx = -1;
    int32_t qIdx = -1;
    if (modeCode == EDLIB_MODE_HW) {
        tIdx = align.endLocations[0];
        for (int32_t i = 0; i < align.alignmentLength; i++) {
            if (align.alignment[i] != EDLIB_EDOP_INSERT) tIdx--;
        }
    }
    std::cerr << std::endl;
    for (int start = 0; start < align.alignmentLength; start += 50) {
      std::cerr << "T: ";
      int32_t startTIdx = -1;
      for (int32_t j = start; ((j < start + 50) && (j < align.alignmentLength)); ++j) {
	if (align.alignment[j] == EDLIB_EDOP_INSERT) std::cerr << "-";
	else std::cerr << target[++tIdx];
	if (j == start) startTIdx = tIdx;
      }
      std::cerr << " (" << std::max(startTIdx, 0) << " - " << tIdx << ")" << std::endl;

      // match / mismatch
      std::cerr << ("   ");
      for (int32_t j = start; j < start + 50 && j < align.alignmentLength; j++) {
	if (align.alignment[j] == EDLIB_EDOP_MATCH) std::cerr <<  "|";
	else std::cerr << " ";
      }
      std::cerr << std::endl;

      // query
      std::cerr << "Q: ";
      int32_t startQIdx = qIdx;
      for (int32_t j = start; j < start + 50 && j < align.alignmentLength; j++) {
	if (align.alignment[j] == EDLIB_EDOP_DELETE) std::cerr << "-";
	else std::cerr << query[++qIdx];
	if (j == start) startQIdx = qIdx;
      }
      std::cerr << " ("<< std::max(startQIdx, 0) << " - " << qIdx << ")" << std::endl;
      std::cerr << std::endl;
    }
  }
  

  inline float
  percentIdentity(std::string const& s1, std::string const& s2, int32_t splitpos, int32_t window) {
    // Get window boundaries
    int32_t ws = std::max(splitpos - window, 0);
    int32_t we = std::min(splitpos + window, (int32_t) s1.size());
    // Find percent identity
    bool varSeen = false;
    bool refSeen = false;
    int32_t refpos = 0;
    uint32_t gapMM = 0;
    uint32_t mm = 0;
    uint32_t ma = 0;
    float leftPerc = -1;
    float rightPerc = -1;
    bool inGap=false;
    for(uint32_t j = 0; j < s1.size(); ++j) {
      if (s2[j] != '-') varSeen = true;
      if (s1[j] != '-') {
	refSeen = true;
	if ((refpos == splitpos) || (refpos == ws) || (refpos == we)) {
	  if (refpos == splitpos) {
	    leftPerc = 0;
	    if (ma + mm > 0) leftPerc = (float) ma / (float) (ma + mm);
	  }
	  if (refpos == we) {
	    rightPerc = 0;
	    if (ma + mm > 0) rightPerc = (float) ma / (float) (ma + mm);
	  }
	  mm = 0;
	  ma = 0;
	  gapMM = 0;
	}
	++refpos;
      }
      if ((refSeen) && (varSeen)) {
	// Internal gap?
	if ((s2[j] == '-') || (s1[j] == '-')) {
	  if (!inGap) {
	    inGap = true;
	    gapMM = 0;
	  }
	  gapMM += 1;
	} else {
	  if (inGap) {
	    mm += gapMM;
	    inGap=false;
	  }
	  if (s2[j] == s1[j]) ma += 1;
	  else mm += 1;
	}
      }
    }
    if (rightPerc == -1) {
      rightPerc = 0;
      if (ma + mm > 0) rightPerc = (float) ma / (float) (ma + mm);
    }
    //std::cerr << ws << ',' << splitpos << ',' << we << ',' << leftPerc << ',' << rightPerc << std::endl;
    return std::min(leftPerc, rightPerc); 
  }


  template<typename TConfig, typename TAlign>
  inline bool
  _trimAlignedSequences(TConfig const& c, TAlign const& align, AlignDescriptor const& ad, Geno& gbp) {
    int32_t s = -1;
    int32_t e = -1;
    int32_t cpos = 0;
    int32_t rpos = 0;
    int32_t startConsPos = 0;
    int32_t startRefPos = 0;
    int32_t endConsPos = 0;
    int32_t endRefPos = 0;
    for(uint32_t j = 0; j<align.shape()[1]; ++j) {
      if (align[0][j] != '-') {
	++cpos;
	if (align[1][j] != '-') {
	  ++rpos;
	  if (s == -1) s = j;
	  e = j + 1;
	}
	if (cpos == ad.cStart) {
	  startConsPos = cpos;
	  startRefPos = rpos;
	}
	if (cpos == ad.cEnd) {
	  endConsPos = cpos;
	  endRefPos = rpos;
	}
      }
    }
    std::string s0;
    std::string s1;
    for(int32_t j = s; j < e; ++j) {
      if (align[0][j] != '-') s0.push_back(align[0][j]);
      if (align[1][j] != '-') s1.push_back(align[1][j]);
    }
    if (gbp.svStart != -1) {
      if (startConsPos < c.minimumFlankSize) return false;
      if ((int) (s0.size()) < startConsPos + c.minimumFlankSize) return false;
      if (startRefPos < c.minimumFlankSize) return false;
      if ((int) (s1.size()) < startRefPos + c.minimumFlankSize) return false;
      gbp.altStart = s0.substr(startConsPos - c.minimumFlankSize, 2 * c.minimumFlankSize);
      gbp.refStart = s1.substr(startRefPos - c.minimumFlankSize, 2 * c.minimumFlankSize);
    }
    if (gbp.svEnd != -1) {
      if (endConsPos < c.minimumFlankSize) return false;
      if ((int) (s0.size()) < endConsPos + c.minimumFlankSize) return false;
      if (endRefPos < c.minimumFlankSize) return false;
      if ((int) (s1.size()) < endRefPos + c.minimumFlankSize) return false;
      gbp.altEnd = s0.substr(endConsPos - c.minimumFlankSize, 2 * c.minimumFlankSize);
      gbp.refEnd = s1.substr(endRefPos - c.minimumFlankSize, 2 * c.minimumFlankSize);
    }
    return true;
  }

  inline int32_t
  _editDistanceHW(std::string const& query, std::string const& target) {
    EdlibAlignResult align = edlibAlign(query.c_str(), query.size(), target.c_str(), target.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
    // Debug: requires EDLIB_TASK_PATH otherwise EDLIB_TASK_DISTANCE
    //printAlignment(query, target, EDLIB_MODE_HW, align);
    int32_t ed = align.editDistance;
    edlibFreeAlignResult(align);
    return ed;
  }

  
  template<typename TConfig, typename TJunctionMap, typename TReadCountMap>
  inline void
  genotypeLR(TConfig& c, std::vector<StructuralVariantRecord>& svs, TJunctionMap& jctMap, TReadCountMap& covMap) {
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
    for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
      hdr[file_c] = sam_hdr_read(samfile[file_c]);
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

    // Coverage distribution
    typedef uint16_t TMaxCoverage;
    uint32_t maxCoverage = std::numeric_limits<TMaxCoverage>::max();

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
	  //std::cerr << "svid:" << itSV->id << "," << itSV->svStart << "," << itSV->svEnd << "," << itSV->svt << ",consensus-to-reference-alignment" << std::endl;
	  //for(uint32_t i = 0; i<align.shape()[0]; ++i) {
	  //if (i == 0) {
	  //int32_t cpos = 0;
	  //for(uint32_t j = 0; j<align.shape()[1]; ++j) {
	  //if (align[i][j] != '-') ++cpos;
	  //if (cpos == ad.cStart) std::cerr << '|';
	  //else if (cpos == ad.cEnd) std::cerr << '|';
	  //else std::cerr << '#';
	  //}
	  //std::cerr << std::endl;
	  //}
	  //for(uint32_t j = 0; j<align.shape()[1]; ++j) std::cerr << align[i][j];
	  //std::cerr << std::endl;
	  //}
	  //std::cerr << std::endl;

	  // Trim aligned sequences
	  gbp[itSV->id].svStart = itSV->svStart;
	  if (itSV->chr2 == refIndex) gbp[itSV->id].svEnd = itSV->svEnd;
	  gbp[itSV->id].svt = itSV->svt;
	  if (!_trimAlignedSequences(c, align, ad, gbp[itSV->id])) continue;
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

	      // Read long enough?
	      if (spHit < c.minimumFlankSize) continue;
	      if (readlen < c.minimumFlankSize + spHit) continue;
	      
	      //std::cerr << ">" << svid << ',' << gbp[svid].svStart << ',' << gbp[svid].svEnd << ',' << gbp[svid].svt << std::endl;
	      std::string refseq;
	      std::string altseq;
	      if (rpHit == gbp[svid].svStart) {
		refseq = gbp[svid].refStart;
		altseq = gbp[svid].altStart;
	      } else if (rpHit == gbp[svid].svEnd) {
		refseq = gbp[svid].refEnd;
		altseq = gbp[svid].altEnd;
	      } else {
		std::cerr << "Error: Unclear reference hit!" << std::endl;
		exit(-1);
	      }

	      // Reverse complement probes
	      std::string revref = refseq;
	      std::string revalt = altseq;
	      reverseComplement(revref);
	      reverseComplement(revalt);

	      // Edit distances
	      int32_t altFwd = _editDistanceHW(altseq, sequence);
	      int32_t altRev = _editDistanceHW(revalt, sequence);
	      int32_t refFwd = _editDistanceHW(refseq, sequence);
	      int32_t refRev = _editDistanceHW(revref, sequence);
	      double scoreAlt = (1.0 - c.flankQuality) * altseq.size();
	      double scoreRef = (1.0 - c.flankQuality) * refseq.size();
	      if (std::min(altFwd, refFwd) < std::min(altRev, refRev)) {
		scoreAlt = scoreAlt / (double) altFwd;
		scoreRef = scoreRef / (double) refFwd;
	      } else {
		scoreAlt = scoreAlt / (double) altRev;
		scoreRef = scoreRef / (double) refRev;
	      }

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

}

#endif
