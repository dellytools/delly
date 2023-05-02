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

  struct Geno {
    int32_t svid;
    int32_t sp;
    int32_t rp;
    
    Geno(int32_t const id, int32_t const s, int32_t const r) : svid(id), sp(s), rp(r) {}
  };
  
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


  template<typename TAlign>
  inline void
  _trimAlignedSequences(TAlign const& align, StructuralVariantRecord& sv, std::vector<std::string>& refseq, std::vector<std::string>& altseq) {
    int32_t s = -1;
    int32_t e = -1;
    for(uint32_t j = 0; j<align.shape()[1]; ++j) {
      if (align[0][j] != '-') {
	if (align[1][j] != '-') {
	  if (s == -1) s = j;
	  e = j + 1;
	}
      }
    }
    std::string s0;
    std::string s1;
    for(int32_t j = s; j < e; ++j) {
      if (align[0][j] != '-') s0.push_back(align[0][j]);
      if (align[1][j] != '-') s1.push_back(align[1][j]);
    }
    int32_t minsize = std::min(s0.size(), s1.size());
    if (sv.svt == 8) {
      // Breakpoint is at the end
      altseq[sv.id] = s0.substr(s0.size() - minsize, minsize);
      refseq[sv.id] = s1.substr(s1.size() - minsize, minsize);
    } else {
      altseq[sv.id] = s0.substr(0, minsize);
      refseq[sv.id] = s1.substr(0, minsize);
    }
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

  template<typename TConfig, typename TSVs>
  inline void
  _generateProbes(TConfig const& c, bam_hdr_t* hdr, TSVs& svs, std::vector<std::string>& refseq, std::vector<std::string>& altseq) {
    typedef boost::multi_array<char, 2> TAlign;
    
    // Preprocess REF and ALT
    boost::posix_time::ptime noww = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(noww) << "] " << "Generate REF and ALT probes" << std::endl;
    
    std::vector<std::string> refProbes(svs.size());
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      char* seq = NULL;

      // Iterate all structural variants
      for(typename TSVs::iterator itSV = svs.begin(); itSV != svs.end(); ++itSV) {
	if ((itSV->chr != refIndex) && (itSV->chr2 != refIndex)) continue;
	
	// Lazy loading of reference sequence
	if (seq == NULL) {
	  int32_t seqlen = -1;
	  std::string tname(hdr->target_name[refIndex]);
	  seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
	}

	// Set tag alleles
	if (itSV->chr == refIndex) {
	  itSV->alleles = _addAlleles(boost::to_upper_copy(std::string(seq + itSV->svStart - 1, seq + itSV->svStart)), std::string(hdr->target_name[itSV->chr2]), *itSV, itSV->svt);
	}
	if (!itSV->precise) continue;

	// Get the reference sequence
	if ((itSV->chr != itSV->chr2) && (itSV->chr2 == refIndex)) {
	  Breakpoint bp(*itSV);
	  _initBreakpoint(hdr, bp, (int32_t) itSV->consensus.size(), itSV->svt);
	  refProbes[itSV->id] = _getSVRef(c, seq, bp, refIndex, itSV->svt);
	}
	if (itSV->chr == refIndex) {
	  Breakpoint bp(*itSV);
	  if (_translocation(itSV->svt)) bp.part1 = refProbes[itSV->id];
	  if (itSV->svt ==4) {
	    int32_t bufferSpace = std::max((int32_t) ((itSV->consensus.size() - itSV->insLen) / 3), c.minimumFlankSize);
	    _initBreakpoint(hdr, bp, bufferSpace, itSV->svt);
	  } else _initBreakpoint(hdr, bp, (int32_t) itSV->consensus.size(), itSV->svt);
	  std::string svRefStr = _getSVRef(c, seq, bp, refIndex, itSV->svt);

	  // Check for Ns
	  if (nContent(svRefStr)) continue;
	  
	  // Find breakpoint to reference
	  TAlign align;
	  if (!_consRefAlignment(itSV->consensus, svRefStr, align, itSV->svt)) continue;

	  AlignDescriptor ad;
	  if (!_findSplit(c, itSV->consensus, svRefStr, align, ad, itSV->svt)) continue;

	  // Get exact alleles for INS and DEL
	  if (itSV->svEnd - itSV->svStart <= c.indelsize) {
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

	  // Generate genotyping probes
	  _trimAlignedSequences(align, *itSV, refseq, altseq);
	}
      }
      if (seq != NULL) free(seq);
    }
    // Clean-up
    fai_destroy(fai);
  }


  template<typename TConfig, typename TSVs, typename TGenoMap, typename TReadCountMap>
  inline void
  _fetchReads(TConfig const& c, std::vector<samFile*>& samfile, std::vector<hts_idx_t*>& idx, std::vector<bam_hdr_t*>& hdr, int32_t const file_c, TSVs& svs, TGenoMap& genoMap, TReadCountMap& covMap, std::vector<bool> const& probe) {
    typedef typename TGenoMap::mapped_type TGeno;
    
    boost::posix_time::ptime noww = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(noww) << "] " << "Select reads" << std::endl;
    
    // Iterate chromosomes
    for(int32_t refIndex=0; refIndex < (int32_t) hdr[file_c]->n_targets; ++refIndex) {
      
      // Fetch breakpoints
      std::vector<int32_t> breakpoint(hdr[file_c]->target_len[refIndex], -1);
      for(typename TSVs::iterator itSV = svs.begin(); itSV != svs.end(); ++itSV) {
	if (itSV->chr != refIndex) continue;
	if (!probe[itSV->id]) continue;
	int32_t pos = itSV->svStart;
	if (itSV->svt == 3) pos = itSV->svEnd; // for duplications the end region is first in the refprobe
	while ((breakpoint[pos] != -1) && (pos + 1 < (int) breakpoint.size())) ++pos;
	breakpoint[pos] = itSV->id;
      }	
      
      // Coverage track
      typedef uint16_t TMaxCoverage;
      uint32_t maxCoverage = std::numeric_limits<TMaxCoverage>::max();
      typedef std::vector<TMaxCoverage> TBpCoverage;
      TBpCoverage covBases(hdr[file_c]->target_len[refIndex], 0);

      // Parse reads
      hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr[file_c]->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
	if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;

	std::size_t seed = hash_lr(rec);
	int32_t readlen = readLength(rec);
	
	// Reference and sequence pointer
	uint32_t rp = rec->core.pos; // reference pointer
	uint32_t sp = 0; // sequence pointer
	  
	// Parse the CIGAR
	std::set<int32_t> processed;
	uint32_t* cigar = bam_get_cigar(rec);
	for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	  if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	    int32_t rpold = rp;
	    for(uint32_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k) {
	      if ((rp < hdr[file_c]->target_len[refIndex]) && (covBases[rp] < maxCoverage - 1)) ++covBases[rp];
	      ++rp;
	      ++sp;
	    }
	    for(uint32_t k = std::max((int) (rpold - c.minRefSep), 0); (k < rp + c.minRefSep) && (k < breakpoint.size()); ++k) {
	      if ((breakpoint[k] != -1) && (processed.find(breakpoint[k]) == processed.end())) {
		// Read long enough?
		if ((int) sp >= c.minimumFlankSize) {
		  if (readlen >= c.minimumFlankSize + (int) sp) {
		    if (genoMap.find(seed) == genoMap.end()) genoMap.insert(std::make_pair(seed, TGeno()));
		    if (rec->core.flag & BAM_FREVERSE) genoMap[seed].push_back(Geno(breakpoint[k], (readlen - sp), rp));
		    else genoMap[seed].push_back(Geno(breakpoint[k], sp, rp));
		    processed.insert(breakpoint[k]);
		  }
		}
	      }
	    }
	  } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	    rp += bam_cigar_oplen(cigar[i]);
	  } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	    sp += bam_cigar_oplen(cigar[i]);
	  } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	    rp += bam_cigar_oplen(cigar[i]);
	  } else if ((bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) {
	    sp += bam_cigar_oplen(cigar[i]);
	  } else {
	    std::cerr << "Unknown Cigar options" << std::endl;
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
  
  
  template<typename TConfig, typename TJunctionMap, typename TReadCountMap>
  inline void
  genotypeLR(TConfig& c, std::vector<StructuralVariantRecord>& svs, TJunctionMap& jctMap, TReadCountMap& covMap) {
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
    
    // Generate probes
    std::vector<std::string> refseq(svs.size());
    std::vector<std::string> altseq(svs.size());
    _generateProbes(c, hdr[0], svs, refseq, altseq);
    std::vector<bool> probeSuccess(svs.size(), true);
    for(uint32_t i = 0; i < svs.size(); ++i) {
      if ((refseq[i].empty()) || (altseq[i].empty())) probeSuccess[i] = false;
    }

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
      dumpOut << "#svid\tbam\tqname\tchr\tpos\tsubsequence\tmapq\ttype" << std::endl;
    }

#pragma omp parallel for default(shared)    
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {

      // Identify reads required for genotyping
      std::map<std::size_t, std::vector<Geno> > genoMap;
      _fetchReads(c, samfile, idx, hdr, file_c, svs, genoMap, covMap, probeSuccess);

      // Genotype SVs
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Align to REF and ALT" << std::endl;

      for(int32_t refIndex=0; refIndex < (int32_t) hdr[file_c]->n_targets; ++refIndex) {
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr[file_c]->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  // Only primary alignments for full sequence
	  if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) continue;
	  std::size_t seed = hash_lr(rec);
	  if (genoMap.find(seed) != genoMap.end()) {
	    std::string sequence;
	    for(uint32_t k = 0; k < genoMap[seed].size(); ++k) {
	      int32_t svid = genoMap[seed][k].svid;
	      if ((refseq[svid].empty()) || (altseq[svid].empty())) continue;
	      int32_t probelen = std::max(refseq[svid].size(), altseq[svid].size());
	      if (genoMap[seed][k].sp < probelen) continue;
	      
	      // Lazy-loading of sequence
	      if (sequence.empty()) {
		sequence.resize(rec->core.l_qseq);
		uint8_t* seqptr = bam_get_seq(rec);
		for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
		if (rec->core.flag & BAM_FREVERSE) reverseComplement(sequence);
	      }
	      if ((int) sequence.size() < genoMap[seed][k].sp + probelen) continue;

	      //std::cerr << svs[svid].svStart << "-" << svs[svid].svEnd << "," << svs[svid].svt << std::endl;
	      //std::cerr << seed << "," << genoMap[seed][k].sp << "," << genoMap[seed][k].rp << std::endl;
	      //std::cerr << sequence.size() << ',' << probelen << std::endl;
	      std::string subseq = sequence.substr(genoMap[seed][k].sp - probelen, 2 * probelen);	      
	      if ((jctMap[file_c][svid].ref.size() + jctMap[file_c][svid].alt.size()) >= c.maxGenoReadCount) continue;
	    
	      // Reverse complement probes
	      std::string revref = refseq[svid];
	      std::string revalt = altseq[svid];
	      reverseComplement(revref);
	      reverseComplement(revalt);
	      
	      // Edit distances
	      int32_t altFwd = _editDistanceHW(altseq[svid], subseq);
	      int32_t altRev = _editDistanceHW(revalt, subseq);
	      int32_t refFwd = _editDistanceHW(refseq[svid], subseq);
	      int32_t refRev = _editDistanceHW(revref, subseq);
	      double scoreAlt = (1.0 - c.flankQuality) * altseq[svid].size();
	      double scoreRef = (1.0 - c.flankQuality) * refseq[svid].size();
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
			dumpOut << svidStr << "\t" << c.files[file_c].string() << "\t" << bam_get_qname(rec) << "\t" << hdr[file_c]->target_name[rec->core.tid] << "\t" << rec->core.pos << "\t" << subseq << "\t" << (int32_t) rec->core.qual << "\tSR" << std::endl;
		      }
		      jctMap[file_c][svid].alt.push_back((uint8_t) std::min(aq, (uint32_t) rec->core.qual));
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
