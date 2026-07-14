#ifndef ASMODE_H
#define ASMODE_H

#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <vector>
#include <cctype>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <math.h>
#include <stdio.h>

#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>

#include "edlib.h"
#include "delly.h"
#include "coverage.h"
#include "genotype.h"
#include "util.h"
#include "junction.h"
#include "cluster.h"
#include "assemble.h"
#include "modvcf.h"
#include "merge.h"

namespace torali {


  struct AsmConfig {
    bool hasVcfFile;
    uint16_t minMapQual;
    uint32_t minClip;
    uint32_t minCliqueSize;
    uint32_t graphPruning;
    uint32_t minRefSep;
    uint32_t maxReadSep;
    int32_t minConsWindow;
    int32_t minimumFlankSize;
    int32_t indelsize;
    int32_t nchr;
    float flankQuality;
    float meiMinFrac;
    float trMinFrac;
    bool diploid;
    uint32_t nsamples;
    std::set<int32_t> svtset;
    boost::filesystem::path outfile;
    std::vector<boost::filesystem::path> files;
    boost::filesystem::path genome;
    std::vector<std::string> sampleName;
    std::vector<uint32_t> fileSample;
    std::vector<uint8_t> fileHap;
  };

  // Strip trailing haplotype
  inline std::string
  _stripHaplotypeSuffix(std::string s) {
    static const std::vector<std::string> suffixes = { ".hap1", ".hap2", ".hapA", ".hapB", ".h1", ".h2", ".mat", ".pat", ".maternal", ".paternal", ".1", ".2" };
    for(uint32_t i = 0; i < suffixes.size(); ++i) {
      if ((s.size() > suffixes[i].size()) && (s.compare(s.size() - suffixes[i].size(), suffixes[i].size(), suffixes[i]) == 0)) {
	s = s.substr(0, s.size() - suffixes[i].size());
	break;
      }
    }
    return s;
  }

  inline std::string
  _commonSampleName(std::string const& a, std::string const& b) {
    std::size_t n = std::min(a.size(), b.size());
    std::size_t k = 0;
    while ((k < n) && (a[k] == b[k])) ++k;
    std::string pre = a.substr(0, k);
    while ((!pre.empty()) && ((pre.back() == '.') || (pre.back() == '_') || (pre.back() == '-'))) pre.pop_back();
    return (!pre.empty()) ? pre : _stripHaplotypeSuffix(a);
  }


  template<typename TConfig, typename TReadBp, typename TSvtSRBamRecord>
  inline void
  findAsmJunctions(TConfig const& c, TReadBp& readBp, boost::unordered_map<std::size_t, uint32_t>& readSample, boost::unordered_map<std::size_t, uint8_t>& readHap, TSvtSRBamRecord& srBR) {
    bool const doDel = ((c.svtset.empty()) || (c.svtset.find(2) != c.svtset.end()));
    bool const doIns = ((c.svtset.empty()) || (c.svtset.find(4) != c.svtset.end()));

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
    
    // Parse genome chr-by-chr
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Split-read scanning" << std::endl;

    // Iterate chromosomes
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {

      // Collect reads from all samples
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	// Read alignments
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {

	  // Keep secondary alignments
	  if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
	  if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;

	  // Fetch sample and haplotype
	  std::size_t seed = hash_lr(rec);
	  boost::hash_combine(seed, c.fileSample[file_c]);
	  readSample[seed] = c.fileSample[file_c];
	  readHap[seed] = c.fileHap[file_c];
	  //std::cerr << bam_get_qname(rec) << '\t' << seed << std::endl;
	  uint32_t rp = rec->core.pos; // reference pointer
	  uint32_t sp = 0; // sequence pointer
	  int32_t readStart = rec->core.pos;
	  if (rec->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) readStart = -1;
	  int32_t const seqlen = readLength(rec);
	  bool const rev = (rec->core.flag & BAM_FREVERSE);

	  // Parse the CIGAR
	  const uint32_t* cigar = bam_get_cigar(rec);
	  for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	    uint32_t const oplen = bam_cigar_oplen(cigar[i]);
	    if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	      sp += oplen;
	      rp += oplen;
	    } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	      // Intra-alignment deletion
	      if ((doDel) && (oplen > c.minRefSep) && ((int32_t) sp <= seqlen)) {
		int32_t ss = rev ? (seqlen - (int32_t) sp) : (int32_t) sp;
		srBR[2].push_back(SRBamRecord(rec->core.tid, rp, rec->core.tid, rp + oplen, readStart, ss, rec->core.qual, 0, seed));
	      }
	      rp += oplen;
	    } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	      // Intra-alignment insertion
	      if ((doIns) && (oplen > c.minRefSep) && ((int32_t) (sp + oplen) <= seqlen)) {
		int32_t ss = rev ? (seqlen - (int32_t) sp - (int32_t) oplen) : (int32_t) sp;
		if (ss < 0) ss = 0;
		srBR[4].push_back(SRBamRecord(rec->core.tid, rp, rec->core.tid, rp + 1, readStart, ss, rec->core.qual, oplen, seed));
	      }
	      sp += oplen;
	    } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	      rp += oplen;
	    } else if ((bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) {
	      int32_t finalsp = sp;
	      bool scleft = false;
	      if (sp == 0) {
		finalsp += bam_cigar_oplen(cigar[i]); // Leading soft-clip / hard-clip
		scleft = true;
	      }
	      sp += bam_cigar_oplen(cigar[i]);
	      if (bam_cigar_oplen(cigar[i]) > c.minClip) _insertJunction(readBp, seed, rec, rp, finalsp, scleft);
	    } else {
	      std::cerr << "Unknown Cigar options" << std::endl;
	    }
	  }
	}
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }
    }
    // Sort junctions
    for(typename TReadBp::iterator it = readBp.begin(); it != readBp.end(); ++it) std::sort(it->second.begin(), it->second.end());

    // Clean-up
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }
 

  template<typename TConfig, typename TSvtSRBamRecord>
  inline void
  _findAsmBreakpoints(TConfig const& c, TSvtSRBamRecord& srBR, boost::unordered_map<std::size_t, uint32_t>& readSample, boost::unordered_map<std::size_t, uint8_t>& readHap) {
    // Breakpoints
    typedef std::vector<Junction> TJunctionVector;
    typedef boost::unordered_map<std::size_t, TJunctionVector> TReadBp;
    TReadBp readBp;
    findAsmJunctions(c, readBp, readSample, readHap, srBR);
    fetchSVs(c, readBp, srBR);
  }


  template<typename TConfig, typename TSRStore>
  inline void
  _findAsmStructuralVariants(TConfig const& c, std::vector<StructuralVariantRecord>& svs, TSRStore& srStore, std::vector<int32_t>& svSample, std::vector<int8_t>& svHap) {
    // Assembly splits
    typedef std::vector<SRBamRecord> TSRBamRecord;
    typedef std::vector<TSRBamRecord> TSvtSRBamRecord;
    TSvtSRBamRecord srBR(2 * DELLY_SVT_TRANS, TSRBamRecord());

    // Sample and haplotype
    boost::unordered_map<std::size_t, uint32_t> readSample;
    boost::unordered_map<std::size_t, uint8_t> readHap;

    // Find SR breakpoints
    _findAsmBreakpoints(c, srBR, readSample, readHap);

    int32_t ci = 10;
    for(uint32_t svt = 0; svt < srBR.size(); ++svt) {
      if (srBR[svt].empty()) continue;
      std::sort(srBR[svt].begin(), srBR[svt].end());
      for(uint32_t i = 0; i < srBR[svt].size(); ++i) {
	int32_t svid = (int32_t) svs.size();
	srBR[svt][i].svid = svid;
	svs.push_back(StructuralVariantRecord(srBR[svt][i].chr, srBR[svt][i].pos, srBR[svt][i].chr2, srBR[svt][i].pos2, -ci, ci, -ci, ci, 1, srBR[svt][i].qual, srBR[svt][i].qual, srBR[svt][i].inslen, svt, svid));
	// Sample and haplotype
	boost::unordered_map<std::size_t, uint32_t>::const_iterator smpIt = readSample.find(srBR[svt][i].id);
	svSample.push_back((smpIt != readSample.end()) ? (int32_t) smpIt->second : -1);
	svHap.push_back((smpIt != readSample.end()) ? (int8_t) readHap[srBR[svt][i].id] : (int8_t) -1);
	// Track split-reads for consensus
	if (srStore.find(srBR[svt][i].id) == srStore.end()) srStore.insert(std::make_pair(srBR[svt][i].id, std::vector<SeqSlice>()));
	srStore[srBR[svt][i].id].push_back(SeqSlice(svid, srBR[svt][i].sstart, srBR[svt][i].inslen, srBR[svt][i].qual));
      }
    }
  }

  // Cigar I or D
  template<typename TConfig>
  inline bool
  _directIndelAllele(TConfig const& c, char const* seq, int32_t const seqlen, StructuralVariantRecord& sv) {
    int32_t bp = sv.svStart;
    if ((bp < 2) || (bp >= seqlen)) return false;
    if (sv.svt == 2) {
      // Deletion
      int32_t dellen = sv.svEnd - sv.svStart;
      if ((dellen <= 0) || (bp + dellen > seqlen)) return false;
      // Left-align
      int32_t homLeft = 0;
      while ((bp > 1) && (_ucBase(seq[bp-1]) == _ucBase(seq[bp+dellen-1]))) {
	--bp;
	++homLeft;
      }
      // Right homology
      int32_t homRight = 0;
      while ((bp + dellen + homRight < seqlen) && (_ucBase(seq[bp+homRight]) == _ucBase(seq[bp+dellen+homRight]))) ++homRight;
      sv.svStart = bp;
      sv.svEnd = bp + dellen;
      // Exact alleles or symbolic
      if (dellen <= c.indelsize) {
	std::string refVCF(dellen + 1, 'N');
	for(int32_t k = 0; k <= dellen; ++k) refVCF[k] = _ucBase(seq[bp-1+k]);
	std::string altVCF(1, _ucBase(seq[bp-1]));
	sv.alleles = _addAlleles(refVCF, altVCF);
      }
      sv.homLen = std::max(0, homLeft + homRight - 1);
      int32_t wig = std::max(homLeft, homRight);
      sv.ciposlow = -wig;
      sv.ciposhigh = wig;
      sv.ciendlow = -wig;
      sv.ciendhigh = wig;
      sv.precise = true;
      sv.srAlignQuality = 1.0;
      return true;
    } else if (sv.svt == 4) {
      // Insertion
      int32_t inslen = sv.insLen;
      int32_t cbp = sv.consBp;
      if ((inslen <= 0) || (cbp < 0) || (cbp + inslen > (int32_t) sv.consensus.size())) return false;
      std::string ins = sv.consensus.substr(cbp, inslen);
      // Left-align
      int32_t homLeft = 0;
      while ((bp > 1) && (_ucBase(seq[bp-1]) == ins[ins.size()-1])) {
	ins = ins[ins.size()-1] + ins.substr(0, ins.size()-1);
	--bp; ++homLeft;
      }
      sv.svStart = bp;
      sv.svEnd = bp;
      std::string altVCF(1, _ucBase(seq[bp-1]));
      altVCF += ins;
      std::string refVCF(1, _ucBase(seq[bp-1]));
      sv.alleles = _addAlleles(refVCF, altVCF);
      sv.homLen = homLeft;
      sv.ciposlow = -homLeft;
      sv.ciposhigh = homLeft;
      sv.ciendlow = -homLeft;
      sv.ciendhigh = homLeft;
      sv.precise = true;
      sv.srAlignQuality = 1.0;
      return true;
    }
    return false;
  }


  template<typename TConfig, typename TSRStore>
  inline void
  _setAsmConsensus(TConfig const& c, std::vector<StructuralVariantRecord>& svs, TSRStore& srStore) {
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

    // Parse BAM
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Parse assembly alleles" << std::endl;

    // Derive consensus sequences from assemblies
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	// Read alignments (full chromosome because primary alignments might be somewhere else)
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  // Only primary alignments with the full sequence information
	  if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;

	  // Seed
	  std::size_t seed = hash_lr(rec);
	  boost::hash_combine(seed, c.fileSample[file_c]);
	  if (srStore.find(seed) != srStore.end()) {
	    // Get sequence
	    std::string sequence;
	    sequence.resize(rec->core.l_qseq);
	    const uint8_t* seqptr = bam_get_seq(rec);
	    for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	    int32_t readlen = sequence.size();

	    // Iterate all spanned SVs
	    for(uint32_t ri = 0; ri < srStore[seed].size(); ++ri) {
	      SeqSlice seqsl = srStore[seed][ri];
	      int32_t svid = seqsl.svid;

	      // Extract subsequence
	      int32_t window = c.minConsWindow;
	      // Add breakpoint uncertainty
	      window += std::max(svs[svid].ciposhigh - svs[svid].ciposlow, svs[svid].ciendhigh - svs[svid].ciendlow);
	      int32_t sPos = seqsl.sstart - window;
	      int32_t ePos = seqsl.sstart + window + seqsl.inslen;
	      if (rec->core.flag & BAM_FREVERSE) {
		sPos = (readlen - seqsl.sstart) - window - seqsl.inslen;
		ePos = (readlen - seqsl.sstart) + window;
	      }
	      if (sPos < 0) sPos = 0;
	      if (ePos > (int32_t) readlen) ePos = readlen;
	      if (((ePos - sPos) > window) && ((ePos - sPos) < 100000)) {
		if ((ePos - sPos) > (int32_t) svs[svid].consensus.size()) {
		  svs[svid].consensus = sequence.substr(sPos, (ePos - sPos));
		  int32_t bpFwdPos = (rec->core.flag & BAM_FREVERSE) ? (readlen - seqsl.sstart - seqsl.inslen) : seqsl.sstart;
		  svs[svid].consBp = bpFwdPos - sPos;
		  //std::cerr << svs[svid].svStart << ',' << bam_get_qname(rec) << ',' << sPos << ',' << ePos << ":" << window << "\t" << svs[svid].consensus << std::endl;
		}
	      }
	    }
	  }
	}
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }
    }

    // Align consensus sequences
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      char* seq = NULL;
      for(int32_t refIndex2 = 0; refIndex2 <= refIndex; ++refIndex2) {
	char* sndSeq = NULL;
	for(uint32_t svid = 0; svid < svs.size(); ++svid) {
	  if (svs[svid].chr == refIndex) {
	    bool computeAlign = false;
	    if (!svs[svid].consensus.empty()) {
	      if (_translocation(svs[svid].svt)) {
		if ((refIndex2 != refIndex) && (svs[svid].chr2 == refIndex2)) {
		  // Lazy loading of references
		  if (sndSeq == NULL) {
		    int32_t seqlen = -1;
		    std::string tname(hdr->target_name[refIndex2]);
		    sndSeq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex2], &seqlen);
		  }
		  computeAlign = true;
		}
	      } else {
		if ((refIndex2 == refIndex) && (svs[svid].chr2 == refIndex2)) computeAlign = true;
	      }
	      if (computeAlign) {
		// Lazy loading of references
		if (seq == NULL) {
		  int32_t seqlen = -1;
		  std::string tname(hdr->target_name[refIndex]);
		  seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
		}
		// Align consensus
		bool success = false;
		if ((svs[svid].svt == 2) || (svs[svid].svt == 4)) success = _directIndelAllele(c, seq, hdr->target_len[refIndex], svs[svid]);
		else success = alignConsensus(c, hdr, seq, sndSeq, svs[svid], true);
		if (!success) {
		  svs[svid].consensus = "";
		  svs[svid].srSupport = 0;
		  svs[svid].srAlignQuality = 0;
		}
	      }
	    } else {
	      svs[svid].consensus = "";
	      svs[svid].srSupport = 0;
	      svs[svid].srAlignQuality = 0;
	    }
	  }
	}
	if (sndSeq != NULL) free(sndSeq);
      }
      // Set tag alleles
      for(uint32_t svid = 0; svid < svs.size(); ++svid) {
	if ((svs[svid].chr == refIndex) && (svs[svid].alleles.empty())) {
	  // Lazy loading of references
	  if (seq == NULL) {
	    int32_t seqlen = -1;
	    std::string tname(hdr->target_name[refIndex]);
	    seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
	  }
	  svs[svid].alleles = _addAlleles(boost::to_upper_copy(std::string(seq + svs[svid].svStart - 1, seq + svs[svid].svStart)), std::string(hdr->target_name[svs[svid].chr2]), svs[svid], svs[svid].svt);
	}
      }
      // Annotate SV subtype
      for(uint32_t svid = 0; svid < svs.size(); ++svid) {
	if ((svs[svid].chr == refIndex) && (!svs[svid].consensus.empty()) && (!_translocation(svs[svid].svt))) {
	  if (seq == NULL) {
	    int32_t seqlen = -1;
	    std::string tname(hdr->target_name[refIndex]);
	    seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
	  }
	  annotateSV(c, hdr, seq, svs[svid]);
	}
      }
      if (seq != NULL) free(seq);
    }

    // Clean-up
    fai_destroy(fai);	
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }
  
  // Cross-sample merge
  template<typename TConfig>
  inline void
  _asmMergeWrite(TConfig& c, std::vector<StructuralVariantRecord>& svs, std::vector<int32_t> const& svSample, std::vector<int8_t> const& svHap) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Merging SVs across samples" << std::endl;

    // Merge config
    MergeConfig mc;
    mc.filterForPass = false;
    mc.filterForPrecise = false;
    mc.cnvMode = false;
    mc.hasGenome = false;
    mc.chunksize = 500;
    mc.svcounter = 1;
    mc.alleleCounter = 1;
    mc.bpoffset = 1000;
    mc.minsize = 0;
    mc.maxsize = 1000000;
    mc.coverage = 0;
    mc.totalSamples = c.nsamples;
    mc.qualthres = 0;
    mc.recoverlap = 0.8f;
    mc.vaf = 0;
    mc.meiOffset = 50;
    mc.meiSizeRatio = 0.85f;
    mc.meiSeqId = 0.8f;
    mc.trOffset = 200;
    mc.trFrac = 0.25f;
    mc.trSeqId = 0.7f;
    mc.normFrac = 0.5f;
    mc.juncSeqId = 0.7f;
    mc.seqCutoff = 10000;
    mc.recurrentSamples = 10;
    mc.repMinAF = 0.005f;
    mc.outfile = c.outfile;
    mc.genome = c.genome;

    // Merged records and per-record haplotype support
    std::vector<StructuralVariantRecord> merged;
    std::vector<std::vector<uint8_t> > mergedHap;

    for(int32_t svt = 0; svt < 2 * DELLY_SVT_TRANS; ++svt) {
      // Build sample,haplotype nodes for this SV type
      std::vector<MergeSV> nodes;
      for(uint32_t i = 0; i < svs.size(); ++i) {
	if (svs[i].svt != svt) continue;
	if (svs[i].consensus.empty()) continue;
	if (svSample[i] < 0) continue;
	MergeSV n;
	n.tid = svs[i].chr;
	n.mtid = svs[i].chr2;
	n.svStart = svs[i].svStart;
	n.size = (svt == 4) ? svs[i].insLen : std::abs(svs[i].svEnd - svs[i].svStart);
	n.svEnd = (svt == 4) ? (svs[i].svStart + n.size) : svs[i].svEnd;
	n.pos2 = svs[i].svEnd;
	n.svt = svt;
	n.homlen = svs[i].homLen;
	n.trperiod = svs[i].anno.trPeriod;
	n.score = svs[i].mapq;
	n.fileIdx = svSample[i];
	n.hap = svHap[i];
	n.supp = 1;
	n.ac = 1;
	n.sr = 10;
	n.srmapq = svs[i].srMapQuality;
	n.altSupport = 10;
	n.subtype = (int8_t) svs[i].anno.seqType;
	n.insStrand = ((svs[i].anno.seqType >= 1) && (svs[i].anno.seqType <= 6)) ? (svs[i].anno.isRC ? 1 : 0) : -1;
	n.precise = true;
	n.fromSiteList = false;
	n.srq = svs[i].srAlignQuality;
	n.ce = 2.0f;
	n.comp = -1;
	n.id = boost::lexical_cast<std::string>(i);
	// Sequence used for allele-similarity clustering
	std::size_t comma = svs[i].alleles.find(',');
	std::string refA = (comma != std::string::npos) ? svs[i].alleles.substr(0, comma) : svs[i].alleles;
	std::string altA = (comma != std::string::npos) ? svs[i].alleles.substr(comma + 1) : std::string();
	std::string seq;
	if (svt == 4) {
	  if ((!altA.empty()) && (altA[0] != '<') && (altA.size() > 1) && (n.size <= mc.seqCutoff)) seq = boost::to_upper_copy(altA.substr(1));
	} else if (svt == 2) {
	  if ((refA.size() > 1) && (n.size <= mc.seqCutoff)) seq = boost::to_upper_copy(refA.substr(1));
	} else {
	  if ((!svs[i].consensus.empty()) && ((int32_t) svs[i].consensus.size() <= mc.seqCutoff)) seq = boost::to_upper_copy(svs[i].consensus);
	}
	n.seq = seq;
	nodes.push_back(n);
      }
      if (nodes.empty()) continue;

      boost::unordered_map<std::string, MergeAgg> selected;
      std::vector<AlleleGroup> groups;
      _clusterAndSelect(mc, nodes, selected, &groups);

      // One merged record
      for(uint32_t g = 0; g < groups.size(); ++g) {
	int32_t svIdx = boost::lexical_cast<int32_t>(nodes[groups[g].rep].id);
	StructuralVariantRecord rec = svs[svIdx];
	rec.id = (int32_t) merged.size();
	rec.alleleid = (int32_t) groups[g].agg.alleleId;
	rec.nallele = groups[g].agg.nAllele;
	rec.ciposlow = groups[g].agg.ciposLo;
	rec.ciposhigh = groups[g].agg.ciposHi;
	rec.ciendlow = groups[g].agg.ciendLo;
	rec.ciendhigh = groups[g].agg.ciendHi;
	rec.srSupport = 10;
	rec.mapq = svs[svIdx].mapq * 10;
	std::vector<uint8_t> hb(c.nsamples, 0);
	for(uint32_t m = 0; m < groups[g].members.size(); ++m) {
	  MergeSV const& mn = nodes[groups[g].members[m]];
	  if ((mn.fileIdx >= 0) && (mn.fileIdx < (int32_t) c.nsamples) && (mn.hap >= 0)) hb[mn.fileIdx] |= (uint8_t) (1 << mn.hap);
	}
	merged.push_back(rec);
	mergedHap.push_back(hb);
      }
    }

    // Sort by coordinate
    {
      std::vector<uint32_t> ord(merged.size());
      for(uint32_t i = 0; i < ord.size(); ++i) ord[i] = i;
      std::sort(ord.begin(), ord.end(), [&](uint32_t a, uint32_t b) { return (merged[a].chr < merged[b].chr) || ((merged[a].chr == merged[b].chr) && (merged[a].svStart < merged[b].svStart)); });
      std::vector<StructuralVariantRecord> ms;
      ms.reserve(merged.size());
      std::vector<std::vector<uint8_t> > mh;
      mh.reserve(merged.size());
      for(uint32_t i = 0; i < ord.size(); ++i) {
	StructuralVariantRecord r = merged[ord[i]];
	r.id = (int32_t) i;
	ms.push_back(r);
	mh.push_back(mergedHap[ord[i]]);
      }
      merged.swap(ms);
      mergedHap.swap(mh);
    }

    // Genotype maps
    typedef std::vector<JunctionCount> TSVJunctionMap;
    std::vector<TSVJunctionMap> jctMap(c.nsamples);
    typedef std::vector<SpanningCount> TSVSpanningMap;
    std::vector<TSVSpanningMap> spanMap(c.nsamples);
    typedef std::vector<ReadCount> TSVReadCount;
    std::vector<TSVReadCount> rcMap(c.nsamples);
    for(uint32_t s = 0; s < c.nsamples; ++s) {
      jctMap[s].resize(merged.size(), JunctionCount());
      spanMap[s].resize(merged.size(), SpanningCount());
      rcMap[s].resize(merged.size());
    }

    // Simple genotyping
    const uint8_t gtQual = 30;
    for(uint32_t nid = 0; nid < merged.size(); ++nid) {
      for(uint32_t s = 0; s < c.nsamples; ++s) {
	uint8_t bits = mergedHap[nid][s];
	if (c.diploid) {
	  bool h0 = (bits & 1);
	  bool h1 = (bits & 2);
	  if ((!h0) && (!h1)) {
	    jctMap[s][nid].ref.push_back(gtQual);
	    jctMap[s][nid].ref.push_back(gtQual);
	  } else if (h0 && h1) {
	    jctMap[s][nid].alt.push_back(gtQual);
	    jctMap[s][nid].hp1alt.push_back(gtQual);
	    jctMap[s][nid].alt.push_back(gtQual);
	    jctMap[s][nid].hp2alt.push_back(gtQual);
	  } else {
	    jctMap[s][nid].ps = 1;
	    if (h0) {
	      jctMap[s][nid].alt.push_back(gtQual);
	      jctMap[s][nid].hp1alt.push_back(gtQual);
	      jctMap[s][nid].ref.push_back(gtQual);
	      jctMap[s][nid].hp2ref.push_back(gtQual);
	    }
	    else {
	      jctMap[s][nid].ref.push_back(gtQual);
	      jctMap[s][nid].hp1ref.push_back(gtQual);
	      jctMap[s][nid].alt.push_back(gtQual);
	      jctMap[s][nid].hp2alt.push_back(gtQual);
	    }
	  }
	} else {
	  if (bits) {
	    jctMap[s][nid].alt.push_back(gtQual);
	    jctMap[s][nid].ref.push_back(gtQual);
	  } else {
	    jctMap[s][nid].ref.push_back(gtQual);
	    jctMap[s][nid].ref.push_back(gtQual);
	  }
	}
      }
    }

    // Write tmp BCF
    boost::filesystem::path finalOut = c.outfile;
    std::string model = (finalOut.string() == "-") ? std::string("./delly_asm_%%%%%%%%.bcf") : (finalOut.string() + ".%%%%%%%%.tmp.bcf");
    boost::filesystem::path tmpOut = boost::filesystem::unique_path(model);
    c.outfile = tmpOut;
    vcfOutput(c, merged, jctMap, rcMap, spanMap);
    c.outfile = finalOut;
    mc.outfile = finalOut;
    mergeBCFs(mc, std::vector<boost::filesystem::path>(1, tmpOut));
    boost::filesystem::remove(tmpOut);
    boost::filesystem::remove(boost::filesystem::path(tmpOut.string() + ".csi"));
  }

 template<typename TConfig>
 inline int32_t
 runAsm(TConfig& c) {
   
#ifdef PROFILE
   ProfilerStart("delly.prof");
#endif
   // Structural Variants
   std::vector<StructuralVariantRecord> svs;
   std::vector<int32_t> svSample;
   std::vector<int8_t> svHap;

   {
     typedef std::vector<SeqSlice> TSeqSliceVector;
     typedef boost::unordered_map<std::size_t, TSeqSliceVector> TReadSeqSlices;
     TReadSeqSlices srStore;

     // Per-(sample, haplotype) candidate SVs
     _findAsmStructuralVariants(c, svs, srStore, svSample, svHap);
     _setAsmConsensus(c, svs, srStore);
   }

   // Cross-sample SV merge
   _asmMergeWrite(c, svs, svSample, svHap);
   
#ifdef PROFILE
   ProfilerStop();
#endif

   return 0;
 }

 int asmode(int argc, char **argv) {
   AsmConfig c;
   c.hasVcfFile = true;
   c.minCliqueSize = 2;
   c.graphPruning = 1000;
   c.meiMinFrac = 0.8;
   c.trMinFrac = 0.85;
   
   // Parameter
   std::string svtype;
   std::string scoring;
   std::string mode;
   boost::program_options::options_description generic("Generic options");
   generic.add_options()
     ("help,?", "show help message")
     ("svtype,t", boost::program_options::value<std::string>(&svtype)->default_value("ALL"), "SV type to compute [DEL, INS, DUP, INV, BND, ALL]")
     ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
     ("mode,y", boost::program_options::value<std::string>(&mode)->default_value("squashed"), "input mode [squashed, diploid]")
     ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "SV output file")
     ;
   
   boost::program_options::options_description disc("Discovery options");
   disc.add_options()
     ("map-qual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(10), "min. mapping quality")
     ("minclip,c", boost::program_options::value<uint32_t>(&c.minClip)->default_value(25), "min. clipping length")
     ("minrefsep,m", boost::program_options::value<uint32_t>(&c.minRefSep)->default_value(30), "min. reference separation")
     ("maxreadsep,n", boost::program_options::value<uint32_t>(&c.maxReadSep)->default_value(100), "max. read separation")
     ;

   boost::program_options::options_description cons("Consensus options");
   cons.add_options()
     ("cons-window,w", boost::program_options::value<int32_t>(&c.minConsWindow)->default_value(1000), "consensus window")
     ("flank-size,f", boost::program_options::value<int32_t>(&c.minimumFlankSize)->default_value(100), "min. flank size")
     ("flank-quality,a", boost::program_options::value<float>(&c.flankQuality)->default_value(0.9), "min. flank quality")
     ("indel-size,i", boost::program_options::value<int32_t>(&c.indelsize)->default_value(10000), "use exact alleles for InDels <10kbp")
     ;

   boost::program_options::options_description hidden("Hidden options");
   hidden.add_options()
     ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
     ;
   
   boost::program_options::positional_options_description pos_args;
   pos_args.add("input-file", -1);
   
   boost::program_options::options_description cmdline_options;
   cmdline_options.add(generic).add(disc).add(cons).add(hidden);
   boost::program_options::options_description visible_options;
   visible_options.add(generic).add(disc).add(cons);
   boost::program_options::variables_map vm;
   boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
   boost::program_options::notify(vm);
   
   // Check command line arguments
   if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
     std::cerr << std::endl;
     std::cerr << "Usage: delly " << argv[0] << " [OPTIONS] -y squashed -g <ref.fa> <asm1.bam> <asm2.bam> ..." << std::endl;
     std::cerr << "       delly " << argv[0] << " [OPTIONS] -y diploid -g <ref.fa> <asm1.h1.bam> <asm1.h2.bam> <asm2.h1.bam> ..." << std::endl;
     std::cerr << visible_options << "\n";
     return 0;
   }

   // SV types to compute?
   if (!_svTypesToCompute(c, svtype)) {
     std::cerr << "Please specify a valid SV type, i.e., -t INV or -t DEL,INV without spaces." << std::endl;
     return 1;
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
   
   // Input mode
   if ((mode != "squashed") && (mode != "diploid")) {
     std::cerr << "Please specify a valid input mode (squashed or diploid)." << std::endl;
     return 1;
   }
   c.diploid = (mode == "diploid");
   if (c.diploid && (c.files.size() % 2 != 0)) {
     std::cerr << "Diploid mode expects an even number of BAM files (hap1 hap2 per sample)." << std::endl;
     return 1;
   }

   // Check input files
   std::vector<std::string> fileSM(c.files.size());
   c.nchr = 0;
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
     if (!c.nchr) c.nchr = hdr->n_targets;
     else {
       if (c.nchr != hdr->n_targets) {
	 std::cerr << "BAM files have different number of chromosomes!" << std::endl;
	 return 1;
       }
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
     std::string sampleName = "unknown";
     getSMTag(std::string(hdr->text), c.files[file_c].stem().string(), sampleName);
     fileSM[file_c] = sampleName;
     bam_hdr_destroy(hdr);
     hts_idx_destroy(idx);
     sam_close(samfile);
   }

   // Map input files to samples/haplotypes
   c.nsamples = (c.diploid) ? (c.files.size() / 2) : c.files.size();
   c.fileSample.resize(c.files.size());
   c.fileHap.resize(c.files.size());
   c.sampleName.resize(c.nsamples);
   for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
     if (c.diploid) {
       c.fileSample[file_c] = file_c / 2;
       c.fileHap[file_c] = (uint8_t) (file_c % 2);
     } else {
       c.fileSample[file_c] = file_c;
       c.fileHap[file_c] = 0;
     }
   }
   for(unsigned int s = 0; s < c.nsamples; ++s) {
     unsigned int f0 = (c.diploid) ? (2 * s) : s;
     c.sampleName[s] = (c.diploid) ? _commonSampleName(fileSM[f0], fileSM[f0 + 1]) : fileSM[f0];
   }
   // De-duplicate sample names
   {
     std::set<std::string> snames;
     uint32_t ucount = 0;
     for(unsigned int s = 0; s < c.nsamples; ++s) {
       while (snames.find(c.sampleName[s]) != snames.end()) {
	 std::cerr << "Warning: Duplicate sample names: " << c.sampleName[s] << std::endl;
	 c.sampleName[s] += "_" + boost::lexical_cast<std::string>(ucount++);
	 std::cerr << "Warning: Changing sample name to " << c.sampleName[s] << std::endl;
       }
       snames.insert(c.sampleName[s]);
     }
   }

   // Check outfile
   if (!vm.count("outfile")) c.outfile = "-";
   else {
     if (c.outfile.string() != "-") {
       if (!_outfileValid(c.outfile)) return 1;
     }
   }

   // Show cmd
   boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
   std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] ";
   std::cerr << "delly ";
   for(int i=0; i<argc; ++i) { std::cerr << argv[i] << ' '; }
   std::cerr << std::endl;
   
   // Run assembly mode
   return runAsm(c);
 }

}

#endif
