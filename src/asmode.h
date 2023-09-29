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
    std::set<int32_t> svtset;
    boost::filesystem::path outfile;
    std::vector<boost::filesystem::path> files;
    boost::filesystem::path genome;
    std::vector<std::string> sampleName;
  };

  template<typename TConfig, typename TReadBp>
  inline void
  findAsmJunctions(TConfig const& c, TReadBp& readBp) {
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

	  std::size_t seed = hash_lr(rec);
	  //std::cerr << bam_get_qname(rec) << '\t' << seed << std::endl;
	  uint32_t rp = rec->core.pos; // reference pointer
	  uint32_t sp = 0; // sequence pointer
	    
	  // Parse the CIGAR
	  uint32_t* cigar = bam_get_cigar(rec);
	  for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	    if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	      sp += bam_cigar_oplen(cigar[i]);
	      rp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	      if (bam_cigar_oplen(cigar[i]) > c.minRefSep) _insertJunction(readBp, seed, rec, rp, sp, false);
	      rp += bam_cigar_oplen(cigar[i]);
	      _insertJunction(readBp, seed, rec, rp, sp, true);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	      if (bam_cigar_oplen(cigar[i]) > c.minRefSep) _insertJunction(readBp, seed, rec, rp, sp, false);
	      sp += bam_cigar_oplen(cigar[i]);
	      _insertJunction(readBp, seed, rec, rp, sp, true);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	      rp += bam_cigar_oplen(cigar[i]);
	    } else if ((bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) {
	      int32_t finalsp = sp;
	      bool scleft = false;
	      if (sp == 0) {
		finalsp += bam_cigar_oplen(cigar[i]); // Leading soft-clip / hard-clip
		scleft = true;
	      }
	      sp += bam_cigar_oplen(cigar[i]);
	      //std::cerr << bam_get_qname(rec) << ',' << rp << ',' << finalsp << ',' << scleft << std::endl;
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
    for(typename TReadBp::iterator it = readBp.begin(); it != readBp.end(); ++it) {
      std::sort(it->second.begin(), it->second.end(), SortJunction<Junction>());
    }

    // Clean-up
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }
 

  template<typename TConfig, typename TSvtSRBamRecord>
  inline void
  _findAsmBreakpoints(TConfig const& c, TSvtSRBamRecord& srBR) {
    // Breakpoints
    typedef std::vector<Junction> TJunctionVector;
    typedef std::map<std::size_t, TJunctionVector> TReadBp;
    TReadBp readBp;
    findAsmJunctions(c, readBp);
    fetchSVs(c, readBp, srBR);
  }


  template<typename TConfig, typename TSRStore>
  inline void
  _findAsmStructuralVariants(TConfig const& c, std::vector<StructuralVariantRecord>& svs, TSRStore& srStore) {
    // Assembly splits
    typedef std::vector<SRBamRecord> TSRBamRecord;
    typedef std::vector<TSRBamRecord> TSvtSRBamRecord;
    TSvtSRBamRecord srBR(2 * DELLY_SVT_TRANS, TSRBamRecord());
    
    // Find SR breakpoints
    _findAsmBreakpoints(c, srBR);
   
    // Cluster SVs
    int32_t ci = 10;
    for(uint32_t svt = 0; svt < srBR.size(); ++svt) {
      if (srBR[svt].empty()) continue;
      // Sort
      std::sort(srBR[svt].begin(), srBR[svt].end(), SortSRBamRecord<SRBamRecord>());
      //outputSRBamRecords(c, srBR, true);

      // Cluster
      cluster(c, srBR[svt], svs, svt);
      // Add singletons
      int32_t svid = svs.size();
      for(uint32_t i = 0; i<srBR[svt].size(); ++i) {
	if (srBR[svt][i].svid == -1) {
	  srBR[svt][i].svid = svid;
	  svs.push_back(StructuralVariantRecord(srBR[svt][i].chr, srBR[svt][i].pos, srBR[svt][i].chr2, srBR[svt][i].pos2, -ci, ci, -ci, ci, 1,srBR[svt][i].qual, srBR[svt][i].qual, srBR[svt][i].inslen, svt, svid));
	  ++svid;
	}
      }
      //outputStructuralVariants(c, svs, srBR, svt, true);

      // Track split-reads
      for(uint32_t i = 0; i < srBR[svt].size(); ++i) {
	if (srBR[svt][i].svid != -1){
	  if (srStore.find(srBR[svt][i].id) == srStore.end()) srStore.insert(std::make_pair(srBR[svt][i].id, std::vector<SeqSlice>()));
	  srStore[srBR[svt][i].id].push_back(SeqSlice(srBR[svt][i].svid, srBR[svt][i].sstart, srBR[svt][i].inslen, srBR[svt][i].qual));
	}
      }
    }
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

	  std::size_t seed = hash_lr(rec);
	  if (srStore.find(seed) != srStore.end()) {
	    // Get sequence
	    std::string sequence;
	    sequence.resize(rec->core.l_qseq);
	    uint8_t* seqptr = bam_get_seq(rec);
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
		bool realign = true;
		if (svs[svid].svt == 4) realign = false;
		if (!alignConsensus(c, hdr, seq, sndSeq, svs[svid], realign)) {
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
  
 template<typename TConfig>
 inline int32_t
 runAsm(TConfig& c) {
   
#ifdef PROFILE
   ProfilerStart("delly.prof");
#endif
   // Structural Variants
   typedef std::vector<StructuralVariantRecord> TVariants;
   TVariants svs;

   {
     TVariants svc;
     
     // Split-read store
     typedef std::vector<SeqSlice> TSeqSliceVector;
     typedef boost::unordered_map<std::size_t, TSeqSliceVector> TReadSeqSlices;
     TReadSeqSlices srStore;

     // Find candidate SVs
     _findAsmStructuralVariants(c, svc, srStore);   
     
     // Set consensus sequence
     _setAsmConsensus(c, svc, srStore);
   
     // Sort SVs
     sort(svc.begin(), svc.end(), SortSVs<StructuralVariantRecord>());

     // Filter SVs
     uint32_t idCount = 0;
     for(uint32_t i = 0; i < svc.size(); ++i) {
       if (svc[i].consensus.empty()) continue;
       svc[i].srSupport = 10;
       svc[i].mapq *= 10;
       svc[i].id = idCount++;
       svs.push_back(svc[i]);
     }
   }

   // Annotate junction reads
   typedef std::vector<JunctionCount> TSVJunctionMap;
   typedef std::vector<TSVJunctionMap> TSampleSVJunctionMap;
   TSampleSVJunctionMap jctMap(c.files.size());

   // Annotate spanning coverage
   typedef std::vector<SpanningCount> TSVSpanningMap;
   typedef std::vector<TSVSpanningMap> TSampleSVSpanningMap;
   TSampleSVSpanningMap spanMap(c.files.size());
   
   // Annotate coverage
   typedef std::vector<ReadCount> TSVReadCount;
   typedef std::vector<TSVReadCount> TSampleSVReadCount;
   TSampleSVReadCount rcMap(c.files.size());

   // Initialize count maps
   for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
     jctMap[file_c].resize(svs.size(), JunctionCount());
     spanMap[file_c].resize(svs.size(), SpanningCount());
     rcMap[file_c].resize(svs.size(), ReadCount());
   }

   // VCF Output
   vcfOutput(c, svs, jctMap, rcMap, spanMap);
   
#ifdef PROFILE
   ProfilerStop();
#endif

   // End
   boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
   std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
  
   return 0;
 }

 int asmode(int argc, char **argv) {
   AsmConfig c;
   c.hasVcfFile = true;
   c.minCliqueSize = 2;
   c.graphPruning = 1000;
   
   // Parameter
   std::string svtype;
   std::string scoring;
   std::string mode;
   boost::program_options::options_description generic("Generic options");
   generic.add_options()
     ("help,?", "show help message")
     ("svtype,t", boost::program_options::value<std::string>(&svtype)->default_value("ALL"), "SV type to compute [DEL, INS, DUP, INV, BND, ALL]")
     ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
     ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "SV output file")
     ;
   
   boost::program_options::options_description disc("Discovery options");
   disc.add_options()
     ("mapqual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(10), "min. mapping quality")
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
     std::cerr << "Usage: delly " << argv[0] << " [OPTIONS] -g <ref.fa> <asm1.bam> <asm2.bam> ..." << std::endl;
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
   
   // Check input files
   c.sampleName.resize(c.files.size());
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
     c.sampleName[file_c] = sampleName;
     bam_hdr_destroy(hdr);
     hts_idx_destroy(idx);
     sam_close(samfile);
   }
   checkSampleNames(c);

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
