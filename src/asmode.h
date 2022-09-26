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
    bool hasDumpFile;
    bool hasVcfFile;
    uint16_t minMapQual;
    uint16_t minGenoQual;
    uint32_t minClip;
    uint32_t minRefSep;
    uint32_t maxReadSep;
    int32_t nchr;
    int32_t minimumFlankSize;
    int32_t indelsize;
    int32_t maxInsertionSize;
    float indelExtension;
    float flankQuality;
    std::set<int32_t> svtset;
    boost::filesystem::path dumpfile;
    boost::filesystem::path outfile;
    boost::filesystem::path vcffile;
    std::vector<boost::filesystem::path> files;
    boost::filesystem::path genome;
    std::vector<std::string> sampleName;
  };

  template<typename TConfig, typename TSRStore>
  inline void
  verifySplit(TConfig const& c, std::vector<StructuralVariantRecord>& svs, TSRStore& srStore) {
    // Sequence store
    typedef std::vector<std::string> TSequences;
    typedef std::vector<TSequences> TSVSequences;
    TSVSequences seqStore(svs.size(), TSequences());

    // SV consensus done
    std::vector<bool> svcons(svs.size(), false);
    
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
    {
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Contig alignment (phase 1)" << std::endl;
      
      faidx_t* fai = fai_load(c.genome.string().c_str());
      for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
	// Load sequence
	int32_t seqlen = -1;
	std::string tname(hdr->target_name[refIndex]);
	char* seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
	
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
		if (!svcons[svid]) {
		  //std::cerr << "Info:" << seqsl.sstart << ',' << seqsl.inslen << ',' << readlen << std::endl;
		  int32_t window = 1000;
		  int32_t sPos = seqsl.sstart - window;
		  int32_t ePos = seqsl.sstart + seqsl.inslen + window;
		  if (rec->core.flag & BAM_FREVERSE) {
		    sPos = (readlen - (seqsl.sstart + seqsl.inslen)) - window;
		    ePos = (readlen - seqsl.sstart) + window;
		  }
		  if (sPos < 0) sPos = 0;
		  if (ePos > (int32_t) readlen) ePos = readlen;
		  // Min. seq length
		  if (((ePos - sPos) > window) && ((ePos - sPos) <= (c.maxInsertionSize + window))) {
		    std::string seqalign = sequence.substr(sPos, (ePos - sPos));
		    // Can we align directly?
		    if ((!_translocation(svs[svid].svt)) && (svs[svid].chr == refIndex)) {
		      //std::cerr << seqalign.size() << ',' << refIndex << ',' << svs[svid].svt << std::endl;
		      svs[svid].consensus = seqalign;
		      if (!alignConsensus(c, hdr, seq, NULL, svs[svid], true)) {
			svs[svid].consensus = "";
			svs[svid].srSupport = 0;
			svs[svid].srAlignQuality = 0;
		      } else {
			seqStore[svid].clear();
			svcons[svid] = true;
		      }
		    } else {
		      seqStore[svid].push_back(seqalign);
		    }
		  }
		}
	      }
	    }
	  }
	  bam_destroy1(rec);
	  hts_itr_destroy(iter);
	}
	// Clean-up
	if (seq != NULL) free(seq);
      }
      // Clean-up
      fai_destroy(fai);
    }

    // Handle left-overs and translocations
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Contig alignment (phase 2)" << std::endl;
    
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      char* seq = NULL;
      for(int32_t refIndex2 = 0; refIndex2 <= refIndex; ++refIndex2) {
	char* sndSeq = NULL;
	for(uint32_t svid = 0; svid < svcons.size(); ++svid) {
	  if (!svcons[svid]) {
	    if (seqStore[svid].size() > 1) {
	      bool computeAlignment = false;
	      if (_translocation(svs[svid].svt)) {
		if ((refIndex2 != refIndex) && (svs[svid].chr == refIndex) && (svs[svid].chr2 == refIndex2)) {
		  computeAlignment = true;
		  // Lazy loading of references
		  if (sndSeq == NULL) {
		    int32_t seqlen = -1;
		    std::string tname(hdr->target_name[refIndex2]);
		    sndSeq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex2], &seqlen);
		  }
		}
	      } else {
		if ((refIndex2 == refIndex) && (svs[svid].chr == refIndex) && (svs[svid].chr2 == refIndex2)) computeAlignment =	true;
	      }
	      if (computeAlignment) {
		// Lazy loading of reference
		if (seq == NULL) {
		  int32_t seqlen = -1;
		  std::string tname(hdr->target_name[refIndex]);
		  sndSeq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
		}
		for(uint32_t k = 0; k < seqStore[svid].size(); ++k) {
		  svs[svid].consensus = seqStore[svid][k];
		  if (!alignConsensus(c, hdr, seq, sndSeq, svs[svid], true)) {
		    svs[svid].consensus = "";
		    svs[svid].srSupport = 0;
		    svs[svid].srAlignQuality = 0;
		  } else {
		    seqStore[svid].clear();
		    svcons[svid] = true;
		  }
		}
	      }
	    }
	  }
	}
	if (sndSeq != NULL) free(sndSeq);
      }
      // Clean-up
      if (seq != NULL) free(seq);
    }
    fai_destroy(fai);
    
    // Clean-up
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }

    // Flag unfinished SVs
    for(uint32_t svid = 0; svid < svcons.size(); ++svid) {
      if (!svcons[svid]) {
	//std::cerr << "Missing: " << svid << ',' << svs[svid].svt << std::endl;
	svs[svid].consensus = "";
	svs[svid].srSupport = 0;
	svs[svid].srAlignQuality = 0;
      }
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

   // Open header
   samFile* samfile = sam_open(c.files[0].string().c_str(), "r");
   bam_hdr_t* hdr = sam_hdr_read(samfile);

   // SV Discovery
   if (!c.hasVcfFile) {
     // Structural Variant Candidates
     typedef std::vector<StructuralVariantRecord> TVariants;
     TVariants svc;
     
     // Assembly splits
     typedef std::vector<SRBamRecord> TSRBamRecord;
     typedef std::vector<TSRBamRecord> TSvtSRBamRecord;
     TSvtSRBamRecord srBR(2 * DELLY_SVT_TRANS, TSRBamRecord());

     typedef boost::icl::interval_set<uint32_t> TChrIntervals;
     typedef std::vector<TChrIntervals> TRegionsGenome;
     typedef typename TChrIntervals::interval_type TIVal;
     TRegionsGenome validRegions(hdr->n_targets);
     for(int32_t i = 0; i < hdr->n_targets; ++i) {
       validRegions[i].insert(TIVal::right_open(0, hdr->target_len[i]));
     }
     _findSRBreakpoints(c, validRegions, srBR);
     //outputSRBamRecords(c, srBR);

     // Split-read store
     typedef std::vector<SeqSlice> TSvPosVector;
     typedef boost::unordered_map<std::size_t, TSvPosVector> TReadSV;
     TReadSV srStore;

     // Track split-reads
     uint32_t svid = 0;
     for(uint32_t svt = 0; svt < srBR.size(); ++svt) {
       for(uint32_t i = 0; i < srBR[svt].size(); ++i) {
	 // Create SV call
	 svc.push_back(StructuralVariantRecord(srBR[svt][i].chr, srBR[svt][i].pos, srBR[svt][i].chr2, srBR[svt][i].pos2, -10, 10, -10, 10, 1, srBR[svt][i].qual, srBR[svt][i].qual, srBR[svt][i].inslen, svt, svid));

	 // Track split-read
	 if (srStore.find(srBR[svt][i].id) == srStore.end()) srStore.insert(std::make_pair(srBR[svt][i].id, std::vector<SeqSlice>()));
	 srStore[srBR[svt][i].id].push_back(SeqSlice(svid, srBR[svt][i].sstart, srBR[svt][i].inslen, srBR[svt][i].qual));
	 ++svid;
       }
     }

     // Verify split
     verifySplit(c, svc, srStore);

     // Sort SVs
     sort(svc.begin(), svc.end(), SortSVs<StructuralVariantRecord>());

     // Keep verified SVs only
     StructuralVariantRecord lastSV;
     for(typename TVariants::iterator svIter = svc.begin(); svIter != svc.end(); ++svIter) {
       if ((svIter->srSupport == 0) && (svIter->peSupport == 0)) continue;
       // Duplicate?
       if (!svs.empty()) {
	 if ((lastSV.chr == svIter->chr) && (lastSV.chr2 == svIter->chr2) && (lastSV.svt == svIter->svt) && (std::abs(svIter->svStart - lastSV.svStart) < c.minRefSep) && (std::abs(svIter->svEnd - lastSV.svEnd) < c.minRefSep)) continue;
       }
       lastSV = *svIter;
       svs.push_back(*svIter);
     }

     // Sort
     sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());
     
     // Re-number SVs
     uint32_t cliqueCount = 0;
     for(typename TVariants::iterator svIt = svs.begin(); svIt != svs.end(); ++svIt, ++cliqueCount) svIt->id = cliqueCount;
     //outputStructuralVariants(c, svs);
   } else vcfParse(c, hdr, svs);
   // Clean-up
   bam_hdr_destroy(hdr);
   sam_close(samfile);

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
      
   // SV Genotyping
   genotypeLR(c, svs, jctMap, rcMap);

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
   
   // Parameter
   std::string svtype;
   std::string scoring;
   std::string mode;
   boost::program_options::options_description generic("Generic options");
   generic.add_options()
     ("help,?", "show help message")
     ("svtype,t", boost::program_options::value<std::string>(&svtype)->default_value("ALL"), "SV type to compute [DEL, INS, DUP, INV, BND, ALL]")
     ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
     ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "SV BCF output file")
     ;
   
   boost::program_options::options_description disc("Discovery options");
   disc.add_options()
     ("mapqual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(10), "min. mapping quality")
     ("minclip,c", boost::program_options::value<uint32_t>(&c.minClip)->default_value(25), "min. clipping length")
     ("minrefsep,m", boost::program_options::value<uint32_t>(&c.minRefSep)->default_value(30), "min. reference separation")
     ("maxreadsep,n", boost::program_options::value<uint32_t>(&c.maxReadSep)->default_value(100), "max. read separation")
     ("extension,e", boost::program_options::value<float>(&c.indelExtension)->default_value(0.5), "enforce indel extension")
     ;

   boost::program_options::options_description cons("Consensus options");
   cons.add_options()
     ("flank-size,f", boost::program_options::value<int32_t>(&c.minimumFlankSize)->default_value(100), "min. flank size")
     ("flank-quality,a", boost::program_options::value<float>(&c.flankQuality)->default_value(0.9), "min. flank quality")
     ("indel-size,i", boost::program_options::value<int32_t>(&c.indelsize)->default_value(10000), "use exact alleles for InDels <10kbp")
     ("max-isize,r", boost::program_options::value<int32_t>(&c.maxInsertionSize)->default_value(10000), "max. insertion size")
     ;     
   
   boost::program_options::options_description geno("Genotyping options");
   geno.add_options()
     ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input VCF/BCF file for genotyping")
     ("geno-qual,u", boost::program_options::value<uint16_t>(&c.minGenoQual)->default_value(5), "min. mapping quality for genotyping")
     ("dump,d", boost::program_options::value<boost::filesystem::path>(&c.dumpfile), "gzipped output file for SV-reads")
     ;

   boost::program_options::options_description hidden("Hidden options");
   hidden.add_options()
     ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
     ;
   
   boost::program_options::positional_options_description pos_args;
   pos_args.add("input-file", -1);
   
   boost::program_options::options_description cmdline_options;
   cmdline_options.add(generic).add(disc).add(cons).add(geno).add(hidden);
   boost::program_options::options_description visible_options;
   visible_options.add(generic).add(disc).add(cons).add(geno);
   boost::program_options::variables_map vm;
   boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
   boost::program_options::notify(vm);
   
   // Check command line arguments
   if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
     std::cerr << std::endl;
     std::cerr << "Usage: delly " << argv[0] << " [OPTIONS] -g <ref.fa> <haplotype1.sort.bam> <haplotype2.sort.bam> ..." << std::endl;
     std::cerr << visible_options << "\n";
     return 0;
   }

   // SV types to compute?
   if (!_svTypesToCompute(c, svtype)) {
     std::cerr << "Please specify a valid SV type, i.e., -t INV or -t DEL,INV without spaces." << std::endl;
     return 1;
   }

   // Dump reads
   if (vm.count("dump")) c.hasDumpFile = true;
   else c.hasDumpFile = false;

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

   // Check input VCF file
   if (vm.count("vcffile")) {
     if (!(boost::filesystem::exists(c.vcffile) && boost::filesystem::is_regular_file(c.vcffile) && boost::filesystem::file_size(c.vcffile))) {
       std::cerr << "Input VCF/BCF file is missing: " << c.vcffile.string() << std::endl;
       return 1;
     }
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
     c.hasVcfFile = true;
   } else c.hasVcfFile = false;
   
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
