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
    bool svtcmd;
    uint16_t minMapQual;
    uint32_t minClip;
    uint32_t minRefSep;
    uint32_t maxReadSep;
    int32_t nchr;
    int32_t indelsize;
    float indelExtension;
    std::set<int32_t> svtset;
    boost::filesystem::path outfile;
    boost::filesystem::path vcffile;
    std::vector<boost::filesystem::path> files;
    boost::filesystem::path genome;
    std::vector<std::string> sampleName;
  };
  
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
     outputSRBamRecords(c, srBR);

     // Process
     for(uint32_t svt = 0; svt < srBR.size(); ++svt) {
       if (srBR[svt].empty()) continue;
     }
   } else vcfParse(c, hdr, svs);
   // Clean-up
   bam_hdr_destroy(hdr);
   sam_close(samfile);

   // VCF Output
   //vcfOutput(c, svs, jctMap, rcMap, spanMap);

#ifdef PROFILE
   ProfilerStop();
#endif

   // End
   boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
   std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
  
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
     ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("sv.bcf"), "SV BCF output file")
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
     ("indel-size,i", boost::program_options::value<int32_t>(&c.indelsize)->default_value(10000), "use exact alleles for InDels <10kbp")
     ;     
   
   boost::program_options::options_description geno("Genotyping options");
   geno.add_options()
     ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input VCF/BCF file for genotyping")
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
     std::cout << std::endl;
     std::cout << "Usage: delly " << argv[0] << " [OPTIONS] -g <ref.fa> <haplotype1.sort.bam> <haplotype2.sort.bam> ..." << std::endl;
     std::cout << visible_options << "\n";
     return 0;
   }

   // SV types to compute?
   _svTypesToCompute(c, svtype, vm.count("svtype"));

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
   
   // Check output directory
   if (!_outfileValid(c.outfile)) return 1;

   // Show cmd
   boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
   std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
   std::cout << "delly ";
   for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
   std::cout << std::endl;
   
   // Run assembly mode
   return runAsm(c);
 }

}

#endif
