#ifndef TEGUA_H
#define TEGUA_H

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

#include "delly.h"
#include "coverage.h"
#include "util.h"
#include "junction.h"
#include "cluster.h"
#include "assemble.h"
#include "modvcf.h"

namespace torali {


  struct TeguaConfig {
    bool hasDumpFile;
    bool hasVcfFile;
    bool isHaplotagged;
    bool svtcmd;
    uint16_t minMapQual;
    uint32_t minClip;
    uint32_t minRefSep;
    uint32_t maxReadSep;
    uint32_t graphPruning;
    int32_t nchr;
    int32_t minimumFlankSize;
    float indelExtension;
    float flankQuality;
    std::set<int32_t> svtset;
    std::vector<int32_t> chrlen;
    std::string svtype;
    DnaScore<int> aliscore;
    boost::filesystem::path dumpfile;
    boost::filesystem::path outfile;
    boost::filesystem::path vcffile;
    std::vector<boost::filesystem::path> files;
    boost::filesystem::path genome;
    std::vector<std::string> sampleName;
  };
  

 template<typename TConfig>
 inline int32_t
 runTegua(TConfig const& c) {

#ifdef PROFILE
   ProfilerStart("delly.prof");
#endif
   
   // Structural Variants
   typedef std::vector<StructuralVariantRecord> TVariants;
   TVariants svs;

   // Annotate junction reads
   typedef std::vector<JunctionCount> TSVJunctionMap;
   typedef std::vector<TSVJunctionMap> TSampleSVJunctionMap;
   TSampleSVJunctionMap junctionCountMap(c.files.size());

   // Annotate spanning coverage
   typedef std::vector<SpanningCount> TSVSpanningMap;
   typedef std::vector<TSVSpanningMap> TSampleSVSpanningMap;
   TSampleSVSpanningMap spanCountMap(c.files.size());
   
   // Annotate coverage
   typedef std::vector<ReadCount> TSVReadCount;
   typedef std::vector<TSVReadCount> TSampleSVReadCount;
   TSampleSVReadCount rcMap(c.files.size());
   
   // Open header
   samFile* samfile = sam_open(c.files[0].string().c_str(), "r");
   bam_hdr_t* hdr = sam_hdr_read(samfile);

   // SV Discovery
   if (!c.hasVcfFile) {
     // SR Store
     typedef std::vector<SeqSlice> TSvPosVector;
     typedef boost::unordered_map<std::size_t, TSvPosVector> TReadSV;
     TReadSV srStore;

     // Find junctions
     {
       // Split-reads
       typedef std::vector<SRBamRecord> TSRBamRecord;
       typedef std::vector<TSRBamRecord> TSvtSRBamRecord;
       TSvtSRBamRecord srBR(2 * DELLY_SVT_TRANS, TSRBamRecord());
       {
	 // Breakpoints
	 typedef std::vector<Junction> TJunctionVector;
	 typedef std::map<std::size_t, TJunctionVector> TReadBp;
	 TReadBp readBp;
	 findJunctions(c, readBp);
	 fetchSVs(c, readBp, srBR);
       }
   
       // Debug
       if (c.hasDumpFile) outputSRBamRecords(c, srBR);

       // Cluster BAM records
       for(uint32_t svt = 0; svt < srBR.size(); ++svt) {
	 if (srBR[svt].empty()) continue;
	 
	 // Sort
	 std::sort(srBR[svt].begin(), srBR[svt].end(), SortSRBamRecord<SRBamRecord>());
	 
	 // Cluster
	 cluster(c, srBR[svt], svs, c.maxReadSep, svt);
	 
	 // Track split-reads
	 for(uint32_t i = 0; i < srBR[svt].size(); ++i) {
	   // Read assigned?
	   if (srBR[svt][i].svid != -1) {
	     if (srStore.find(srBR[svt][i].id) == srStore.end()) srStore.insert(std::make_pair(srBR[svt][i].id, TSvPosVector()));
	     srStore[srBR[svt][i].id].push_back(SeqSlice(srBR[svt][i].svid, srBR[svt][i].sstart, srBR[svt][i].inslen, srBR[svt][i].qual));
	   }
	 }
       }
     }
     
     // Assemble
     assemble(c, srStore, svs);

     // Re-sort SVs
     typedef std::map<int32_t, int32_t> TIdMap;
     TIdMap idMap;
     std::sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());
     uint32_t cliqueCount = 0;
     for(typename TVariants::iterator svIt = svs.begin(); svIt != svs.end(); ++svIt, ++cliqueCount) {
       idMap[svIt->id] = cliqueCount;
       svIt->id = cliqueCount;       
     }

     // Initialize count maps
     for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
       junctionCountMap[file_c].resize(svs.size(), JunctionCount());
       spanCountMap[file_c].resize(svs.size(), SpanningCount());
       rcMap[file_c].resize(svs.size(), ReadCount());
     }
      
     // Annotate ALT support
     for(typename TReadSV::const_iterator it = srStore.begin(); it != srStore.end(); ++it) {
       for(uint32_t i = 0; i < it->second.size(); ++i) {
	 junctionCountMap[0][idMap[it->second[i].svid]].alt.push_back(it->second[i].qual);
       }
     }
   } else {
     std::cerr << "Only de novo discovery is supported for long-reads!" << std::endl;
     return 1;
   }
   
   // Reference SV Genotyping
   trackRef(c, svs, junctionCountMap, rcMap);

   // VCF Output
   //outputVcf(c, svs);
   vcfOutput(c, svs, junctionCountMap, rcMap, spanCountMap);

   // Clean-up
   bam_hdr_destroy(hdr);
   sam_close(samfile);


#ifdef PROFILE
   ProfilerStop();
#endif

   // End
   boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
   std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
  
   return 0;
 }

 int teguaMain(int argc, char **argv) {
   TeguaConfig c;
   c.isHaplotagged = false;
   
   // Parameter
   boost::program_options::options_description generic("Generic options");
   generic.add_options()
     ("help,?", "show help message")
     ("svtype,t", boost::program_options::value<std::string>(&c.svtype)->default_value("ALL"), "SV type to compute [DEL, INS, DUP, INV, BND, ALL]")
     ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
     ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("sv.bcf"), "SV BCF output file")
     ;
   
   boost::program_options::options_description disc("Discovery options");
   disc.add_options()
     ("extension,e", boost::program_options::value<float>(&c.indelExtension)->default_value(0.6), "enforce indel extension")
     ("mapqual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(10), "min. mapping quality")
     ("minclip,c", boost::program_options::value<uint32_t>(&c.minClip)->default_value(50), "min. clipping length")
     ("minrefsep,m", boost::program_options::value<uint32_t>(&c.minRefSep)->default_value(50), "min. reference separation")
     ("maxreadsep,n", boost::program_options::value<uint32_t>(&c.maxReadSep)->default_value(75), "max. read separation")
     ;

   boost::program_options::options_description geno("Genotyping options");
   geno.add_options()
     ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input VCF/BCF file for genotyping")
     ("dump,d", boost::program_options::value<boost::filesystem::path>(&c.dumpfile), "gzipped output file for SV-reads (optional)")
     ;

   
   
   boost::program_options::options_description hidden("Hidden options");
   hidden.add_options()
     ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
     ("pruning,j", boost::program_options::value<uint32_t>(&c.graphPruning)->default_value(1000), "graph pruning cutoff")
     ;
   
   boost::program_options::positional_options_description pos_args;
   pos_args.add("input-file", -1);
   
   boost::program_options::options_description cmdline_options;
   cmdline_options.add(generic).add(disc).add(geno).add(hidden);
   boost::program_options::options_description visible_options;
   visible_options.add(generic).add(disc).add(geno);
   boost::program_options::variables_map vm;
   boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
   boost::program_options::notify(vm);
   
   // Check command line arguments
   if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
     std::cout << std::endl;
     std::cout << "Usage: dellyLR " << argv[0] << " [OPTIONS] -g <ref.fa> <sample1.sort.bam> <sample2.sort.bam> ..." << std::endl;
     std::cout << visible_options << "\n";
     return 0;
   }

   // SV types to compute?
   _svTypesToCompute(c, svtype, vm.count("svtype"));

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
     if (!c.nchr) {
       c.nchr = hdr->n_targets;
       c.chrlen.resize(hdr->n_targets);
       for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) c.chrlen[refIndex] = hdr->target_len[refIndex];
     } else {
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
   std::cout << "dellyLR ";
   for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
   std::cout << std::endl;
   
   // Run Tegua
   c.aliscore = DnaScore<int>(3, -2, -3, -1);
   c.flankQuality = 0.8;
   c.minimumFlankSize = 50;
   return runTegua(c);
 }

}

#endif
