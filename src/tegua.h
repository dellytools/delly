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


  struct TeguaConfig {
    bool hasDumpFile;
    bool hasExcludeFile;
    bool hasVcfFile;
    uint16_t minMapQual;
    uint16_t minGenoQual;
    uint32_t minClip;
    uint32_t minRefSep;
    uint32_t maxReadSep;
    uint32_t graphPruning;
    uint32_t minCliqueSize;
    uint32_t maxReadPerSV;
    uint32_t maxGenoReadCount;
    int32_t minConsWindow;
    int32_t nchr;
    int32_t minimumFlankSize;
    int32_t indelsize;
    float indelExtension;
    float flankQuality;
    std::set<int32_t> svtset;
    DnaScore<int> aliscore;
    boost::filesystem::path dumpfile;
    boost::filesystem::path outfile;
    boost::filesystem::path vcffile;
    std::vector<boost::filesystem::path> files;
    boost::filesystem::path genome;
    boost::filesystem::path exclude;
    std::vector<std::string> sampleName;
  };
  

 template<typename TConfig>
 inline int32_t
 runTegua(TConfig& c) {

#ifdef PROFILE
   ProfilerStart("delly.prof");
#endif

   // Structural Variants
   typedef std::vector<StructuralVariantRecord> TVariants;
   TVariants svs;

   // Open header
   samFile* samfile = sam_open(c.files[0].string().c_str(), "r");
   bam_hdr_t* hdr = sam_hdr_read(samfile);

   // Exclude intervals
   typedef boost::icl::interval_set<uint32_t> TChrIntervals;
   typedef std::vector<TChrIntervals> TRegionsGenome;
   TRegionsGenome validRegions;
   if (!_parseExcludeIntervals(c, hdr, validRegions)) {
     std::cerr << "Delly couldn't parse exclude intervals!" << std::endl;
     bam_hdr_destroy(hdr);
     sam_close(samfile);
     return 1;
   }
     
   // SV Discovery
   if (!c.hasVcfFile) {
       
     // Structural Variant Candidates
     typedef std::vector<StructuralVariantRecord> TVariants;
     TVariants svc;

     // Split-read store
     typedef std::pair<int32_t, std::size_t> TPosRead;
     typedef std::vector<SeqSlice> TSvPosVector;
     typedef boost::unordered_map<TPosRead, TSvPosVector> TPosReadSV;
     typedef std::vector<TPosReadSV> TGenomicPosReadSV;
     TGenomicPosReadSV srStore(c.nchr, TPosReadSV());

     // SV Discovery
     _clusterSRReads(c, validRegions, svc, srStore);

     // Assemble
     assemble(c, validRegions, svc, srStore);

     // Sort SVs
     sort(svc.begin(), svc.end(), SortSVs<StructuralVariantRecord>());
      
     // Remove duplicates
     StructuralVariantRecord lastSV;
     for(typename TVariants::iterator svIter = svc.begin(); svIter != svc.end(); ++svIter) {
       if ((svIter->srSupport == 0) && (svIter->peSupport == 0)) continue;
       // Duplicate?
       if (!svs.empty()) {
	 if ((lastSV.chr == svIter->chr) && (lastSV.chr2 == svIter->chr2) && (lastSV.svt == svIter->svt) && (std::abs(svIter->svStart - lastSV.svStart) < c.minRefSep) && (std::abs(svIter->svEnd - lastSV.svEnd) < c.minRefSep)) {
	   // Check length
	   int32_t lengthvar = std::min(0.1 * (svIter->svEnd - svIter->svStart), 0.1 * (lastSV.svEnd - lastSV.svStart));
	   int32_t lengthdiff = std::min(svIter->svEnd - svIter->svStart, lastSV.svEnd - lastSV.svStart);
	   if (lengthvar < 15) lengthvar = 15;
	   if (lengthdiff < lengthvar) continue;
	 }
       }
       lastSV = *svIter;
       svs.push_back(*svIter);
     }

     // Sort
     sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());
     
     // Re-number SVs
     uint32_t cliqueCount = 0;
     for(typename TVariants::iterator svIt = svs.begin(); svIt != svs.end(); ++svIt, ++cliqueCount) svIt->id = cliqueCount;
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

 int tegua(int argc, char **argv) {
   TeguaConfig c;
   
   // Parameter
   std::string svtype;
   std::string scoring;
   std::string mode;
   boost::program_options::options_description generic("Generic options");
   generic.add_options()
     ("help,?", "show help message")
     ("svtype,t", boost::program_options::value<std::string>(&svtype)->default_value("ALL"), "SV type to compute [DEL, INS, DUP, INV, BND, ALL]")
     ("technology,y", boost::program_options::value<std::string>(&mode)->default_value("ont"), "seq. technology [pb, ont]")
     ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
     ("exclude,x", boost::program_options::value<boost::filesystem::path>(&c.exclude), "file with regions to exclude")
     ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "BCF output file")
     ;
   
   boost::program_options::options_description disc("Discovery options");
   disc.add_options()
     ("mapqual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(1), "min. mapping quality")
     ("minclip,c", boost::program_options::value<uint32_t>(&c.minClip)->default_value(25), "min. clipping length")
     ("min-clique-size,z", boost::program_options::value<uint32_t>(&c.minCliqueSize)->default_value(3), "min. clique size")     
     ("minrefsep,m", boost::program_options::value<uint32_t>(&c.minRefSep)->default_value(30), "min. reference separation")
     ("maxreadsep,n", boost::program_options::value<uint32_t>(&c.maxReadSep)->default_value(100), "max. read separation")
     ;

   boost::program_options::options_description cons("Consensus options");
   cons.add_options()
     ("max-reads,p", boost::program_options::value<uint32_t>(&c.maxReadPerSV)->default_value(15), "max. reads for consensus computation")
     ("cons-window,w", boost::program_options::value<int32_t>(&c.minConsWindow)->default_value(1000), "consensus window")
     ("flank-size,f", boost::program_options::value<int32_t>(&c.minimumFlankSize)->default_value(100), "min. flank size")
     ("flank-quality,a", boost::program_options::value<float>(&c.flankQuality)->default_value(0.9), "min. flank quality")
     ("indel-size,i", boost::program_options::value<int32_t>(&c.indelsize)->default_value(10000), "use exact alleles for InDels <10kbp")
     ;
   
   boost::program_options::options_description geno("Genotyping options");
   geno.add_options()
     ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input VCF/BCF file for genotyping")
     ("geno-qual,u", boost::program_options::value<uint16_t>(&c.minGenoQual)->default_value(5), "min. mapping quality for genotyping")
     ("max-geno-count,b", boost::program_options::value<uint32_t>(&c.maxGenoReadCount)->default_value(250), "max. reads aligned for SR genotyping")
     ("dump,d", boost::program_options::value<boost::filesystem::path>(&c.dumpfile), "gzipped output file for SV-reads")
     ;

   boost::program_options::options_description hidden("Hidden options");
   hidden.add_options()
     ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
     ("pruning,j", boost::program_options::value<uint32_t>(&c.graphPruning)->default_value(1000), "graph pruning cutoff")
     ("extension,e", boost::program_options::value<float>(&c.indelExtension)->default_value(0.5), "enforce indel extension")
     ("scoring,s", boost::program_options::value<std::string>(&scoring)->default_value("3,-2,-3,-1"), "alignment scoring")
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
     std::cerr << "Usage: delly " << argv[0] << " [OPTIONS] -g <ref.fa> <sample1.sort.bam> <sample2.sort.bam> ..." << std::endl;
     std::cerr << visible_options << "\n";
     return 0;
   }

   // Set alignment score
   _alignmentScore(c, scoring);
   
   // SV types to compute?
   if (!_svTypesToCompute(c, svtype)) {
     std::cerr << "Please specify a valid SV type, i.e., -t INV or -t DEL,INV without spaces." << std::endl;
     return 1;
   }

   // Dump reads
   if (vm.count("dump")) c.hasDumpFile = true;
   else c.hasDumpFile = false;

   // Clique size
   if (c.minCliqueSize < 2) c.minCliqueSize = 2;

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

   // Check exclude file
   if (vm.count("exclude")) {
     if (!(boost::filesystem::exists(c.exclude) && boost::filesystem::is_regular_file(c.exclude) && boost::filesystem::file_size(c.exclude))) {
       std::cerr << "Exclude file is missing: " << c.exclude.string() << std::endl;
       return 1;
     }
     c.hasExcludeFile = true;
   } else c.hasExcludeFile = false;

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
   
   // Run Tegua
   if (mode == "pb") c.indelExtension = 0.7;
   else if (mode == "ont") c.indelExtension = 0.5;
   return runTegua(c);
 }

}

#endif
