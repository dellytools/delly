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
#include <htslib/faidx.h>

#include "util.h"
#include "junction.h"
#include "cluster.h"
#include "modvcf.h"

namespace torali {


  struct TeguaConfig {
    bool hasSvInput;
    bool hasDumpFile;
    uint16_t minMapQual;
    uint32_t minClip;
    uint32_t minRefSep;
    uint32_t maxReadSep;
    uint32_t graphPruning;
    int32_t nchr;
    float indelExtension;
    std::vector<int32_t> chrlen;
    std::string svtype;
    boost::filesystem::path svinput;
    boost::filesystem::path inputfile;
    boost::filesystem::path dumpfile;
    boost::filesystem::path outfile;
    boost::filesystem::path genome;
  };
  

 template<typename TConfig>
 inline int32_t
 loadSV(TConfig const& c, std::vector<std::vector<SRBamRecord> >& br) {
   samFile* samfile = sam_open(c.inputfile.string().c_str(), "r");
   hts_set_fai_filename(samfile, c.genome.string().c_str());
   bam_hdr_t* hdr = sam_hdr_read(samfile);
   
   // Annotate intervals
   std::ifstream svFile(c.svinput.string().c_str(), std::ifstream::in);
   if (svFile.is_open()) {
     uint32_t svcount = 0;
     while (svFile.good()) {
       std::string svFromFile;
       getline(svFile, svFromFile);
       typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
       boost::char_separator<char> sep(" \t,;");
       Tokenizer tokens(svFromFile, sep);
       Tokenizer::iterator tokIter = tokens.begin();
       if (tokIter!=tokens.end()) {
	 std::string chr = *tokIter++;
	 int32_t tid = bam_name2id(hdr, chr.c_str());
	 if ((tid < 0) || (tid >= (int32_t) hdr->n_targets)) continue;
	 int32_t pos = boost::lexical_cast<int32_t>(*tokIter++);
	 if (tokIter!=tokens.end()) {
	   std::string chr2 = *tokIter++;
	   int32_t mid = bam_name2id(hdr, chr2.c_str());
	   if ((mid < 0) || (mid >= (int32_t) hdr->n_targets)) continue;
	   int32_t pos2 = boost::lexical_cast<int32_t>(*tokIter++);
	   // Parse ID
	   if (tokIter!=tokens.end()) {
	     std::string idname = *tokIter++;
	     if (idname.size() >= 3) {
	       if (idname.substr(0,3) == "DEL") br[2].push_back(SRBamRecord(tid, pos, mid, pos2, 2, 0, 0, svcount++));
	       else if (idname.substr(0,3) == "DUP") br[3].push_back(SRBamRecord(tid, pos, mid, pos2, 3, 0, 0, svcount++));
	       else if (idname.substr(0,3) == "INS") br[4].push_back(SRBamRecord(tid, pos, mid, pos2, 4, 0, 0, svcount++));
	       else if (idname.substr(0,3) == "INV") {
		 // Both inversion types
		 br[0].push_back(SRBamRecord(tid, pos, mid, pos2, 0, 0, 0, svcount++));
		 br[1].push_back(SRBamRecord(tid, pos, mid, pos2, 1, 0, 0, svcount++));
	       }
	     }
	   }
	 }
       }
     }
     svFile.close();
   }
   
   // Clean-up
   bam_hdr_destroy(hdr);
   sam_close(samfile);
   
   return br.size();
 }

 template<typename TConfig>
 inline int32_t
 runTegua(TConfig const& c) {
   // Load SVs
   typedef std::vector<SRBamRecord> TSRBamRecord;
   typedef std::vector<TSRBamRecord> TSvtSRBamRecord;
   TSvtSRBamRecord insv;
   if (c.hasSvInput) {
     int32_t insvlen = loadSV(c, insv);
     if (!insvlen) {
       std::cerr << "Error in loading SVs!" << std::endl;
       return -1;
     } else {
       std::cout << "Loaded " << insvlen << " SVs" << std::endl;
     }
   }

   // SR Store
   typedef std::pair<int32_t, std::size_t> TPosRead;
   typedef boost::unordered_map<TPosRead, int32_t> TPosReadSV;
   typedef std::vector<TPosReadSV> TGenomicPosReadSV;
   TGenomicPosReadSV srStore(c.chrlen.size(), TPosReadSV());

   // Structural Variants
   std::vector<StructuralVariantRecord> srSVs;

   
   // Find junctions
   {
     TSvtSRBamRecord srBR(2 * DELLY_SVT_TRANS, TSRBamRecord());
     findJunctions(c, srBR);
   
     // Debug
     if (c.hasDumpFile) outputSRBamRecords(c, srBR);

     // Cluster BAM records
     for(uint32_t svt = 0; svt < srBR.size(); ++svt) {
       if (srBR[svt].empty()) continue;
       
       // Sort
       std::sort(srBR[svt].begin(), srBR[svt].end(), SortSRBamRecord<SRBamRecord>());
       
       // Cluster
       cluster(c, srBR[svt], srSVs, c.maxReadSep, svt);

       // Track split-reads
       for(uint32_t i = 0; i < srBR[svt].size(); ++i) {
	 // Read assigned?
	 if ((srBR[svt][i].svid != -1) && (srBR[svt][i].rstart != -1)) {
	   if (srBR[svt][i].rstart < (int32_t) c.chrlen[srBR[svt][i].chr]) srStore[srBR[svt][i].chr].insert(std::make_pair(std::make_pair(srBR[svt][i].rstart, srBR[svt][i].id), srBR[svt][i].svid));
	   if (srBR[svt][i].chr != srBR[svt][i].chr2) {
	     // Unclear which chr was primary alignment so insert both if and only if rstart < reference length
	     if (srBR[svt][i].rstart < (int32_t) c.chrlen[srBR[svt][i].chr2]) srStore[srBR[svt][i].chr2].insert(std::make_pair(std::make_pair(srBR[svt][i].rstart, srBR[svt][i].id), srBR[svt][i].svid));
	   }
	 }
       }
     }
   }

   // Sort SVs
   std::sort(srSVs.begin(), srSVs.end(), SortSVs<StructuralVariantRecord>());
   
   // VCF Output
   outputVcf(c, srSVs);
   
   return 0;
 }

 int teguaMain(int argc, char **argv) {
   TeguaConfig c;
   
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
   
   boost::program_options::options_description hidden("Hidden options");
   hidden.add_options()
     ("input-file", boost::program_options::value<boost::filesystem::path>(&c.inputfile), "input bam file")
     ("sv,s", boost::program_options::value<boost::filesystem::path>(&c.svinput), "input SV file, format: chr, pos, chr2, pos2 (optional)")
     ("dump,d", boost::program_options::value<boost::filesystem::path>(&c.dumpfile), "gzipped read output file")
     ("pruning,j", boost::program_options::value<uint32_t>(&c.graphPruning)->default_value(1000), "graph pruning cutoff")
     ;
   
   boost::program_options::positional_options_description pos_args;
   pos_args.add("input-file", -1);
   
   boost::program_options::options_description cmdline_options;
   cmdline_options.add(generic).add(disc).add(hidden);
   boost::program_options::options_description visible_options;
   visible_options.add(generic).add(disc);
   boost::program_options::variables_map vm;
   boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
   boost::program_options::notify(vm);
   
   // Check command line arguments
   if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
     std::cout << "Usage: tegua " << argv[0] << " [OPTIONS] -g <ref.fa> <input.bam>" << std::endl;
     std::cout << visible_options << "\n";
    return 1;
   }

   // SV input
   if (vm.count("svinput")) c.hasSvInput = true;
   else c.hasSvInput = false;

   // Dump reads
   if (vm.count("dump")) c.hasDumpFile = true;
   else c.hasDumpFile = false;
   
   // Check input BAM
   {
     if (!(boost::filesystem::exists(c.inputfile) && boost::filesystem::is_regular_file(c.inputfile) && boost::filesystem::file_size(c.inputfile))) {
       std::cerr << "Alignment file is missing: " << c.inputfile.string() << std::endl;
       return 1;
     }
     samFile* samfile = sam_open(c.inputfile.string().c_str(), "r");
     if (samfile == NULL) {
       std::cerr << "Fail to open file " << c.inputfile.string() << std::endl;
       return 1;
     }
     hts_idx_t* idx = sam_index_load(samfile, c.inputfile.string().c_str());
     if (idx == NULL) {
       std::cerr << "Fail to open index for " << c.inputfile.string() << std::endl;
       return 1;
     }
     bam_hdr_t* hdr = sam_hdr_read(samfile);
     if (hdr == NULL) {
       std::cerr << "Fail to open header for " << c.inputfile.string() << std::endl;
       return 1;
     }
     c.chrlen.resize(hdr->n_targets);
     for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) c.chrlen[refIndex] = hdr->target_len[refIndex];
     bam_hdr_destroy(hdr);
     hts_idx_destroy(idx);
     sam_close(samfile);
   }
   c.nchr = c.chrlen.size();
   
   // Run Tegua
   return runTegua(c);
 }

}

#endif
