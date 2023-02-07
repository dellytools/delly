#ifndef PANGENOME_H
#define PANGENOME_H

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
#include "gfa.h"
#include "gaf.h"

namespace torali {


  struct GraphConfig {
    bool hasDumpFile;
    bool hasVcfFile;
    uint16_t minMapQual;
    uint16_t minGenoQual;
    uint32_t minClip;
    uint32_t minRefSep;
    uint32_t maxReadSep;
    uint32_t graphPruning;
    uint32_t minCliqueSize;
    uint32_t maxReadPerSV;
    int32_t nchr;
    int32_t minimumFlankSize;
    int32_t indelsize;
    int32_t maxInsertionSize;
    float indelExtension;
    float flankQuality;
    std::set<int32_t> svtset;
    DnaScore<int> aliscore;
    boost::filesystem::path dumpfile;
    boost::filesystem::path outfile;
    boost::filesystem::path seqfile;
    boost::filesystem::path vcffile;
    std::vector<boost::filesystem::path> files;
    boost::filesystem::path genome;
    boost::filesystem::path exclude;
    std::vector<std::string> sampleName;
  };

  template<typename TReadBp>
  inline void
  _insertGraphJunction(TReadBp& readBp, std::size_t const seed, AlignRecord const& ar, uint32_t const pathidx, int32_t const rp, int32_t const sp, bool const scleft) {
    bool fw = ar.path[pathidx].first;
    int32_t readStart = ar.pstart; // Not needed I think
    typedef typename TReadBp::mapped_type TJunctionVector;
    typename TReadBp::iterator it = readBp.find(seed);
    if (sp <= ar.qlen) {
      if (!fw) {
	if (it != readBp.end()) it->second.push_back(Junction(fw, scleft, ar.path[pathidx].second, readStart, rp, ar.qlen - sp, ar.mapq));
	else readBp.insert(std::make_pair(seed, TJunctionVector(1, Junction(fw, scleft, ar.path[pathidx].second, readStart, rp, ar.qlen - sp, ar.mapq))));
      } else {
	if (it != readBp.end()) it->second.push_back(Junction(fw, scleft, ar.path[pathidx].second, readStart, rp, sp, ar.mapq));
	else readBp.insert(std::make_pair(seed, TJunctionVector(1, Junction(fw, scleft, ar.path[pathidx].second, readStart, rp, sp, ar.mapq))));
      }
    }
  }


  
  template<typename TConfig, typename TReadBp>
  inline bool
  findGraphJunctions(TConfig& c, Graph const& g, TReadBp& readBp) {
    // Vertex map
    std::vector<std::string> idSegment(g.smap.size());
    for(typename Graph::TSegmentIdMap::const_iterator it = g.smap.begin(); it != g.smap.end(); ++it) idSegment[it->second] = it->first;

    // Iterate graph alignments
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      
      // Open GAF
      std::ifstream gafFile;
      boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
      if (is_gz(c.files[file_c])) {
	gafFile.open(c.files[file_c].string().c_str(), std::ios_base::in | std::ios_base::binary);
	dataIn.push(boost::iostreams::gzip_decompressor(), 16*1024);
      } else gafFile.open(c.files[file_c].string().c_str(), std::ios_base::in);
      dataIn.push(gafFile);

      // Parse GAF
      std::istream instream(&dataIn);
      bool parseAR = true;
      while (parseAR) {
	AlignRecord ar;
	std::string qname;
	if (parseAlignRecord(instream, g, ar, qname)) {
	  if (ar.mapq < c.minMapQual) continue;

	  // Iterate path alignments
	  uint32_t refstart = 0;
	  for(uint32_t pi = 0; pi < ar.path.size(); ++pi) {
	    // Vertex coordinates
	    std::string seqname = idSegment[ar.path[pi].second];
	    uint32_t seqlen = g.segments[ar.path[pi].second].len;
	    uint32_t pstart = 0;
	    uint32_t plen = seqlen;
	    if (pi == 0) {
	      plen -= ar.pstart;
	      if (ar.path[pi].first) pstart = ar.pstart;
	    }
	    if (pi + 1 == ar.path.size()) {
	      plen = ar.pend - ar.pstart - refstart;
	      if (!ar.path[pi].first) {
		if (pi == 0) pstart = seqlen - ar.pend;
		else pstart = ar.pstart + refstart + seqlen - ar.pend;
	      }
	    }

	    // Compute local alignment end
	    uint32_t refend = refstart + plen;
	    uint32_t rp = 0;  // Reference pointer
	    uint32_t srpend = 0; // Segment reference pointer
	    for (uint32_t i = 0; i < ar.cigarop.size(); ++i) {
	      if ((ar.cigarop[i] == BAM_CMATCH) || (ar.cigarop[i] == BAM_CEQUAL) || (ar.cigarop[i] == BAM_CDIFF)) {
		for(uint32_t k = 0; k < ar.cigarlen[i]; ++k, ++rp) {
		  if ((rp >= refstart) && (rp < refend)) ++srpend;
		}
	      }
	      else if (ar.cigarop[i] == BAM_CDEL) {
		for(uint32_t k = 0; k < ar.cigarlen[i]; ++k, ++rp) {
		  if ((rp >= refstart) && (rp < refend)) ++srpend;
		}
	      }
	    }
	    
	    // Parse CIGAR
	    rp = 0;  // Reference pointer
	    uint32_t srp = 0;
	    uint32_t sp = ar.qstart;
	    for (uint32_t i = 0; i < ar.cigarop.size(); ++i) {
	      if ((ar.cigarop[i] == BAM_CMATCH) || (ar.cigarop[i] == BAM_CEQUAL) || (ar.cigarop[i] == BAM_CDIFF)) {
		for(uint32_t k = 0; k < ar.cigarlen[i]; ++k, ++sp, ++rp) {
		  if ((rp >= refstart) && (rp < refend)) ++srp;
		}
	      }
	      else if (ar.cigarop[i] == BAM_CDEL) {
		// Insert start-junction
		if (ar.cigarlen[i] > c.minRefSep) {
		  if ((rp >= refstart) && (rp < refend)) {
		    int32_t locbeg = pstart + 1 + srp;
		    if (!ar.path[pi].first) locbeg = pstart + 1 + (srpend - srp - ar.cigarlen[i]);
		    if ((locbeg > 0) && (locbeg < (int32_t) seqlen)) {
		      //std::cerr << seqname << '\t' << locbeg << "\tRead\t" << qname << '\t' << ar.seed << '\t' << ar.qlen << "\tPath\t" << pi << '\t' << (int32_t) ar.path[pi].first << '\t' << ar.pstart << "\tReadBp\t" << sp << '\t' << ar.cigarlen[i] << std::endl;
		      _insertGraphJunction(readBp, ar.seed, ar, pi, locbeg, sp, false);
		    }
		  }
		}
		// Adjust segment reference pointer for deletion
		for(uint32_t k = 0; k < ar.cigarlen[i]; ++k, ++rp) {
		  if ((rp >= refstart) && (rp < refend)) ++srp;
		}
		// Insert end-junction
		if (ar.cigarlen[i] > c.minRefSep) {
		  if ((rp >= refstart) && (rp < refend)) {
		    int32_t locbeg = pstart + 1 + srp;
		    if (!ar.path[pi].first) locbeg = pstart + 1 + (srpend - srp) + ar.cigarlen[i];
		    if ((locbeg > 0) && (locbeg < (int32_t) seqlen)) {
		      //std::cerr << seqname << '\t' << locbeg << "\tRead\t" << qname << '\t' << ar.seed << '\t' << ar.qlen << "\tPath\t" << pi << '\t' << (int32_t) ar.path[pi].first << '\t' << ar.pstart << "\tReadBp\t" << sp << '\t' << ar.cigarlen[i] << std::endl;
		      _insertGraphJunction(readBp, ar.seed, ar, pi, locbeg, sp, true);
		    }
		  }
		}
	      }
	      else if (ar.cigarop[i] == BAM_CINS) {
		sp += ar.cigarlen[i];
	      }
	      else {
		std::cerr << "Warning: Unknown Cigar option " << ar.cigarop[i] << std::endl;
		return false;
	      }
	    }
	
	    // Next segment
	    refstart = refend;
	  }
	  
	  //std::cerr << ar.seed << ',' << ar.qlen << ',' << ar.qstart << ',' << ar.qend << ',' << ar.strand << ',' << ar.plen << ',' << ar.pstart << ',' << ar.pend << ',' << ar.matches << ',' << ar.alignlen << ',' << ar.mapq << std::endl;
	} else parseAR = false;
      }

      // Sort junctions
      for(typename TReadBp::iterator it = readBp.begin(); it != readBp.end(); ++it) {
	std::sort(it->second.begin(), it->second.end(), SortJunction<Junction>());
      }
      
      // Close file
      dataIn.pop();
      if (is_gz(c.files[file_c])) dataIn.pop();
      gafFile.close();
    }
      
    return true;
  }

  template<typename TConfig, typename TSvtSRBamRecord>
  inline void
  _findGraphSRBreakpoints(TConfig const& c, Graph const& g, TSvtSRBamRecord& srBR) {
    // Breakpoints
    typedef std::vector<Junction> TJunctionVector;
    typedef std::map<std::size_t, TJunctionVector> TReadBp;
    TReadBp readBp;
    findGraphJunctions(c, g, readBp);
    fetchSVs(c, readBp, srBR);
  }


  template<typename TConfig>
  inline void
  outputGraphSRBamRecords(TConfig const& c, Graph const& g, std::vector<std::vector<SRBamRecord> > const& br) {
    // Vertex map
    std::vector<std::string> idSegment(g.smap.size());
    for(typename Graph::TSegmentIdMap::const_iterator it = g.smap.begin(); it != g.smap.end(); ++it) idSegment[it->second] = it->first;
    
    // Header
    std::cerr << "id\tsegment1\tpos1\tsegment2\tpos2\tsvtype\tct\tqual\tinslen" << std::endl;
    
    // SVs
    for(uint32_t svt = 0; svt < br.size(); ++svt) {
      for(uint32_t i = 0; i < br[svt].size(); ++i) {
	std::cerr << br[svt][i].id << '\t' << idSegment[br[svt][i].chr] << '\t' << br[svt][i].pos << '\t' << idSegment[br[svt][i].chr2] << '\t' << br[svt][i].pos2 << '\t' << _addID(svt) << '\t' << _addOrientation(svt) << '\t' << br[svt][i].qual << '\t' << br[svt][i].inslen << std::endl;
      }
    }
  }




  template<typename TConfig, typename TSvtSRBamRecord>  
  inline void
  outputGraphStructuralVariants(TConfig const& c, Graph const& g, std::vector<StructuralVariantRecord> const& svs, TSvtSRBamRecord const& srBR, int32_t const svt) {
    // Header
    std::cerr << "segment1\tpos1\tsegment2\tpos2\tsvtype\tct\tpeSupport\tsrSupport" << std::endl;
    
    // Hash reads
    typedef std::map<std::size_t, std::string> THashMap;
    THashMap hm;
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {

      // Open GAF
      std::ifstream gafFile;
      boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
      if (is_gz(c.files[file_c])) {
	gafFile.open(c.files[file_c].string().c_str(), std::ios_base::in | std::ios_base::binary);
	dataIn.push(boost::iostreams::gzip_decompressor(), 16*1024);
      } else gafFile.open(c.files[file_c].string().c_str(), std::ios_base::in);
      dataIn.push(gafFile);

      // Parse GAF
      std::istream instream(&dataIn);
      bool parseAR = true;
      while (parseAR) {
	AlignRecord ar;
	std::string qname;
	if (parseAlignRecord(instream, g, ar, qname)) {
	  if (hm.find(ar.seed) == hm.end()) hm.insert(std::make_pair(ar.seed, qname));
	  else {
	    if (hm[ar.seed] != qname) {
	      std::cerr << "Warning: Hash collision! " << ar.seed << ',' << hm[ar.seed] << ',' << qname << std::endl;
	    }
	  }
	} else parseAR = false;
      }

      // Close file
      dataIn.pop();
      if (is_gz(c.files[file_c])) dataIn.pop();
      gafFile.close();

      // Track split-reads
      typedef std::vector<std::string> TReadNameVector;
      typedef std::vector<TReadNameVector> TSVReadNames;
      TSVReadNames svReadNames(svs.size(), TReadNameVector());
      for(uint32_t i = 0; i < srBR[svt].size(); ++i) {
	if (srBR[svt][i].svid != -1) {
	  svReadNames[srBR[svt][i].svid].push_back(hm[srBR[svt][i].id]);
	}
      }
      hm.clear();

      // Vertex map
      std::vector<std::string> idSegment(g.smap.size());
      for(typename Graph::TSegmentIdMap::const_iterator it = g.smap.begin(); it != g.smap.end(); ++it) idSegment[it->second] = it->first;
  
      // SVs
      for(uint32_t i = 0; i < svs.size(); ++i) {
	if (svs[i].svt != svt) continue;
	std::cerr << idSegment[svs[i].chr] << '\t' << svs[i].svStart << '\t' << idSegment[svs[i].chr2] << '\t' << svs[i].svEnd << '\t' << _addID(svs[i].svt) << '\t' << _addOrientation(svs[i].svt) << '\t' << svs[i].peSupport << '\t' << svs[i].srSupport << '\t';
	for(uint32_t k = 0; k < svReadNames[svs[i].id].size(); ++k) std::cerr << svReadNames[svs[i].id][k] << ',';
	std::cerr << std::endl;
      }
    }
  }

  
 template<typename TConfig>
 inline int32_t
 runGraph(TConfig& c) {

#ifdef PROFILE
   ProfilerStart("delly.prof");
#endif

   // Structural Variants
   typedef std::vector<StructuralVariantRecord> TVariants;
   TVariants svs;

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

     // Load pan-genome graph
     std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Load pan-genome graph" << std::endl;
     Graph g;
     parseGfa(c, g, false);
     c.nchr = g.smap.size();
     
     // Split-reads
     std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] SV discovery" << std::endl;
     typedef std::vector<SRBamRecord> TSRBamRecord;
     typedef std::vector<TSRBamRecord> TSvtSRBamRecord;
     TSvtSRBamRecord srBR(2 * DELLY_SVT_TRANS, TSRBamRecord());
     _findGraphSRBreakpoints(c, g, srBR);

     // Debug
     //outputGraphSRBamRecords(c, g, srBR);

     // Cluster BAM records
     for(uint32_t svt = 0; svt < srBR.size(); ++svt) {
       if (srBR[svt].empty()) continue;

       // Sort
       std::sort(srBR[svt].begin(), srBR[svt].end(), SortSRBamRecord<SRBamRecord>());

       // Cluster
       cluster(c, srBR[svt], svc, c.maxReadSep, svt);

       // Debug
       outputGraphStructuralVariants(c, g, svc, srBR, svt);
     }
       

     /*
     // Assemble
     assemble(c, validRegions, svc, srStore);

     // Sort SVs
     sort(svc.begin(), svc.end(), SortSVs<StructuralVariantRecord>());
      
     // Keep assembled SVs only
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
     */
   } else {
     //vcfParse(c, hdr, svs);
   }

   /*
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

   */
   
#ifdef PROFILE
   ProfilerStop();
#endif

   // End
   boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
   std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
  
   return 0;
 }

 int pg(int argc, char **argv) {
   GraphConfig c;
   
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
     ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "BCF output file")
     ;
   
   boost::program_options::options_description disc("Discovery options");
   disc.add_options()
     ("mapqual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(10), "min. mapping quality")
     ("minclip,c", boost::program_options::value<uint32_t>(&c.minClip)->default_value(25), "min. clipping length")
     ("min-clique-size,z", boost::program_options::value<uint32_t>(&c.minCliqueSize)->default_value(3), "min. clique size")
     ("extension,e", boost::program_options::value<float>(&c.indelExtension)->default_value(0.5), "indel extension penalty")
     ("minrefsep,m", boost::program_options::value<uint32_t>(&c.minRefSep)->default_value(30), "min. graph separation")
     ("maxreadsep,n", boost::program_options::value<uint32_t>(&c.maxReadSep)->default_value(100), "max. read separation")
     ;

   boost::program_options::options_description cons("Consensus options");
   cons.add_options()
     ("max-reads,p", boost::program_options::value<uint32_t>(&c.maxReadPerSV)->default_value(15), "max. reads for consensus computation")
     ("flank-size,f", boost::program_options::value<int32_t>(&c.minimumFlankSize)->default_value(100), "min. flank size")
     ("flank-quality,a", boost::program_options::value<float>(&c.flankQuality)->default_value(0.9), "min. flank quality")
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
     ("pruning,j", boost::program_options::value<uint32_t>(&c.graphPruning)->default_value(1000), "graph pruning cutoff")
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
     std::cerr << "Usage: delly " << argv[0] << " [OPTIONS] -g <pan-genome.gfa.gz> <sample1.gaf.gz> <sample2.gaf.gz> ..." << std::endl;
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
     std::cerr << "Pan-genome graph is missing: " << c.genome.string() << std::endl;
     return 1;
   }
   
   // Check input files
   c.sampleName.resize(c.files.size());
   for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
     if (!(boost::filesystem::exists(c.files[file_c]) && boost::filesystem::is_regular_file(c.files[file_c]) && boost::filesystem::file_size(c.files[file_c]))) {
       std::cerr << "Graph alignment file is missing: " << c.files[file_c].string() << std::endl;
       return 1;
     }
     c.sampleName[file_c] = c.files[file_c].stem().string();
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
   
   // Run pan-genome SV discovery
   return runGraph(c);
 }

}

#endif
