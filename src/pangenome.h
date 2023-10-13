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
    uint16_t minMapQual;
    uint32_t minClip;
    uint32_t minRefSep;
    uint32_t maxReadSep;
    uint32_t graphPruning;
    uint32_t minCliqueSize;
    uint32_t maxReadPerSV;
    int32_t nchr;
    int32_t minimumFlankSize;
    float flankQuality;
    std::set<int32_t> svtset;
    DnaScore<int> aliscore;
    boost::filesystem::path outfile;
    boost::filesystem::path fastqfile;
    std::vector<boost::filesystem::path> files;
    boost::filesystem::path genome;
    boost::filesystem::path exclude;
    std::vector<std::string> sampleName;
  };

  template<typename TReadBp>
  inline void
  _insertGraphJunction(TReadBp& readBp, std::size_t const seed, AlignRecord const& ar, uint32_t const pathidx, int32_t const rp, int32_t const sp, bool const scleft) {
    typedef typename TReadBp::mapped_type TJunctionVector;

    bool fw = ar.path[pathidx].first;
    int32_t readStart = ar.qstart; // Query start (not needed)
    if (sp <= ar.qlen) {
      typename TReadBp::iterator it = readBp.find(seed);
      if (it != readBp.end()) it->second.push_back(Junction(fw, scleft, ar.path[pathidx].second, readStart, rp, sp, ar.mapq));
      else readBp.insert(std::make_pair(seed, TJunctionVector(1, Junction(fw, scleft, ar.path[pathidx].second, readStart, rp, sp, ar.mapq))));
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
	    uint32_t seqlen = g.nodelen(ar.path[pi].second);
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

	    // Leading junction
	    if ((pi == 0) && (sp > c.minRefSep)) {
	      int32_t locbeg = pstart + 1 + srp;
	      if (!ar.path[pi].first) locbeg = pstart + 1 + (srpend - srp);
	      if ((locbeg > 0) && (locbeg < (int32_t) seqlen)) {
		//std::cerr << seqname << '\t' << locbeg << "\tRead\t" << qname << '\t' << ar.seed << '\t' << ar.qlen << "\tPath\t" << pi << '\t' << (int32_t) ar.path[pi].first << '\t' << ar.pstart << "\tReadBp\t" << sp << std::endl;
		_insertGraphJunction(readBp, ar.seed, ar, pi, locbeg, sp, ar.path[pi].first);
	      }
	    }
	    // Internal junctions
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
		// Insert start-junction
		if (ar.cigarlen[i] > c.minRefSep) {
		  if ((rp >= refstart) && (rp < refend)) {
		    int32_t locbeg = pstart + 1 + srp;
		    if (!ar.path[pi].first) locbeg = pstart + 1 + (srpend - srp);
		    if ((locbeg > 0) && (locbeg < (int32_t) seqlen)) {
		      //std::cerr << seqname << '\t' << locbeg << "\tRead\t" << qname << '\t' << ar.seed << '\t' << ar.qlen << "\tPath\t" << pi << '\t' << (int32_t) ar.path[pi].first << '\t' << ar.pstart << "\tReadBp\t" << sp << '\t' << ar.cigarlen[i] << std::endl;
		      _insertGraphJunction(readBp, ar.seed, ar, pi, locbeg, sp, !ar.path[pi].first);
		    }
		  }
		}
		sp += ar.cigarlen[i];
		// Insert end-junction
		if (ar.cigarlen[i] > c.minRefSep) {
		  if ((rp >= refstart) && (rp < refend)) {
		    int32_t locbeg = pstart + 1 + srp;
		    if (!ar.path[pi].first) locbeg = pstart + 1 + (srpend - srp);
		    if ((locbeg > 0) && (locbeg < (int32_t) seqlen)) {
		      //std::cerr << seqname << '\t' << locbeg << "\tRead\t" << qname << '\t' << ar.seed << '\t' << ar.qlen << "\tPath\t" << pi << '\t' << (int32_t) ar.path[pi].first << '\t' << ar.pstart << "\tReadBp\t" << sp << '\t' << ar.cigarlen[i] << std::endl;
		      _insertGraphJunction(readBp, ar.seed, ar, pi, locbeg, sp, ar.path[pi].first);
		    }
		  }
		}
	      }
	      else {
		std::cerr << "Warning: Unknown Cigar option " << ar.cigarop[i] << std::endl;
		return false;
	      }
	    }
	    // Trailing junction
	    if ((pi + 1 == ar.path.size()) && ((int32_t) (sp + c.minRefSep) < ar.qlen)) {
	      int32_t locbeg = pstart + 1 + srp;
	      if (!ar.path[pi].first) locbeg = pstart + 1 + (srpend - srp);
	      if ((locbeg > 0) && (locbeg < (int32_t) seqlen)) {
		//std::cerr << seqname << '\t' << locbeg << "\tRead\t" << qname << '\t' << ar.seed << '\t' << ar.qlen << "\tPath\t" << pi << '\t' << (int32_t) ar.path[pi].first << '\t' << ar.pstart << "\tReadBp\t" << sp << std::endl;
		_insertGraphJunction(readBp, ar.seed, ar, pi, locbeg, sp, !ar.path[pi].first);
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
    }
    
    // Header
    std::cerr << "id\tsegment1\tpos1\tsegment2\tpos2\tsvtype\tct\tqual\tinslen" << std::endl;
    
    // SVs
    for(uint32_t svt = 0; svt < br.size(); ++svt) {
      for(uint32_t i = 0; i < br[svt].size(); ++i) {
	std::cerr << hm[br[svt][i].id] << '\t' << br[svt][i].id << '\t' << idSegment[br[svt][i].chr] << '\t' << br[svt][i].pos << '\t' << idSegment[br[svt][i].chr2] << '\t' << br[svt][i].pos2 << '\t' << _addID(svt) << '\t' << _addOrientation(svt) << '\t' << br[svt][i].qual << '\t' << br[svt][i].inslen << std::endl;
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

  inline void
  outputGraphStructuralVariants(Graph const& g, std::vector<StructuralVariantRecord> const& svs) {
    // Header
    std::cerr << "segment1\tpos1\tsegment2\tpos2\tid\tsvtype\tct\tpeSupport\tsrSupport\tconsensus" << std::endl;

    // Vertex map
    std::vector<std::string> idSegment(g.smap.size());
    for(typename Graph::TSegmentIdMap::const_iterator it = g.smap.begin(); it != g.smap.end(); ++it) idSegment[it->second] = it->first;
    
    // SVs
    for(uint32_t i = 0; i < svs.size(); ++i) {
      std::string idname(_addID(svs[i].svt));
      std::string padNumber = boost::lexical_cast<std::string>(i);
      padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
      idname += padNumber;
      std::cerr << idSegment[svs[i].chr] << '\t' << svs[i].svStart << '\t' << idSegment[svs[i].chr2] << '\t' << svs[i].svEnd << '\t' << idname << '\t' << _addID(svs[i].svt) << '\t' << _addOrientation(svs[i].svt) << '\t' << svs[i].peSupport << '\t' << svs[i].srSupport << '\t' << svs[i].consensus << std::endl;
    }
  }

  template<typename TConfig, typename TSRStore>
  inline void
  _clusterGraphSRReads(TConfig c, Graph const& g, std::vector<StructuralVariantRecord>& svc, TSRStore& srStore) {
    // Split-reads
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] SV discovery" << std::endl;
    typedef std::vector<SRBamRecord> TSRBamRecord;
    typedef std::vector<TSRBamRecord> TSvtSRBamRecord;
    TSvtSRBamRecord srBR(2 * DELLY_SVT_TRANS, TSRBamRecord());
    _findGraphSRBreakpoints(c, g, srBR);

    // Debug
    //outputGraphSRBamRecords(c, g, srBR);

    // Cluster BAM records
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Cluster SVs" << std::endl;
    for(uint32_t svt = 0; svt < srBR.size(); ++svt) {
      if (srBR[svt].empty()) continue;

      // Sort
      std::sort(srBR[svt].begin(), srBR[svt].end(), SortSRBamRecord<SRBamRecord>());

      // Cluster
      cluster(c, srBR[svt], svc, svt);

      // Debug
      //outputGraphStructuralVariants(c, g, svc, srBR, svt);
      
      // Track split-reads
      for(uint32_t i = 0; i < srBR[svt].size(); ++i) {
	// Read assigned?
	if ((srBR[svt][i].svid != -1) && (srBR[svt][i].rstart != -1)) {
	  if (srStore.find(srBR[svt][i].id) == srStore.end()) srStore.insert(std::make_pair(srBR[svt][i].id, std::vector<SeqSlice>()));
	  srStore[srBR[svt][i].id].push_back(SeqSlice(srBR[svt][i].svid, srBR[svt][i].sstart, srBR[svt][i].inslen, srBR[svt][i].qual));
	}
      }
    }
  }


  template<typename TConfig, typename TSRStore>
  inline void
  assembleGraph(TConfig const& c, Graph const& g, std::vector<StructuralVariantRecord>& svs, TSRStore& srStore) {
    // Assembly
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Split-read assembly" << std::endl;
    
    // Sequence store
    typedef std::vector<std::string> TSequences;
    typedef std::vector<TSequences> TSVSequences;
    TSVSequences seqStore(svs.size(), TSequences());

    // SV consensus done
    std::vector<bool> svcons(svs.size(), false);

    // Vertex map
    std::vector<std::string> idSegment(g.smap.size());
    for(typename Graph::TSegmentIdMap::const_iterator it = g.smap.begin(); it != g.smap.end(); ++it) idSegment[it->second] = it->first;

    // Load FASTQ
    std::ifstream fqfile;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    if (is_gz(c.fastqfile)) {
      fqfile.open(c.fastqfile.string().c_str(), std::ios_base::in | std::ios_base::binary);
      dataIn.push(boost::iostreams::gzip_decompressor(), 16*1024);
    } else fqfile.open(c.fastqfile.string().c_str(), std::ios_base::in);
    dataIn.push(fqfile);
    std::istream instream(&dataIn);
    std::string gline;
    uint64_t lnum = 0;
    std::string qname;
    std::string sequence;
    bool validRec = true;
    while(std::getline(instream, gline)) {
      if (lnum % 2 == 0) {
	// FASTA or FASTQ
	if ((gline[0] == '>') || (gline[0] == '@')) {
	  validRec = true;
	  qname = gline.substr(1);
	  qname = qname.substr(0, qname.find(' '));
	  qname = qname.substr(0, qname.find('\t'));
	} else validRec = false;
      } else if (lnum % 2 == 1) {
	if (validRec) {
	  std::size_t seed = hash_lr(qname);
	  if (srStore.find(seed) != srStore.end()) {
	    sequence = gline;
	    int32_t readlen = sequence.size();

	    // Iterate all spanned SVs
	    for(uint32_t ri = 0; ri < srStore[seed].size(); ++ri) {
	      SeqSlice seqsl = srStore[seed][ri];
	      int32_t svid = seqsl.svid;

	      // Debug SV read
	      //std::cerr << "SV:" << svid << '\t' << idSegment[svs[svid].chr] << '\t' << svs[svid].svStart << '\t' << idSegment[svs[svid].chr2] << '\t' << svs[svid].svEnd << '\t' << _addID(svs[svid].svt) << '\t' << _addOrientation(svs[svid].svt) << '\t' << svs[svid].srSupport << '\t' << qname << '\t' << seqsl.sstart << '\t' << sequence.substr(std::max(0, seqsl.sstart - 20), (std::min(seqsl.sstart + seqsl.inslen + 20, readlen) - std::max(0, seqsl.sstart - 20))) << std::endl;

	      
	      if ((!svcons[svid]) && (seqStore[svid].size() < c.maxReadPerSV)) {
		// Extract subsequence (otherwise MSA takes forever)
		int32_t window = 1000; // MSA should be larger
		int32_t sPos = seqsl.sstart - window;
		int32_t ePos = seqsl.sstart + seqsl.inslen + window;
		if (sPos < 0) sPos = 0;
		if (ePos > (int32_t) readlen) ePos = readlen;
		// Min. seq length and max insertion size, 10kbp?
		if ((ePos - sPos) > window) {
		  seqStore[svid].push_back(sequence.substr(sPos, (ePos - sPos)));

		  // Enough split-reads?
		  if ((seqStore[svid].size() == c.maxReadPerSV) || ((int32_t) seqStore[svid].size() == svs[svid].srSupport)) {
		    if (seqStore[svid].size() > 1) {
		      if (svs[svid].svt != 4) msaEdlib(c, seqStore[svid], svs[svid].consensus);
		      else msaWfa(c, seqStore[svid], svs[svid].consensus);
		      // Check  svs[svid].consensus.size() < svs[svid].insLen + 4 * c.minConsWindow to see if it worked

		      // Debug
		      //std::string idname(_addID(svs[svid].svt));
		      //std::string padNumber = boost::lexical_cast<std::string>(svid);
		      //padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
		      //idname += padNumber;
		      //std::cerr << "SV\t" << idname << '\t' << idSegment[svs[svid].chr] << '\t' << svs[svid].svStart << '\t' << idSegment[svs[svid].chr2] << '\t' << svs[svid].svEnd << '\t' << _addID(svs[svid].svt) << '\t' << _addOrientation(svs[svid].svt) << '\t' << svs[svid].srSupport << '\t' << svs[svid].consensus << std::endl;
		    }
		    seqStore[svid].clear();
		    svcons[svid] = true;
		  }
		}
	      }
	    }
	  }
	}
      }
      ++lnum;
    }
    // Clean-up
    dataIn.pop();
    if (is_gz(c.fastqfile)) dataIn.pop();
    fqfile.close();

    // Handle left-overs
    for(uint32_t svid = 0; svid < svcons.size(); ++svid) {
      if (!svcons[svid]) {
	if (seqStore[svid].size() > 1) {
	  if (svs[svid].svt != 4) msaEdlib(c, seqStore[svid], svs[svid].consensus);
	  else msaWfa(c, seqStore[svid], svs[svid].consensus);
	}
	seqStore[svid].clear();
	svcons[svid] = true;

	// Debug
	//std::string idname(_addID(svs[svid].svt));
	//std::string padNumber = boost::lexical_cast<std::string>(svid);
	//padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
	//idname += padNumber;
	//std::cerr << "SV\t" << idname << '\t' << idSegment[svs[svid].chr] << '\t' << svs[svid].svStart << '\t' << idSegment[svs[svid].chr2] << '\t' << svs[svid].svEnd << '\t' << _addID(svs[svid].svt) << '\t' << _addOrientation(svs[svid].svt) << '\t' << svs[svid].srSupport << '\t' << svs[svid].consensus << std::endl;
      }
    }
    
    // Clean-up unfinished SVs
    for(uint32_t svid = 0; svid < svcons.size(); ++svid) {
      if (!svcons[svid]) {
	//std::cerr << "Missing: " << svid << ',' << svs[svid].svt << std::endl;
	svs[svid].consensus = "";
	svs[svid].srSupport = 0;
	svs[svid].srAlignQuality = 0;
      }
    }
  }

  inline void
  _fillPrefix(Graph& g, uint32_t const nodeid, std::string prefix, std::vector<std::string>& all, int32_t const reqlen) {
    for(uint32_t i = 0; i < g.links.size(); ++i) {
      if ((g.links[i].to == nodeid) && (g.links[i].tofwd)) {
	if (g.links[i].fromfwd) {
	  std::cerr << "PrefixLink\t" << g.links[i].from << ',' << (int) g.links[i].fromfwd << ':' << g.links[i].to << ',' << (int) g.links[i].tofwd << ':' << g.nodelen(g.links[i].from) << ',' << reqlen << std::endl;
	  int32_t localstart = g.nodelen(g.links[i].from)-reqlen;
	  if (localstart < 0) {
	    // Recurse
	    _fillPrefix(g, g.links[i].from, g.nodeseq(g.links[i].from) + prefix, all, reqlen - g.nodelen(g.links[i].from));
	  } else {
	    all.push_back(g.nodeseq(g.links[i].from).substr(localstart) + prefix);
	  }
	}
      }
    }
  }


  inline void
  _fillSuffix(Graph& g, uint32_t const nodeid, std::string suffix, std::vector<std::string>& all, int32_t const reqlen) {
    for(uint32_t i = 0; i < g.links.size(); ++i) {
      if ((g.links[i].from == nodeid) && (g.links[i].fromfwd)) {
	if (g.links[i].tofwd) {
	  std::cerr << "SuffixLink\t" << g.links[i].from << ',' << (int) g.links[i].fromfwd << ':' << g.links[i].to << ',' << (int) g.links[i].tofwd << ':' << g.nodelen(g.links[i].to) << ',' << reqlen << std::endl;
	  int32_t localend = reqlen;
	  if (localend > (int) g.nodelen(g.links[i].to)) {
	    // Recurse
	    _fillSuffix(g, g.links[i].to, suffix + g.nodeseq(g.links[i].to), all, reqlen - g.nodelen(g.links[i].to));
	  } else {
	    all.push_back(suffix + g.nodeseq(g.links[i].to).substr(0, localend));
	  }	  
	}
      }
    }
  }

  template<typename TConfig>
  inline void
  alignToGraph(TConfig const& c, Graph& g, std::vector<StructuralVariantRecord>& svs) {
    // Generate pairwise alignment
    for(uint32_t svid = 0; svid < svs.size(); ++svid) {
      std::cerr << "SV:" << svid << ',' << svs[svid].svStart << ',' << svs[svid].svEnd << '\t' << svs[svid].svt << ',' << svs[svid].consensus.size() << std::endl;
      bool validConsensusAlignment = false;
      AlignDescriptor ad;
      if ( (int32_t) svs[svid].consensus.size() >= (2 * c.minimumFlankSize + svs[svid].insLen)) {
	if (svs[svid].svt == 2) { 	// Deletion
	  uint32_t svsize = svs[svid].svEnd - svs[svid].svStart;
	  if (svsize < svs[svid].consensus.size()) {
	    int32_t startpos = svs[svid].svStart - svs[svid].consensus.size();
	    std::vector<std::string> prefix;
	    if (startpos < 0) {
	      startpos = 0;
	      _fillPrefix(g, svs[svid].chr, "", prefix, svs[svid].consensus.size() - svs[svid].svStart);
	    }
	    for(uint32_t i = 0; i < prefix.size(); ++i) std::cerr << "Prefix: " << prefix[i] << std::endl;
	    int32_t endpos = svs[svid].svEnd + svs[svid].consensus.size();
	    std::vector<std::string> suffix;
	    if (endpos > (int) g.nodelen(svs[svid].chr)) {
	      endpos = g.nodelen(svs[svid].chr);
	      _fillSuffix(g, svs[svid].chr, "", suffix, svs[svid].consensus.size() - (g.nodelen(svs[svid].chr) - svs[svid].svEnd));
	    }
	    for(uint32_t i = 0; i < prefix.size(); ++i) std::cerr << "Suffix: " << suffix[i] << std::endl;
	    std::string svRefStr = g.nodeseq(svs[svid].chr).substr(startpos, (endpos - startpos));
	    std::cerr << "SV:" << svid << ',' << svRefStr << std::endl;
	    std::cerr << "SV:" << svid << ',' << svs[svid].consensus << std::endl;
	    if (_alignConsensus(c, svs[svid].consensus, svRefStr, svs[svid].svt, ad, true)) {
	      validConsensusAlignment = true;
	      std::cerr << "Valid alignment" << std::endl;
	    }
	  } else {
	    int32_t startpos1 = svs[svid].svStart - svs[svid].consensus.size();
	    if (startpos1 < 0) startpos1 = 0;
	    int32_t endpos1 = svs[svid].svStart + svs[svid].consensus.size();
	    if (endpos1 > (int) g.nodelen(svs[svid].chr)) endpos1 = g.nodelen(svs[svid].chr);
	    std::string svRefStr = g.nodeseq(svs[svid].chr).substr(startpos1, (endpos1 - startpos1));
	    int32_t startpos2 = svs[svid].svEnd - svs[svid].consensus.size();
	    if (startpos2 < 0) startpos2 = 0;
	    int32_t endpos2 = svs[svid].svEnd + svs[svid].consensus.size();
	    if (endpos2 > (int) g.nodelen(svs[svid].chr)) endpos2 = g.nodelen(svs[svid].chr);
	    svRefStr += g.nodeseq(svs[svid].chr).substr(startpos2, (endpos2 - startpos2));
	    if (_alignConsensus(c, svs[svid].consensus, svRefStr, svs[svid].svt, ad, true)) {
	      validConsensusAlignment = true;
	      std::cerr << "Valid alignment" << std::endl;
	    }	    
	  }
	}
      }
      if (!validConsensusAlignment) {
	// Requires path alignments
	//svs[svid].consensus = "";
	//svs[svid].srSupport = 0;
	//svs[svid].srAlignQuality = 0;
      } else {
	svs[svid].precise = true;
	//sv.svStart=finalGapStart;
	//sv.svEnd=finalGapEnd;
	svs[svid].srAlignQuality = ad.percId;
	svs[svid].insLen=ad.cEnd - ad.cStart - 1;
	svs[svid].homLen=std::max(0, ad.homLeft + ad.homRight - 2);
	int32_t ci_wiggle = std::max(ad.homLeft, ad.homRight);
	svs[svid].ciposlow = -ci_wiggle;
	svs[svid].ciposhigh = ci_wiggle;
	svs[svid].ciendlow = -ci_wiggle;
	svs[svid].ciendhigh = ci_wiggle;
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
   if (svs.empty()) {
       
     // Structural Variant Candidates
     typedef std::vector<StructuralVariantRecord> TVariants;
     TVariants svc;

     // Load pan-genome graph
     std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Load pan-genome graph" << std::endl;
     Graph g;
     parseGfa(c, g);
     c.nchr = g.smap.size();

     // Split-read store
     typedef std::vector<SeqSlice> TSvPosVector;
     typedef boost::unordered_map<std::size_t, TSvPosVector> TReadSV;
     TReadSV srStore;

     // SV Discovery
     _clusterGraphSRReads(c, g, svc, srStore);

     // Assemble
     assembleGraph(c, g, svc, srStore);
     srStore.clear();

     // Sort SVs
     sort(svc.begin(), svc.end(), SortSVs<StructuralVariantRecord>());
           
     // Consensus alignment
     //alignToGraph(c, g, svc);
     
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
     outputGraphStructuralVariants(g, svs);
   }

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
     ("fastq,x", boost::program_options::value<boost::filesystem::path>(&c.fastqfile), "input FASTA/FASTQ file")
     ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "BCF output file")
     ;
   
   boost::program_options::options_description disc("Discovery options");
   disc.add_options()
     ("mapqual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(10), "min. mapping quality")
     ("minclip,c", boost::program_options::value<uint32_t>(&c.minClip)->default_value(25), "min. clipping length")
     ("min-clique-size,z", boost::program_options::value<uint32_t>(&c.minCliqueSize)->default_value(3), "min. clique size")
     ("minrefsep,m", boost::program_options::value<uint32_t>(&c.minRefSep)->default_value(30), "min. graph separation")
     ("maxreadsep,n", boost::program_options::value<uint32_t>(&c.maxReadSep)->default_value(100), "max. read separation")
     ;

   boost::program_options::options_description cons("Consensus options");
   cons.add_options()
     ("max-reads,p", boost::program_options::value<uint32_t>(&c.maxReadPerSV)->default_value(15), "max. reads for consensus computation")
     ("flank-size,f", boost::program_options::value<int32_t>(&c.minimumFlankSize)->default_value(100), "min. flank size")
     ("flank-quality,a", boost::program_options::value<float>(&c.flankQuality)->default_value(0.9), "min. flank quality")
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
   cmdline_options.add(generic).add(disc).add(cons).add(hidden);
   boost::program_options::options_description visible_options;
   visible_options.add(generic).add(disc).add(cons);
   boost::program_options::variables_map vm;
   boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
   boost::program_options::notify(vm);
   
   // Check command line arguments
   if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome")) || (!vm.count("fastq"))) {
     std::cerr << std::endl;
     std::cerr << "Usage: delly " << argv[0] << " [OPTIONS] -g <pan-genome.gfa.gz> -x <sample.fq.gz> <sample1.gaf.gz> <sample2.gaf.gz> ..." << std::endl;
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
