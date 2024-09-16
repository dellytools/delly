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

#include "gfa.h"
#include "gaf.h"
#include "util.h"

namespace torali {


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
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Split-read scanning" << std::endl;
    
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
      for(typename TReadBp::iterator it = readBp.begin(); it != readBp.end(); ++it) std::sort(it->second.begin(), it->second.end());
      
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

}

#endif
