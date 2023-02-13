#ifndef GFA_H
#define GFA_H


#include <iostream>
#include <fstream>

#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>

#include "util.h"

namespace torali
{

  struct Link {
    bool fromfwd;
    bool tofwd;
    uint32_t from;
    uint32_t to;

    Link() {}
    Link(bool const fv, bool const tv, uint32_t const fr, uint32_t tos) : fromfwd(fv), tofwd(tv), from(fr), to(tos) {}
  };

  struct LinkCargo {
    bool fromfwd;
    bool tofwd;
    uint32_t from;
    uint32_t to;
    uint32_t support;
    uint32_t mapq;

    LinkCargo() {}
    LinkCargo(Link const lk) : fromfwd(lk.fromfwd), tofwd(lk.tofwd), from(lk.from), to(lk.to), support(0), mapq(0) {}
    LinkCargo(bool const fv, bool const tv, uint32_t const fr, uint32_t tos) : fromfwd(fv), tofwd(tv), from(fr), to(tos), support(0) {}
  };


  template<typename TLink>
  struct SortLinks : public std::binary_function<TLink, TLink, bool>
  {
    inline bool operator()(TLink const& l1, TLink const& l2) {
      return ((l1.from < l2.from) || ((l1.from==l2.from) && (l1.to < l2.to)));
    }
  };

  

  struct Graph {
    typedef std::map<std::string, uint32_t> TSegmentIdMap;
    
    std::vector<uint32_t> offset;
    std::vector<Link> links;
    TSegmentIdMap smap;
    std::string sequence;

    bool empty() const { return sequence.empty(); }
    uint32_t nodelen(uint32_t const id) const {
      // Node length
      if (id + 1 == offset.size()) return sequence.size() - offset[id];
      else if (id + 1 < offset.size()) return (offset[id+1] - offset[id]);
      else return 0;
    }
    std::string nodeseq(uint32_t const id) const {
      return sequence.substr(offset[id], nodelen(id));
    }
  };


  template<typename TConfig>
  inline bool
  parseGfa(TConfig const& c, Graph& g) {
    // Open GFA
    std::ifstream gfaFile;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    if (is_gz(c.genome)) {
      gfaFile.open(c.genome.string().c_str(), std::ios_base::in | std::ios_base::binary);
      dataIn.push(boost::iostreams::gzip_decompressor(), 16*1024);
    } else gfaFile.open(c.genome.string().c_str(), std::ios_base::in);
    dataIn.push(gfaFile);

    // Parse GFA
    uint32_t id_counter = 0;
    std::istream instream(&dataIn);
    std::string gline;
    while(std::getline(instream, gline)) {
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep("\t");
      Tokenizer tokens(gline, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      if (tokIter != tokens.end()) {
	// What element
	if (*tokIter == "#") continue;
	else if (*tokIter == "S") {
	  // Segment
	  ++tokIter;
	  if (tokIter != tokens.end()) {
	    // Name
	    std::string segname = *tokIter;
	    ++tokIter;
	    if (tokIter != tokens.end()) {
	      // Sequence
	      std::string sequence = *tokIter;
	      // New segment
	      g.offset.push_back(g.sequence.size());
	      g.sequence += sequence;
	      // Keep segment name <-> id relationship
	      g.smap.insert(std::make_pair(segname, id_counter));
	      ++id_counter;
	    } else {
	      std::cerr << "S segment lacks sequence information!" << std::endl;
	      return false;
	    }
	  } else {
	    std::cerr << "S line lacks segment name!" << std::endl;
	    return false;
	  }
	}
	else if (*tokIter == "L") {
	  // Link
	  ++tokIter;
	  if (tokIter != tokens.end()) {
	    // From
	    if (g.smap.find(*tokIter) == g.smap.end()) {
	      std::cerr << "Link with unknown from segment! " << *tokIter << std::endl;
	      return false;
	    }
	    uint32_t fromId = g.smap[*tokIter];
	    ++tokIter;
	    if (tokIter != tokens.end()) {
	      // FromOrient
	      bool fromfwd = true;
	      if (*tokIter == "-") fromfwd = false;
	      ++tokIter;
	      if (tokIter != tokens.end()) {
		// To
		if (g.smap.find(*tokIter) == g.smap.end()) {
		  std::cerr << "Link with unknown to segment! " << *tokIter << std::endl;
		  return false;
		}
		uint32_t toId = g.smap[*tokIter];
		++tokIter;
		if (tokIter != tokens.end()) {
		  // ToOrient
		  bool tofwd = true;
		  if (*tokIter == "-") tofwd = false;
		  ++tokIter;
		  if (tokIter != tokens.end()) {
		    // Overlap CIGAR
		    if (*tokIter != "0M") {
		      std::cerr << "Currently only 0M links are supported!" << std::endl;
		      return false;
		    }
		    g.links.push_back(Link(fromfwd, tofwd, fromId, toId));
		  }
		}
	      }
	    }
	  }
	} else {
	  // Todo
	  std::cerr << "Warning: Unknown line " << *tokIter << std::endl;
	  continue;
	}
      }
    }
    dataIn.pop();
    if (is_gz(c.genome)) dataIn.pop();
    gfaFile.close();

    // Graph statistics
    std::cerr << "Parsed: " << g.offset.size() << " segments, " << g.links.size() << " links" << std::endl;
    std::cerr << "Total sequence size: " << g.sequence.size() << std::endl;

    return true;
  }

  inline void
  writeGfa(Graph const& g) {
    // Vertex map
    std::vector<std::string> idSegment(g.smap.size());
    for(typename Graph::TSegmentIdMap::const_iterator it = g.smap.begin(); it != g.smap.end(); ++it) idSegment[it->second] = it->first;
    
    // Temporary output file
    std::string filename = "test.out.gfa";
    
    // Output rGFA
    std::ofstream sfile;
    sfile.open(filename.c_str());

    // Output segments
    for(uint32_t i = 0; i < g.offset.size(); ++i) {
      std::string seqid = idSegment[i];
      sfile << "S\t" << seqid;
      sfile << "\t" << g.nodeseq(i);
      sfile << std::endl;
    }

    // Output links
    for(uint32_t i = 0; i < g.links.size(); ++i) {
      sfile << "L\t" << idSegment[g.links[i].from];
      if (g.links[i].fromfwd) sfile << "\t+";
      else sfile << "\t-";
      sfile << "\t" << idSegment[g.links[i].to];
      if (g.links[i].tofwd) sfile << "\t+";
      else sfile << "\t-";
      sfile << "\t0M" << std::endl;
    }
    sfile.close();
  }

}

#endif
