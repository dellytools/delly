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

  #ifndef POS_UNDEF
  #define POS_UNDEF std::numeric_limits<uint32_t>::max()
  #endif
  

  struct Segment {
    uint32_t tid;
    uint32_t pos;
    uint32_t len;

    Segment(uint32_t const t, uint32_t const p, uint32_t const l) : tid(t), pos(p), len(l) {}
  };

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
    
    std::vector<Segment> segments;
    std::vector<Link> links;
    std::vector<std::string> chrnames;
    std::vector<uint32_t> ranks;
    TSegmentIdMap smap;

    bool empty() const { return segments.empty(); }
  };


  template<typename TConfig>
  inline bool
  parseGfa(TConfig const& c, Graph& g, bool const storeseq) {
    typedef std::map<uint32_t, uint32_t> TRankMap;
    TRankMap rmap;
    
    // Chromosome map
    typedef std::pair<uint32_t, uint32_t> TRankIdx;
    typedef std::map<std::string, TRankIdx> TChrMap;
    TChrMap chrmap; // Each chromosome is mapped to a rank and globally unique integer id
    
    // Segment FASTA sequences
    std::ofstream sfile;
    if (storeseq) {
      // Vertex coordinates, store node sequences
      boost::filesystem::remove(c.seqfile.string());
      boost::filesystem::remove(c.seqfile.string() + ".fai");
      sfile.open(c.seqfile.string().c_str());
    }
    uint64_t seqsize = 0;

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
	      std::string chrn;
	      uint32_t rank = POS_UNDEF;
	      uint32_t tid = POS_UNDEF;
	      uint32_t pos = POS_UNDEF;
	      for(;tokIter != tokens.end(); ++tokIter) {
		// Optional fields
		boost::char_separator<char> kvsep(":");
		Tokenizer tokens(*tokIter, kvsep);
		Tokenizer::iterator tikv = tokens.begin();
		if (*tikv == "SN") {
		  ++tikv; ++tikv;
		  chrn = *tikv;
		} else if (*tikv == "SO") {
		  ++tikv; ++tikv;
		  pos = boost::lexical_cast<uint32_t>(*tikv);
		} else if (*tikv == "SR") {
		  ++tikv; ++tikv;
		  rank = boost::lexical_cast<uint32_t>(*tikv);
		}
	      }
	      // Remap ranks
	      if (rmap.find(rank) == rmap.end()) {
		// Create new rank level
		uint32_t newrank = rmap.size();
		rmap.insert(std::make_pair(rank, newrank));
	      }
	      rank = rmap[rank];
	      // Set chromosome names to NA and pos to 0 if not present
	      if (chrn.empty()) chrn = "NA";
	      if (pos == POS_UNDEF) pos = 0;

	      // New chromosome?
	      if (chrmap.find(chrn) != chrmap.end()) {
		if (rank != chrmap[chrn].first) {
		  std::cerr << "Identical chromosome names in different ranks!" << std::endl;
		  return false;
		}
		tid = chrmap[chrn].second;
	      } else {
		uint32_t chrid = chrmap.size();
		tid = chrid;
		chrmap.insert(std::make_pair(chrn, std::make_pair(rank, chrid)));
	      }
	      
	      // New segment
	      g.segments.push_back(Segment(tid, pos, sequence.size()));
	      if (storeseq) {
		// Store sequence
		sfile << ">" << segname << " " << chrn << ":" << (pos + 1) << "-" << (pos + sequence.size()) << ":" << rank << std::endl;
		sfile << sequence << std::endl;
	      }
	      seqsize += sequence.size();
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

    // Store chr names and ranks in the graph
    g.chrnames.resize(chrmap.size());
    g.ranks.resize(chrmap.size());
    for(typename TChrMap::const_iterator it = chrmap.begin(); it != chrmap.end(); ++it) {
      g.chrnames[it->second.second] = it->first;
      g.ranks[it->second.second] = it->second.first;
    }
    
    // Graph statistics
    std::cerr << "Parsed: " << g.segments.size() << " segments, " << g.links.size() << " links" << std::endl;
    std::cerr << "Total sequence size: " << seqsize << std::endl;

    if (storeseq) {
      // Close FASTA file
      sfile.close();

      // Build index
      if (fai_build(c.seqfile.string().c_str())) {
	std::cerr << "Could not build FASTA index!" << std::endl;
	return false;
      }
    }
    
    return true;
  }

  template<typename TConfig>
  inline void
  writeGfa(TConfig const& c, Graph const& g) {
    // Vertex map
    std::vector<std::string> idSegment(g.smap.size());
    for(typename Graph::TSegmentIdMap::const_iterator it = g.smap.begin(); it != g.smap.end(); ++it) idSegment[it->second] = it->first;
    
    // Temporary output file
    std::string filename = "test.out.gfa";
    
    // Output rGFA
    std::ofstream sfile;
    sfile.open(filename.c_str());

    // Output segments
    faidx_t* fai = fai_load(c.seqfile.string().c_str());
    for(uint32_t i = 0; i < g.segments.size(); ++i) {
      std::string seqid;
      if (c.seqCoords) seqid = g.chrnames[g.segments[i].tid];
      else seqid = idSegment[i];
      int32_t seqlen;
      char* seq = faidx_fetch_seq(fai, seqid.c_str(), 0, faidx_seq_len(fai, seqid.c_str()), &seqlen);
      sfile << "S\t" << seqid;
      if (c.seqCoords) sfile << "\t" << std::string(seq + g.segments[i].pos, seq + g.segments[i].pos + g.segments[i].len);
      else sfile << "\t" << seq;
      sfile << "\tLN:i:" << g.segments[i].len;
      sfile << "\tSN:Z:" << g.chrnames[g.segments[i].tid];
      sfile << "\tSO:i:" << g.segments[i].pos;
      sfile << "\tSR:i:" << g.ranks[g.segments[i].tid];
      sfile << std::endl;
      free(seq);
    }
    fai_destroy(fai);

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
