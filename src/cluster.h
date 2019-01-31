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

#ifndef CLUSTER_H
#define CLUSTER_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/math/distributions/binomial.hpp>

#include <htslib/sam.h>

#include "util.h"
#include "junction.h"

namespace torali
{

  // Edge struct
  template<typename TWeight, typename TVertex>
  struct SREdgeRecord {
    typedef TVertex TVertexType;
    TVertex source;
    TVertex target;
    TWeight weight;

    SREdgeRecord(TVertex s, TVertex t, TWeight w) : source(s), target(t), weight(w) {}
  };

  // Sort edge records
  template<typename TRecord>
  struct SortSREdgeRecords : public std::binary_function<TRecord, TRecord, bool>
  {
    inline bool operator()(TRecord const& e1, TRecord const& e2) const {
      return ((e1.weight < e2.weight) || ((e1.weight == e2.weight) && (e1.source < e2.source)) || ((e1.weight == e2.weight) && (e1.source == e2.source) && (e1.target < e2.target)));
    }
  };


  template<typename TConfig, typename TCompEdgeList>
  inline void
  _searchCliques(TConfig const& c, TCompEdgeList& compEdge, std::vector<SRBamRecord>& br, std::vector<StructuralVariantRecord>& sv) {
    typedef typename TCompEdgeList::mapped_type TEdgeList;
    typedef typename TEdgeList::value_type TEdgeRecord;
    typedef typename TEdgeRecord::TVertexType TVertex;

    // Iterate all components
    for(typename TCompEdgeList::iterator compIt = compEdge.begin(); compIt != compEdge.end(); ++compIt) {
      // Sort edges by weight
      std::sort(compIt->second.begin(), compIt->second.end(), SortSREdgeRecords<TEdgeRecord>());

      // Find a large clique
      typename TEdgeList::const_iterator itWEdge = compIt->second.begin();
      typename TEdgeList::const_iterator itWEdgeEnd = compIt->second.end();
      typedef std::set<TVertex> TCliqueMembers;
      TCliqueMembers clique;
      TCliqueMembers incompatible;
      
      // Initialize clique
      clique.insert(itWEdge->source);
      int32_t chr = br[itWEdge->source].chr;
      int32_t chr2 = br[itWEdge->source].chr2;
      int32_t ciposlow = br[itWEdge->source].pos;
      uint64_t pos = br[itWEdge->source].pos;
      int32_t ciposhigh = br[itWEdge->source].pos; 
      int32_t ciendlow = br[itWEdge->source].pos2;
      uint64_t pos2 = br[itWEdge->source].pos2;
      int32_t ciendhigh = br[itWEdge->source].pos2;
      int32_t mapq = br[itWEdge->source].qual;
      int32_t inslen = br[itWEdge->source].inslen;

      // Grow clique
      bool cliqueGrow = true;
      while (cliqueGrow) {
	itWEdge = compIt->second.begin();
	cliqueGrow = false;
	// Find next best edge for extension
	for(;(!cliqueGrow) && (itWEdge != itWEdgeEnd);++itWEdge) {
	  TVertex v;
	  if ((clique.find(itWEdge->source) == clique.end()) && (clique.find(itWEdge->target) != clique.end())) v = itWEdge->source;
	  else if ((clique.find(itWEdge->source) != clique.end()) && (clique.find(itWEdge->target) == clique.end())) v = itWEdge->target;
	  else continue;
	  if (incompatible.find(v) != incompatible.end()) continue;
	  // Try to update clique with this vertex
	  int32_t newCiPosLow = std::min(br[v].pos, ciposlow);
	  int32_t newCiPosHigh = std::max(br[v].pos, ciposhigh);
	  int32_t newCiEndLow = std::min(br[v].pos2, ciendlow);
	  int32_t newCiEndHigh = std::max(br[v].pos2, ciendhigh);
	  if (((newCiPosHigh - newCiPosLow) < (int32_t) c.maxReadSep) && ((newCiEndHigh - newCiEndLow) < (int32_t) c.maxReadSep)) cliqueGrow = true;
	  if (cliqueGrow) {
	    // Accept new vertex
	    clique.insert(v);
	    ciposlow = newCiPosLow;
	    pos += br[v].pos;
	    ciposhigh = newCiPosHigh;
	    ciendlow = newCiEndLow;
	    pos2 += br[v].pos2;
	    ciendhigh = newCiEndHigh;
	    mapq += br[v].qual;
	    inslen += br[v].inslen;
	  } else incompatible.insert(v);
	}
      }

      // At least 2 split reads?
      if (clique.size()>1) {
	int32_t svStart = (int32_t) (pos / (uint64_t) clique.size());
	int32_t svEnd = (int32_t) (pos2 / (uint64_t) clique.size());
	int32_t svInsLen = (int32_t) (inslen / (int32_t) clique.size());
	if ((ciposlow > svStart) || (ciposhigh < svStart) || (ciendlow > svEnd) || (ciendhigh < svEnd)) {
	  std::cerr << "Warning: Confidence intervals out of bounds: " << ciposlow << ',' << svStart << ',' << ciposhigh << ':' << ciendlow << ',' << svEnd << ',' << ciendhigh << std::endl;
	}
	sv.push_back(StructuralVariantRecord(chr, svStart, chr2, svEnd, (ciposlow - svStart), (ciposhigh - svStart), (ciendlow - svEnd), (ciendhigh - svEnd), clique.size(), mapq, svInsLen));
      }
    }
  }
  

  template<typename TConfig>
  inline void
  cluster(TConfig const& c, std::vector<SRBamRecord>& br, std::vector<StructuralVariantRecord>& sv) {
    uint32_t count = 0;
    for(int32_t refIdx = 0; refIdx < c.nchr; ++refIdx) {
      
      // Components
      typedef std::vector<uint32_t> TComponent;
      TComponent comp;
      comp.resize(br.size(), 0);
      uint32_t numComp = 0;

      // Edge lists for each component
      typedef uint32_t TWeightType;
      typedef uint32_t TVertex;
      typedef SREdgeRecord<TWeightType, TVertex> TEdgeRecord;
      typedef std::vector<TEdgeRecord> TEdgeList;
      typedef std::map<uint32_t, TEdgeList> TCompEdgeList;
      TCompEdgeList compEdge;

	
      std::size_t lastConnectedNode = 0;
      std::size_t lastConnectedNodeStart = 0;
      for(uint32_t i = 0; i<br.size(); ++i) {
	if (br[i].chr == refIdx) {
	  ++count;
	  // Safe to clean the graph?
	  if (i > lastConnectedNode) {
	    // Clean edge lists
	    if (!compEdge.empty()) {
	      // Search cliques
	      _searchCliques(c, compEdge, br, sv);
	      lastConnectedNodeStart = lastConnectedNode;
	      compEdge.clear();
	    }
	  }
	  
	  
	  for(uint32_t j = i + 1; j<br.size(); ++j) {
	    if (br[j].chr == refIdx) {
	      if ( (uint32_t) (br[j].pos - br[i].pos) > c.maxReadSep) break;
	      if ( (uint32_t) std::abs(br[j].pos2 - br[i].pos2) < c.maxReadSep) {
		// Update last connected node
		if (j > lastConnectedNode) lastConnectedNode = j;
		
		// Assign components
		uint32_t compIndex = 0;
		if (!comp[i]) {
		  if (!comp[j]) {
		    // Both vertices have no component
		    compIndex = ++numComp;
		    comp[i] = compIndex;
		    comp[j] = compIndex;
		    compEdge.insert(std::make_pair(compIndex, TEdgeList()));
		  } else {
		    compIndex = comp[j];
		    comp[i] = compIndex;
		  }	
		} else {
		  if (!comp[j]) {
		    compIndex = comp[i];
		    comp[j] = compIndex;
		  } else {
		    // Both vertices have a component
		    if (comp[j] == comp[i]) {
		      compIndex = comp[j];
		    } else {
		      // Merge components
		      compIndex = comp[i];
		      uint32_t otherIndex = comp[j];
		      if (otherIndex < compIndex) {
			compIndex = comp[j];
			otherIndex = comp[i];
		      }
		      // Re-label other index
		      for(uint32_t k = lastConnectedNodeStart; k <= lastConnectedNode; ++k) {
			if (otherIndex == comp[k]) comp[k] = compIndex;
		      }
		      // Merge edge lists
		      TCompEdgeList::iterator compEdgeIt = compEdge.find(compIndex);
		      TCompEdgeList::iterator compEdgeOtherIt = compEdge.find(otherIndex);
		      compEdgeIt->second.insert(compEdgeIt->second.end(), compEdgeOtherIt->second.begin(), compEdgeOtherIt->second.end());
		      compEdge.erase(compEdgeOtherIt);
		    }
		  }
		}
		
		// Append new edge
		TCompEdgeList::iterator compEdgeIt = compEdge.find(compIndex);
		if (compEdgeIt->second.size() < c.graphPruning) {
		  // Breakpoint distance
		  TWeightType weight = std::abs(br[j].pos2 - br[i].pos2) + std::abs(br[j].pos - br[i].pos);
		  compEdgeIt->second.push_back(TEdgeRecord(i, j, weight));
		}
	      }
	    }
	  }
	}
      }
      // Search cliques
      if (!compEdge.empty()) {
	_searchCliques(c, compEdge, br, sv);
	compEdge.clear();
      }
    }

    // Sort SVs
    std::sort(sv.begin(), sv.end(), SortSVs<StructuralVariantRecord>());
  }

}

#endif
