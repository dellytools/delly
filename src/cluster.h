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


  // Reduced bam alignment record data structure
  struct BamAlignRecord {
    int32_t tid;         
    int32_t pos;
    int32_t mtid; 
    int32_t mpos;
    int32_t alen;
    int32_t malen;
    int32_t Median;
    int32_t Mad;
    int32_t maxNormalISize;
    uint32_t flag;
    uint8_t MapQuality;
  
    BamAlignRecord(bam1_t* rec, uint8_t pairQuality, uint16_t a, uint16_t ma, int32_t median, int32_t mad, int32_t maxISize) : tid(rec->core.tid), pos(rec->core.pos), mtid(rec->core.mtid), mpos(rec->core.mpos), alen(a), malen(ma), Median(median), Mad(mad), maxNormalISize(maxISize), flag(rec->core.flag), MapQuality(pairQuality) {}
  };

  // Sort reduced bam alignment records
  template<typename TRecord>
  struct SortBamRecords : public std::binary_function<TRecord, TRecord, bool>
  {
    inline bool operator()(TRecord const& s1, TRecord const& s2) const {
      if (s1.tid==s1.mtid) {
	return ((std::min(s1.pos, s1.mpos) < std::min(s2.pos, s2.mpos)) || 
		((std::min(s1.pos, s1.mpos) == std::min(s2.pos, s2.mpos)) && (std::max(s1.pos, s1.mpos) < std::max(s2.pos, s2.mpos))) ||
		((std::min(s1.pos, s1.mpos) == std::min(s2.pos, s2.mpos)) && (std::max(s1.pos, s1.mpos) == std::max(s2.pos, s2.mpos)) && (s1.maxNormalISize < s2.maxNormalISize)));
      } else {
	return ((s1.pos < s2.pos) ||
		((s1.pos == s2.pos) && (s1.mpos < s2.mpos)) ||
		((s1.pos == s2.pos) && (s1.mpos == s2.mpos) && (s1.maxNormalISize < s2.maxNormalISize)));
      }
    }
  };
  

  // Edge struct
  template<typename TWeight, typename TVertex>
  struct EdgeRecord {
    typedef TVertex TVertexType;
    TVertex source;
    TVertex target;
    TWeight weight;
    
    EdgeRecord(TVertex s, TVertex t, TWeight w) : source(s), target(t), weight(w) {}
  };

  // Sort edge records
  template<typename TRecord>
  struct SortEdgeRecords : public std::binary_function<TRecord, TRecord, bool>
  {
    inline bool operator()(TRecord const& e1, TRecord const& e2) const {
      return ((e1.weight < e2.weight) || ((e1.weight == e2.weight) && (e1.source < e2.source)) || ((e1.weight == e2.weight) && (e1.source == e2.source) && (e1.target < e2.target)));
    }
  };

  // Initialize clique, deletions
  template<typename TBamRecord, typename TSize>
  inline void
  _initClique(TBamRecord const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, int32_t const svt) {
    if (_translocation(svt)) {
      uint8_t ct = _getSpanOrientation(svt);
      if (ct%2==0) {
	svStart = el.pos + el.alen;
	if (ct>=2) svEnd = el.mpos;
	else svEnd = el.mpos + el.malen;
      } else {
	svStart = el.pos;
	if (ct>=2) svEnd = el.mpos + el.malen;
	else svEnd = el.mpos;
      }
      wiggle=el.maxNormalISize;
    } else {
      if (svt == 0) {
	svStart = el.mpos + el.malen;
	svEnd = el.pos + el.alen;
	wiggle = el.maxNormalISize - std::max(el.alen, el.malen);
      } else if (svt == 1) {
	svStart = el.mpos;
	svEnd = el.pos;
	wiggle = el.maxNormalISize - std::max(el.alen, el.malen);
      } else if (svt == 2) {
	svStart = el.mpos + el.malen;
	svEnd = el.pos;
	wiggle =  -el.maxNormalISize;
      } else if (svt == 3) {
	svStart = el.mpos;
	svEnd = el.pos + el.alen;
	wiggle = el.maxNormalISize;
      }
    } 
  }

  // Update clique, deletions
  template<typename TBamRecord, typename TSize>
  inline bool 
  _updateClique(TBamRecord const& el, TSize& svStart, TSize& svEnd, TSize& wiggle, int32_t const svt) 
  {
    if (_translocation(svt)) {
      int ct = _getSpanOrientation(svt);
      TSize newSvStart;
      TSize newSvEnd;
      TSize newWiggle = wiggle;
      if (ct%2==0) {
	newSvStart = std::max(svStart, el.pos + el.alen);
	newWiggle -= (newSvStart - svStart);
	if (ct>=2) {
	  newSvEnd = std::min(svEnd, el.mpos);
	  newWiggle -= (svEnd - newSvEnd);
	} else  {
	  newSvEnd = std::max(svEnd, el.mpos + el.malen);
	  newWiggle -= (newSvEnd - svEnd);
	}
      } else {
	newSvStart = std::min(svStart, el.pos);
	newWiggle -= (svStart - newSvStart);
	if (ct>=2) {
	  newSvEnd = std::max(svEnd, el.mpos + el.malen);
	  newWiggle -= (newSvEnd - svEnd);
	} else {
	  newSvEnd = std::min(svEnd, el.mpos);
	  newWiggle -= (svEnd - newSvEnd);
	}
      }
      // Is this still a valid translocation cluster?
      if (newWiggle>0) {
	svStart = newSvStart;
	svEnd = newSvEnd;
	wiggle = newWiggle;
	return true;
      }
      return false;
    } else {
      if ((svt == 0) || (svt == 1)) { 
	int ct = _getSpanOrientation(svt);
	TSize newSvStart;
	TSize newSvEnd;
	TSize newWiggle;
	TSize wiggleChange;
	if (!ct) {
	  newSvStart = std::max(svStart, el.mpos + el.malen);
	  newSvEnd = std::max(svEnd, el.pos + el.alen);
	  newWiggle = std::min(el.maxNormalISize - (newSvStart - el.mpos), el.maxNormalISize - (newSvEnd - el.pos));
	  wiggleChange = wiggle - std::max(newSvStart - svStart, newSvEnd - svEnd);
	} else {
	  newSvStart = std::min(svStart, el.mpos);
	  newSvEnd = std::min(svEnd, el.pos);
	  newWiggle = std::min(el.maxNormalISize - (el.mpos + el.malen - newSvStart), el.maxNormalISize - (el.pos + el.alen - newSvEnd));
	  wiggleChange = wiggle - std::max(svStart - newSvStart, svEnd - newSvEnd);
	}
	if (wiggleChange < newWiggle) newWiggle=wiggleChange;
	
	// Does the new inversion size agree with all pairs
	if ((newSvStart < newSvEnd) && (newWiggle>=0)) {
	  svStart = newSvStart;
	  svEnd = newSvEnd;
	  wiggle = newWiggle;
	  return true;
	}
	return false;
      } else if (svt == 2) {
	TSize newSvStart = std::max(svStart, el.mpos + el.malen);
	TSize newSvEnd = std::min(svEnd, el.pos);
	TSize newWiggle = el.pos + el.alen - el.mpos - el.maxNormalISize - (newSvEnd - newSvStart);
	TSize wiggleChange = wiggle + (svEnd-svStart) - (newSvEnd - newSvStart);
	if (wiggleChange > newWiggle) newWiggle=wiggleChange;
	
	// Does the new deletion size agree with all pairs
	if ((newSvStart < newSvEnd) && (newWiggle<=0)) {
	  svStart = newSvStart;
	  svEnd = newSvEnd;
	  wiggle = newWiggle;
	  return true;
	}
	return false;
      } else if (svt == 3) {
	TSize newSvStart = std::min(svStart, el.mpos);
	TSize newSvEnd = std::max(svEnd, el.pos + el.alen);
	TSize newWiggle = el.pos - (el.mpos + el.malen) + el.maxNormalISize - (newSvEnd - newSvStart);
	TSize wiggleChange = wiggle - ((newSvEnd - newSvStart) - (svEnd-svStart));
	if (wiggleChange < newWiggle) newWiggle = wiggleChange;
	
	// Does the new duplication size agree with all pairs
	if ((newSvStart < newSvEnd) && (newWiggle>=0)) {
	  svStart = newSvStart;
	  svEnd = newSvEnd;
	  wiggle = newWiggle;
	  return true;
	}
	return false;
      }
    }
    return false;
  }


  template<typename TConfig, typename TCompEdgeList>
  inline void
  _searchCliques(TConfig const& c, TCompEdgeList& compEdge, std::vector<SRBamRecord>& br, std::vector<StructuralVariantRecord>& sv, int32_t const svt) {
    typedef typename TCompEdgeList::mapped_type TEdgeList;
    typedef typename TEdgeList::value_type TEdgeRecord;
    typedef typename TEdgeRecord::TVertexType TVertex;

    // Iterate all components
    for(typename TCompEdgeList::iterator compIt = compEdge.begin(); compIt != compEdge.end(); ++compIt) {
      // Sort edges by weight
      std::sort(compIt->second.begin(), compIt->second.end(), SortEdgeRecords<TEdgeRecord>());

      // Find a large clique
      typename TEdgeList::const_iterator itWEdge = compIt->second.begin();
      typename TEdgeList::const_iterator itWEdgeEnd = compIt->second.end();
      typedef std::set<TVertex> TCliqueMembers;
      typedef std::set<std::size_t> TSeeds;
      TCliqueMembers clique;
      TCliqueMembers incompatible;
      TSeeds seeds;
      
      // Initialize clique
      clique.insert(itWEdge->source);
      seeds.insert(br[itWEdge->source].id);
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

      // Initialize wiggle
      uint32_t wiggle = c.maxReadSep;
      if (_translocation(svt)) wiggle = 2 * c.maxReadSep;
      else {
	// At most 1000bp breakpoint offset
	uint32_t svvar = std::abs(0.1 * (pos2 - pos));
	if (svt == 4) svvar = std::abs(0.1 * inslen);
	if (wiggle < svvar) wiggle = svvar;
	if (wiggle > 1000) wiggle = 1000;
      }

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
	  if (seeds.find(br[v].id) != seeds.end()) continue;
	  // Try to update clique with this vertex
	  int32_t newCiPosLow = std::min(br[v].pos, ciposlow);
	  int32_t newCiPosHigh = std::max(br[v].pos, ciposhigh);
	  int32_t newCiEndLow = std::min(br[v].pos2, ciendlow);
	  int32_t newCiEndHigh = std::max(br[v].pos2, ciendhigh);
	  if (((newCiPosHigh - newCiPosLow) < (int32_t) wiggle) && ((newCiEndHigh - newCiEndLow) < (int32_t) wiggle)) cliqueGrow = true;
	  if (cliqueGrow) {
	    // Accept new vertex
	    clique.insert(v);
	    seeds.insert(br[v].id);
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

      // Count paired-end fragments only once
      std::size_t prevSeed = 0;
      uint32_t cliqSize = 0;
      TSeeds selectedSeeds;
      for(typename TSeeds::const_iterator itS = seeds.begin(); itS != seeds.end(); ++itS) {
	if (prevSeed + 1 != *itS) {
	  selectedSeeds.insert(*itS);
	  ++cliqSize;
	}
	prevSeed = *itS;
      }
	    
      // Enough split reads?
      if (cliqSize >= c.minCliqueSize) {
	int32_t svStart = (int32_t) (pos / (uint64_t) clique.size());
	int32_t svEnd = (int32_t) (pos2 / (uint64_t) clique.size());
	int32_t svInsLen = (int32_t) (inslen / (int32_t) clique.size());
	if (_svSizeCheck(svStart, svEnd, svt, svInsLen)) {
	  if ((ciposlow > svStart) || (ciposhigh < svStart) || (ciendlow > svEnd) || (ciendhigh < svEnd)) {
	    std::cerr << "Warning: Confidence intervals out of bounds: " << ciposlow << ',' << svStart << ',' << ciposhigh << ':' << ciendlow << ',' << svEnd << ',' << ciendhigh << std::endl;
	  }
	  int32_t svid = sv.size();
	  sv.push_back(StructuralVariantRecord(chr, svStart, chr2, svEnd, (ciposlow - svStart), (ciposhigh - svStart), (ciendlow - svEnd), (ciendhigh - svEnd), cliqSize, mapq / clique.size(), mapq, svInsLen, svt, svid));
	  // Reads assigned
	  for(typename TCliqueMembers::iterator itC = clique.begin(); itC != clique.end(); ++itC) {
	    if (selectedSeeds.find(br[*itC].id) != selectedSeeds.end()) {
	      //std::cerr << svid << ',' << br[*itC].id << std::endl;
	      br[*itC].svid = svid;
	    }
	  }
	}
      }
    }
  }
  

  template<typename TConfig>
  inline void
  cluster(TConfig const& c, std::vector<SRBamRecord>& br, std::vector<StructuralVariantRecord>& sv, int32_t const svt) {
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
      typedef EdgeRecord<TWeightType, TVertex> TEdgeRecord;
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
	      _searchCliques(c, compEdge, br, sv, svt);
	      lastConnectedNodeStart = lastConnectedNode;
	      compEdge.clear();
	    }
	  }

	  // Get size variability
	  uint32_t varisize = c.maxReadSep;
	  if (_translocation(svt)) varisize = 2 * c.maxReadSep;
	  else {
	    // At most 1000bp breakpoint offset
	    uint32_t svvar = std::abs(0.1 * (br[i].pos2 - br[i].pos));
	    if (svt == 4) svvar = std::abs(0.1 * br[i].inslen);
	    if (varisize < svvar) varisize = svvar;
	    if (varisize > 1000) varisize = 1000;
	  }
	  for(uint32_t j = i + 1; j<br.size(); ++j) {
	    if (br[j].chr == refIdx) {
	      if ( (uint32_t) (br[j].pos - br[i].pos) > varisize) break;
	      if ((svt == 4) && (std::abs(br[j].inslen - br[i].inslen) > varisize)) continue;
	      if ( (uint32_t) std::abs(br[j].pos2 - br[i].pos2) < varisize) {
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
	_searchCliques(c, compEdge, br, sv, svt);
	compEdge.clear();
      }
    }
  }


  template<typename TConfig, typename TCompEdgeList, typename TBamRecord, typename TSVs>
  inline void
  _searchCliques(TConfig const& c, TCompEdgeList& compEdge, TBamRecord const& bamRecord, TSVs& svs, int32_t const svt) {
    typedef typename TCompEdgeList::mapped_type TEdgeList;
    typedef typename TEdgeList::value_type TEdgeRecord;

    // Iterate all components
    for(typename TCompEdgeList::iterator compIt = compEdge.begin(); compIt != compEdge.end(); ++compIt) {
      // Sort edges by weight
      std::sort(compIt->second.begin(), compIt->second.end(), SortEdgeRecords<TEdgeRecord>());
      
      // Find a large clique
      typename TEdgeList::const_iterator itWEdge = compIt->second.begin();
      typename TEdgeList::const_iterator itWEdgeEnd = compIt->second.end();
      typedef std::set<std::size_t> TCliqueMembers;
      
      TCliqueMembers clique;
      TCliqueMembers incompatible;
      int32_t svStart = -1;
      int32_t svEnd = -1;
      int32_t wiggle = 0;
      int32_t clusterRefID=bamRecord[itWEdge->source].tid;
      int32_t clusterMateRefID=bamRecord[itWEdge->source].mtid;
      _initClique(bamRecord[itWEdge->source], svStart, svEnd, wiggle, svt);
      if ((clusterRefID==clusterMateRefID) && (svStart >= svEnd))  continue;
      clique.insert(itWEdge->source);
      
      // Grow the clique from the seeding edge
      bool cliqueGrow=true;
      while (cliqueGrow) {
	itWEdge = compIt->second.begin();
	cliqueGrow = false;
	for(;(!cliqueGrow) && (itWEdge != itWEdgeEnd);++itWEdge) {
	  std::size_t v;
	  if ((clique.find(itWEdge->source) == clique.end()) && (clique.find(itWEdge->target) != clique.end())) v = itWEdge->source;
	  else if ((clique.find(itWEdge->source) != clique.end()) && (clique.find(itWEdge->target) == clique.end())) v = itWEdge->target;
	  else continue;
	  if (incompatible.find(v) != incompatible.end()) continue;
	  cliqueGrow = _updateClique(bamRecord[v], svStart, svEnd, wiggle, svt);
	  if (cliqueGrow) clique.insert(v);
	  else incompatible.insert(v);
	}
      }

      // Enough paired-ends
      if ((clique.size() >= c.minCliqueSize) && (_svSizeCheck(svStart, svEnd, svt))) {
	StructuralVariantRecord svRec;
	svRec.chr = clusterRefID;
	svRec.chr2 = clusterMateRefID;
	svRec.svStart = (uint32_t) svStart + 1;
	svRec.svEnd = (uint32_t) svEnd + 1;
	svRec.peSupport = clique.size();
	int32_t ci_wiggle = std::max(abs(wiggle), 50);
	svRec.ciposlow = -ci_wiggle;
	svRec.ciposhigh = ci_wiggle;
	svRec.ciendlow = -ci_wiggle;
	svRec.ciendhigh = ci_wiggle;
	svRec.mapq = 0;
	std::vector<uint8_t> mapQV;
	for(typename TCliqueMembers::const_iterator itC = clique.begin(); itC!=clique.end(); ++itC) {
	  mapQV.push_back(bamRecord[*itC].MapQuality);
	  svRec.mapq += bamRecord[*itC].MapQuality;
	}
	std::sort(mapQV.begin(), mapQV.end());
	svRec.peMapQuality = mapQV[mapQV.size()/2];
	svRec.srSupport=0;
	svRec.srAlignQuality=0;
	svRec.precise=false;
	svRec.svt = svt;
	svRec.insLen = 0;
	svRec.homLen = 0;
	svs.push_back(svRec);
      }
    }
  }
  
  

  template<typename TConfig>
  inline void
  cluster(TConfig const& c, std::vector<BamAlignRecord>& bamRecord, std::vector<StructuralVariantRecord>& svs, uint32_t const varisize, int32_t const svt) {
    typedef typename std::vector<BamAlignRecord> TBamRecord;
    // Components
    typedef std::vector<uint32_t> TComponent;
    TComponent comp;
    comp.resize(bamRecord.size(), 0);
    uint32_t numComp = 0;
      
    // Edge lists for each component
    typedef uint32_t TWeightType;
    typedef uint32_t TVertex;
    typedef EdgeRecord<TWeightType, TVertex> TEdgeRecord;
    typedef std::vector<TEdgeRecord> TEdgeList;
    typedef std::map<uint32_t, TEdgeList> TCompEdgeList;
    TCompEdgeList compEdge;
    
    // Iterate the chromosome range
    std::size_t lastConnectedNode = 0;
    std::size_t lastConnectedNodeStart = 0;
    std::size_t bamItIndex = 0;
    for(TBamRecord::const_iterator bamIt = bamRecord.begin(); bamIt != bamRecord.end(); ++bamIt, ++bamItIndex) {
      // Safe to clean the graph?
      if (bamItIndex > lastConnectedNode) {
	// Clean edge lists
	if (!compEdge.empty()) {
	  _searchCliques(c, compEdge, bamRecord, svs, svt);
	  lastConnectedNodeStart = lastConnectedNode;
	  compEdge.clear();
	}
      }
      int32_t const minCoord = _minCoord(bamIt->pos, bamIt->mpos, svt);
      int32_t const maxCoord = _maxCoord(bamIt->pos, bamIt->mpos, svt);
      TBamRecord::const_iterator bamItNext = bamIt;
      ++bamItNext;
      std::size_t bamItIndexNext = bamItIndex + 1;
      for(; ((bamItNext != bamRecord.end()) && ((uint32_t) std::abs(_minCoord(bamItNext->pos, bamItNext->mpos, svt) + bamItNext->alen - minCoord) <= varisize)) ; ++bamItNext, ++bamItIndexNext) {
	  // Check that mate chr agree (only for translocations)
	if (bamIt->mtid != bamItNext->mtid) continue;
	
	// Check combinability of pairs
	if (_pairsDisagree(minCoord, maxCoord, bamIt->alen, bamIt->maxNormalISize, _minCoord(bamItNext->pos, bamItNext->mpos, svt), _maxCoord(bamItNext->pos, bamItNext->mpos, svt), bamItNext->alen, bamItNext->maxNormalISize, svt)) continue;
	
	// Update last connected node
	if (bamItIndexNext > lastConnectedNode ) lastConnectedNode = bamItIndexNext;
	
	// Assign components
	uint32_t compIndex = 0;
	if (!comp[bamItIndex]) {
	  if (!comp[bamItIndexNext]) {
	    // Both vertices have no component
	    compIndex = ++numComp;
	    comp[bamItIndex] = compIndex;
	    comp[bamItIndexNext] = compIndex;
	    compEdge.insert(std::make_pair(compIndex, TEdgeList()));
	  } else {
	    compIndex = comp[bamItIndexNext];
	    comp[bamItIndex] = compIndex;
	  }
	} else {
	  if (!comp[bamItIndexNext]) {
	    compIndex = comp[bamItIndex];
	    comp[bamItIndexNext] = compIndex;
	  } else {
	    // Both vertices have a component
	    if (comp[bamItIndexNext] == comp[bamItIndex]) {
	      compIndex = comp[bamItIndexNext];
	    } else {
	      // Merge components
	      compIndex = comp[bamItIndex];
	      uint32_t otherIndex = comp[bamItIndexNext];
	      if (otherIndex < compIndex) {
		compIndex = comp[bamItIndexNext];
		otherIndex = comp[bamItIndex];
	      }
	      // Re-label other index
	      for(std::size_t i = lastConnectedNodeStart; i <= lastConnectedNode; ++i) {
		if (otherIndex == comp[i]) comp[i] = compIndex;
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
	  TWeightType weight = (TWeightType) ( std::log((double) abs( abs( (_minCoord(bamItNext->pos, bamItNext->mpos, svt) - minCoord) - (_maxCoord(bamItNext->pos, bamItNext->mpos, svt) - maxCoord) ) - abs(bamIt->Median - bamItNext->Median)) + 1) / std::log(2) );
	  compEdgeIt->second.push_back(TEdgeRecord(bamItIndex, bamItIndexNext, weight));
	}
      }
    }
    if (!compEdge.empty()) {
      _searchCliques(c, compEdge, bamRecord, svs, svt);
      compEdge.clear();
    }
  }
  

    
}

#endif
