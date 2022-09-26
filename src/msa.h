#ifndef MSA_H
#define MSA_H

#include <boost/multi_array.hpp>
#include "needle.h"
#include "gotoh.h"

namespace torali {

  inline int32_t
  lcs(std::string const& s1, std::string const& s2) {
    uint32_t m = s1.size();
    uint32_t n = s2.size();
    int32_t prevdiag = 0;
    int32_t prevprevdiag = 0;
    std::vector<int32_t> onecol(n+1, 0);
    for(uint32_t i = 0; i <= m; ++i) {
      for(uint32_t j = 0; j <= n; ++j) {
	if ((i==0) || (j==0)) {
	  onecol[j] = 0;
	  prevprevdiag = 0;
	  prevdiag = 0;
	} else {
	  prevprevdiag = prevdiag;
	  prevdiag = onecol[j];
	  if (s1[i-1] == s2[j-1]) onecol[j] = prevprevdiag + 1;
	  else onecol[j] = (onecol[j] > onecol[j-1]) ? onecol[j] : onecol[j-1];
	}
      }
    }
    return onecol[n];
  }

  template<typename TSplitReadSet, typename TDistArray>
  inline void
  distanceMatrix(TSplitReadSet const& sps, TDistArray& d) {
    typedef typename TDistArray::index TDIndex;
    typename TSplitReadSet::const_iterator sIt1 = sps.begin();
    for (TDIndex i = 0; sIt1 != sps.end(); ++sIt1, ++i) {
      typename TSplitReadSet::const_iterator sIt2 = sIt1;
      ++sIt2;
      for (TDIndex j = i+1; sIt2 != sps.end(); ++sIt2, ++j) {
	d[i][j] = (lcs(*sIt1, *sIt2) * 100) / std::min(sIt1->size(), sIt2->size());
      }
    }
  }

  template<typename TDistArray, typename TDIndex>
  inline int
  closestPair(TDistArray const& d, TDIndex num, TDIndex& dI, TDIndex& dJ) {
    int dMax = -1;
    for (TDIndex i = 0; i<num; ++i) {
      for (TDIndex j = i+1; j<num; ++j) {
	if (d[i][j]>dMax) {
	  dMax = d[i][j];
	  dI = i;
	  dJ = j;
	}
      }
    }
    return dMax;
  }

  template<typename TDistArray, typename TPhylogeny, typename TDIndex>
  inline void
  updateDistanceMatrix(TDistArray& d, TPhylogeny const& p, TDIndex num, TDIndex& dI, TDIndex& dJ) {
    for (TDIndex i = 0; i < num; ++i) 
      if (p[i][0] == -1) 
	d[i][num] = (((dI < i) ? d[dI][i] : d[i][dI]) + ((dJ < i) ? d[dJ][i] : d[i][dJ])) / 2;
    for (TDIndex i = 0; i<dI; ++i) d[i][dI] = -1;
    for (TDIndex i = dI+1; i<num+1; ++i) d[dI][i] = -1;
    for (TDIndex i = 0; i<dJ; ++i) d[i][dJ] = -1;
    for (TDIndex i = dJ+1; i<num+1; ++i) d[dJ][i] = -1;
  }

  template<typename TDistArray, typename TPhylogeny, typename TDIndex>
  inline TDIndex
  upgma(TDistArray& d, TPhylogeny& p, TDIndex num) {
    TDIndex nn = num;
    for(;nn<2*num+1; ++nn) {
      TDIndex dI = 0;
      TDIndex dJ = 0;
      if (closestPair(d, nn, dI, dJ) == -1) break;
      p[dI][0] = nn;
      p[dJ][0] = nn;
      p[nn][1] = dI;
      p[nn][2] = dJ;
      updateDistanceMatrix(d, p, nn, dI, dJ);
    }
    return (nn > 0) ? (nn - 1) : 0;
  }

  template<typename TConfig, typename TSplitReadSet, typename TPhylogeny, typename TDIndex, typename TAlign>
  inline void
  palign(TConfig const& c, TSplitReadSet const& sps, TPhylogeny const& p, TDIndex root, TAlign& align) {
    typedef typename TAlign::index TAIndex;
    if ((p[root][1] == -1) && (p[root][2] == -1)) {
      typename TSplitReadSet::const_iterator sIt = sps.begin();
      if (root) std::advance(sIt, root);
      align.resize(boost::extents[1][sIt->size()]);
      TAIndex ind = 0;
      for(typename std::string::const_iterator str = sIt->begin(); str != sIt->end(); ++str) align[0][ind++] = *str;
    } else {
      TAlign align1;
      palign(c, sps, p, p[root][1], align1);
      TAlign align2;
      palign(c, sps, p, p[root][2], align2);
      AlignConfig<true, true> endFreeAlign;
      gotoh(align1, align2, align, endFreeAlign, c.aliscore);
    }
  }

  template<typename TAlign>
  inline void
  sprealign(TAlign& align) {
    typedef typename TAlign::index TAIndex;
    for(TAIndex i = 0; i<align.shape()[0]; ++i) {
      // Precompute sub-alignment shapes
      int seqLength = 0;
      std::set<TAIndex> gaps;
      for(TAIndex j = 0; j<align.shape()[1]; ++j) {	
	if (align[i][j] != '-') {
	  ++seqLength;
	  bool gapOnly = true;
	  for(TAIndex k = 0; k<align.shape()[0]; ++k) {
	    if ((k!=i) && (align[k][j] != '-')) {
	      gapOnly = false;
	      break;
	    }
	  }
	  if (gapOnly) gaps.insert(j);
	}
      }

      // Create sub-alignments
      TAlign align1;
      TAIndex aind1 = 0;
      align1.resize(boost::extents[1][seqLength]);
      TAlign align2;
      TAIndex aind2 = 0;
      align2.resize(boost::extents[align.shape()[0] - 1][align.shape()[1] - gaps.size()]);
      for(TAIndex j = 0; j<align.shape()[1]; ++j) {
	if (align[i][j] != '-') align1[0][aind1++] = align[i][j];
	if (gaps.find(j) == gaps.end()) {
	  TAIndex kr = 0;
	  for(TAIndex k = 0; k<align.shape()[0]; ++k)
	    if (k!=i) align2[kr++][aind2] = align[k][j];
	  ++aind2;
	}
      }
      
      // Re-align sequence to profile
      AlignConfig<true, true> endFreeAlign;
      gotoh(align1, align2, align, endFreeAlign);
    }
  }


  template<typename TConfig, typename TAlign>
  inline void
  consensus(TConfig const& c, TAlign const& align, std::string& gapped, std::string& cs) {
    typedef typename TAlign::index TAIndex;

    // Calculate coverage
    typedef boost::multi_array<bool, 2> TFlag;
    TFlag fl(boost::extents[align.shape()[0]][align.shape()[1]]);
    typedef std::vector<int> TCoverage;
    TCoverage cov;
    cov.resize(align.shape()[1], 0);
    for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
      int start = 0;
      int end = -1;
      for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
	fl[i][j] = false;
	if (align[i][j] != '-') end = j;
	else if (end == -1) start = j + 1;
      }
      for(TAIndex j = start; j<=end; ++j) {
	++cov[j];
	fl[i][j] = true;
      }
    }
    
    int covThreshold = c.minCliqueSize;
    TAIndex j = 0;
    std::vector<char> cons(align.shape()[1], '-');
    for(typename TCoverage::const_iterator itCov = cov.begin(); itCov != cov.end(); ++itCov, ++j) {
      int32_t maxIdx = 4;  // Leading/trailing gaps until min. coverage is reached
      if (*itCov >= covThreshold) {
	// Get consensus letter
	std::vector<int32_t> count(5, 0); // ACGT-
	for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
	  if (fl[i][j]) {
	    if ((align[i][j] == 'A') || (align[i][j] == 'a')) ++count[0];
	    else if ((align[i][j] == 'C') || (align[i][j] == 'c')) ++count[1];
	    else if ((align[i][j] == 'G') || (align[i][j] == 'g')) ++count[2];
	    else if ((align[i][j] == 'T') || (align[i][j] == 't')) ++count[3];
	    else ++count[4];
	  }
	}
	maxIdx = 0;
	int32_t maxCount = count[0];
	for(uint32_t i = 1; i<5; ++i) {
	  if (count[i] > maxCount) {
	    maxCount = count[i];
	    maxIdx = i;
	  }
	}
      }
      switch (maxIdx) {
      case 0: cons[j] = 'A'; break;
      case 1: cons[j] = 'C'; break;
      case 2: cons[j] = 'G'; break;
      case 3: cons[j] = 'T'; break;
      default: break;
      }
    }
    gapped = std::string(cons.begin(), cons.end());
    for(uint32_t i = 0; i<cons.size(); ++i) {
      if (cons[i] != '-') cs.push_back(cons[i]);
    }
  }

  template<typename TConfig, typename TAlign>
  inline void
  consensus(TConfig const& c, TAlign const& align, std::string& cs) {
    std::string gapped;
    consensus(c, align, gapped, cs);
    //std::cerr << "Consensus:" << std::endl;
    //std::cerr << gapped << std::endl;
    //std::cerr << cs << std::endl;
  }

  template<typename TConfig, typename TSplitReadSet>
  inline int
  msa(TConfig const& c, TSplitReadSet const& sps, std::string& cs) {
    // Compute distance matrix
    typedef boost::multi_array<int, 2> TDistArray;
    typedef typename TDistArray::index TDIndex;
    TDIndex num = sps.size();
    TDistArray d(boost::extents[2*num+1][2*num+1]);
    for (TDIndex i = 0; i<(2*num+1); ++i) 
      for (TDIndex j = i+1; j<(2*num+1); ++j) 
	d[i][j]=-1;
    distanceMatrix(sps, d);

    // UPGMA
    typedef boost::multi_array<int, 2> TPhylogeny;
    TPhylogeny p(boost::extents[2*num+1][3]);
    for(TDIndex i = 0; i<(2*num+1); ++i) 
      for (TDIndex j = 0; j<3; ++j) p[i][j] = -1;
    TDIndex root = upgma(d, p, num);

    // Debug guide tree
    //std::cerr << "Phylogeny" << std::endl;
    //std::cerr << "#Sequences: " << sps.size() << std::endl;
    //std::cerr << "Root: " << root << std::endl;
    //std::cerr << "Node:Parent\tLeftChild\tRightChild" << std::endl;
    //for(TDIndex i = 0; i<(2*num+1); ++i) {
    //std::cerr << i << ':' << '\t';
    //for (TDIndex j = 0; j<3; ++j) {
    //std::cerr << p[i][j] << '\t';
    //}
    //std::cerr << std::endl;
    //}
    
    // Progressive Alignment
    typedef boost::multi_array<char, 2> TAlign;
    TAlign align;
    palign(c, sps, p, root, align);

    // Debug MSA
    //for(uint32_t i = 0; i<align.shape()[0]; ++i) {
    //for(uint32_t j = 0; j<align.shape()[1]; ++j) {
    //std::cerr << align[i][j];
    //}
    //std::cerr << std::endl;
    //}

    // Sequence to profile re-alignment
    //sprealign(align);

    // Consensus calling
    consensus(c, align, cs);
    
    // Return split-read support
    return align.shape()[0];
  }

  template<typename TStructuralVariant>
  inline void
  outputConsensus(bam_hdr_t* hdr, TStructuralVariant const& sv, std::string const& cons) {
    std::cerr << ">" << hdr->target_name[sv.chr] << ':' << sv.svStart << ',' << hdr->target_name[sv.chr2] << ':' << sv.svEnd << " SVT:" << sv.svt << " SR:" << sv.srSupport << " PE:" << sv.peSupport << std::endl;
    std::cerr << cons << std::endl;
  }

}

#endif
