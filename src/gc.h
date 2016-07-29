/*
============================================================================
DELLY: Structural variant discovery by integrated PE mapping and SR analysis
============================================================================
Copyright (C) 2012 Tobias Rausch

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

#ifndef GC_H
#define GC_H

#include <boost/dynamic_bitset.hpp>

#include "util.h"


namespace torali
{

  struct GCBias {
    int32_t median;
    int32_t minNormalISize;
    int32_t maxNormalISize;
    uint8_t defaultOrient;
    double gc[101];

    GCBias() {
      for(uint32_t i = 0; i<101; ++i) gc[i] = 1.0;
    }
  };


  template<typename TConfig, typename TInvariantRegion, typename TSampleLibrary, typename TGCBias>
  inline void
  gcBiasPerRG(TConfig const& c, TInvariantRegion const& invariantRegions, TSampleLibrary const& sampleLib, TGCBias& gcBias) {
    typedef typename TSampleLibrary::value_type TLibraryMap;
    typedef typename TInvariantRegion::value_type TChrIntervals;
    typedef typename TGCBias::value_type TGCMap;

    // Defaults
    uint32_t minFCObservations = 100;
    uint32_t maxFCObservations = 100000;
    uint32_t fragPerBin = 500;
    
    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile[0]);

    // Initialize GC maps
    gcBias.resize(c.files.size(), TGCMap());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      for(typename TLibraryMap::const_iterator libIt = sampleLib[file_c].begin(); libIt != sampleLib[file_c].end(); ++libIt) {
	typename TGCMap::iterator itGC = gcBias[file_c].insert(std::make_pair(libIt->first, GCBias())).first;
	itGC->second.median = libIt->second.median;
	itGC->second.minNormalISize = libIt->second.minNormalISize;
	itGC->second.maxNormalISize = libIt->second.maxNormalISize;
	itGC->second.defaultOrient = libIt->second.defaultOrient;
      }
    }
    
    // Process chromosome by chromosome
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "GC bias calculation" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );

    // Parse genome
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      ++show_progress;
      if (invariantRegions[refIndex].empty()) continue;
      bool nodata = true;
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	std::string suffix("cram");
	std::string str(c.files[file_c].string());
	if ((str.size() >= suffix.size()) && (str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0)) {
	  nodata = false;
	  break;
	}
	uint64_t mapped = 0;
	uint64_t unmapped = 0;
	hts_idx_get_stat(idx[file_c], refIndex, &mapped, &unmapped);
	if (mapped) {
	  nodata = false;
	  break;
	}
      }
      if (nodata) continue;

      // Load chromosome
      char* seq = NULL;
      int32_t seqlen = -1;
      std::string tname(hdr->target_name[refIndex]);
      seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);

      // Preprocess GC-content
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet gc(hdr->target_len[refIndex]);
      TBitSet nrun(hdr->target_len[refIndex]);
      for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	if ((seq[i] == 'c') || (seq[i] == 'C') || (seq[i] == 'g') || (seq[i] == 'G')) gc[i] = 1;
	else if ((seq[i] == 'n') || (seq[i] == 'N')) nrun[i] = 1;
      }

      // Debug: Chromosome GC content
      //std::cerr << refIndex << ',' << hdr->target_len[refIndex] << ',' << nrun.count() << ',' << gc.count() << ',' << (double) gc.count() / (double) (hdr->target_len[refIndex] - nrun.count()) << std::endl;


      
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	typedef std::vector<int32_t> TBinFragmentCount;
	typedef std::vector<TBinFragmentCount> TGCFragmentCountVector;
	typedef boost::unordered_map<std::string, TGCFragmentCountVector> TRgBinFragmentCount;
	TRgBinFragmentCount rgGCP;
	for(typename TLibraryMap::const_iterator libIt = sampleLib[file_c].begin(); libIt != sampleLib[file_c].end(); ++libIt) {
	  TRgBinFragmentCount::iterator it = rgGCP.insert(std::make_pair(libIt->first, TGCFragmentCountVector())).first;
	  it->second.resize(101, TBinFragmentCount());
	}
	for(typename TChrIntervals::const_iterator vRIt = invariantRegions[refIndex].begin(); vRIt != invariantRegions[refIndex].end(); ++vRIt) {
	  uint32_t rsize = vRIt->upper() - vRIt->lower();
	  typedef std::vector<uint8_t> TFragmentCount;
	  typedef boost::unordered_map<std::string, TFragmentCount> TRgFragmentCount;
	  TRgFragmentCount rcounts;
	  for(typename TLibraryMap::const_iterator libIt = sampleLib[file_c].begin(); libIt != sampleLib[file_c].end(); ++libIt) {
	    TRgFragmentCount::iterator it = rcounts.insert(std::make_pair(libIt->first, TFragmentCount())).first;
	    it->second.resize(rsize, 0);
	  }
	  hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, vRIt->lower(), vRIt->upper());
	  bam1_t* rec = bam_init1();
	  while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	    if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	    if (rec->core.qual < c.minGenoQual) continue;

	    uint32_t midPoint = rec->core.pos + halfAlignmentLength(rec);
	    if ((midPoint >= vRIt->lower()) && (midPoint < vRIt->upper())) {
	      // Get the library information
	      std::string rG = "DefaultLib";
	      uint8_t *rgptr = bam_aux_get(rec, "RG");
	      if (rgptr) {
		char* rg = (char*) (rgptr + 1);
		rG = std::string(rg);
	      }
	      TRgFragmentCount::iterator it = rcounts.find(rG);
	      ++it->second[midPoint - vRIt->lower()];
	    }
	  }
	  
	  bam_destroy1(rec);
	  hts_itr_destroy(iter);

	  // Process counts of this region
	  for(typename TRgFragmentCount::const_iterator itRg = rcounts.begin(); itRg != rcounts.end(); ++itRg) {
	    int32_t totalcount = 0;
	    for(TFragmentCount::const_iterator itFC = itRg->second.begin(); itFC != itRg->second.end(); ++itFC) totalcount += *itFC;
	    if (totalcount) {
	      typename TRgBinFragmentCount::iterator itBinGCP = rgGCP.find(itRg->first);
	      typename TLibraryMap::const_iterator libIt = sampleLib[file_c].find(itRg->first);
	      uint32_t binsize = libIt->second.avgdist * fragPerBin;
	      if (rsize < binsize) continue;
	      int32_t rsum = 0;
	      int32_t gcsum = 0;
	      int32_t nsum = 0;
	      for(uint32_t i = 0; i<binsize; ++i) {
		rsum += itRg->second[i];
		gcsum += gc[vRIt->lower() + i];
		nsum += nrun[vRIt->lower() + i];
	      }
	      for(uint32_t i = binsize; i<rsize; ++i) {
		if (!nsum) {
		  int32_t bin = (gcsum * 100) / binsize;
		  if (itBinGCP->second[bin].size() < maxFCObservations) itBinGCP->second[bin].push_back(rsum);
		}
		rsum -= itRg->second[i - binsize];
		gcsum -= gc[vRIt->lower() + i - binsize];
		nsum -= nrun[vRIt->lower() + i - binsize];
		rsum += itRg->second[i];
		gcsum += gc[vRIt->lower() + i];
		nsum += nrun[vRIt->lower() + i];
	      }
	      if (!nsum) {
		int32_t bin = (gcsum * 100) / binsize;
		if (itBinGCP->second[bin].size() < maxFCObservations) itBinGCP->second[bin].push_back(rsum);
	      }
	    }
	  }
	}

	// Calculate GC bias
	for(TRgBinFragmentCount::iterator itBinGCP = rgGCP.begin(); itBinGCP != rgGCP.end(); ++itBinGCP) {
	  typename TGCMap::iterator itGCMap = gcBias[file_c].find(itBinGCP->first);
	  std::vector<int32_t> meds(101, -1);
	  std::vector<int32_t> validMedCollection;
	  for(int32_t gci = 0; gci < 101; ++gci) {
	    if (itBinGCP->second[gci].size() > minFCObservations) {
	      int32_t med = 0;
	      getMedian(itBinGCP->second[gci].begin(), itBinGCP->second[gci].end(), med);
	      validMedCollection.push_back(med);
	      meds[gci] = med;
	    }
	  }
	  int32_t globalMed = 0;
	  getMedian(validMedCollection.begin(), validMedCollection.end(), globalMed);
	  for(int32_t gci = 0; gci < 101; ++gci) {
	    if (meds[gci] != -1) itGCMap->second.gc[gci] = (double) globalMed / (double) meds[gci];
	  }
	}
      }
      if (seq != NULL) free(seq);
    }
    fai_destroy(fai);
      
    // Clean-up
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }

  template<typename TConfig, typename TGCBias, typename TSVs, typename TCountMap, typename TBpLevelType>
  inline void
  annotateCoverageGCAware(TConfig const& c, TGCBias const& gcBias, TSVs& svs, TCountMap& countMap, TBpLevelType bpLevel)
  {
    typedef typename TSVs::value_type TStructuralVariant;
    typedef typename TGCBias::value_type TGCMap;
    
    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile[0]);

    // Sort Structural Variants
    sort(svs.begin(), svs.end(), SortSVs<TStructuralVariant>());

    // Initialize count maps
    countMap.resize(c.files.size());
    for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) countMap[file_c].resize(svs.size());

    // Process chromosome by chromosome
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Read-depth annotation" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );
    
    // Parse genome
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      ++show_progress;

      // Fetch all SV boundaries
      typename TSVs::const_iterator itSV = std::lower_bound(svs.begin(), svs.end(), TStructuralVariant(refIndex, 0, hdr->target_len[refIndex]), SortSVs<TStructuralVariant>());
      typename TSVs::const_iterator itSVEnd = std::upper_bound(svs.begin(), svs.end(), TStructuralVariant(refIndex, hdr->target_len[refIndex], hdr->target_len[refIndex]), SortSVs<TStructuralVariant>());
      if (itSV != itSVEnd) {
	typedef std::vector< std::pair<int32_t, int32_t> > TInterval;
	TInterval intervals;
	typedef std::vector<int32_t> TBreaks;
	TBreaks iBounds;
	iBounds.push_back(itSV->svStart);
	iBounds.push_back(itSV->svEnd);
	int32_t maxPos = itSV->svEnd;
	typename TSVs::const_iterator itSVFor = itSV;
	for(++itSVFor;itSVFor!=itSVEnd; ++itSVFor) {
	  if (itSVFor->svStart>maxPos) {
	    std::sort(iBounds.begin(), iBounds.end());
	    typename TBreaks::const_iterator itBreak = iBounds.begin();
	    typename TBreaks::const_iterator itBreakNext = iBounds.begin();
	    for(++itBreakNext; itBreakNext!=iBounds.end(); ++itBreakNext, ++itBreak)
	      if (*itBreak < *itBreakNext) intervals.push_back(std::make_pair(*itBreak, *itBreakNext));
	    maxPos=itSVFor->svEnd;
	    iBounds.clear();
	  }
	  iBounds.push_back(itSVFor->svStart);
	  iBounds.push_back(itSVFor->svEnd);
	  if (itSVFor->svEnd > maxPos) maxPos=itSVFor->svEnd;
	}
	std::sort(iBounds.begin(), iBounds.end());
	typename TBreaks::const_iterator itBreak = iBounds.begin();
	typename TBreaks::const_iterator itBreakNext = iBounds.begin();
	for(++itBreakNext; itBreakNext!=iBounds.end(); ++itBreakNext, ++itBreak) 
	  if (*itBreak < *itBreakNext) intervals.push_back(std::make_pair(*itBreak, *itBreakNext));
	iBounds.clear();

	// Load chromosome
	char* seq = NULL;
	int32_t seqlen = -1;
	std::string tname(hdr->target_name[refIndex]);
	seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);

	// Preprocess GC-content
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet gc(hdr->target_len[refIndex]);
	TBitSet nrun(hdr->target_len[refIndex]);
	for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	  if ((seq[i] == 'c') || (seq[i] == 'C') || (seq[i] == 'g') || (seq[i] == 'G')) gc[i] = 1;
	  else if ((seq[i] == 'n') || (seq[i] == 'N')) nrun[i] = 1;
	}
    
	// Iterate all samples
#pragma omp parallel for default(shared)
	for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	  // Read alignments
	  typedef std::vector<double> TIntervalSum;
	  TIntervalSum readSum;
	  TIntervalSum bpSum;
      	  
	  // Process sub-intervals
	  for(typename TInterval::const_iterator itInt = intervals.begin(); itInt!=intervals.end(); ++itInt) {
	    // Count reads / aligned base-pairs
	    double bp_sum = 0;
	    double read_sum = 0;
	    
	    hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, itInt->first, itInt->second);
	    bam1_t* rec = bam_init1();
	    while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	      if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	      if (rec->core.qual < c.minGenoQual) continue;

	      // Get the library information
	      std::string rG = "DefaultLib";
	      uint8_t *rgptr = bam_aux_get(rec, "RG");
	      if (rgptr) {
		char* rg = (char*) (rgptr + 1);
		rG = std::string(rg);
	      }
	      typename TGCMap::const_iterator gcIt = gcBias[file_c].find(rG);
	      int32_t midPoint = rec->core.pos;
	      int32_t regStart = 0;
	      int32_t regEnd = 0;
	      if (gcIt->second.median == 0) {
		// Single-end library
		midPoint = rec->core.pos + halfAlignmentLength(rec);
		regStart = std::max(0, midPoint - 200);
		regEnd = std::min(midPoint + 200, (int32_t) hdr->target_len[refIndex]);
	      } else {
		// Paired-end library
		if (!(rec->core.flag & BAM_FPAIRED)) continue;
		if (rec->core.flag & BAM_FMUNMAP) continue;
		if (rec->core.pos < rec->core.mpos) continue;

		int32_t outerISize = rec->core.pos - rec->core.mpos + rec->core.l_qseq;
		if ((rec->core.tid!=rec->core.mtid) || (getStrandIndependentOrientation(rec->core) != gcIt->second.defaultOrient) || (outerISize < gcIt->second.minNormalISize) || (outerISize > gcIt->second.maxNormalISize)) continue;
		midPoint = rec->core.pos + outerISize / 2;
		regStart = rec->core.pos;
		regEnd = std::min(rec->core.mpos + rec->core.l_qseq, (int32_t) hdr->target_len[refIndex]);
	      }
	      int32_t nsum = 0;
	      int32_t gcsum = 0;
	      for(int32_t i = regStart; i<regEnd; ++i) {
		nsum += nrun[i];
		gcsum += gc[i];
	      }
	      int32_t size = (regEnd - regStart) - nsum;
	      int32_t gcbin = 50;
	      if (size > 0) gcbin = (gcsum * 100) / size;
	      //std::cerr << itInt->first << '-' << itInt->second << ':' << rec->core.pos << ',' << midPoint << ',' << regStart << ',' << regEnd << ',' << gcbin << std::endl;
	      if ((midPoint >= itInt->first) && (midPoint < itInt->second)) read_sum += gcIt->second.gc[gcbin];
	      _addBpCounts(rec, itInt->first, itInt->second, bp_sum, bpLevel);
	    }
	    bam_destroy1(rec);
	    hts_itr_destroy(iter);
	    bpSum.push_back(bp_sum);
	    readSum.push_back(read_sum);
	  }

	  typename TSVs::const_iterator itFileSV = itSV;
	  for(; itFileSV != itSVEnd; ++itFileSV) {
	    // Sum up sub-intervals
	    double cumBpSum = 0;
	    double cumReadSum = 0;
	    typename TInterval::const_iterator itInt = std::lower_bound(intervals.begin(), intervals.end(), std::make_pair(itFileSV->svStart, 0));
	    typename TIntervalSum::const_iterator itRead = readSum.begin();
	    itRead += (itInt - intervals.begin());
	    typename TIntervalSum::const_iterator itBp = bpSum.begin();
	    itBp += (itInt -intervals.begin());
	    for(;itInt!=intervals.end(); ++itInt, ++itRead, ++itBp) {
	      if (itInt->first >= itFileSV->svEnd) break;
	      cumBpSum+=(*itBp);
	      cumReadSum+=(*itRead);
	    }
	
	    // Store counts
	    countMap[file_c][itFileSV->id].first = (int32_t) cumBpSum;
	    countMap[file_c][itFileSV->id].second = (int32_t) cumReadSum;
	  }
	}
	if (seq != NULL) free(seq);
      }
    }
    fai_destroy(fai);
      
    // Clean-up
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }

}

#endif
