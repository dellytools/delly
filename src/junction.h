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

#ifndef JUNCTION_H
#define JUNCTION_H

#include "tags.h"
#include "api/BamReader.h"
#include "api/BamIndex.h"
#include "api/BamAlignment.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

namespace torali {

unsigned int const MAXKMERLENGTH=64;

inline std::string
  _reverseComplement(std::string const& ref) {
  std::string rev=ref;
  std::string::const_reverse_iterator itR = ref.rbegin();
  std::string::const_reverse_iterator itREnd = ref.rend();
  for(unsigned int i = 0; itR!=itREnd; ++itR, ++i) {
    switch (*itR) {
    case 'A': rev[i]='T'; break;
    case 'C': rev[i]='G'; break;
    case 'G': rev[i]='C'; break;
    case 'T': rev[i]='A'; break;
    case 'N': rev[i]='N'; break;
    default: break;
    }
  }
  return rev;
}


  template<typename TReference, typename TKmer>
inline unsigned int
    _getMinHammingDistance(TReference const& ref, TKmer const& kmer, unsigned int minReqDist) {
    unsigned int minHammingDist=kmer.size();
    if (ref.size()>=kmer.size()) {
      for(unsigned int seqIndex=0;((seqIndex<=(ref.size()-kmer.size())) && (minHammingDist>minReqDist));++seqIndex) {
	unsigned int hammingDist=0;
	for(unsigned int k=0;k<kmer.size();++k) 
	  if (ref[seqIndex+k]!=kmer[k]) ++hammingDist;
	if (hammingDist<minHammingDist) minHammingDist=hammingDist;
      }
    }
    return minHammingDist;
  }


  template<typename TReference, typename TKmerSet>
inline void
    __getKmers(TReference const& ref, TKmerSet& kmerSet, unsigned int kmerLength, unsigned int alphsize) {
    typedef typename TKmerSet::key_type TUSize;
    
    char currentk[MAXKMERLENGTH];
    for(unsigned int i = 0; i<MAXKMERLENGTH; ++i) currentk[i] = 0;
    unsigned int kmerlen=0;
    TUSize bucket = 0;
    std::string::const_iterator sIt= ref.begin();
    std::string::const_iterator sItEnd= ref.end();
    for(unsigned int seqIndex=0;sIt!=sItEnd;++sIt,++seqIndex) {
      if (*sIt!='N') {
	if (kmerlen==kmerLength) {
	  kmerSet.insert(bucket);
	  bucket -= ((TUSize) currentk[seqIndex % kmerLength] * (TUSize) std::pow((long double) alphsize, (long double) kmerLength - 1));
	  bucket *= (TUSize) alphsize;
	  bucket += (TUSize) dna_encode[int(*sIt)];
	} else bucket += ((TUSize) dna_encode[int(*sIt)] * (TUSize) std::pow((long double) alphsize, (long double) kmerLength - (++kmerlen)));
	currentk[seqIndex % kmerLength] = dna_encode[int(*sIt)];
      } else {
	if (kmerlen == kmerLength) kmerSet.insert(bucket);
	kmerlen = 0;
	bucket = 0;
      }
    }
    if (kmerlen == kmerLength) kmerSet.insert(bucket);
  }

  template<typename TReference, typename TKmerSet>
inline void
    _getKmers(TReference const& ref, TKmerSet& kmerSet, unsigned int kmerLength, unsigned int alphsize) {
    __getKmers(ref,kmerSet, kmerLength, alphsize);
    __getKmers(_reverseComplement(ref),kmerSet, kmerLength, alphsize);
  }

  template<typename TReference, typename TKmerSet, typename TUniqueKmer>
inline void
    _getUniqueKmers(TReference const& ref, TKmerSet const& kmerSet, TUniqueKmer& uniqueKmer, unsigned int kmerLength, unsigned int alphsize) {
    typedef typename TKmerSet::key_type TUSize;

    char currentk[MAXKMERLENGTH];
    for(unsigned int i = 0; i<MAXKMERLENGTH; ++i) currentk[i] = 0;
    unsigned int kmerlen=0;
    TUSize bucket = 0;
    std::string::const_iterator sIt= ref.begin();
    std::string::const_iterator sItEnd= ref.end();
    for(unsigned int seqIndex=0;sIt!=sItEnd;++sIt,++seqIndex) {
      if (*sIt!='N') {
	if (kmerlen==kmerLength) {
	  if (kmerSet.find(bucket)==kmerSet.end()) uniqueKmer.insert(ref.substr(seqIndex-kmerlen, kmerlen));
	  bucket -= ((TUSize) currentk[seqIndex % kmerLength] * (TUSize) std::pow((long double) alphsize, (long double) kmerLength - 1));
	  bucket *= (TUSize) alphsize;
	  bucket += (TUSize) dna_encode[int(*sIt)];
	} else bucket += ((TUSize) dna_encode[int(*sIt)] * (TUSize) std::pow((long double) alphsize, (long double) kmerLength - (++kmerlen)));
	currentk[seqIndex % kmerLength] = dna_encode[int(*sIt)];
      } else {
	if ((kmerlen == kmerLength) && (kmerSet.find(bucket)==kmerSet.end())) uniqueKmer.insert(ref.substr(seqIndex-kmerlen, kmerlen));
	kmerlen = 0;
	bucket = 0;
      }
    }
    if ((kmerlen == kmerLength) && (kmerSet.find(bucket)==kmerSet.end())) uniqueKmer.insert(ref.substr(ref.size()-kmerlen, kmerlen));
  }


  template<typename TFiles, typename TGenome, typename TSampleLibrary, typename TSVs, typename TCountMap>
inline void
    annotateJunctionReads(TFiles const& files, TGenome const& genome, uint16_t const minMapQual, TSampleLibrary& sampleLib, TSVs& svs, TCountMap& countMap)
{
  typedef typename TCountMap::key_type TSampleSVPair;
  typedef typename TCountMap::mapped_type TCountPair;
  typedef typename TCountPair::first_type TQualVector;

  // Get the references
  BamTools::BamReader readerRef;
  if ( ! readerRef.Open(files[0].string()) ) return;
  BamTools::RefVector references = readerRef.GetReferenceData();

  // Sort Structural Variants
  sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());

  // Initialize count map
  for(typename TSampleLibrary::iterator sIt = sampleLib.begin(); sIt!=sampleLib.end();++sIt) 
    for(typename TSVs::const_iterator itSV = svs.begin(); itSV!=svs.end(); ++itSV) countMap.insert(std::make_pair(std::make_pair(sIt->first, itSV->id), TCountPair()));

  // Process chromosome by chromosome
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Junction read annotation" << std::endl;
  boost::progress_display show_progress( references.size() );

  // Get junction probes for all SVs on this chromosome
  typedef boost::unordered_map<unsigned int, std::string> TProbes;
  TProbes altProbes;
  TProbes refProbes;

  // Parse genome
  kseq_t *seq;
  int l;
  gzFile fp = gzopen(genome.string().c_str(), "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    // Find reference index
    BamTools::RefVector::const_iterator itRef = references.begin();
    for(int32_t refIndex=0;itRef!=references.end();++itRef, ++refIndex) {
      if (seq->name.s == references[refIndex].RefName) {
        ++show_progress;
	
	// Iterate SVs
	typename TSVs::const_iterator itSV = svs.begin();
	typename TSVs::const_iterator itSVEnd = svs.end();
	for(;itSV!=itSVEnd;++itSV) {
	  if (itSV->consensus.empty()) continue;
	  if (itSV->chr2 == refIndex) {
	    int consLen = itSV->consensus.size();
            // For translocations temporarily store the first reference part in the probe strings
	    refProbes[itSV->id]=boost::to_upper_copy(std::string(seq->seq.s + std::max(itSV->svEnd - consLen, 0), seq->seq.s + itSV->svEnd + consLen));
	  }
          if (itSV->chr == refIndex) {
	    int consLen = itSV->consensus.size();
	    std::string sourceRight=refProbes[itSV->id];
	    refProbes[itSV->id]="";
	    std::string sourceLeft=boost::to_upper_copy(std::string(seq->seq.s + std::max(itSV->svStart - consLen, 0), seq->seq.s + itSV->svStart + consLen));
	    std::string cons=boost::to_upper_copy(itSV->consensus);
	    for (unsigned int kmerLength=11; kmerLength<MAXKMERLENGTH; kmerLength+=2) {
	      unsigned int minReqHammingDist = (unsigned int) ((std::log(kmerLength) * std::log(kmerLength)) / std::log(4));
	      typedef std::set<boost::multiprecision::uint128_t> TKmerSet;
	      TKmerSet refKmerSet;
	      std::string refLeft=sourceLeft;
	      std::string refRight=sourceRight;
	      _getKmers(refLeft, refKmerSet, kmerLength, 4);
	      _getKmers(refRight, refKmerSet, kmerLength, 4);
	      typedef boost::unordered_set<std::string> TUniqueKmers;
	      TUniqueKmers uniqueKmers;
	      _getUniqueKmers(cons, refKmerSet, uniqueKmers, kmerLength, 4);
	      unsigned int maxHammingDistance=0;
	      std::string maxHammingKmer="";
	      if (uniqueKmers.size()) {
		TUniqueKmers::iterator itK = uniqueKmers.begin();
		for(;itK!=uniqueKmers.end();++itK) {
		  unsigned int hammingDist = std::min(_getMinHammingDistance(refLeft, *itK, minReqHammingDist), _getMinHammingDistance(refRight, *itK, minReqHammingDist));
		  if ((hammingDist>minReqHammingDist) && (hammingDist>maxHammingDistance)) {
		    maxHammingDistance=hammingDist;
		    maxHammingKmer=*itK;
		  }
		}
	      }
	      if (maxHammingDistance) {
		// Found suitable kmer
		TKmerSet consKmerSet;
		_getKmers(cons, consKmerSet, kmerLength, 4);
		TUniqueKmers uniqueRefKmers;
		refLeft=refLeft.substr(std::max((int) refLeft.size()/2 - (int) kmerLength, 0), 2 * kmerLength);
		refRight=refRight.substr(std::max((int) refRight.size()/2 - (int) kmerLength, 0), 2 * kmerLength);
		_getUniqueKmers(refLeft, consKmerSet, uniqueRefKmers, kmerLength, 4);
		_getUniqueKmers(refRight, consKmerSet, uniqueRefKmers, kmerLength, 4);
		unsigned int maxRefHammingDistance=0;
		std::string maxRefHammingKmer="";
		if (uniqueRefKmers.size()) {
		  TUniqueKmers::iterator itK = uniqueRefKmers.begin();
		  for(;itK!=uniqueRefKmers.end();++itK) {
		    unsigned int hammingDist = _getMinHammingDistance(cons, *itK, minReqHammingDist);
		    if ((hammingDist>minReqHammingDist) && (hammingDist>maxRefHammingDistance)) {
		      maxRefHammingDistance=hammingDist;
		      maxRefHammingKmer=*itK;
		    }
		  }
		}
		if (maxRefHammingDistance) {
		  //std::cerr << references[itSV->chr].RefName << ',' << itSV->svStart << ',' << itSV->svEnd << std::endl;
		  //std::cerr << "RefLeft: " << sourceLeft << std::endl;
		  //std::cerr << "RefRight: " << sourceRight << std::endl;
		  //std::cerr << "Contig: " << cons << std::endl;
		  //std::cerr << "KmerLength: " << kmerLength << std::endl;
		  //std::cerr << maxRefHammingKmer << ',' << maxRefHammingDistance << ',' << maxHammingKmer << ',' << maxHammingDistance << std::endl;
		  altProbes[itSV->id]=maxHammingKmer;
		  refProbes[itSV->id]=maxRefHammingKmer;
		  break;
		}
	      }
	    }
	  }
	}
	
	// Iterate all samples
	for(unsigned int file_c = 0; file_c < files.size(); ++file_c) {
	  // Get a sample name
	  std::string sampleName(files[file_c].stem().string());
	  typename TSampleLibrary::iterator sampleIt=sampleLib.find(sampleName);

	  // Initialize bam file
	  BamTools::BamReader reader;
	  reader.Open(files[file_c].string());
	  reader.LocateIndex();

	  // Read alignments
	  itSV = svs.begin();
	  itSVEnd = svs.end();
	  for(;itSV!=itSVEnd;++itSV) {
	    if ((itSV->chr == refIndex) && (!itSV->consensus.empty()) && (altProbes.find(itSV->id)!=altProbes.end()) && (refProbes.find(itSV->id)!=refProbes.end())) {
	      // Consensus length
	      unsigned int consLen = itSV->consensus.size();
	      std::string altKmer = altProbes[itSV->id];
	      std::string rAltKmer = _reverseComplement(altKmer);
	      std::string refKmer = refProbes[itSV->id];
	      std::string rRefKmer = _reverseComplement(refKmer);
	      TQualVector altKmerCount;
	      TQualVector refKmerCount;
	      

	      // Unique reads
	      typedef std::set<int> TUniqueReads;
	      TUniqueReads unique_reads;

	      // Scan left and right breakpoint
	      BamTools::BamAlignment al;
	      for (unsigned int bpPoint = 0; bpPoint<2; ++bpPoint) {
		int32_t regionChr = itSV->chr;
		int regionStart = std::max(0, (int) itSV->svStart - (int) consLen/2);
		int regionEnd = itSV->svStart + consLen/2;
		if (bpPoint==1) {
		  regionChr = itSV->chr2;
		  regionStart = std::max(0, (int) itSV->svEnd - (int) consLen/2);
		  regionEnd = itSV->svEnd + consLen/2;
		}
		if (reader.SetRegion(regionChr, regionStart, regionChr, regionEnd)) {
		  while( reader.GetNextAlignmentCore(al) ) {
		    if (al.RefID != refIndex) break;
		    if ((al.AlignmentFlag & 0x0004) || (al.AlignmentFlag & 0x0100) || (al.AlignmentFlag & 0x0200) || (al.AlignmentFlag & 0x0400) || (al.AlignmentFlag & 0x0800) || (al.MapQuality < minMapQual)) continue;
		    
		    // Get the sequence 
		    al.BuildCharData();

		    // Is it a unique pair
		    boost::hash<std::string> string_hash;
		    typename TUniqueReads::const_iterator pos = unique_reads.begin();
		    bool inserted;
		    boost::tie(pos, inserted) = unique_reads.insert(string_hash(al.QueryBases));
		    if (inserted) {
		      unsigned int hamCutoff=2;
		      unsigned int altHam=_getMinHammingDistance(al.QueryBases, altKmer, 0);
		      if (altHam>hamCutoff) altHam=_getMinHammingDistance(al.QueryBases, rAltKmer, 0);
		      unsigned int refHam=_getMinHammingDistance(al.QueryBases, refKmer, 0);
		      if (refHam>hamCutoff) refHam=_getMinHammingDistance(al.QueryBases, rRefKmer, 0);
		      if ((altHam<=hamCutoff) && (refHam>hamCutoff)) {
			altKmerCount.push_back(al.MapQuality);
		      } else if ((altHam>hamCutoff) && (refHam<=hamCutoff)) {
			refKmerCount.push_back(al.MapQuality);
		      }
		    }
		  }
		}
	      }

	      // Insert counts
	      TSampleSVPair svSample = std::make_pair(sampleName, itSV->id);
	      typename TCountMap::iterator countMapIt=countMap.find(svSample);
	      countMapIt->second.first=refKmerCount;
	      countMapIt->second.second=altKmerCount;
	    }
	  }
	}
      }
    }
  }
  kseq_destroy(seq);
  gzclose(fp);
}




}

#endif


