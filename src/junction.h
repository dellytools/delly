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


#include <boost/multiprecision/cpp_int.hpp>
#include <htslib/kseq.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <stdio.h>
#include "fasta_reader.h"
#include "tags.h"

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
    unsigned int seqIndex=0;
    for(;sIt!=sItEnd;++sIt,++seqIndex) {
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

template<typename TSequence, typename TBreaks>
inline void
  _getReferenceBreakpoint(TSequence const& ref, TSequence const& var, TSequence& outVar, TBreaks& bpcoord) 
  {
    typedef std::vector<char> TAlignSeq;
    typedef FastaRecord<std::string, unsigned int, Dna5GapAlphabet, TAlignSeq, void> TFastaRecord;
    std::vector<TFastaRecord> alignFwd;
    std::vector<TFastaRecord> alignRev;
    std::vector<TFastaRecord> align;

    // Convert to alignment data structure
    TAlignSeq r;
    for(typename TSequence::const_iterator itS = ref.begin(); itS!=ref.end(); ++itS) r.push_back(dna5_encode[(int) *itS]);
    TAlignSeq v;
    for(typename TSequence::const_iterator itS = var.begin(); itS!=var.end(); ++itS) v.push_back(dna5_encode[(int) *itS]);

    // Decide forward or reverse alignment
    DnaScore<int> scr(1, -2, -5, 0);
    AlignConfig<true, false, false, true> alc;
    int fwdScore=globalGotohAlignment(alignFwd, r, v, scr, alc);
    reverseComplement(v);
    int revScore=globalGotohAlignment(alignRev, r, v, scr, alc);
    if (fwdScore > revScore) reverseComplement(v);

    // Compute final alignment
    globalGotohAlignment(align, r, v, scr, alc);
      
    // Find longest internal gap
    int refIndex=-1;
    int varIndex=-1;
    int gapStartRefIndex=0;
    int gapEndRefIndex=0;
    int gapStartVarIndex=0;
    int gapEndVarIndex=0;
    int longestGap=0;
    TAlignSeq::const_iterator alRef = align[0].seq.begin();
    TAlignSeq::const_iterator alVar = align[1].seq.begin();
    bool inGap=false;
    for (; (alRef != align[0].seq.end()) && (alVar != align[1].seq.end()); ++alRef, ++alVar) {
      if (*alRef != dna5gap_encode[(int) '-']) ++refIndex;
      if (*alVar != dna5gap_encode[(int) '-']) ++varIndex;
      // Internal gap?
      if (((*alVar == dna5gap_encode[(int) '-']) || (*alRef == dna5gap_encode[(int) '-'])) && (refIndex>0) && (varIndex>0)) {
	if (!inGap) {
	  gapStartRefIndex=refIndex;
	  gapEndRefIndex=refIndex;
	  gapStartVarIndex=varIndex;
	  gapEndVarIndex=varIndex;
	} else {
	  gapEndRefIndex=refIndex;
	  gapEndVarIndex=varIndex;
	}
	inGap=true;
      } else {
	if ((inGap) && ((gapEndRefIndex - gapStartRefIndex) > longestGap)) {
	  longestGap=gapEndRefIndex - gapStartRefIndex;
	  bpcoord[0]=gapStartRefIndex;
	  bpcoord[1]=gapEndRefIndex;
	  bpcoord[2]=gapStartVarIndex;
	  bpcoord[3]=gapEndVarIndex;
	}
	inGap=false;
      }
    }

    // Output sequence might be reverse complemented
    outVar.clear();
    for(TAlignSeq::const_iterator itA = v.begin(); itA!=v.end(); ++itA) outVar.push_back(dna5_decode[(int) *itA]);

    // Debug code
    /*
    std::cout << "FwdScore: " << fwdScore << ", RevScore: " << revScore << std::endl;
    std::cout << "Alignment: " << fwd << std::endl;
    std::cout << "Ref: ";
    int ind=0;
    for(TAlignSeq::iterator alItTmp = r.begin(); alItTmp != r.end(); ++alItTmp, ++ind) {
      std::cout << dna5_decode[(int) *alItTmp];
      if ((ind==bpcoord[0]) || (ind==bpcoord[1])) std::cout << " | ";
    }
    std::cout << std::endl;
    std::cout << "Var: ";
    ind=0;
    for(TAlignSeq::iterator alItTmp = v.begin(); alItTmp != v.end(); ++alItTmp, ++ind) {
      std::cout << dna5_decode[(int) *alItTmp];
      if ((ind==bpcoord[2]) || (ind==bpcoord[3])) std::cout << " | ";
    }
    std::cout << std::endl;
    std::cout << ">Ref" << std::endl;
    for(TAlignSeq::iterator alItTmp = align[0].seq.begin(); alItTmp != align[0].seq.end(); ++alItTmp) std::cout << dna5gap_decode[(int) *alItTmp];
    std::cout << std::endl;
    std::cout << ">Contig" << std::endl;
    for(TAlignSeq::iterator alItTmp = align[1].seq.begin(); alItTmp != align[1].seq.end(); ++alItTmp) std::cout << dna5gap_decode[(int) *alItTmp];
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    */
  }



  template<typename TFiles, typename TGenome, typename TSampleLibrary, typename TSVs, typename TCountMap, typename TTag>
inline void
    annotateJunctionReads(TFiles const& files, TGenome const& genome, uint16_t const minMapQual, TSampleLibrary& sampleLib, TSVs& svs, TCountMap& countMap, SVType<TTag> svType)
{
  typedef typename TCountMap::key_type TSampleSVPair;
  typedef typename TCountMap::mapped_type TCountPair;
  typedef typename TCountPair::first_type TQualVector;

  // Open file handles
  typedef std::vector<std::string> TRefNames;
  typedef std::vector<uint32_t> TRefLength;
  TRefNames refnames;
  TRefLength reflen;
  typedef std::vector<samFile*> TSamFile;
  typedef std::vector<hts_idx_t*> TIndex;
  TSamFile samfile;
  TIndex idx;
  samfile.resize(files.size());
  idx.resize(files.size());
  for(unsigned int file_c = 0; file_c < files.size(); ++file_c) {
    samfile[file_c] = sam_open(files[file_c].string().c_str(), "r");
    idx[file_c] = sam_index_load(samfile[file_c], files[file_c].string().c_str());
    if (!file_c) {
      bam_hdr_t* hdr = sam_hdr_read(samfile[file_c]);
      for (int i = 0; i<hdr->n_targets; ++i) {
	refnames.push_back(hdr->target_name[i]);
	reflen.push_back(hdr->target_len[i]);
      }
      bam_hdr_destroy(hdr);
    }
  }

  // Sort Structural Variants
  sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());

  // Initialize count map
  for(typename TSampleLibrary::iterator sIt = sampleLib.begin(); sIt!=sampleLib.end();++sIt) 
    for(typename TSVs::const_iterator itSV = svs.begin(); itSV!=svs.end(); ++itSV) countMap.insert(std::make_pair(std::make_pair(sIt->first, itSV->id), TCountPair()));

  // Process chromosome by chromosome
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Junction read annotation" << std::endl;
  boost::progress_display show_progress( refnames.size() );

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
    for(int32_t refIndex=0; refIndex < refnames.size(); ++refIndex) {
      if (seq->name.s == refnames[refIndex]) {
        ++show_progress;
	
	// Iterate SVs
	typename TSVs::const_iterator itSV = svs.begin();
	typename TSVs::const_iterator itSVEnd = svs.end();
	for(;itSV!=itSVEnd;++itSV) {
	  if (itSV->consensus.empty()) continue;
	  int consLen = itSV->consensus.size();

	  // Create a pseudo structural variant record
	  StructuralVariantRecord svRec;
	  svRec.chr = itSV->chr;
	  svRec.chr2 = itSV->chr2;
	  svRec.svStartBeg = std::max(itSV->svStart - consLen, 0);
	  svRec.svStart = itSV->svStart;
	  svRec.svStartEnd = std::min((uint32_t) itSV->svStart + consLen, reflen[itSV->chr]);
	  svRec.svEndBeg = std::max(itSV->svEnd - consLen, 0);
	  svRec.svEnd = itSV->svEnd;
	  svRec.svEndEnd = std::min((uint32_t) itSV->svEnd + consLen, reflen[itSV->chr2]);
	  svRec.ct = itSV->ct;
	  if ((itSV->chr != itSV->chr2) && (itSV->chr2 == refIndex)) {
            refProbes[itSV->id] = _getSVRef(seq->seq.s, svRec, refIndex, svType);
          }
          if (itSV->chr == refIndex) {
	    // Get the reference string
	    if (itSV->chr != itSV->chr2) svRec.consensus=refProbes[itSV->id];
	    refProbes[itSV->id]="";
	    std::string svRefStr = _getSVRef(seq->seq.s, svRec, refIndex, svType);
	    std::string consensus;
	    std::vector<int> bpcoord(4);
	    std::fill(bpcoord.begin(), bpcoord.end(), 0);
	    _getReferenceBreakpoint(svRefStr, itSV->consensus, consensus, bpcoord);
	    // Has a breakpoint been found
	    if ((bpcoord[2]==0) || (bpcoord[3]==0)) continue;
	    // Find tagging kmers
	    for (unsigned int kmerLength=11; kmerLength<MAXKMERLENGTH; kmerLength+=2) {
	      unsigned int minReqHammingDist = (unsigned int) ((std::log(kmerLength) * std::log(kmerLength)) / std::log(4));
	      typedef std::set<boost::multiprecision::uint128_t> TKmerSet;
	      TKmerSet refKmerSet;
	      _getKmers(svRefStr, refKmerSet, kmerLength, 4);
	      typedef boost::unordered_set<std::string> TUniqueKmers;
	      TUniqueKmers uniqueKmers;
	      std::string consPiece=consensus.substr(std::max(bpcoord[2] - (int) kmerLength, 0), 2 * kmerLength);
	      _getUniqueKmers(consPiece, refKmerSet, uniqueKmers, kmerLength, 4);
	      unsigned int maxHammingDistance=0;
	      std::string maxHammingKmer="";
	      if (uniqueKmers.size()) {
		TUniqueKmers::iterator itK = uniqueKmers.begin();
		for(;itK!=uniqueKmers.end();++itK) {
		  unsigned int hammingDist = _getMinHammingDistance(svRefStr, *itK, minReqHammingDist);
		  if ((hammingDist>minReqHammingDist) && (hammingDist>maxHammingDistance)) {
		    maxHammingDistance=hammingDist;
		    maxHammingKmer=*itK;
		  }
		}
	      }
	      if (maxHammingDistance) {
		// Found suitable kmer, find reference tagging kmer
		TKmerSet consKmerSet;
		_getKmers(consensus, consKmerSet, kmerLength, 4);
		TUniqueKmers uniqueRefKmers;
		std::string refLeft=svRefStr.substr(std::max(bpcoord[0] - (int) kmerLength, 0), 2 * kmerLength);
		std::string refRight=svRefStr.substr(std::max(bpcoord[1] - (int) kmerLength, 0), 2 * kmerLength);
		_getUniqueKmers(refLeft, consKmerSet, uniqueRefKmers, kmerLength, 4);
		_getUniqueKmers(refRight, consKmerSet, uniqueRefKmers, kmerLength, 4);
		unsigned int maxRefHammingDistance=0;
		std::string maxRefHammingKmer="";
		if (uniqueRefKmers.size()) {
		  TUniqueKmers::iterator itK = uniqueRefKmers.begin();
		  for(;itK!=uniqueRefKmers.end();++itK) {
		    unsigned int hammingDist = _getMinHammingDistance(consensus, *itK, minReqHammingDist);
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
	      for (unsigned int bpPoint = 0; bpPoint<2; ++bpPoint) {
		int32_t regionChr = itSV->chr;
		int regionStart = std::max(0, (int) itSV->svStart - (int) consLen/2);
		int regionEnd = itSV->svStart + consLen/2;
		if (bpPoint==1) {
		  regionChr = itSV->chr2;
		  regionStart = std::max(0, (int) itSV->svEnd - (int) consLen/2);
		  regionEnd = itSV->svEnd + consLen/2;
		}
		hts_itr_t* iter = sam_itr_queryi(idx[file_c], regionChr, regionStart, regionEnd);
		bam1_t* rec = bam_init1();
		while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
		  if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
		  if (rec->core.qual < minMapQual) continue;

		  // Get the sequence
		  std::string sequence; 
		  sequence.resize(rec->core.l_qseq);
		  uint8_t* seqptr = bam_get_seq(rec);
		  for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

		  // Is it a unique pair
		  boost::hash<std::string> string_hash;
		  typename TUniqueReads::const_iterator pos = unique_reads.begin();
		  bool inserted;
		  boost::tie(pos, inserted) = unique_reads.insert(string_hash(sequence));
		  if (inserted) {
		    unsigned int hamCutoff=2;
		    unsigned int altHam=_getMinHammingDistance(sequence, altKmer, 0);
		    if (altHam>hamCutoff) altHam=_getMinHammingDistance(sequence, rAltKmer, 0);
		    unsigned int refHam=_getMinHammingDistance(sequence, refKmer, 0);
		    if (refHam>hamCutoff) refHam=_getMinHammingDistance(sequence, rRefKmer, 0);
		    if ((altHam<=hamCutoff) && (refHam>hamCutoff)) {
		      altKmerCount.push_back(rec->core.qual);
		    } else if ((altHam>hamCutoff) && (refHam<=hamCutoff)) {
		      refKmerCount.push_back(rec->core.qual);
		    }
		  }
		}
		bam_destroy1(rec);
		hts_itr_destroy(iter);
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

  // Clean-up
  for(unsigned int file_c = 0; file_c < files.size(); ++file_c) {
    hts_idx_destroy(idx[file_c]);
    sam_close(samfile[file_c]);
  }
}




}

#endif


