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

#ifndef SPANNING_H
#define SPANNING_H

#include <boost/container/flat_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <htslib/sam.h>
#include "tags.h"

namespace torali {

  struct SpanningCount {
    int32_t refh1;
    int32_t refh2;
    int32_t alth1;
    int32_t alth2;
    std::vector<uint8_t> ref;
    std::vector<uint8_t> alt;

    SpanningCount() : refh1(0), refh2(0), alth1(0), alth2(0) {}
  };

  template<typename TDefaultOrientation>
    inline bool
    _mateIsUpstream(TDefaultOrientation defOrient, bool firstRead, bool reverse) {
    if (firstRead) {
      if (reverse) {
	if (defOrient % 2 == 0) return false;
	else return true;
      } else {
	if (defOrient % 2 == 0) return true;
	else return false;
      }
    } else {
      if (reverse) {
	if ((defOrient==1) || (defOrient==2)) return false;
	else return true;
      } else {
	if ((defOrient==1) || (defOrient==2)) return true;
	else return false;
      }
    }
  }   

  template<typename TConfig, typename TSampleLibrary, typename TSVs, typename TCountMap>
  inline void
  annotateSpanningCoverage(TConfig& c, TSampleLibrary& sampleLib, TSVs& svs, TCountMap& spanCountMap)
  {
    typedef typename TCountMap::value_type::value_type TCountPair;
    typedef typename TSampleLibrary::value_type TLibraryMap;

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
    std::sort(svs.begin(), svs.end(), SortSVs<StructuralVariantRecord>());

    // Initialize count map
    spanCountMap.resize(c.files.size());
    uint32_t lastId = svs.size();
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) spanCountMap[file_c].resize(2 * lastId, TCountPair());

    // Dump file
    boost::iostreams::filtering_ostream dumpOut;
    if (c.pedumpflag) {
      dumpOut.push(boost::iostreams::gzip_compressor());
      dumpOut.push(boost::iostreams::file_sink(c.pedump.string().c_str(), std::ios_base::out | std::ios_base::binary));
      dumpOut << "#svid\tbam\tqname\tchr\tpos\tmatechr\tmatepos\tmapq" << std::endl;
    }

    // Iterate all samples
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Breakpoint spanning coverage annotation" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );

    // Qualities
    typedef boost::unordered_map<std::size_t, uint8_t> TQualities;
    std::vector<TQualities> qualitiestra;
    qualitiestra.resize(c.files.size());
    typedef boost::unordered_map<std::size_t, int32_t> TAlignmentLength;
    std::vector<TAlignmentLength> alentra;
    alentra.resize(c.files.size());

    // Iterate chromosomes
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      ++show_progress;
      
#pragma omp parallel for default(shared)
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	// Get variability
	int32_t varisize = getVariability(c, sampleLib[file_c]);
	
	// Qualities and alignment length
	TQualities qualities;
	TAlignmentLength alen;

	// Read alignments
	for(typename TSVs::const_iterator itSV = svs.begin(); itSV!=svs.end(); ++itSV) {
	  if (itSV->peSupport == 0) continue;
	  if ((itSV->chr != refIndex) && (itSV->chr2 != refIndex)) continue;
	  
	  // Scan left and right breakpoint
	  Breakpoint bp(*itSV);
	  _initBreakpoint(hdr, bp, varisize, itSV->svt);
	  unsigned int maxBp = 2;
	  if ((bp.chr == bp.chr2) && (bp.svStart + 3 * varisize >= bp.svEnd)) maxBp = 1;
	  for (unsigned int bpPoint = 0; bpPoint < maxBp; ++bpPoint) {
	    int32_t regionChr = 0;
	    int32_t regionStart = 0;
	    int32_t regionEnd = 0;
	    if (bpPoint == (unsigned int)(bp.chr == bp.chr2)) {
	      regionChr = bp.chr2;
	      regionStart = bp.svEndBeg;
	      regionEnd = bp.svEndEnd;
	    } else {
	      regionChr = bp.chr;
	      regionStart = bp.svStartBeg;
	      if (maxBp == 2) regionEnd = bp.svStartEnd;
	      else regionEnd = bp.svEndEnd;
	    }
	    hts_itr_t* iter = sam_itr_queryi(idx[file_c], regionChr, regionStart, regionEnd);
	    bam1_t* rec = bam_init1();
	    int32_t lastAlignedPos = 0;
	    std::set<std::size_t> lastAlignedPosReads;
	    while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	      if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP | BAM_FMUNMAP)) continue;
	      if (!(rec->core.flag & BAM_FPAIRED) || (rec->core.qual < c.minGenoQual)) continue;

	      // Clean-up the read store for identical alignment positions
	      if (rec->core.pos > lastAlignedPos) {
		lastAlignedPosReads.clear();
		lastAlignedPos = rec->core.pos;
	      }

	      // Mapping positions valid?
	      if (_mappingPosGeno(rec->core.tid, rec->core.mtid, rec->core.pos, rec->core.mpos, itSV->svt)) continue;

	      // Get the library information
	      typename TLibraryMap::iterator libIt = _findLib(rec, sampleLib[file_c]);
	      if (libIt->second.median == 0) continue; // Single-end library
		    
	      // Get or store the mapping quality for the partner
	      if (_firstPairObs(rec, lastAlignedPosReads)) {
		// First read
		lastAlignedPosReads.insert(hash_string(bam_get_qname(rec)));
		std::size_t hv = hash_pair(rec);
		if (_translocation(itSV->svt)) {
		  qualitiestra[file_c][hv]= rec->core.qual;
		  alentra[file_c][hv]= alignmentLength(rec);
		} else {
		  qualities[hv]= rec->core.qual;
		  alen[hv]= alignmentLength(rec);
		}
	      } else {
		// Second read
		std::size_t hv=hash_pair_mate(rec);
		int32_t alenmate = 0;
		uint8_t pairQuality = 0;
		if (_translocation(itSV->svt)) {
		  if (qualitiestra[file_c].find(hv) == qualitiestra[file_c].end()) continue; // Mate discarded
		  pairQuality = std::min((uint8_t) qualitiestra[file_c][hv], (uint8_t) rec->core.qual);
		  alenmate = alentra[file_c][hv];
		  qualitiestra[file_c][hv] = 0;
		} else {
		  if (qualities.find(hv) == qualities.end()) continue; // Mate discarded
		  pairQuality = std::min((uint8_t) qualities[hv], (uint8_t) rec->core.qual);
		  alenmate = alen[hv];
		  qualities[hv] = 0;
		}

		// Pair quality
		if (pairQuality < c.minGenoQual) continue;
		
		// Outer insert size
		int outerISize = 0;
		if (rec->core.pos < rec->core.mpos) outerISize = rec->core.mpos + alenmate - rec->core.pos;
		else outerISize = rec->core.pos + alignmentLength(rec) - rec->core.mpos;
		
		// Insert the interval
		if ((getSVType(rec->core) == libIt->second.defaultOrient) && (outerISize >= libIt->second.minNormalISize) && (outerISize <= libIt->second.maxNormalISize) && (rec->core.tid==rec->core.mtid)) {
		  // Normal spanning coverage
		  int32_t sPos = 0;
		  int32_t ePos = 0;
		  if (rec->core.pos < rec->core.mpos) {
		    sPos = rec->core.pos;
		    ePos = rec->core.mpos + alenmate;
		  } else {
		    sPos = rec->core.mpos;
		    ePos = rec->core.pos + alignmentLength(rec);
		  }
		  int32_t midPoint = sPos+(ePos-sPos)/2;
		  int32_t svidnum = -1;
		  if (std::abs(midPoint - itSV->svStart) < std::abs(itSV->svEnd - midPoint)) {
		    if ((itSV->chr == rec->core.tid) && (itSV->svStart > sPos + c.minimumFlankSize) && (itSV->svStart + c.minimumFlankSize < ePos)) svidnum = itSV->id;
		  } else {
		    if ((itSV->chr2 == rec->core.tid) && (itSV->svEnd > sPos + c.minimumFlankSize) && (itSV->svEnd + c.minimumFlankSize < ePos)) svidnum = itSV->id + lastId;
		  }
		  if (svidnum >= 0) {
		    uint8_t* hpptr = bam_aux_get(rec, "HP");
#pragma omp critical
		    {
		      spanCountMap[file_c][svidnum].ref.push_back(pairQuality);
		      if (hpptr) {
			c.isHaplotagged = true;
			int hap = bam_aux2i(hpptr);
			if (hap == 1) ++spanCountMap[file_c][svidnum].refh1;
			else ++spanCountMap[file_c][svidnum].refh2;
		      }
		    }
		  }
		} else if ((getSVType(rec->core) != libIt->second.defaultOrient) || (outerISize < libIt->second.minNormalISize) || (outerISize > libIt->second.maxNormalISize) || (rec->core.tid!=rec->core.mtid)) {

		  if (_isizeMappingPos(rec, varisize) != itSV->svt) continue;
		  if (!(((itSV->chr == rec->core.tid) && (itSV->chr2 == rec->core.mtid)) || ((itSV->chr == rec->core.mtid) && (itSV->chr2 == rec->core.tid)))) continue;
		  
		  // Does the pair confirm the SV
		  int32_t const minPos = _minCoord(rec->core.pos, rec->core.mpos, itSV->svt);
		  int32_t const maxPos = _maxCoord(rec->core.pos, rec->core.mpos, itSV->svt);
		  
		  if (rec->core.tid==itSV->chr) {
		    if (minPos < itSV->svStart) {
		      if (_pairsDisagree(minPos, maxPos, rec->core.l_qseq, libIt->second.maxNormalISize, itSV->svStart, itSV->svEnd, rec->core.l_qseq, libIt->second.maxNormalISize, itSV->svt)) continue;
		    } else {
		      if (_pairsDisagree(itSV->svStart, itSV->svEnd, rec->core.l_qseq, libIt->second.maxNormalISize, minPos, maxPos, rec->core.l_qseq, libIt->second.maxNormalISize, itSV->svt)) continue;
		    }
		  }
		  
		  int32_t svidnumLeft = -1;
		  int32_t svidnumRight = -1;
		  // Missing spanning coverage
		  if (_mateIsUpstream(libIt->second.defaultOrient, (rec->core.flag & BAM_FREAD1), (rec->core.flag & BAM_FREVERSE))) {
		    if ((itSV->chr==rec->core.tid) && (itSV->svStart>=rec->core.pos) && (itSV->svStart<=(rec->core.pos + libIt->second.maxNormalISize))) svidnumLeft = itSV->id;
		    if ((itSV->chr2==rec->core.tid) && (itSV->svEnd>=rec->core.pos) && (itSV->svEnd<=(rec->core.pos + libIt->second.maxNormalISize))) svidnumRight = itSV->id + lastId;
		  } else {
		    if ((itSV->chr==rec->core.tid) && (itSV->svStart>=std::max(0, rec->core.pos + (int) alignmentLength(rec) - libIt->second.maxNormalISize)) && (itSV->svStart<=(rec->core.pos + (int) alignmentLength(rec)))) svidnumLeft = itSV->id;
		    if ((itSV->chr2==rec->core.tid) && (itSV->svEnd>=std::max(0, rec->core.pos + (int) alignmentLength(rec) - libIt->second.maxNormalISize)) && (itSV->svEnd<=(rec->core.pos + (int) alignmentLength(rec)))) svidnumRight = itSV->id + lastId;
		  }
		  if (_mateIsUpstream(libIt->second.defaultOrient, !(rec->core.flag & BAM_FREAD1), (rec->core.flag & BAM_FMREVERSE))) {
		    if ((itSV->chr==rec->core.mtid) && (itSV->svStart>=rec->core.mpos) && (itSV->svStart<=(rec->core.mpos + libIt->second.maxNormalISize))) svidnumLeft = itSV->id;
		    if ((itSV->chr2==rec->core.mtid) && (itSV->svEnd>=rec->core.mpos) && (itSV->svEnd<=(rec->core.mpos + libIt->second.maxNormalISize)))  svidnumRight = itSV->id + lastId;
		  } else {
		    if ((itSV->chr==rec->core.mtid) && (itSV->svStart>=std::max(0, rec->core.mpos + (int) alenmate - libIt->second.maxNormalISize)) && (itSV->svStart<=(rec->core.mpos + (int) alenmate))) svidnumLeft = itSV->id;
		    if ((itSV->chr2==rec->core.mtid) && (itSV->svEnd>=std::max(0,rec->core.mpos + (int) alenmate - libIt->second.maxNormalISize)) && (itSV->svEnd<=(rec->core.mpos + (int) alenmate))) svidnumRight = itSV->id + lastId;
		  }
		  if ((svidnumLeft >= 0) && (svidnumRight >= 0)) {
		    uint8_t* hpptr = bam_aux_get(rec, "HP");
#pragma omp critical
		    {
		      if (c.pedumpflag) {
			std::string svid(_addID(itSV->svt));
			std::string padNumber = boost::lexical_cast<std::string>(itSV->id);
			padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
			svid += padNumber;
			dumpOut << svid << "\t" << c.files[file_c].string() << "\t" << bam_get_qname(rec) << "\t" << hdr->target_name[rec->core.tid] << "\t" << rec->core.pos << "\t" << hdr->target_name[rec->core.mtid] << "\t" << rec->core.mpos << "\t" << rec->core.qual << std::endl;
		      }
		      spanCountMap[file_c][svidnumLeft].alt.push_back(pairQuality);
		      if (hpptr) {
			c.isHaplotagged = true;
			int hap = bam_aux2i(hpptr);
			if (hap == 1) ++spanCountMap[file_c][svidnumLeft].alth1;
			else ++spanCountMap[file_c][svidnumLeft].alth2;
		      }
		      spanCountMap[file_c][svidnumRight].alt.push_back(pairQuality);
		      if (hpptr) {
			c.isHaplotagged = true;
			int hap = bam_aux2i(hpptr);
			if (hap == 1) ++spanCountMap[file_c][svidnumRight].alth1;
			else ++spanCountMap[file_c][svidnumRight].alth2;
		      }
		    }
		  }
		}
	      }
	    }
	    bam_destroy1(rec);
	    hts_itr_destroy(iter);
	  }
	}
      }
    }
    // Clean-up
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }

}

#endif
