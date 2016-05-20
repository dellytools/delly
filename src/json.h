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

#ifndef JSON_H
#define JSON_H

#include <iostream>
#include <fstream>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>


namespace torali
{

  template<typename TConfig, typename TStructuralVariantRecord, typename TJunctionCountMap, typename TReadCountMap, typename TCountMap, typename TTag>
    inline void
    jsonOutput(TConfig const& c, std::vector<TStructuralVariantRecord> const& svs, TJunctionCountMap const& jctCountMap, TReadCountMap const& readCountMap, TCountMap const& spanCountMap, SVType<TTag> svType) 
  {
    // Typedefs
    typedef typename TCountMap::key_type TSampleSVPair;

    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Genotyping" << std::endl;
    boost::progress_display show_progress( svs.size() );

    // Open one bam file header
    samFile* samfile = sam_open(c.files[0].string().c_str(), "r");
    bam_hdr_t* bamhd = sam_hdr_read(samfile);
    
    // Output all structural variants
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    
    // Genotype arrays
    int32_t *gts = (int*) malloc(c.files.size() * 2 * sizeof(int));
    float *gls = (float*) malloc(c.files.size() * 3 * sizeof(float));
    int32_t *rcl = (int*) malloc(c.files.size() * sizeof(int));
    int32_t *rc = (int*) malloc(c.files.size() * sizeof(int));
    int32_t *rcr = (int*) malloc(c.files.size() * sizeof(int));
    int32_t *cnest = (int*) malloc(c.files.size() * sizeof(int));
    int32_t *drcount = (int*) malloc(c.files.size() * sizeof(int));
    int32_t *dvcount = (int*) malloc(c.files.size() * sizeof(int));
    int32_t *rrcount = (int*) malloc(c.files.size() * sizeof(int));
    int32_t *rvcount = (int*) malloc(c.files.size() * sizeof(int));
    int32_t *gqval = (int*) malloc(c.files.size() * sizeof(int));
    std::vector<std::string> ftarr;
    ftarr.resize(c.files.size());

    // Iterate all structural variants
    typedef std::vector<TStructuralVariantRecord> TSVs;
    for(typename TSVs::const_iterator svIter = svs.begin(); svIter!=svs.end(); ++svIter) {
      ++show_progress;
      
      // Set values
      std::string filterVal = "PASS";
      if ((svIter->precise) && (svIter->chr == svIter->chr2)) {
	if ((svIter->srSupport < 3) || (svIter->peMapQuality < 20)) filterVal = "LowQual";
      } else {
	if ((svIter->peSupport < 3) || (svIter->peMapQuality < 20) || ( (svIter->chr != svIter->chr2) && (svIter->peSupport < 5) ) ) filterVal = "LowQual";
      }
      std::string preciseVal = "false";
      if (svIter->precise) preciseVal = "true";
      std::string dellyVersion("EMBL.DELLYv");
      dellyVersion += dellyVersionNumber;

      // Output main fields
      dataOut << '{';
      dataOut << "\"chrom\": \"" << bamhd->target_name[svIter->chr] << "\", ";
      dataOut << "\"start\": " << svIter->svStart - 1 << ", ";
      if (_addID(svType) != "TRA") dataOut << "\"end\": " << svIter->svEnd << ", ";
      else dataOut << "\"end\": " << svIter->svStart << ", ";
      dataOut << "\"id\": " << svIter->id << ", ";
      if (_addID(svType) == "TRA") {
	dataOut << "\"$children\": [{";
	dataOut << "\"chrom\": \"" << bamhd->target_name[svIter->chr2] << "\", ";
	dataOut << "\"start\": " << svIter->svEnd << ", ";
	dataOut << "\"end\": " << svIter->svEnd + 1 << ", ";
	dataOut << "\"id\": " << (svIter->id + svs.size()) << "";
	dataOut << "}] ,";
      }
      dataOut << "\"ref\": \"N\", \"alt\": \"" << _addID(svType) << "\", ";
      dataOut << "\"qual\": 0, ";
      dataOut << "\"filter\": \"" << filterVal << "\", ";
      dataOut << "\"info\": {";
      dataOut << "\"precise\": " << preciseVal << ", ";
      dataOut << "\"svtype\": \"" << _addID(svType) << "\", ";
      dataOut << "\"svmethod\": \"" << dellyVersion << "\", ";
      dataOut << "\"inslen\": " << (int) _addInsertionLength(svIter->insLen, svType) << ", ";
      dataOut << "\"pe\": " << svIter->peSupport << ", ";
      dataOut << "\"mapq\": " << (int) (svIter->peMapQuality) << ", ";
      dataOut << "\"cipos\": [" << -svIter->wiggle << ", " << svIter->wiggle << "], ";
      dataOut << "\"ciend\": [" << -svIter->wiggle << ", " << svIter->wiggle << "], ";
      dataOut << "\"ct\": \"" << _addOrientation(svIter->ct) << "\"";
      if (svIter->precise)  {
	dataOut << ", ";
	dataOut << "\"sr\": " << svIter->srSupport << ", ";
	dataOut << "\"srq\": " << (float) (svIter->srAlignQuality) << ", ";
	dataOut << "\"consensus\": \"" << svIter->consensus << "\", ";
	dataOut << "\"ce\": " << (float) (entropy(svIter->consensus));
      }
      dataOut << "}, ";


      // Add genotype columns
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	TSampleSVPair sampleSVPairLeft = std::make_pair(file_c, svIter->id);
	TSampleSVPair sampleSVPairRight = std::make_pair(file_c, -svIter->id);
	typename TJunctionCountMap::const_iterator jctCountMapIt=jctCountMap.find(sampleSVPairLeft);
	typename TCountMap::const_iterator spanLeftIt=spanCountMap.find(sampleSVPairLeft);
	typename TCountMap::const_iterator spanRightIt=spanCountMap.find(sampleSVPairRight);

	// Counters
	rcl[file_c] = 0;
	rc[file_c] = 0;
	rcr[file_c] = 0;
	cnest[file_c] = 0;
	drcount[file_c] = 0;
	dvcount[file_c] = 0;
	rrcount[file_c] = 0;
	rvcount[file_c] = 0;
	if ((spanLeftIt!=spanCountMap.end()) && (spanRightIt!=spanCountMap.end())) {
	  if (spanLeftIt->second.first.size()<spanRightIt->second.first.size()) {
	    drcount[file_c] = spanLeftIt->second.first.size();
	    dvcount[file_c] = spanLeftIt->second.second.size();
	  } else {
	    drcount[file_c] = spanRightIt->second.first.size();
	    dvcount[file_c] = spanRightIt->second.second.size();
	  }
	}
	if (jctCountMapIt!=jctCountMap.end()) {
	  rrcount[file_c] = jctCountMapIt->second.first.size();
	  rvcount[file_c] = jctCountMapIt->second.second.size();
	}
	
	// Compute GLs
	if (svIter->precise) {
	  if (jctCountMapIt!=jctCountMap.end()) _computeGLs(jctCountMapIt->second.first, jctCountMapIt->second.second, gls, gqval, gts, file_c);
	  else {
	    gls[file_c * 3 + 2] = 0;
	    gls[file_c * 3 + 1] = 0;
	    gls[file_c * 3] = 0;
	    gqval[file_c] = 0;
	    gts[file_c * 2] = bcf_gt_missing;
	    gts[file_c * 2 + 1] = bcf_gt_missing;
	  }
	} else {  // Imprecise SVs
	  if ((spanLeftIt!=spanCountMap.end()) && (spanRightIt!=spanCountMap.end())) {
	    if (spanLeftIt->second.first.size()<spanRightIt->second.first.size()) _computeGLs(spanLeftIt->second.first, spanLeftIt->second.second, gls, gqval, gts, file_c);
	    else _computeGLs(spanRightIt->second.first, spanRightIt->second.second, gls, gqval, gts, file_c);
	  } else {
	    gls[file_c * 3 + 2] = 0;
	    gls[file_c * 3 + 1] = 0;
	    gls[file_c * 3] = 0;
	    gqval[file_c] = 0;
	    gts[file_c * 2] = bcf_gt_missing;
	    gts[file_c * 2 + 1] = bcf_gt_missing;
	  }
	}
	
	// Compute RCs
	typename TReadCountMap::const_iterator readCountMapIt=readCountMap.find(sampleSVPairLeft);
	if (readCountMapIt!=readCountMap.end()) {
	  rcl[file_c] = readCountMapIt->second.leftRC;
	  rc[file_c] = readCountMapIt->second.rc;
	  rcr[file_c] = readCountMapIt->second.rightRC;
	}	
	cnest[file_c] = -1;
	if ((rcl[file_c] + rcr[file_c]) > 0) cnest[file_c] = boost::math::iround( 2.0 * (double) rc[file_c] / (double) (rcl[file_c] + rcr[file_c]) );
      
	// Genotype filter
	if (gqval[file_c] < 15) ftarr[file_c] = "LowQual";
	else ftarr[file_c] = "PASS";
      }

      // Output genotypes
      dataOut << "\"genotypes\": {";
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	dataOut << "\"" << c.sampleName[file_c] << "\": {";
	dataOut << "\"gt\": [" << bcf_gt_allele(gts[file_c*2]) << ", " << bcf_gt_allele(gts[file_c*2+1]) << "], ";
	dataOut << "\"gl\": [" << gls[file_c*3] << ", " << gls[file_c*3+1] << ", " << gls[file_c*3+2] << "], ";
	dataOut << "\"gq\": " << gqval[file_c] << ", ";
	dataOut << "\"ft\": \"" << ftarr[file_c] << "\", ";
	dataOut << "\"rcl\": " << rcl[file_c] << ", ";
	dataOut << "\"rc\": " << rc[file_c] << ", ";
	dataOut << "\"rcr\": " << rcr[file_c] << ", ";
	dataOut << "\"cn\": " << cnest[file_c] << ", ";
	dataOut << "\"dr\": " << drcount[file_c] << ", ";
	dataOut << "\"dv\": " << dvcount[file_c] << ", ";
	dataOut << "\"rr\": " << rrcount[file_c] << ", ";
	dataOut << "\"rv\": " << rvcount[file_c] << "";
	dataOut << "}";
      }
      dataOut << "}";
      dataOut << '}' << std::endl;

      // Add children
      if (_addID(svType) == "TRA") {
	dataOut << '{';
	dataOut << "\"chrom\": \"" << bamhd->target_name[svIter->chr2] << "\", ";
	dataOut << "\"start\": " << svIter->svEnd << ", ";
	dataOut << "\"end\": " << svIter->svEnd + 1 << ", ";
	dataOut << "\"id\": " << (svIter->id + svs.size()) << ", ";
	dataOut << "\"$parent\": {";
	dataOut << "\"chrom\": \"" << bamhd->target_name[svIter->chr] << "\", ";
	dataOut << "\"start\": " << svIter->svStart - 1 << ", ";
	dataOut << "\"end\": " << svIter->svStart << ", ";
	dataOut << "\"id\": " << svIter->id << "";
	dataOut << "}";
	dataOut << '}' << std::endl;
      }
    }
    
    // Clean-up
    free(gts);
    free(gls);
    free(rcl);
    free(rc);
    free(rcr);
    free(cnest);
    free(drcount);
    free(dvcount);
    free(rrcount);
    free(rvcount);
    free(gqval);
    
    // Close BAM file
    bam_hdr_destroy(bamhd);
    sam_close(samfile);
  }  

}

#endif
