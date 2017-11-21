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

#ifndef MERGE_H
#define MERGE_H

#include <iostream>
#include <fstream>
#include <boost/unordered_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/icl/interval_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/progress.hpp>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "tags.h"
#include "version.h"
#include "util.h"
#include "modvcf.h"


namespace torali
{


struct MergeConfig {
  bool filterForPass;
  bool filterForPrecise;
  uint32_t svcounter;
  uint32_t bpoffset;
  uint32_t minsize;
  uint32_t maxsize;
  float recoverlap;
  boost::filesystem::path outfile;
  std::vector<boost::filesystem::path> files;
};

struct IntervalScore {
  uint32_t start;
  uint32_t end;
  uint32_t score;
  
  IntervalScore(uint32_t s, uint32_t e, uint32_t c) : start(s), end(e), score(c) {}
};

template<typename TRecord>
struct SortIScores : public std::binary_function<TRecord, TRecord, bool>
{
  inline bool operator()(TRecord const& s1, TRecord const& s2) const {
    return ((s1.start < s2.start) || ((s1.start == s2.start) && (s1.end < s2.end)));
  }

};

template<typename TPos>
double recOverlap(TPos const s1, TPos const e1, TPos const s2, TPos const e2) {
  if ((e1 < s2) || (s1 > e2)) return 0;
  double lenA = (double) (e1-s1);
  if (lenA <= 0) return 0;
  double lenB = (double) (e2-s2);
  if (lenB <= 0) return 0;
  double overlapLen = double(std::min(e1, e2) - std::max(s1, s2));
  if (overlapLen <= 0) return 0;
  return (overlapLen / std::max(lenA, lenB));
}


template<typename TGenomeIntervals, typename TContigMap>
void _fillIntervalMap(MergeConfig const& c, TGenomeIntervals& iScore, TContigMap& cMap, int32_t const svtin) {
  typedef typename TGenomeIntervals::value_type TIntervalScores;
  typedef typename TIntervalScores::value_type IntervalScore;

  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Reading input VCF/BCF files" << std::endl;
  boost::progress_display show_progress( c.files.size() );


  boost::unordered_map<int32_t, std::string> refmap;
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    ++show_progress;
    htsFile* ifile = bcf_open(c.files[file_c].string().c_str(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);
    bcf1_t* rec = bcf_init();

    int32_t nsvend = 0;
    int32_t* svend = NULL;
    int32_t ninslen = 0;
    int32_t* inslen = NULL;
    int32_t npe = 0;
    int32_t* pe = NULL;
    int32_t nsr = 0;
    int32_t* sr = NULL;
    int32_t nmapq = 0;
    int32_t* mapq = NULL;
    int32_t nct = 0;
    char* ct = NULL;
    int32_t nsrq = 0;
    float* srq = NULL;
    int32_t nsvt = 0;
    char* svt = NULL;
    int32_t nchr2 = 0;
    char* chr2 = NULL;
    while (bcf_read(ifile, hdr, rec) == 0) {
      bcf_unpack(rec, BCF_UN_INFO);
      // Check PASS
      bool pass = true;
      if (c.filterForPass) pass = (bcf_has_filter(hdr, rec, const_cast<char*>("PASS"))==1);
      if (!pass) continue;

      // Correct SV type
      int32_t recsvt = -1;
      if ((bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt) > 0) && (bcf_get_info_string(hdr, rec, "CT", &ct, &nct) > 0)) recsvt = _decodeOrientation(std::string(ct), std::string(svt));
      if (recsvt != svtin) continue;

      // Correct size
      std::string chrName(bcf_hdr_id2name(hdr, rec->rid));
      uint32_t tid = cMap[chrName];
      uint32_t svStart = rec->pos;
      uint32_t svEnd = svStart + 1;
      if (bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend) > 0) svEnd = *svend;

      // SV size check 
      if (std::string(svt) != "BND") {
	if (std::string(svt) == "INS") {
	  unsigned int inslenVal = 0;
	  if (bcf_get_info_int32(hdr, rec, "INSLEN", &inslen, &ninslen) > 0) inslenVal = *inslen;
	  if ((inslenVal < c.minsize) || (inslenVal > c.maxsize)) continue;
	} else {
	  if ((svEnd - svStart < c.minsize) || (svEnd - svStart > c.maxsize)) continue;
	}
      }
      bool precise = false;
      if (bcf_get_info_flag(hdr, rec, "PRECISE", 0, 0) > 0) precise=true;
      if ((c.filterForPrecise) && (!precise)) continue;

      unsigned int peSupport = 0;
      if (bcf_get_info_int32(hdr, rec, "PE", &pe, &npe) > 0) peSupport = *pe;
      unsigned int srSupport = 0;
      if (bcf_get_info_int32(hdr, rec, "SR", &sr, &nsr) > 0) srSupport = *sr;
      // Remove this line
      //if (srSupport > 0) precise = true;

      uint8_t peMapQuality = 0;
      if (bcf_get_info_int32(hdr, rec, "MAPQ", &mapq, &nmapq) > 0) peMapQuality = (uint8_t) *mapq;
      float srAlignQuality = 0;
      if (bcf_get_info_float(hdr, rec, "SRQ", &srq, &nsrq) > 0) srAlignQuality = *srq;
      if (bcf_get_info_string(hdr, rec, "CHR2", &chr2, &nchr2) > 0) {
	std::string chr2Name(chr2);
	//mtid = cMap[chr2Name];
      }

      // Quality score for the SV
      uint32_t score = 0;
      if (_isKeyPresent(hdr, "SCORE")) {
	int32_t nvcfscore = 0;
	if (_getInfoType(hdr, "SCORE") == BCF_HT_INT) {
	  int32_t* vcfscore = NULL;
	  bcf_get_info_int32(hdr, rec, "SCORE", &vcfscore, &nvcfscore);
	  score = *vcfscore;
	  free(vcfscore);
	} else if (_getInfoType(hdr, "SCORE") == BCF_HT_REAL) {
	  float* vcfscore = NULL;
	  bcf_get_info_float(hdr, rec, "SCORE", &vcfscore, &nvcfscore);
	  score = (uint32_t) (*vcfscore * 10000); // for scores in [0,1] 
	  free(vcfscore);
	}
      } else {
	if (precise) score = srSupport * (uint32_t) (100 * srAlignQuality);
	else score = peSupport * (uint32_t) peMapQuality;
      }

      // Store the interval
      iScore[tid].push_back(IntervalScore(svStart, svEnd, score));
    }
    if (svend != NULL) free(svend);
    if (inslen != NULL) free(inslen);
    if (pe != NULL) free(pe);
    if (sr != NULL) free(sr);
    if (mapq != NULL) free(mapq);
    if (ct != NULL) free(ct);
    if (srq != NULL) free(srq);
    if (svt != NULL) free(svt);
    if (chr2 != NULL) free(chr2);
    bcf_hdr_destroy(hdr);
    bcf_close(ifile);
    bcf_destroy(rec);
  }
}

template<typename TGenomeIntervals>
void _processIntervalMap(MergeConfig const& c, TGenomeIntervals const& iScore, TGenomeIntervals& iSelected, int32_t const svtin) {
  typedef typename TGenomeIntervals::value_type TIntervalScores;
  typedef typename TIntervalScores::value_type IntervalScore;

  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Merging SVs" << std::endl;
  boost::progress_display show_progress( iScore.size() );

  unsigned int seqId = 0;
  for(typename TGenomeIntervals::const_iterator iG = iScore.begin(); iG != iScore.end(); ++iG, ++seqId) {
    ++show_progress;
    typedef std::vector<bool> TIntervalSelector;
    TIntervalSelector keepInterval;
    keepInterval.resize(iG->size(), true);
    typename TIntervalSelector::iterator iK = keepInterval.begin();
    for(typename TIntervalScores::const_iterator iS = iG->begin(); iS != iG->end(); ++iS, ++iK) {
      typename TIntervalScores::const_iterator iSNext = iS;
      typename TIntervalSelector::iterator iKNext = iK;
      ++iSNext; ++iKNext;
      for(; iSNext != iG->end(); ++iSNext, ++iKNext) {
	if (iSNext->start - iS->start > c.bpoffset) break;
	else {
	  if (((iSNext->end > iS->end) && (iSNext->end - iS->end < c.bpoffset)) || ((iSNext->end <= iS->end) &&(iS->end - iSNext->end < c.bpoffset))) {
	    if ((_translocation(svtin)) || (recOverlap(iS->start, iS->end, iSNext->start, iSNext->end) >= c.recoverlap)) {
	      if (iS->score < iSNext->score) *iK = false;
	      else if (iSNext ->score < iS->score) *iKNext = false;
	      else {
		if (iS->start < iSNext->start) *iKNext = false;
		else if (iS->end < iSNext->end) *iKNext = false;
		else *iK = false;
	      }
	    }
	  }
	}
      }
      if (*iK) iSelected[seqId].push_back(IntervalScore(iS->start, iS->end, iS->score));
    }
  }
}

template<typename TGenomeIntervals, typename TContigMap>
void _outputSelectedIntervals(MergeConfig& c, TGenomeIntervals const& iSelected, TContigMap& cMap, int32_t const svtin) {
  typedef typename TGenomeIntervals::value_type TIntervalScores;
  typedef typename TIntervalScores::value_type IntervalScore;

  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Filtering SVs" << std::endl;

  // Open output VCF file
  htsFile *fp = hts_open(c.outfile.string().c_str(), "wb");
  bcf_hdr_t *hdr_out = bcf_hdr_init("w");

  // Write VCF header
  boost::gregorian::date today = now.date();
  std::string datestr("##fileDate=");
  datestr += boost::gregorian::to_iso_string(today);
  bcf_hdr_append(hdr_out, datestr.c_str());
  bcf_hdr_append(hdr_out, "##ALT=<ID=DEL,Description=\"Deletion\">");
  bcf_hdr_append(hdr_out, "##ALT=<ID=DUP,Description=\"Duplication\">");
  bcf_hdr_append(hdr_out, "##ALT=<ID=INV,Description=\"Inversion\">");
  bcf_hdr_append(hdr_out, "##ALT=<ID=BND,Description=\"Translocation\">");
  bcf_hdr_append(hdr_out, "##ALT=<ID=INS,Description=\"Insertion\">");
  bcf_hdr_append(hdr_out, "##FILTER=<ID=LowQual,Description=\"PE/SR support below 3 or mapping quality below 20.\">");
  bcf_hdr_append(hdr_out, "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"PE confidence interval around END\">");
  bcf_hdr_append(hdr_out, "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"PE confidence interval around POS\">");
  bcf_hdr_append(hdr_out, "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">");
  bcf_hdr_append(hdr_out, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">");
  bcf_hdr_append(hdr_out, "##INFO=<ID=PE,Number=1,Type=Integer,Description=\"Paired-end support of the structural variant\">");
  bcf_hdr_append(hdr_out, "##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of paired-ends\">");
  bcf_hdr_append(hdr_out, "##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Split-read support\">");
  bcf_hdr_append(hdr_out, "##INFO=<ID=SRQ,Number=1,Type=Float,Description=\"Split-read consensus alignment quality\">");
  bcf_hdr_append(hdr_out, "##INFO=<ID=CONSENSUS,Number=1,Type=String,Description=\"Split-read consensus sequence\">");
  bcf_hdr_append(hdr_out, "##INFO=<ID=CE,Number=1,Type=Float,Description=\"Consensus sequence entropy\">");
  bcf_hdr_append(hdr_out, "##INFO=<ID=CT,Number=1,Type=String,Description=\"Paired-end signature induced connection type\">");
  bcf_hdr_append(hdr_out, "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">");
  bcf_hdr_append(hdr_out, "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">");
  bcf_hdr_append(hdr_out, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">");
  bcf_hdr_append(hdr_out, "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">");
  bcf_hdr_append(hdr_out, "##INFO=<ID=INSLEN,Number=1,Type=Integer,Description=\"Predicted length of the insertion\">");
  bcf_hdr_append(hdr_out, "##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description=\"Predicted microhomology length using a max. edit distance of 2\">");
  // Add reference contigs
  uint32_t numseq = 0;
  typedef std::map<uint32_t, std::string> TReverseMap;
  TReverseMap rMap;
  for(typename TContigMap::iterator cIt = cMap.begin(); cIt != cMap.end(); ++cIt, ++numseq) rMap[cIt->second] = cIt->first;
  for(typename TReverseMap::iterator rIt = rMap.begin(); rIt != rMap.end(); ++rIt) {
    std::string refname("##contig=<ID=");
    refname += rIt->second + ">";
    bcf_hdr_append(hdr_out, refname.c_str());
  }
  bcf_hdr_add_sample(hdr_out, NULL);
  bcf_hdr_write(fp, hdr_out);

  // Duplicate filter (identical start, end, score values)
  typedef std::pair<uint32_t, uint32_t> TStartEnd;
  typedef std::set<TStartEnd> TIntervalSet;
  typedef std::vector<TIntervalSet> TGenomicIntervalSet;
  TGenomicIntervalSet gis(numseq);

  // Parse input VCF files
  bcf1_t *rout = bcf_init();
  typedef std::vector<htsFile*> THtsFile;
  typedef std::vector<bcf_hdr_t*> TBcfHeader;
  typedef std::vector<bcf1_t*> TBcfRecord;
  typedef std::vector<bool> TEof;
  THtsFile ifile(c.files.size());
  TBcfHeader hdr(c.files.size());
  TBcfRecord rec(c.files.size());
  TEof eof(c.files.size());
  uint32_t allEOF = 0;
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    ifile[file_c] = bcf_open(c.files[file_c].string().c_str(), "r");
    hdr[file_c] = bcf_hdr_read(ifile[file_c]);
    bcf_hdr_set_samples(hdr[file_c], NULL, false); // Do not read the sample information
    rec[file_c] = bcf_init();
    if (bcf_read(ifile[file_c], hdr[file_c], rec[file_c]) == 0) {
      bcf_unpack(rec[file_c], BCF_UN_INFO);
      eof[file_c] = false;
    } else {
      ++allEOF;
      eof[file_c] = true;
    }
  }

  int32_t nsvend = 0;
  int32_t* svend = NULL;
  int32_t npe = 0;
  int32_t* pe = NULL;
  int32_t nsr = 0;
  int32_t* sr = NULL;
  int32_t ninslen = 0;
  int32_t* inslen = NULL;
  int32_t nhomlen = 0;
  int32_t* homlen = NULL;
  int32_t nmapq = 0;
  int32_t* mapq = NULL;
  int32_t nct = 0;
  char* ct = NULL;
  int32_t nsrq = 0;
  float* srq = NULL;
  int32_t nsvt = 0;
  char* svt = NULL;
  int32_t nchr2 = 0;
  char* chr2 = NULL;
  int32_t ncipos = 0;
  int32_t* cipos = NULL;
  int32_t nciend = 0;
  int32_t* ciend = NULL;
  int32_t nce = 0;
  float* ce = NULL;
  int32_t ncons = 0;
  char* cons = NULL;
  while (allEOF < c.files.size()) {
    // Find next sorted record
    int32_t idx = -1;
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      if (!eof[file_c]) {
	if ((idx < 0) || (rec[idx]->rid > rec[file_c]->rid) || ((rec[idx]->rid == rec[file_c]->rid) && (rec[idx]->pos > rec[file_c]->pos))) idx = file_c;
      }
    }

    // Correct SV type
    int32_t recsvt = -1;
    if ((bcf_get_info_string(hdr[idx], rec[idx], "SVTYPE", &svt, &nsvt) > 0) && (bcf_get_info_string(hdr[idx], rec[idx], "CT", &ct, &nct) > 0)) recsvt = _decodeOrientation(std::string(ct), std::string(svt));
    if (recsvt == svtin) {
      // Check PASS
      bool pass = true;
      if (c.filterForPass) pass = (bcf_has_filter(hdr[idx], rec[idx], const_cast<char*>("PASS"))==1);

      // Check PRECISE
      bool precise = false;
      bool passPrecise = true;
      if (bcf_get_info_flag(hdr[idx], rec[idx], "PRECISE", 0, 0) > 0) precise=true;
      if ((c.filterForPrecise) && (!precise)) passPrecise = false;

      // Check PASS and precise
      if ((passPrecise) && (pass)) {
	// Correct size
	std::string chrName(bcf_hdr_id2name(hdr[idx], rec[idx]->rid));
	uint32_t tid = cMap[chrName];
	uint32_t svStart = rec[idx]->pos;
	uint32_t svEnd = svStart + 1;
	if (bcf_get_info_int32(hdr[idx], rec[idx], "END", &svend, &nsvend) > 0) svEnd = *svend;
	unsigned int inslenVal = 0;
	if (bcf_get_info_int32(hdr[idx], rec[idx], "INSLEN", &inslen, &ninslen) > 0) inslenVal = *inslen;

	// Parse INFO fields
	if ((std::string(svt) == "BND") || ((std::string(svt) == "INS") && (inslenVal >= c.minsize) && (inslenVal <= c.maxsize)) || ((std::string(svt) != "BND") && (std::string(svt) != "INS") && (svEnd - svStart >= c.minsize) && (svEnd - svStart <= c.maxsize))) {
	  unsigned int peSupport = 0;
	  if (bcf_get_info_int32(hdr[idx], rec[idx], "PE", &pe, &npe) > 0) peSupport = *pe;
	  unsigned int srSupport = 0;
	  if (bcf_get_info_int32(hdr[idx], rec[idx], "SR", &sr, &nsr) > 0) srSupport = *sr;
	  // Remove this line
	  //if (srSupport > 0) precise = true;
	  
	  uint8_t peMapQuality = 0;
	  if (bcf_get_info_int32(hdr[idx], rec[idx], "MAPQ", &mapq, &nmapq) > 0) peMapQuality = (uint8_t) *mapq;
	  float srAlignQuality = 0;
	  if (bcf_get_info_float(hdr[idx], rec[idx], "SRQ", &srq, &nsrq) > 0) srAlignQuality = *srq;
	  std::string chr2Name = chrName;
	  if (bcf_get_info_string(hdr[idx], rec[idx], "CHR2", &chr2, &nchr2) > 0) {
	    chr2Name = std::string(chr2);
	    //mtid = cMap[chr2Name];
	  }

	  // Proxy quality score for the SV
	  uint32_t score = 0;
	  if (_isKeyPresent(hdr[idx], "SCORE")) {
	    int32_t nvcfscore = 0;
	    if (_getInfoType(hdr[idx], "SCORE") == BCF_HT_INT) {
	      int32_t* vcfscore = NULL;
	      bcf_get_info_int32(hdr[idx], rec[idx], "SCORE", &vcfscore, &nvcfscore);
	      score = *vcfscore;
	      free(vcfscore);
	    } else if (_getInfoType(hdr[idx], "SCORE") == BCF_HT_REAL) {
	      float* vcfscore = NULL;
	      bcf_get_info_float(hdr[idx], rec[idx], "SCORE", &vcfscore, &nvcfscore);
	      score = (uint32_t) (*vcfscore * 10000); // for scores in [0,1] 
	      free(vcfscore);
	    }
	  } else {
	    if (precise) score = srSupport * (uint32_t) (100 * srAlignQuality);
	    else score = peSupport * (uint32_t) peMapQuality;
	  }
	  
	  // Is this a selected interval
	  typename TIntervalScores::const_iterator iter = std::lower_bound(iSelected[tid].begin(), iSelected[tid].end(), IntervalScore(svStart, svEnd, score), SortIScores<IntervalScore>());
	  bool foundInterval = false;
	  for(; (iter != iSelected[tid].end()) && (iter->start == svStart); ++iter) {
	    if ((iter->start == svStart) && (iter->end == svEnd) && (iter->score == score)) {
	      // Duplicate?
	      if (gis[tid].find(std::make_pair(svStart, svEnd)) == gis[tid].end()) {
		foundInterval = true;
		gis[tid].insert(std::make_pair(svStart, svEnd));
	      }
	      break;
	    }
	  }
	  if (foundInterval) {
	    // Fetch missing INFO fields
	    unsigned int homlenVal = 0;
	    if (bcf_get_info_int32(hdr[idx], rec[idx], "HOMLEN", &homlen, &nhomlen) > 0) homlenVal = *homlen;
	    bcf_get_info_int32(hdr[idx], rec[idx], "CIPOS", &cipos, &ncipos);
	    bcf_get_info_int32(hdr[idx], rec[idx], "CIEND", &ciend, &nciend);
	    std::string consensus;
	    float ceVal = 0;
	    if (precise) {
	      if (bcf_get_info_float(hdr[idx], rec[idx], "CE", &ce, &nce) > 0) ceVal = *ce;
	      bcf_get_info_string(hdr[idx], rec[idx], "CONSENSUS", &cons, &ncons);
	      consensus = boost::to_upper_copy(std::string(cons));
	    }
	    
	    // Create new record
	    rout->rid = bcf_hdr_name2id(hdr_out, chrName.c_str());
	    rout->pos = rec[idx]->pos;
	    rout->qual = 0;
	    std::string id;
	    if (c.files.size() == 1) id = std::string(rec[idx]->d.id); // Within one VCF file IDs are unique
	    else {
	      id += _addID(svtin);
	      std::string padNumber = boost::lexical_cast<std::string>(c.svcounter++);
	      padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
	      id += padNumber;
	    }
	    bcf_update_id(hdr_out, rout, id.c_str());
	    std::string refAllele = rec[idx]->d.allele[0];
	    std::string altAllele = rec[idx]->d.allele[1];
	    std::string alleles = refAllele + "," + altAllele;
	    bcf_update_alleles_str(hdr_out, rout, alleles.c_str());
	    int32_t tmppass = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "PASS");
	    bcf_update_filter(hdr_out, rout, &tmppass, 1);
	    
	    // Add INFO fields
	    if (precise) bcf_update_info_flag(hdr_out, rout, "PRECISE", NULL, 1);
	    else bcf_update_info_flag(hdr_out, rout, "IMPRECISE", NULL, 1);
	    bcf_update_info_string(hdr_out, rout, "SVTYPE", _addID(svtin).c_str());
	    std::string dellyVersion("EMBL.DELLYv");
	    dellyVersion += dellyVersionNumber;
	    bcf_update_info_string(hdr_out,rout, "SVMETHOD", dellyVersion.c_str());
	    bcf_update_info_string(hdr_out,rout, "CHR2", chr2Name.c_str());
	    bcf_update_info_int32(hdr_out, rout, "END", &svEnd, 1);
	    bcf_update_info_int32(hdr_out, rout, "PE", &peSupport, 1);
	    int32_t tmpi = peMapQuality;
	    bcf_update_info_int32(hdr_out, rout, "MAPQ", &tmpi, 1);
	    bcf_update_info_string(hdr_out, rout, "CT", _addOrientation(svtin).c_str());
	    bcf_update_info_int32(hdr_out, rout, "CIPOS", cipos, 2);
	    bcf_update_info_int32(hdr_out, rout, "CIEND", ciend, 2);
	    if (precise) {
	      bcf_update_info_int32(hdr_out, rout, "INSLEN", &inslenVal, 1);
	      bcf_update_info_int32(hdr_out, rout, "HOMLEN", &homlenVal, 1);
	      bcf_update_info_int32(hdr_out, rout, "SR", &srSupport, 1);
	      bcf_update_info_float(hdr_out, rout, "SRQ", &srAlignQuality, 1);	
	      bcf_update_info_string(hdr_out, rout, "CONSENSUS", consensus.c_str());
	      bcf_update_info_float(hdr_out, rout, "CE", &ceVal, 1);
	    }
	
	    // Write record
	    bcf_write1(fp, hdr_out, rout);
	    bcf_clear1(rout);	  
	    //std::cerr << bcf_hdr_id2name(hdr[idx], tid) << '\t' << svStart << '\t' << svEnd << std::endl;
	  }
	}
      }
    }

    // Fetch next record
    if (bcf_read(ifile[idx], hdr[idx], rec[idx]) == 0) bcf_unpack(rec[idx], BCF_UN_INFO);
    else {
      ++allEOF;
      eof[idx] = true;
    }
  }
  if (svend != NULL) free(svend);
  if (pe != NULL) free(pe);
  if (sr != NULL) free(sr);
  if (homlen != NULL) free(homlen);
  if (inslen != NULL) free(inslen);
  if (mapq != NULL) free(mapq);
  if (ct != NULL) free(ct);
  if (srq != NULL) free(srq);
  if (svt != NULL) free(svt);
  if (chr2 != NULL) free(chr2);
  if (cipos != NULL) free(cipos);
  if (ciend != NULL) free(ciend);
  if (ce != NULL) free(ce);
  if (cons != NULL) free(cons);

  // Clean-up
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    bcf_hdr_destroy(hdr[file_c]);
    bcf_close(ifile[file_c]);
    bcf_destroy(rec[file_c]);
  }

  // Close VCF file
  bcf_destroy(rout);
  bcf_hdr_destroy(hdr_out);
  hts_close(fp);

  // Build index
  bcf_index_build(c.outfile.string().c_str(), 14);
}

inline void
mergeBCFs(MergeConfig& c, std::vector<boost::filesystem::path> const& cts) {
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Merging SV types" << std::endl;
  boost::progress_display show_progress( 1 );

  // Parse temporary input VCF files
  typedef std::vector<htsFile*> THtsFile;
  typedef std::vector<bcf_hdr_t*> TBcfHeader;
  typedef std::vector<bcf1_t*> TBcfRecord;
  typedef std::vector<bool> TEof;
  THtsFile ifile(cts.size());
  TBcfHeader hdr(cts.size());
  TBcfRecord rec(cts.size());
  TEof eof(cts.size());
  uint32_t allEOF = 0;
  for(unsigned int file_c = 0; file_c < cts.size(); ++file_c) {
    ifile[file_c] = bcf_open(cts[file_c].string().c_str(), "r");
    hdr[file_c] = bcf_hdr_read(ifile[file_c]);
    rec[file_c] = bcf_init();
    if (bcf_read(ifile[file_c], hdr[file_c], rec[file_c]) == 0) {
      bcf_unpack(rec[file_c], BCF_UN_INFO);
      eof[file_c] = false;
    } else {
      ++allEOF;
      eof[file_c] = true;
    }
  }

  // Open output VCF file
  htsFile *fp = hts_open(c.outfile.string().c_str(), "wb");
  bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr[0]);
  bcf_hdr_write(fp, hdr_out);

  // Merge files
  while (allEOF < cts.size()) {
    // Find next sorted record
    int32_t idx = -1;
    for(unsigned int file_c = 0; file_c < cts.size(); ++file_c) {
      if (!eof[file_c]) {
	if ((idx < 0) || (rec[idx]->rid > rec[file_c]->rid) || ((rec[idx]->rid == rec[file_c]->rid) && (rec[idx]->pos > rec[file_c]->pos))) idx = file_c;
      }
    }

    // Write record
    bcf_write1(fp, hdr_out, rec[idx]);

    // Fetch next record
    if (bcf_read(ifile[idx], hdr[idx], rec[idx]) == 0) bcf_unpack(rec[idx], BCF_UN_INFO);
    else {
      ++allEOF;
      eof[idx] = true;
    }
  }
  ++show_progress;
  
  // Clean-up
  for(unsigned int file_c = 0; file_c < cts.size(); ++file_c) {
    bcf_hdr_destroy(hdr[file_c]);
    bcf_close(ifile[file_c]);
    bcf_destroy(rec[file_c]);
  }

  // Close VCF file
  bcf_hdr_destroy(hdr_out);
  hts_close(fp);

  // Build index
  bcf_index_build(c.outfile.string().c_str(), 14);

  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
}

inline int
mergeRun(MergeConfig& c, int32_t const svt) {

  // All files may use a different set of chromosomes
  typedef std::map<std::string, uint32_t> TContigMap;
  TContigMap contigMap;
  uint32_t numseq = 0;
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    htsFile* ifile = bcf_open(c.files[file_c].string().c_str(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);
    const char** seqnames = NULL;
    int nseq=0;
    seqnames = bcf_hdr_seqnames(hdr, &nseq);
    for(int32_t i = 0; i<nseq;++i) {
      std::string chrName(bcf_hdr_id2name(hdr, i));
      if (contigMap.find(chrName) == contigMap.end()) contigMap[chrName] = numseq++;
    }
    if (seqnames!=NULL) free(seqnames);
    bcf_hdr_destroy(hdr);
    bcf_close(ifile);
  }

  // Interval maps
  typedef std::vector<IntervalScore> TIntervalScores;
  typedef std::vector<TIntervalScores> TGenomeIntervals;
  TGenomeIntervals iScore;
  iScore.resize(numseq, TIntervalScores());
  _fillIntervalMap(c, iScore, contigMap, svt);
  for(uint32_t i = 0; i<numseq; ++i) std::sort(iScore[i].begin(), iScore[i].end(), SortIScores<IntervalScore>());

  // Filter intervals
  TGenomeIntervals iSelected;
  iSelected.resize(numseq, TIntervalScores());
  _processIntervalMap(c, iScore, iSelected, svt);
  iScore.clear();
  for(uint32_t i = 0; i<numseq; ++i) std::sort(iSelected[i].begin(), iSelected[i].end(), SortIScores<IntervalScore>());

  // Output best intervals
  _outputSelectedIntervals(c, iSelected, contigMap, svt);

  // End
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;

  return 0;
}

int merge(int argc, char **argv) {
  MergeConfig c;
  c.svcounter = 1;

  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("sv.bcf"), "Merged SV BCF output file")
    ("minsize,m", boost::program_options::value<uint32_t>(&c.minsize)->default_value(0), "min. SV size")
    ("maxsize,n", boost::program_options::value<uint32_t>(&c.maxsize)->default_value(1000000), "max. SV size")
    ("precise,c", "Filter sites for PRECISE")
    ("pass,p", "Filter sites for PASS")
    ;

  // Define overlap options
  boost::program_options::options_description overlap("Overlap options");
  overlap.add_options()
    ("bp-offset,b", boost::program_options::value<uint32_t>(&c.bpoffset)->default_value(1000), "max. breakpoint offset")
    ("rec-overlap,r", boost::program_options::value<float>(&c.recoverlap)->default_value(0.8), "min. reciprocal overlap")
    ;

  // Define hidden options
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
    ;
  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  // Set the visibility
  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(overlap).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic).add(overlap);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);


  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file"))) { 
    std::cout << std::endl;
    std::cout << "Usage: delly " << argv[0] << " [OPTIONS] <sample1.bcf> <sample2.bcf> ..." << std::endl;
    std::cout << visible_options << "\n"; 
    return 0; 
  }

  // Filter for PASS
  if (vm.count("pass")) c.filterForPass = true;
  else c.filterForPass = false;

  // Filter for PRECISE
  if (vm.count("precise")) c.filterForPrecise = true;
  else c.filterForPrecise = false;

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  std::cout << "delly ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Check input BCF files
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    htsFile* ifile = bcf_open(c.files[file_c].string().c_str(), "r");
    if (!ifile) {
      std::cerr << "Fail to load " << c.files[file_c].string() << "!" << std::endl;
      return 1;
    }
    bcf_close(ifile);
  }
  
  // Run merging
  boost::filesystem::path oldPath = c.outfile;
  int rVal = 0;
  int32_t maxSvt = 15;
  std::vector<boost::filesystem::path> svtCollect(maxSvt);
  for(int32_t svt = 0; svt < maxSvt; ++svt) {
    svtCollect[svt] = boost::filesystem::unique_path();
    c.outfile = svtCollect[svt];
    rVal += mergeRun(c, svt);
  }
  // Merge temporary files
  c.outfile = oldPath;
  mergeBCFs(c, svtCollect);
  for(int32_t svt = 0; svt < maxSvt; ++svt) {
    boost::filesystem::remove(svtCollect[svt]);
    boost::filesystem::remove(boost::filesystem::path(svtCollect[svt].string() + ".csi"));
  }
  return rVal;
}

}

#endif

