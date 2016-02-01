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

#ifdef OPENMP
#include <omp.h>
#endif

#include "tags.h"
#include "version.h"
#include "util.h"

using namespace torali;

struct Config {
  uint8_t reqCT;
  uint32_t bpoffset;
  uint32_t minsize;
  uint32_t maxsize;
  float recoverlap;
  std::string svType;
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

void _remove_info_tag(bcf_hdr_t* hdr, bcf1_t* rec, std::string const& tag) {
  bcf_update_info(hdr, rec, tag.c_str(), NULL, 0, BCF_HT_INT);  // Type does not matter for n = 0
}

void _remove_format_tag(bcf_hdr_t* hdr, bcf1_t* rec, std::string const& tag) {
  bcf_update_format(hdr, rec, tag.c_str(), NULL, 0, BCF_HT_INT);  // Type does not matter for n = 0
}

void _remove_info(bcf_hdr_t* hdr, bcf1_t* rec) {
  std::string tmp[] = {"CT", "PRECISE", "IMPRECISE", "SVTYPE", "SVMETHOD", "CIEND", "CIPOS", "CHR2", "END", "PE", "MAPQ", "SR", "SRQ", "CONSENSUS"};
  std::set<std::string> keepInfo(tmp, tmp + sizeof(tmp)/sizeof(tmp[0]));

  if (!(rec->unpacked & BCF_UN_INFO)) bcf_unpack(rec, BCF_UN_INFO);
  
  for (int i = 0; i < rec->n_info; ++i){
    bcf_info_t* inf = &rec->d.info[i];
    const char* key = bcf_hdr_int2id(hdr, BCF_DT_ID, inf->key);
    if (keepInfo.find(std::string(key)) != keepInfo.end()) continue;
    
    if (inf->vptr_free) {
      free(inf->vptr - inf->vptr_off);
      inf->vptr_free = 0;
    }
    rec->d.shared_dirty |= BCF1_DIRTY_INF;
    inf->vptr = NULL;
  }
}

void _remove_format(bcf_hdr_t* hdr, bcf1_t* rec) {
  if (!(rec->unpacked & BCF_UN_FMT)) bcf_unpack(rec, BCF_UN_FMT);
  
  for(int i = 0; i<rec->n_fmt; ++i) {
    bcf_fmt_t* fmt = &rec->d.fmt[i];
    const char* key = bcf_hdr_int2id(hdr, BCF_DT_ID, fmt->id);
    bcf_update_format(hdr, rec, key, NULL, 0, BCF_HT_INT); // the type is irrelevant for n = 0
    // Keep GT
    //if ((key[0]=='G') && key[1]=='T' && (!key[2])) continue;

    if (fmt->p_free) {
      free(fmt->p - fmt->p_off);
      fmt->p_free = 0;
    }
    rec->d.indiv_dirty = 1;
    fmt->p = NULL;
  }
}

template<typename TGenomeIntervals, typename TContigMap, typename TSVType>
void _fillIntervalMap(Config const& c, TGenomeIntervals& iScore, TContigMap& cMap, TSVType svType) {
  typedef typename TGenomeIntervals::value_type TIntervalScores;
  typedef typename TIntervalScores::value_type IntervalScore;

  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Reading input VCF files" << std::endl;
  boost::progress_display show_progress( c.files.size() );


  boost::unordered_map<int32_t, std::string> refmap;
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    ++show_progress;
    htsFile* ifile = bcf_open(c.files[file_c].string().c_str(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);
    bcf1_t* rec = bcf_init();

    int32_t nsvend = 0;
    int32_t* svend = NULL;
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
      // Correct SV type
      bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt);
      if (std::string(svt) != _addID(svType)) continue;

      // Correct CT
      uint8_t ict = 0;
      if (bcf_get_info_string(hdr, rec, "CT", &ct, &nct) > 0) ict = _decodeOrientation(std::string(ct));
      if (ict != c.reqCT) continue;

      // Correct size
      std::string chrName(bcf_hdr_id2name(hdr, rec->rid));
      uint32_t tid = cMap[chrName];
      uint32_t svStart = rec->pos;
      uint32_t svEnd = svStart + 1;
      if (bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend) > 0) svEnd = *svend;

      // Parse INFO fields
      if ((svEnd - svStart < c.minsize) || (svEnd - svStart > c.maxsize)) continue;
      bool precise = false;
      if (bcf_get_info_flag(hdr, rec, "PRECISE", 0, 0) > 0) precise=true;
      unsigned int peSupport = 0;
      if (bcf_get_info_int32(hdr, rec, "PE", &pe, &npe) > 0) peSupport = *pe;
      unsigned int srSupport = 0;
      if (bcf_get_info_int32(hdr, rec, "SR", &sr, &nsr) > 0) srSupport = *sr;
      uint8_t peMapQuality = 0;
      if (bcf_get_info_int32(hdr, rec, "MAPQ", &mapq, &nmapq) > 0) peMapQuality = (uint8_t) *mapq;
      double srAlignQuality = 0;
      if (bcf_get_info_float(hdr, rec, "SRQ", &srq, &nsrq) > 0) srAlignQuality = (double) *srq;
      uint32_t mtid = tid;
      if (bcf_get_info_string(hdr, rec, "CHR2", &chr2, &nchr2) > 0) {
	std::string chr2Name(chr2);
	mtid = cMap[chr2Name];
      }
      if (mtid != tid) continue;

      // Proxy quality score for the SV
      uint32_t score = 0;
      if (precise) score = srSupport * (100 * srAlignQuality);
      else score = peSupport * (uint32_t) peMapQuality;
      
      // Store the interval
      iScore[tid].push_back(IntervalScore(svStart, svEnd, score));
    }
    free(svend);
    free(pe);
    free(sr);
    free(mapq);
    free(ct);
    free(srq);
    free(svt);
    free(chr2);
    bcf_hdr_destroy(hdr);
    bcf_close(ifile);
    bcf_destroy(rec);
  }
}

template<typename TGenomeIntervals>
void _processIntervalMap(Config const& c, TGenomeIntervals const& iScore, TGenomeIntervals& iSelected) {
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
      if (*iK) iSelected[seqId].push_back(IntervalScore(iS->start, iS->end, iS->score));
    }
  }
}

template<typename TGenomeIntervals, typename TContigMap, typename TSVType>
void _outputSelectedIntervals(Config const& c, TGenomeIntervals const& iSelected, TContigMap& cMap, TSVType svType) {
  typedef typename TGenomeIntervals::value_type TIntervalScores;
  typedef typename TIntervalScores::value_type IntervalScore;

  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Filtering SVs" << std::endl;
  boost::progress_display show_progress( c.files.size() );

  // Copy VCF header from first file and remove samples
  htsFile* ifile0 = bcf_open(c.files[0].string().c_str(), "r");
  bcf_hdr_t* hdr0 = bcf_hdr_read(ifile0);
  htsFile *fp = hts_open(c.outfile.string().c_str(), "wg");
  bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr0);
  // Combine the remaining headers
  for(unsigned int file_c = 1; file_c < c.files.size(); ++file_c) {
    htsFile* ifile = bcf_open(c.files[file_c].string().c_str(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);
    hdr_out = bcf_hdr_merge(hdr_out, hdr);
    bcf_hdr_destroy(hdr);
    bcf_close(ifile);
  }
  bcf_hdr_set_samples(hdr_out, NULL, false);
  bcf_hdr_write(fp, hdr_out);
  bcf_hdr_destroy(hdr0);
  bcf_close(ifile0);

  // Write VCF header
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    ++show_progress;
    htsFile* ifile = bcf_open(c.files[file_c].string().c_str(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);
    bcf_hdr_set_samples(hdr, NULL, false); // Do not read the sample information
    bcf1_t* rec = bcf_init();

    int32_t nsvend = 0;
    int32_t* svend = NULL;
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
      // Correct SV type
      bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt);
      if (std::string(svt) != _addID(svType)) continue;

      // Correct CT
      uint8_t ict = 0;
      if (bcf_get_info_string(hdr, rec, "CT", &ct, &nct) > 0) ict = _decodeOrientation(std::string(ct));
      if (ict != c.reqCT) continue;

      // Correct size
      std::string chrName(bcf_hdr_id2name(hdr, rec->rid));
      uint32_t tid = cMap[chrName];
      uint32_t svStart = rec->pos;
      uint32_t svEnd = svStart + 1;
      if (bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend) > 0) svEnd = *svend;

      // Parse INFO fields
      if ((svEnd - svStart < c.minsize) || (svEnd - svStart > c.maxsize)) continue;
      bool precise = false;
      if (bcf_get_info_flag(hdr, rec, "PRECISE", 0, 0) > 0) precise=true;
      unsigned int peSupport = 0;
      if (bcf_get_info_int32(hdr, rec, "PE", &pe, &npe) > 0) peSupport = *pe;
      unsigned int srSupport = 0;
      if (bcf_get_info_int32(hdr, rec, "SR", &sr, &nsr) > 0) srSupport = *sr;
      uint8_t peMapQuality = 0;
      if (bcf_get_info_int32(hdr, rec, "MAPQ", &mapq, &nmapq) > 0) peMapQuality = (uint8_t) *mapq;
      double srAlignQuality = 0;
      if (bcf_get_info_float(hdr, rec, "SRQ", &srq, &nsrq) > 0) srAlignQuality = (double) *srq;
      uint32_t mtid = tid;
      if (bcf_get_info_string(hdr, rec, "CHR2", &chr2, &nchr2) > 0) {
	std::string chr2Name(chr2);
	mtid = cMap[chr2Name];
      }
      if (mtid != tid) continue;

      // Proxy quality score for the SV
      uint32_t score = 0;
      if (precise) score = srSupport * (100 * srAlignQuality);
      else score = peSupport * (uint32_t) peMapQuality;
      
      typename TIntervalScores::const_iterator iter = std::lower_bound(iSelected[tid].begin(), iSelected[tid].end(), IntervalScore(svStart, svEnd, score), SortIScores<IntervalScore>());
      if ((iter != iSelected[tid].end()) && (iter->start == svStart) && (iter->end == svEnd) && (iter->score == score)) {
	rec->rid = bcf_hdr_name2id(hdr_out, chrName.c_str());

	// Remove spurious FORMAT and INFO tags
	_remove_format(hdr, rec);
	_remove_format(hdr_out, rec);
	_remove_info(hdr, rec);
	_remove_info(hdr_out, rec);

	// Write record
	bcf_write1(fp, hdr_out, rec);
	//std::cerr << bcf_hdr_id2name(hdr, tid) << '\t' << svStart << '\t' << svEnd << std::endl;
      }
    }

    free(svend);
    free(pe);
    free(sr);
    free(mapq);
    free(ct);
    free(srq);
    free(svt);
    free(chr2);
    bcf_hdr_destroy(hdr);
    bcf_close(ifile);
    bcf_destroy(rec);
  }

  // Close VCF file
  bcf_hdr_destroy(hdr_out);
  hts_close(fp);
}


template<typename TSVType>
inline int run(Config const& c, TSVType svType) {

  // All files may use a different set of chromosomes
  typedef boost::unordered_map<std::string, uint32_t> TContigMap;
  TContigMap contigMap;
  uint32_t numseq = 0;
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    const char **seqnames = NULL;
    int nseq = 0;

    htsFile* ifile = bcf_open(c.files[file_c].string().c_str(), "r");
    if (!ifile) {
      std::cerr << "Fail to load " << c.files[file_c].string() << "!" << std::endl;
      return 1;
    }
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);
    seqnames = bcf_hdr_seqnames(hdr, &nseq);
    if (nseq<=0) {
      std::cerr << "Fail to load chromosome names. Please index files with tabix." << std::endl;
      return 1;
    }
    for(int i = 0; i<nseq; ++i) {
      std::string chrName(seqnames[i]);
      if (contigMap.find(chrName) == contigMap.end()) contigMap[chrName] = numseq++;
    }
    // Clean-up
    if (seqnames != NULL) free(seqnames);

    bcf_hdr_destroy(hdr);
    bcf_close(ifile);
  }

  // Interval maps
  typedef std::vector<IntervalScore> TIntervalScores;
  typedef std::vector<TIntervalScores> TGenomeIntervals;
  TGenomeIntervals iScore;
  iScore.resize(numseq, TIntervalScores());
  _fillIntervalMap(c, iScore, contigMap, svType);
  for(int i = 0; i<numseq; ++i) std::sort(iScore[i].begin(), iScore[i].end(), SortIScores<IntervalScore>());

  // Filter intervals
  TGenomeIntervals iSelected;
  iSelected.resize(numseq, TIntervalScores());
  _processIntervalMap(c, iScore, iSelected);
  iScore.clear();
  for(int i = 0; i<numseq; ++i) std::sort(iSelected[i].begin(), iSelected[i].end(), SortIScores<IntervalScore>());

  // Output best intervals
  _outputSelectedIntervals(c, iSelected, contigMap, svType);

  // End
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;

  return 0;
}

int main(int argc, char **argv) {
  Config c;

  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("type,t", boost::program_options::value<std::string>(&c.svType)->default_value("DEL"), "SV type (DEL, DUP, INV, TRA, INS)")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("sv.vcf.gz"), "Merged SV VCF output file")
    ("minsize,m", boost::program_options::value<uint32_t>(&c.minsize)->default_value(0), "min. SV size")
    ("maxsize,n", boost::program_options::value<uint32_t>(&c.maxsize)->default_value(1000000), "max. SV size")
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
    ("license,l", "show license")
    ("warranty,w", "show warranty")
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
    printTitle("Delly SV VCF site merging");
    if (vm.count("warranty")) {
      displayWarranty();
    } else if (vm.count("license")) {
      gplV3();
    } else {
      std::cout << "Usage: " << argv[0] << " [OPTIONS] <sample1.vcf> <sample2.vcf> ..." << std::endl;
      std::cout << visible_options << "\n"; 
    }
    return 1; 
  }

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;
 
  // Run merging
  if (c.svType == "DEL") {
    c.reqCT = 2;
    return run(c, SVType<DeletionTag>());
  } else if (c.svType == "DUP") {
    c.reqCT = 3;
    return run(c, SVType<DuplicationTag>());
  }
  //else if (c.svType == "INV") return run(c, SVType<InversionTag>());
  //else if (c.svType == "TRA") return run(c, SVType<TranslocationTag>());
  //else if (c.svType == "INS") return run(c, SVType<InsertionTag>());
  else {
    std::cerr << "SV analysis type not supported by Delly: " << c.svType << std::endl;
    return 1;
  }
}
