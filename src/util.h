#ifndef UTIL_H
#define UTIL_H

#include <boost/multi_array.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <htslib/sam.h>
#include <sstream>
#include <math.h>

#include "edlib.h"
#include "tags.h"


namespace torali
{

  #ifndef LAST_BIN
  #define LAST_BIN 65535
  #endif

  #ifndef MAX_CN
  #define MAX_CN 10
  #endif
  
  struct LibraryInfo {
    int32_t rs;
    int32_t median;
    int32_t mad;
    int32_t minNormalISize;
    int32_t minISizeCutoff;
    int32_t maxNormalISize;
    int32_t maxISizeCutoff;
    uint32_t abnormal_pairs;

    LibraryInfo() : rs(0), median(0), mad(0), minNormalISize(0), minISizeCutoff(0), maxNormalISize(0), maxISizeCutoff(0), abnormal_pairs(0) {}
  };

  struct CNV {
    int32_t chr;
    int32_t start;
    int32_t end;
    int32_t ciposlow;
    int32_t ciposhigh;
    int32_t ciendlow;
    int32_t ciendhigh;
    int32_t qval;
    double cn;
    double mappable;
    double sd;


    CNV() : chr(0), start(0), end(0), ciposlow(0), ciposhigh(0), ciendlow(0), ciendhigh(0), qval(0), cn(-1), mappable(0), sd(1) {}
    CNV(int32_t const c, int32_t const s, int32_t const e, int32_t const cil, int32_t const cih, int32_t const cel, int32_t ceh, double const estcn, double const mp) : chr(c), start(s), end(e), ciposlow(cil), ciposhigh(cih), ciendlow(cel), ciendhigh(ceh), qval(0), cn(estcn), mappable(mp), sd(1) {}
  };

  template<typename TCNV>
  struct SortCNVs : public std::binary_function<TCNV, TCNV, bool>
  {
    inline bool operator()(TCNV const& sv1, TCNV const& sv2) {
      return ((sv1.chr<sv2.chr) || ((sv1.chr==sv2.chr) && (sv1.start<sv2.start)) || ((sv1.chr==sv2.chr) && (sv1.start==sv2.start) && (sv1.end<sv2.end)) || ((sv1.chr==sv2.chr) && (sv1.start==sv2.start) && (sv1.end==sv2.end) && (sv1.cn < sv2.cn)));
    }
  };




  // Read count struct
  struct ReadCount {
    int32_t leftRC;
    int32_t rc;
    int32_t rightRC;

    ReadCount() {}
    ReadCount(int32_t l, int32_t m, int32_t r) : leftRC(l), rc(m), rightRC(r) {}
  };


  inline uint32_t
  infixStart(EdlibAlignResult& cigar) {
    int32_t tIdx = cigar.endLocations[0];
    for (int32_t i = 0; i < cigar.alignmentLength; i++) {
      if (cigar.alignment[i] != EDLIB_EDOP_INSERT) tIdx--;
    }
    if (tIdx >= 0) return tIdx + 1;
    else return 0;
  }

  inline uint32_t
  infixEnd(EdlibAlignResult& cigar) {
    return cigar.endLocations[0];
  }
  
  inline void
  printAlignmentPretty(std::string const& query, std::string const& target, EdlibAlignMode const modeCode, EdlibAlignResult& align) {
    int32_t tIdx = -1;
    int32_t qIdx = -1;
    if (modeCode == EDLIB_MODE_HW) {
        tIdx = align.endLocations[0];
        for (int32_t i = 0; i < align.alignmentLength; i++) {
            if (align.alignment[i] != EDLIB_EDOP_INSERT) tIdx--;
        }
    }
    std::cerr << std::endl;
    for (int start = 0; start < align.alignmentLength; start += 50) {
      std::cerr << "T: ";
      int32_t startTIdx = -1;
      for (int32_t j = start; ((j < start + 50) && (j < align.alignmentLength)); ++j) {
	if (align.alignment[j] == EDLIB_EDOP_INSERT) std::cerr << "-";
	else std::cerr << target[++tIdx];
	if (j == start) startTIdx = tIdx;
      }
      std::cerr << " (" << std::max(startTIdx, 0) << " - " << tIdx << ")" << std::endl;

      // match / mismatch
      std::cerr << ("   ");
      for (int32_t j = start; j < start + 50 && j < align.alignmentLength; j++) {
	if (align.alignment[j] == EDLIB_EDOP_MATCH) std::cerr <<  "|";
	else std::cerr << " ";
      }
      std::cerr << std::endl;

      // query
      std::cerr << "Q: ";
      int32_t startQIdx = qIdx;
      for (int32_t j = start; j < start + 50 && j < align.alignmentLength; j++) {
	if (align.alignment[j] == EDLIB_EDOP_DELETE) std::cerr << "-";
	else std::cerr << query[++qIdx];
	if (j == start) startQIdx = qIdx;
      }
      std::cerr << " ("<< std::max(startQIdx, 0) << " - " << qIdx << ")" << std::endl;
      std::cerr << std::endl;
    }
  }

  inline void
  printAlignment(std::string const& seqI, std::string const& seqJ, EdlibAlignMode const modeCode, EdlibAlignResult& cigar) {
    int32_t tIdx = -1;
    int32_t qIdx = -1;
    uint32_t missingEnd = 0;
    uint32_t missingStart = 0;
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      tIdx = cigar.endLocations[0];
      if (tIdx < (int32_t) seqJ.size()) missingEnd = seqJ.size() - tIdx - 1;
      for (int32_t i = 0; i < cigar.alignmentLength; i++) {
	if (cigar.alignment[i] != EDLIB_EDOP_INSERT) tIdx--;
      }
      if (tIdx >= 0) missingStart = tIdx + 1;
    }
    // infix alignment, fix start
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      if (missingStart) {
	for (uint32_t j = 0; j < missingStart; ++j) std::cerr << '-';
      }
    }
    // seqI
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_DELETE) std::cerr << '-';
      else std::cerr << seqI[++qIdx];
    }
    // infix alignment, fix end
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      if (missingEnd) {
	for (uint32_t j = 0; j < missingEnd; ++j) std::cerr << '-';
      }
    }
    std::cerr << std::endl;
    // infix alignment, fix start
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      if (missingStart) {
	for (uint32_t j = 0; j < missingStart; ++j) std::cerr << seqJ[j];
      }
    }
    // seqJ
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_INSERT) std::cerr << '-';
      else std::cerr << seqJ[++tIdx];
    }
    // infix alignment, fix end
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      if (missingEnd) {
	for (uint32_t j = 0; j < missingEnd; ++j) std::cerr << seqJ[++tIdx];
      }
    }
    std::cerr << std::endl;
  }

  
  template<typename TConfig>
  inline void
  checkSampleNames(TConfig& c) {
    uint32_t ucount = 0;
    std::set<std::string> snames;
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      while (snames.find(c.sampleName[file_c]) != snames.end()) {
	std::cerr << "Warning: Duplicate sample names: " << c.sampleName[file_c] << std::endl;
	c.sampleName[file_c] += "_" + boost::lexical_cast<std::string>(ucount++);
	std::cerr << "Warning: Changing sample name to " << c.sampleName[file_c] << std::endl;
      }
      snames.insert(c.sampleName[file_c]);
    }
  }

  inline bool
  nContent(std::string const& s) {
    for(uint32_t i = 0; i < s.size(); ++i) {
      if ((s[i] == 'N') || (s[i] == 'n')) return true;
    }
    return false;
  }
  
  // Decode Orientation
  inline int32_t
    _decodeOrientation(std::string const& value) {
    if (value=="3to3") return 0;
    else if (value=="5to5") return 1;
    else if (value=="3to5") return 2;
    else if (value=="5to3") return 3;
    else return 4;
  }
  
  // Decode Orientation
  inline int32_t
  _decodeOrientation(std::string const& value, std::string const& svt) {
    if (svt == "BND") {
      if (value=="3to3") return DELLY_SVT_TRANS + 0;
      else if (value=="5to5") return DELLY_SVT_TRANS + 1;
      else if (value=="3to5") return DELLY_SVT_TRANS + 2;
      else if (value=="5to3") return DELLY_SVT_TRANS + 3;
      else return -1;
    } else if (svt == "CNV") {
      return 9;
    } else {
      if (value=="3to3") return 0;
      else if (value=="5to5") return 1;
      else if (value=="3to5") return 2;
      else if (value=="5to3") return 3;
      else return 4;
    }
  }
  
  
  // Deletions
  inline std::string
  _addID(int32_t const svt) {
    if (svt == 0) return "INV";
    else if (svt == 1) return "INV";
    else if (svt == 2) return "DEL";
    else if (svt == 3) return "DUP";
    else if (svt == 4) return "INS";
    else if (svt == 9) return "CNV";
    else return "BND";
  }

  inline std::string
  _addAlleles(std::string const& ref, std::string const& alt) {
    return ref + "," + alt;
  }
      
  inline std::string
  _addAlleles(std::string const& ref, std::string const& chr2, StructuralVariantRecord const& sv, int32_t const svt) {
    if (_translocation(svt)) {
      uint8_t ct = _getSpanOrientation(svt);
      if (ct == 0) {
	return ref + "," + ref + "]" + chr2 + ":" + boost::lexical_cast<std::string>(sv.svEnd) + "]";
      } else if (ct == 1) {
	return ref + "," + "[" + chr2 + ":" + boost::lexical_cast<std::string>(sv.svEnd) + "[" + ref;
      } else if (ct == 2) {
	return ref + "," + ref + "[" + chr2 + ":" + boost::lexical_cast<std::string>(sv.svEnd) + "[";
      } else if (ct == 3) {
	return ref + "," + "]" + chr2 + ":" + boost::lexical_cast<std::string>(sv.svEnd) + "]" + ref;
      } else {
	return ref + ",<" + _addID(svt) + ">";
      }
    } else return ref + ",<" + _addID(svt) + ">";
  }

  
  // Add Orientation
  inline std::string
  _addOrientation(int32_t const svt) {
    uint8_t ct = _getSpanOrientation(svt);
    if (ct==0) return "3to3";
    else if (ct==1) return "5to5";
    else if (ct==2) return "3to5";
    else if (ct==3) return "5to3";
    else return "NtoN";
  }

  inline bool
  is_gz(boost::filesystem::path const& f) {
    std::ifstream bfile(f.string().c_str(), std::ios_base::binary | std::ios::ate);
    bfile.seekg(0, std::ios::beg);
    char byte1;
    bfile.read(&byte1, 1);
    char byte2;
    bfile.read(&byte2, 1);
    bfile.close();
    if ((byte1 == '\x1F') && (byte2 == '\x8B')) return true;
    else return false;
  }

  // Output directory/file checks
  inline bool
  _outfileValid(boost::filesystem::path const& outfile) {
    try {
      boost::filesystem::path outdir;
      if (outfile.has_parent_path()) outdir = outfile.parent_path();
      else outdir = boost::filesystem::current_path();
      if (!boost::filesystem::exists(outdir)) {
	std::cerr << "Output directory does not exist: " << outdir << std::endl;
	return false;
      } else {
	boost::filesystem::file_status s = boost::filesystem::status(outdir);
	boost::filesystem::ofstream file(outfile.string());
	file.close();
	if (!(boost::filesystem::exists(outfile) && boost::filesystem::is_regular_file(outfile))) {
	  std::cerr << "Fail to open output file " << outfile.string() << std::endl;
	  std::cerr << "Output directory permissions: " << s.permissions() << std::endl;
	  return false;
	} else {
	  boost::filesystem::remove(outfile.string());
	}
      }
    } catch (boost::filesystem::filesystem_error const& e) {
      std::cerr << e.what() << std::endl;
      return false;
    }
    return true;
  }

  template<typename TConfig>
  inline bool
    _svTypesToCompute(TConfig& c, std::string const& svtype) {
    if (svtype == "ALL") return true;
    else {
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep(",");
      Tokenizer tokens(svtype, sep);
      for(Tokenizer::iterator tokIter = tokens.begin(); tokIter!=tokens.end(); ++tokIter) {
	if (*tokIter == "DEL") {
	  c.svtset.insert(2);
	} else if (*tokIter == "INS") {
	  c.svtset.insert(4);
	} else if (*tokIter == "DUP") {
	  c.svtset.insert(3);
	} else if (*tokIter == "INV") {
	  c.svtset.insert(0);
	  c.svtset.insert(1);
	} else if (*tokIter == "INV_3to3") {
	  c.svtset.insert(0);
	} else if (*tokIter == "INV_5to5") {
	  c.svtset.insert(1);
	} else if (*tokIter == "BND") {
	  c.svtset.insert(DELLY_SVT_TRANS + 0);
	  c.svtset.insert(DELLY_SVT_TRANS + 1);
	  c.svtset.insert(DELLY_SVT_TRANS + 2);
	  c.svtset.insert(DELLY_SVT_TRANS + 3);
	} else if (*tokIter == "BND_3to3") {
	  c.svtset.insert(DELLY_SVT_TRANS + 0);
	} else if (*tokIter == "BND_5to5") {
	  c.svtset.insert(DELLY_SVT_TRANS + 1);
	} else if (*tokIter == "BND_3to5") {
	  c.svtset.insert(DELLY_SVT_TRANS + 2);
	} else if (*tokIter == "BND_5to3") {
	  c.svtset.insert(DELLY_SVT_TRANS + 3);
	} else {
	  return false;
	}
      }
    }
    return true;    
  }

  inline uint32_t sequenceLength(bam1_t const* rec) {
    uint32_t* cigar = bam_get_cigar(rec);
    uint32_t slen = 0;
    for (uint32_t i = 0; i < rec->core.n_cigar; ++i)
      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF) || (bam_cigar_op(cigar[i]) == BAM_CINS) || (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) slen += bam_cigar_oplen(cigar[i]);
    return slen;
  }

  inline int32_t
  readLength(bam1_t const* rec) {
    //int32_t slen = rec->core.l_qseq;  # Incorrect for seq. with hard-clips
    return sequenceLength(rec);
  }
    
  inline uint32_t alignmentLength(bam1_t const* rec) {
    uint32_t* cigar = bam_get_cigar(rec);
    uint32_t alen = 0;
    for (uint32_t i = 0; i < rec->core.n_cigar; ++i)
      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF) || (bam_cigar_op(cigar[i]) == BAM_CDEL) || (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP)) alen += bam_cigar_oplen(cigar[i]);
    return alen;
  }

  inline double
  _sharedCarriers(std::vector<int32_t> const& gt1, std::vector<int32_t> const& gt2) {
    // Percentage of shared carriers
    uint32_t carnum = 0;
    uint32_t carshared = 0;
    for(uint32_t k = 0; k < gt1.size(); ++k) {
      if ((gt1[k] != 0) || (gt2[k] != 0)) {
	++carnum;
	if ((gt1[k] != 0) && (gt2[k] != 0)) ++carshared;
      }
    }
    return (double) (carshared) / (double) (carnum);
  }

  inline double
  nonrefGtConc(std::vector<int32_t> const& gt1, std::vector<int32_t> const& gt2) {
    uint32_t totgt = 0;
    uint32_t matchgt = 0;
    for(uint32_t k = 0; k < gt1.size(); ++k) {
      if ((gt1[k] != 0) || (gt2[k] != 0)) {
	++totgt;
	if (gt1[k] == gt2[k]) ++matchgt;
      }
    }
    return (double) (matchgt) / (double) (totgt);
  }
  
  inline double
  gtConc(std::vector<int32_t> const& gt1, std::vector<int32_t> const& gt2) {
    uint32_t totgt = 0;
    uint32_t matchgt = 0;
    for(uint32_t k = 0; k < gt1.size(); ++k) {
      ++totgt;
      if (gt1[k] == gt2[k]) ++matchgt;
    }
    return (double) (matchgt) / (double) (totgt);
  }
  
  inline uint32_t halfAlignmentLength(bam1_t const* rec) {
    return (alignmentLength(rec) / 2);
  }

  inline uint32_t
  lastAlignedPosition(bam1_t const* rec) {
    return rec->core.pos + alignmentLength(rec);
  }
  
  inline std::size_t hash_pair(bam1_t* rec) {
    std::size_t seed = hash_string(bam_get_qname(rec));
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    boost::hash_combine(seed, rec->core.mtid);
    boost::hash_combine(seed, rec->core.mpos);
    return seed;
  }

  inline std::size_t hash_pair_mate(bam1_t* rec) {
    std::size_t seed = hash_string(bam_get_qname(rec));
    boost::hash_combine(seed, rec->core.mtid);
    boost::hash_combine(seed, rec->core.mpos);
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    return seed;
  }

  inline std::size_t hash_sr(bam1_t* rec) {
    boost::hash<std::string> string_hash;
    std::string qname = bam_get_qname(rec);
    std::size_t seed = hash_string(qname.c_str());
    boost::hash_combine(seed, string_hash(qname));
    if ((rec->core.flag & BAM_FREAD1) && (seed > 0))  --seed;
    return seed;
  }

  inline std::size_t hash_lr(bam1_t* rec) {
    boost::hash<std::string> string_hash;
    std::string qname = bam_get_qname(rec);
    std::size_t seed = hash_string(qname.c_str());
    boost::hash_combine(seed, string_hash(qname));
    return seed;
  }

  inline std::size_t hash_lr(std::string const& qname) {
    boost::hash<std::string> string_hash;
    std::size_t seed = hash_string(qname.c_str());
    boost::hash_combine(seed, string_hash(qname));
    return seed;
  }
  
  inline void
  reverseComplement(std::string& sequence) {
    std::string rev = boost::to_upper_copy(std::string(sequence.rbegin(), sequence.rend()));
    std::size_t i = 0;
    for(std::string::iterator revIt = rev.begin(); revIt != rev.end(); ++revIt, ++i) {
      switch (*revIt) {
      case 'A': sequence[i]='T'; break;
      case 'C': sequence[i]='G'; break;
      case 'G': sequence[i]='C'; break;
      case 'T': sequence[i]='A'; break;
      case 'N': sequence[i]='N'; break;
      default: break;
      }
    }
  }

  inline std::string
  compressStr(std::string const& data) {
    std::stringstream compressed;
    std::stringstream origin(data);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> out;
    out.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_speed)));
    out.push(origin);
    boost::iostreams::copy(out, compressed);
    return compressed.str();
  }

  inline std::string
  decompressStr(std::string const& data) {
    std::stringstream compressed(data);
    std::stringstream decompressed;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> out;
    out.push(boost::iostreams::gzip_decompressor());
    out.push(compressed);
    boost::iostreams::copy(out, decompressed);
    return decompressed.str();
  }

  inline double
  entropy(std::string const& st) {
    typedef double TPrecision;
    std::vector<char> stvec(st.begin(), st.end());
    std::set<char> alphabet(stvec.begin(), stvec.end());
    TPrecision ent = 0;
    for(std::set<char>::const_iterator c = alphabet.begin(); c != alphabet.end(); ++c) {
      int ctr = 0;
      for (std::vector<char>::const_iterator s = stvec.begin(); s != stvec.end(); ++s)
	if (*s == *c) ++ctr;
      TPrecision freq = (TPrecision) ctr / (TPrecision) stvec.size();
      ent += (freq) * log(freq)/log(2);
    }
    return -ent;
  }


  inline uint32_t
  setMinChrLen(bam_hdr_t const* hdr, double const xx) {
    uint32_t minChrLen = 0;
    std::vector<uint32_t> chrlen(hdr->n_targets, 0);
    uint64_t genomelen = 0;
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      chrlen[refIndex] = hdr->target_len[refIndex];
      genomelen += hdr->target_len[refIndex];
    }
    std::sort(chrlen.begin(), chrlen.end(), std::greater<uint32_t>());
    uint64_t cumsum = 0;
    for(uint32_t i = 0; i < chrlen.size(); ++i) {
      cumsum += chrlen[i];
      minChrLen = chrlen[i];
      if (cumsum > genomelen * xx) break;
    }
    return minChrLen;
  }
  
  
  template<typename TConfig>
  inline bool
  chrNoData(TConfig const& c, uint32_t const refIndex, hts_idx_t const* idx) {
    // Check we have mapped reads on this chromosome
    std::string suffix("cram");
    std::string str(c.bamFile.string());
    if ((str.size() >= suffix.size()) && (str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0)) return false;
    uint64_t mapped = 0;
    uint64_t unmapped = 0;
    hts_idx_get_stat(idx, refIndex, &mapped, &unmapped);
    if (mapped) return false;
    else return true;
  }    
  
  inline std::size_t hash_se(bam1_t* rec) {
    std::size_t seed = hash_string(bam_get_qname(rec));
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    return seed;
  }
  
  inline void
  getSMTag(std::string const& header, std::string const& fileName, std::string& sampleName) {
    std::set<std::string> smIdentifiers;
    std::string delimiters("\n");
    typedef std::vector<std::string> TStrParts;
    TStrParts lines;
    boost::split(lines, header, boost::is_any_of(delimiters));
    TStrParts::const_iterator itH = lines.begin();
    TStrParts::const_iterator itHEnd = lines.end();
    bool rgPresent = false;
    for(;itH!=itHEnd; ++itH) {
      if (itH->find("@RG")==0) {
	std::string delim("\t");
	TStrParts keyval;
	boost::split(keyval, *itH, boost::is_any_of(delim));
	TStrParts::const_iterator itKV = keyval.begin();
	TStrParts::const_iterator itKVEnd = keyval.end();
	for(;itKV != itKVEnd; ++itKV) {
	  size_t sp = itKV->find(":");
	  if (sp != std::string::npos) {
	    std::string field = itKV->substr(0, sp);
	    if (field == "SM") {
	      rgPresent = true;
	      std::string rgSM = itKV->substr(sp+1);
	      smIdentifiers.insert(rgSM);
	    }
	  }
	}
      }
    }
    if (!rgPresent) {
      sampleName = fileName;
    } else if (smIdentifiers.size() == 1) {
      sampleName = *(smIdentifiers.begin());
    } else if (smIdentifiers.size() > 1) {
      sampleName = *(smIdentifiers.begin());
      std::cerr << "Warning: Multiple sample names (@RG:SM) present in the BAM file!" << std::endl;
    }
  }

  template<typename TConfig, typename TRegionsGenome>
  inline int32_t
  _parseExcludeIntervals(TConfig const& c, bam_hdr_t* hdr, TRegionsGenome& validRegions) {
    typedef typename TRegionsGenome::value_type TChrIntervals;
    typedef typename TChrIntervals::interval_type TIVal;
    
    validRegions.resize(hdr->n_targets);
    TRegionsGenome exclg;
    exclg.resize(hdr->n_targets);
    std::vector<bool> validChr;
    validChr.resize(hdr->n_targets, true);
    if (c.hasExcludeFile) {
      std::ifstream chrFile(c.exclude.string().c_str(), std::ifstream::in);
      if (chrFile.is_open()) {
	while (chrFile.good()) {
	  std::string chrFromFile;
	  getline(chrFile, chrFromFile);
	  typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	  boost::char_separator<char> sep(" \t,;");
	  Tokenizer tokens(chrFromFile, sep);
	  Tokenizer::iterator tokIter = tokens.begin();
	  if (tokIter!=tokens.end()) {
	    std::string chrName = *tokIter++;
	    int32_t tid = bam_name2id(hdr, chrName.c_str());
	    if (tid >= 0) {
	      if (tokIter!=tokens.end()) {
		int32_t start = 0;
		try {
		  start = boost::lexical_cast<int32_t>(*tokIter++);
		} catch (boost::bad_lexical_cast&) {
		  std::cerr << "Exclude file needs to be in tab-delimited format: chr, start, end" << std::endl;
		  std::cerr << "Offending line: " << chrFromFile << std::endl;
		  return false;
		}
		if (tokIter!=tokens.end()) {
		  int32_t end = start + 1;
		  try {
		    end = boost::lexical_cast<int32_t>(*tokIter++);
		  } catch (boost::bad_lexical_cast&) {
		    std::cerr << "Exclude file needs to be in tab-delimited format: chr, start, end" << std::endl;
		    std::cerr << "Offending line: " << chrFromFile << std::endl;
		    return false;
		  }
		  if (start < end) {
		    exclg[tid].insert(TIVal::right_open(start, end));
		  } else {
		    std::cerr << "Exclude file needs to be in tab-delimited format (chr, start, end) and start < end." << std::endl;
		    std::cerr << "Offending line: " << chrFromFile << std::endl;
		    return false;
		  }
		} else {
		  std::cerr << "Exclude file needs to be in tab-delimited format: chr, start, end" << std::endl;
		  std::cerr << "Offending line: " << chrFromFile << std::endl;
		  return false;
		}
	      } else validChr[tid] = false; // Exclude entire chromosome
	    }
	  }
	}
	chrFile.close();
      }
    }
    // Create the valid regions
    for (int32_t i = 0; i<hdr->n_targets; ++i) {
      if (!validChr[i]) continue;
      uint32_t istart = 0;
      for(typename TChrIntervals::iterator it = exclg[i].begin(); it != exclg[i].end(); ++it) {
	if (istart + 1 < it->lower()) validRegions[i].insert(TIVal::right_open(istart, it->lower() - 1));
	istart = it->upper();
      }
      if (istart + 1 < hdr->target_len[i]) validRegions[i].insert(TIVal::right_open(istart, hdr->target_len[i]));
    }
    exclg.clear();
    return true;
  }

  template<typename TIterator, typename TValue>
  inline void
  getMedian(TIterator begin, TIterator end, TValue& median) 
  {
    std::nth_element(begin, begin + (end - begin) / 2, end);
    median = *(begin + (end - begin) / 2);
  }
  
  template<typename TVector, typename TPercentile, typename TValue>
  inline void
  getPercentile(TVector& vec, TPercentile p, TValue& percentile) 
  {
    std::nth_element(vec.begin(), vec.begin() + int((vec.size() * p)), vec.end());
    percentile = *(vec.begin() + int(vec.size() * p));
  }


  template<typename TConfig>
  inline int32_t
  getVariability(TConfig const&, std::vector<LibraryInfo> const& lib) {
    int32_t overallVariability = 0;
    for(uint32_t libIdx = 0; libIdx < lib.size(); ++libIdx) {
      if (lib[libIdx].maxNormalISize > overallVariability) overallVariability = lib[libIdx].maxNormalISize;
      if (lib[libIdx].rs > overallVariability) overallVariability = lib[libIdx].rs;
    }
    return overallVariability;
  }

  
  template<typename TIterator, typename TValue>
  inline void
  getMAD(TIterator begin, TIterator end, TValue median, TValue& mad) 
  {
    std::vector<TValue> absDev;
    for(;begin<end;++begin) 
      absDev.push_back(std::abs((TValue)*begin - median));
    getMedian(absDev.begin(), absDev.end(), mad);
  }

  template<typename TIterator, typename TValue>
  inline void
  getMean(TIterator begin, TIterator end, TValue& mean) 
  {
    mean = 0;
    unsigned int count = 0;
    for(; begin<end; ++begin,++count) mean += *begin;
    mean /= count;
  }

  template<typename TIterator, typename TValue>
  inline void
  getStdDev(TIterator begin, TIterator end, TValue mean, TValue& stdDev) 
  {
    stdDev = 0;
    unsigned int count = 0;
    for(;begin<end;++begin,++count) stdDev += ((TValue)*begin - mean) * ((TValue)*begin - mean);
    stdDev = sqrt(stdDev / (TValue) count);
  }

  template<typename TConfig, typename TValidRegion, typename TSampleLibrary>
  inline void
  getLibraryParams(TConfig const& c, TValidRegion const& validRegions, TSampleLibrary& sampleLib) {
    typedef typename TValidRegion::value_type TChrIntervals;

    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    typedef std::vector<bam_hdr_t*> TSamHeader;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    TSamHeader hdr(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
      hdr[file_c] = sam_hdr_read(samfile[file_c]);
    }

    // Iterate all samples
    for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
      uint32_t maxAlignmentsScreened=10000000;
      uint32_t maxNumAlignments=1000000;
      uint32_t minNumAlignments=1000;
      uint32_t alignmentCount=0;
      uint32_t processedNumPairs = 0;
      uint32_t processedNumReads = 0;
      uint32_t rplus = 0;
      uint32_t nonrplus = 0;
      typedef std::vector<uint32_t> TSizeVector;
      TSizeVector vecISize;
      TSizeVector readSize;

      // Collect insert sizes
      bool libCharacterized = false;
      for(uint32_t refIndex=0; refIndex < (uint32_t) hdr[0]->n_targets; ++refIndex) {
	if (validRegions[refIndex].empty()) continue;
	for(typename TChrIntervals::const_iterator vRIt = validRegions[refIndex].begin(); ((vRIt != validRegions[refIndex].end()) && (!libCharacterized)); ++vRIt) {
	  hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, vRIt->lower(), vRIt->upper());
	  bam1_t* rec = bam_init1();
	  while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	    if (!(rec->core.flag & BAM_FREAD2) && (rec->core.l_qseq < 65000)) {
	      if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	      if ((alignmentCount > maxAlignmentsScreened) || ((processedNumReads >= maxNumAlignments) && (processedNumPairs == 0)) || (processedNumPairs >= maxNumAlignments)) {
		  // Paired-end library with enough pairs
		  libCharacterized = true;
		  break;
	      }
	      ++alignmentCount;
	      
	      // Single-end
	      if (processedNumReads < maxNumAlignments) {
		readSize.push_back(rec->core.l_qseq);
		++processedNumReads;
	      }
	      
	      // Paired-end
	      if ((rec->core.flag & BAM_FPAIRED) && !(rec->core.flag & BAM_FMUNMAP) && (rec->core.tid==rec->core.mtid)) {
		if (processedNumPairs < maxNumAlignments) {
		  vecISize.push_back(abs(rec->core.isize));
		  if (getSVType(rec) == 2) ++rplus;
		  else ++nonrplus;
		  ++processedNumPairs;
		}
	      }
	    }
	  }
	  bam_destroy1(rec);
	  hts_itr_destroy(iter);
	  if (libCharacterized) break;
	}
	if (libCharacterized) break;
      }
    
      // Get library parameters
      if (processedNumReads >= minNumAlignments) {
	std::sort(readSize.begin(), readSize.end());
	sampleLib[file_c].rs = readSize[readSize.size() / 2];
      }
      if (processedNumPairs >= minNumAlignments) {
	std::sort(vecISize.begin(), vecISize.end());
	int32_t median = vecISize[vecISize.size() / 2];
	std::vector<uint32_t> absDev;
	for(uint32_t i = 0; i < vecISize.size(); ++i) absDev.push_back(std::abs((int32_t) vecISize[i] - median));
	std::sort(absDev.begin(), absDev.end());
	int32_t mad = absDev[absDev.size() / 2];

	// Get default library orientation
	if ((median >= 50) && (median<=100000)) {
	  if (rplus < nonrplus) {
	    std::cerr << "Warning: Sample has a non-default paired-end layout! File: " << c.files[file_c].string() << std::endl;
	    std::cerr << "The expected paired-end orientation is   ---Read1--->      <---Read2---  which is the default illumina paired-end layout." << std::endl;
	    
	  } else {
	    sampleLib[file_c].median = median;
	    sampleLib[file_c].mad = mad;
	    sampleLib[file_c].maxNormalISize = median + (c.madNormalCutoff * mad);
	    sampleLib[file_c].minNormalISize = median - (c.madNormalCutoff * mad);
	    if (sampleLib[file_c].minNormalISize < 0) sampleLib[file_c].minNormalISize=0;
	    sampleLib[file_c].maxISizeCutoff = median + (c.madCutoff * mad);
	    sampleLib[file_c].minISizeCutoff = median - (c.madCutoff * mad);

	    // Deletion insert-size sanity checks
	    sampleLib[file_c].maxISizeCutoff = std::max(sampleLib[file_c].maxISizeCutoff, 2*sampleLib[file_c].rs);
	    sampleLib[file_c].maxISizeCutoff = std::max(sampleLib[file_c].maxISizeCutoff, 500);

	    if (sampleLib[file_c].minISizeCutoff < 0) sampleLib[file_c].minISizeCutoff=0;
	  }
	}
      }
    }

    // Clean-up
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      bam_hdr_destroy(hdr[file_c]);
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }


}

#endif
