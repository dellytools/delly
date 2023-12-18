#ifndef COMPVCF_H
#define COMPVCF_H

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
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>


#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <stdio.h>

#include "edlib.h"
#include "tags.h"
#include "version.h"
#include "util.h"
#include "modvcf.h"

namespace torali
{

  struct CompSVRecord {
    int32_t match;
    int32_t tid;
    int32_t mtid;
    int32_t svStart;
    int32_t svEnd;
    int32_t svLen;
    int32_t svt;
    int32_t qual;
    int32_t consBp;
    int32_t score;
    int32_t bestMatchId;
    double gtConc;
    double nonrefGtConc;
    std::string id;
    std::string allele;
    std::vector<int32_t> gt;

    CompSVRecord() : match(0), tid(0), svStart(0), svEnd(0), svLen(0), svt(0), qual(0), consBp(0), score(0), bestMatchId(0), gtConc(0), nonrefGtConc(0), id(""), allele("") {}
    CompSVRecord(int32_t const t, int32_t const svS) : match(0), tid(t), svStart(svS), svEnd(0), svLen(0), svt(0), qual(0), consBp(0), score(0), bestMatchId(0), gtConc(0), nonrefGtConc(0), id(""), allele("") {}
  };

  template<typename TSV>
  struct SortCompSVRecord : public std::binary_function<TSV, TSV, bool>
  {
    inline bool operator()(TSV const& sv1, TSV const& sv2) {
      return ((sv1.tid<sv2.tid) || ((sv1.tid==sv2.tid) && (sv1.svStart<sv2.svStart)) || ((sv1.tid==sv2.tid) && (sv1.svStart==sv2.svStart) && (sv1.svEnd<sv2.svEnd)));
    }
  };
  
  
  struct CompvcfConfig {
    typedef std::map<std::string, uint32_t> TChrMap;
    bool filterForPass;
    bool checkID;
    bool checkCT;
    int32_t qualthres;
    int32_t bpdiff;
    int32_t minsize;
    int32_t maxsize;
    int32_t minac;
    int32_t maxac;
    float sizeratio;
    float divergence;
    boost::filesystem::path vcffile;
    boost::filesystem::path base;
    std::string outprefix;
    std::vector<std::string> samples;
    TChrMap chrmap;
  };

  inline void
  compareSVs(CompvcfConfig const& c, std::vector<CompSVRecord>& basesv, std::vector<CompSVRecord>& compsv) {
    typedef std::vector<CompSVRecord> TCompSVType;
    
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Comparing " << compsv.size() << " SVs with " << basesv.size() << " SVs in the base VCF/BCF file " << std::endl;
    for(uint32_t i = 0; i < basesv.size(); ++i) {
      int32_t earliestStart = std::max(basesv[i].svStart - (c.bpdiff + 1), 0);
      typename TCompSVType::const_iterator itsv = std::lower_bound(compsv.begin(), compsv.end(), CompSVRecord(basesv[i].tid, earliestStart), SortCompSVRecord<CompSVRecord>());
      for(uint32_t j = (itsv - compsv.begin()); j < compsv.size(); ++j) {
	if (basesv[i].tid < compsv[j].tid) break;  // Sorted by tid
	if (c.checkCT) {
	  if (basesv[i].svt != compsv[j].svt) continue;
	} else {
	  if (_addID(basesv[i].svt) != _addID(compsv[j].svt)) continue;  // Compare SV types (BND, INV,...) but not CT
	}
	if (basesv[i].mtid != compsv[j].mtid) continue;
	if (std::abs(basesv[i].svStart - compsv[j].svStart) > c.bpdiff) {
	  if ((compsv[j].svStart > basesv[i].svStart) && ((compsv[j].svStart - basesv[i].svStart) > c.bpdiff)) break; // Sorted by tid & svStart
	  continue;
	}
	if (std::abs(basesv[i].svEnd - compsv[j].svEnd) > c.bpdiff) continue;
	float bprat = 1 - (float) std::max(std::abs(basesv[i].svStart - compsv[j].svStart), std::abs(basesv[i].svEnd - compsv[j].svEnd)) / (float) (c.bpdiff);
	float sizerat = 1;
	if ((basesv[i].svLen) && (compsv[j].svLen)) {
	  sizerat = (float) basesv[i].svLen / (float) compsv[j].svLen;
	  if (basesv[i].svLen > compsv[j].svLen) sizerat = (float) compsv[j].svLen / (float) basesv[i].svLen;
	}
	if (sizerat < c.sizeratio) continue;
	// Check SV similarity
	float scorerat = 0;
	if ((!basesv[i].allele.empty()) && (!compsv[j].allele.empty())) {
	  int32_t leftOffset = std::min(basesv[i].consBp, compsv[j].consBp);
	  int32_t rightOffset = std::min(basesv[i].allele.size() - basesv[i].consBp, compsv[j].allele.size() - compsv[j].consBp);
	  std::string seqI = basesv[i].allele.substr(basesv[i].consBp - leftOffset, leftOffset+rightOffset);
	  std::string seqJ = compsv[j].allele.substr(compsv[j].consBp - leftOffset, leftOffset+rightOffset);
	  EdlibAlignResult cigar = edlibAlign(seqI.c_str(), seqI.size(), seqJ.c_str(), seqJ.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
	  //printAlignment(seqI, seqJ, EDLIB_MODE_NW, cigar);
	  double score = (double) cigar.editDistance / (double) (leftOffset + rightOffset);
	  edlibFreeAlignResult(cigar);
	  if (score > c.divergence) continue;
	  scorerat = 1 - score / c.divergence;
	}
	// Match
	++basesv[i].match;
	++compsv[j].match;
	double gtconc = gtConc(basesv[i].gt, compsv[j].gt);
	if (gtconc > basesv[i].gtConc) basesv[i].gtConc = gtconc;
	if (gtconc > compsv[j].gtConc) compsv[j].gtConc = gtconc;
	double nonrefgtconc = nonrefGtConc(basesv[i].gt, compsv[j].gt);
	if (nonrefgtconc > basesv[i].nonrefGtConc) basesv[i].nonrefGtConc = nonrefgtconc;
	if (nonrefgtconc > compsv[j].nonrefGtConc) compsv[j].nonrefGtConc = nonrefgtconc;

	// Update best match
	int32_t matchScore = (int32_t) ((bprat + sizerat + scorerat + gtconc) * 100);
	if (compsv[j].score < matchScore) {
	  compsv[j].score = matchScore;
	  compsv[j].bestMatchId = i;
	}
	if (basesv[i].score < matchScore) {
	  basesv[i].score = (int32_t) ((bprat + sizerat + scorerat + gtconc) * 100);
	  basesv[i].bestMatchId = j;
	}
	//std::cerr << basesv[i].tid << ',' << basesv[i].svStart << ',' << basesv[i].svEnd << ',' << basesv[i].id << ',' << basesv[i].svLen << std::endl;
	//std::cerr << compsv[j].tid << ',' << compsv[j].svStart << ',' << compsv[j].svEnd << ',' << compsv[j].id << ',' << compsv[j].svLen << std::endl;
      }
    }
  }
  
  inline bool
  _loadCompSVs(CompvcfConfig& c, std::string const& filename, std::vector<CompSVRecord>& allsv) {
    bool success = true;
    std::set<std::string> allIds;

    // Sample map
    std::map<std::string, uint32_t> smap;
    for(uint32_t i = 0; i < c.samples.size(); ++i) smap.insert(std::make_pair(c.samples[i], i));
    
    // Load bcf file
    htsFile* ifile = hts_open(filename.c_str(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);

    // VCF fields
    int32_t nsvend = 0;
    int32_t* svend = NULL;
    int32_t ninslen = 0;
    int32_t* inslen = NULL;
    int32_t nconsbp = 0;
    int32_t* consbp = NULL;
    int32_t svEndVal = -1;
    int32_t nsvt = 0;
    char* svt = NULL;
    std::string svtVal;
    int32_t ncons = 0;
    char* cons = NULL;
    int32_t nchr2 = 0;
    char* chr2 = NULL;
    int32_t npos2 = 0;
    int32_t* pos2 = NULL;
    int32_t nct = 0;
    char* ct = NULL;
    std::string ctVal;
    int ngt = 0;
    int32_t* gt = NULL;
    int32_t inslenVal = -1;
    int32_t svLenVal = 0;
    std::string chr2Name;
    
    // Parse BCF
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Parsing VCF/BCF file " << filename << std::endl;
    uint32_t svcounter = 0;
    bcf1_t* rec = bcf_init1();
    while (bcf_read(ifile, hdr, rec) == 0) {
      bcf_unpack(rec, BCF_UN_INFO);

      // Check SV type
      svtVal = "NA";
      ctVal = "NA";
      svEndVal = -1;
      inslenVal = -1;
      if (_isKeyPresent(hdr, "SVTYPE")) {
	bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt);
	svtVal = std::string(svt);
	if (_isKeyPresent(hdr, "CT")) {
	  bcf_get_info_string(hdr, rec, "CT", &ct, &nct);
	  ctVal = std::string(ct);
	} else {
	  if (svtVal == "INS") ctVal = "NtoN";
	  else if (svtVal == "DEL") ctVal = "3to5";
	  else if (svtVal == "DUP") ctVal = "5to3";
	  else if (svtVal == "INV") ctVal = "3to3"; // or 5to5
	}
      }

      // Symbolic ALT?
      std::string refAllele = rec->d.allele[0];
      std::string altAllele = rec->d.allele[1];
      std::string altSymbol = "NA";
      if ((!altAllele.empty()) && (altAllele[0] == '<') && (altAllele[altAllele.size() - 1] == '>')) {
	// Symbolic ALT
	altSymbol = altAllele.substr(1, altAllele.size() - 2);
	if (svtVal == "NA") svtVal = altSymbol;
	else {
	  if ((altSymbol == "CN0") || (altSymbol == "CN1")) altSymbol = "DEL";
	  if (altSymbol.substr(0,3) == "INS") altSymbol = "INS";
	  if (altSymbol.substr(0,3) == "DEL") altSymbol = "DEL";
	  if (svtVal != altSymbol) {
	    success=false;
	    std::cerr << "Error: SV type " << svtVal << " disagrees with symbolic ALT." << std::endl;
	  }
	}
      }

      // Still unknown SV type?
      if (svtVal == "NA") {
	if (refAllele.size() > altAllele.size()) {
	  svtVal = "DEL";
	  ctVal = "3to5";
	} else {
	  svtVal = "INS";
	  ctVal = "NtoN";
	}
      }

      // SV end
      if (_isKeyPresent(hdr, "END")) {
	if (bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend) > 0) svEndVal = *svend;
      }

      // Insertion length
      if (_isKeyPresent(hdr, "INSLEN")) {
	if (bcf_get_info_int32(hdr, rec, "INSLEN", &inslen, &ninslen) > 0) inslenVal = *inslen;
      }

      // BNDs
      if (svtVal == "BND") {
	int32_t pos2val = -1;
	if ((_isKeyPresent(hdr, "CHR2")) && (_isKeyPresent(hdr, "POS2"))) {
	  if (bcf_get_info_int32(hdr, rec, "POS2", &pos2, &npos2) > 0) pos2val = *pos2;
	  if (bcf_get_info_string(hdr, rec, "CHR2", &chr2, &nchr2) > 0) chr2Name = std::string(chr2);
	} else {
	  // Parse ALT
	  std::vector<std::size_t> posBND;
	  std::size_t found = altAllele.find('[');
	  while (found!=std::string::npos) {
	    posBND.push_back(found);
	    found = altAllele.find('[', found+1);
	  }
	  found = altAllele.find(']');
	  while (found!=std::string::npos) {
	    posBND.push_back(found);
	    found = altAllele.find(']', found+1);
	  }
	  std::sort(posBND.begin(), posBND.end());
	  altAllele = altAllele.substr(posBND[0] + 1, posBND[1] - posBND[0] - 1);
	  found = altAllele.find(':');
	  if (found!=std::string::npos) {
	    chr2Name = altAllele.substr(0, found);
	    pos2val = boost::lexical_cast<int32_t>(altAllele.substr(found + 1));
	  }
	}
	if ((!chr2Name.empty()) && (pos2val != -1)) {
	  svEndVal = pos2val;
	  svLenVal = 0;
	} else {
	  success=false;
	  std::cerr << "Error: Could not parse BND ALT allele " << rec->d.allele[1] << std::endl;
	}
      }
      
      // Assign SV end position
      if (svEndVal == -1) {
	if (altSymbol == "NA") {
	  if (svtVal == "DEL") svEndVal = rec->pos + (refAllele.size() - altAllele.size());
	  else if (svtVal == "INS") svEndVal = rec->pos + 1;
	  else {
	    success=false;
	    std::cerr << "Error: SV end position is unknown for " << svtVal << "!" << std::endl;
	  }
	} else {
	  success=false;
	  std::cerr << "Error: Missing SV end position!" << std::endl;
	}
      }

      // Assign insertion length
      if (inslenVal == -1) {
	if (svtVal == "INS") inslenVal = altAllele.size() - refAllele.size();
	else inslenVal = 0;
      }

      // Assign pass state
      bool pass = true;
      if (c.filterForPass) pass = (bcf_has_filter(hdr, rec, const_cast<char*>("PASS"))==1);

      // Assign site quality
      int32_t qualVal = 0;
      if (!boost::math::isnan(rec->qual)) qualVal = rec->qual;

      // Assign SV length
      if (svtVal != "INS") {
	if (svtVal != "BND") svLenVal = svEndVal - rec->pos;
      } else {
	svEndVal = rec->pos + 1;
	svLenVal = inslenVal;
      }

      // Filter sites
      if ((qualVal >= c.qualthres) && (pass) && (svLenVal >= c.minsize) && (svLenVal < c.maxsize)) {
	// Define SV event
	CompSVRecord sv;
	sv.match = 0;
	std::string chrname = std::string(bcf_hdr_id2name(hdr, rec->rid));
	if (c.chrmap.find(chrname) == c.chrmap.end()) c.chrmap.insert(std::make_pair(chrname, c.chrmap.size()));
	sv.tid = c.chrmap[chrname];
	sv.mtid = sv.tid;
	if (svtVal == "BND") {
	  if (c.chrmap.find(chr2Name) == c.chrmap.end()) c.chrmap.insert(std::make_pair(chr2Name, c.chrmap.size()));
	  sv.mtid = c.chrmap[chr2Name];
	}
	sv.svStart = rec->pos;
	sv.svEnd = svEndVal;
	// Use some canonical ordering, swap if necessary
	if ((svtVal == "BND") && (sv.tid > sv.mtid)) {
	  int32_t tmp = sv.tid;
	  sv.tid = sv.mtid;
	  sv.mtid = tmp;
	  tmp = sv.svStart;
	  sv.svStart = sv.svEnd;
	  sv.svEnd = tmp;
	}
	sv.svLen = svLenVal;
	sv.svt = _decodeOrientation(ctVal, svtVal);
	sv.qual = qualVal;

	// Check genotypes
	bcf_unpack(rec, BCF_UN_ALL);
	bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
	sv.gt.resize(c.samples.size(), -1); // Missing GT initialization
	int32_t gtsum = 0;
	for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
	  if ((bcf_gt_allele(gt[i*2]) != -1) && (bcf_gt_allele(gt[i*2 + 1]) != -1)) {
	    std::string sname = hdr->samples[i];
	    if (smap.find(sname) != smap.end()) {
	      sv.gt[smap[sname]] = bcf_gt_allele(gt[i*2]) + bcf_gt_allele(gt[i*2 + 1]);
	      gtsum += bcf_gt_allele(gt[i*2]) + bcf_gt_allele(gt[i*2 + 1]);
	    }
	  }
	}
	if (gtsum < 0) gtsum = 0;
	
	// Min. and max. allele count
	if ((gtsum >= c.minac) && (gtsum < c.maxac)) {
	  sv.allele = "";
	  if (_isKeyPresent(hdr, "CONSENSUS")) {
	    if (bcf_get_info_string(hdr, rec, "CONSENSUS", &cons, &ncons) > 0) {
	      if (bcf_get_info_int32(hdr, rec, "CONSBP", &consbp, &nconsbp) > 0) {
		sv.allele = std::string(cons);
		sv.consBp = *consbp;
	      }
	    }
	  }
	  sv.id = std::string(rec->d.id);
	  if (sv.id == ".") {
	    sv.id = std::string(bcf_hdr_id2name(hdr, rec->rid)) + "-" + boost::lexical_cast<std::string>(rec->pos) + "-ID" + boost::lexical_cast<std::string>(++svcounter);
	  }
	  if ((c.checkID) && (allIds.find(sv.id) != allIds.end())) {
	    success=false;
	    std::cerr << "Error: Duplicate IDs " << sv.id << std::endl;
	  } else {
	    //std::cerr << sv.tid << ',' << sv.svStart << ',' << sv.svEnd << ',' << sv.id << ',' << sv.svLen << ',' << sv.svt << std::endl;
	    allIds.insert(sv.id);
	    allsv.push_back(sv);
	  }
	}
      }
    }
    bcf_destroy(rec);

    // Clean-up
    if (svend != NULL) free(svend);
    if (inslen != NULL) free(inslen);
    if (consbp != NULL) free(consbp);
    if (svt != NULL) free(svt);
    if (gt != NULL) free(gt);
    if (ct != NULL) free(ct);
    if (cons != NULL) free(cons);
    if (chr2 != NULL) free(chr2);
    if (pos2 != NULL) free(pos2);

    // Close VCF
    bcf_hdr_destroy(hdr);
    bcf_close(ifile);

    return success;
  }

  inline int
  compvcfRun(CompvcfConfig& c) {

    // Load SVs
    std::vector<CompSVRecord> basesv;
    if (!_loadCompSVs(c, c.base.string(), basesv)) return -1;
    
    std::vector<CompSVRecord> compsv;
    if (!_loadCompSVs(c, c.vcffile.string(), compsv)) return -1;

    // Sort SVs
    sort(basesv.begin(), basesv.end(), SortCompSVRecord<CompSVRecord>());
    sort(compsv.begin(), compsv.end(), SortCompSVRecord<CompSVRecord>());

    // Recall, precission, GT concordance
    compareSVs(c, basesv, compsv);

    // Metrics
    uint32_t tp_base = 0;
    uint32_t tp_comp = 0;
    uint32_t redundant = 0;
    uint32_t fn = 0;
    uint32_t fp = 0;
    double gtconc = 0;
    double nonrefgtconc = 0;
    for(uint32_t i = 0; i < basesv.size(); ++i) {
      if (basesv[i].match) {
	++tp_base;
	redundant += basesv[i].match;
	gtconc += basesv[i].gtConc;
	nonrefgtconc += basesv[i].nonrefGtConc;
      } else ++fn;
    }
    for(uint32_t j = 0; j < compsv.size(); ++j) {
      if (compsv[j].match) ++tp_comp;
      else ++fp;
    }
    double recall = (double) (tp_base) / (double) (basesv.size());
    double precision = (double) (tp_base) / (double) (tp_base + fp);
    double redundancyRation = (double) redundant / (double) (tp_base);
    double f1 = 2 * recall * precision / (recall + precision);
    gtconc /= (double) (tp_base);
    nonrefgtconc /= (double) (tp_base);

    std::string filename = c.outprefix + ".tsv";
    std::ofstream ofile(filename.c_str());
    ofile << "Size\tAC\tTP_Base\tFN\tTP_Comp\tFP\tRecall\tPrecision\tF1\tRedundancyRatio\tGTConc\tNonRefGTConc" << std::endl;
    ofile << '[' << boost::lexical_cast<std::string>(c.minsize) << ',' << boost::lexical_cast<std::string>(c.maxsize) << '[' << '\t';
    ofile << '[' << boost::lexical_cast<std::string>(c.minac) << ',' << boost::lexical_cast<std::string>(c.maxac) << '[' << '\t';
    ofile << tp_base << '\t' << fn << '\t' << tp_comp << '\t' << fp << '\t' << recall << '\t' << precision << '\t' << f1 << '\t' << redundancyRation << '\t' << gtconc << '\t' << nonrefgtconc << std::endl;
    ofile.close();

    // Output classification
    std::map<uint32_t, std::string> idxchr;
    for(typename CompvcfConfig::TChrMap::const_iterator it = c.chrmap.begin(); it != c.chrmap.end(); ++it) idxchr.insert(std::make_pair(it->second, it->first));
    filename = c.outprefix + ".sv.classification";
    std::ofstream svfile(filename.c_str());
    svfile << "ID\tClassification\tScore\tBestMatchId\tChrom\tStart\tEnd\tLength\tSVType\tCT\tQuality\tGTConc\tNonRefGTConc\tMatchCountBase\tAlignmentAllele" << std::endl;
    for(uint32_t j = 0; j < compsv.size(); ++j) {
      std::string label;
      std::string baseSVLabel = "NA";
      if (compsv[j].match) {
	label="TP";
	baseSVLabel = basesv[compsv[j].bestMatchId].id;
      } else label="FP";
      svfile << compsv[j].id << '\t' << label << '\t' << compsv[j].score << "\tbase" << baseSVLabel << '\t' << idxchr[compsv[j].tid] << '\t' << compsv[j].svStart << '\t' << compsv[j].svEnd << '\t' << compsv[j].svLen << '\t' << _addID(compsv[j].svt) << '\t' << _addOrientation(compsv[j].svt) << '\t' << compsv[j].qual << '\t' << compsv[j].gtConc << '\t' << compsv[j].nonrefGtConc << '\t' << compsv[j].match << '\t' << compsv[j].allele << std::endl;
    }
    svfile.close();
    
    // Done
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    return 0;
  }


  int compvcf(int argc, char **argv) {
    CompvcfConfig c;
    
    // Define generic options
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("base,a", boost::program_options::value<boost::filesystem::path>(&c.base), "base VCF/BCF file")
      ("quality,y", boost::program_options::value<int32_t>(&c.qualthres)->default_value(0), "min. SV site quality")
      ("minsize,m", boost::program_options::value<int32_t>(&c.minsize)->default_value(50), "min. SV size")
      ("maxsize,n", boost::program_options::value<int32_t>(&c.maxsize)->default_value(100000), "max. SV size")
      ("minac,e", boost::program_options::value<int32_t>(&c.minac)->default_value(1), "min. allele count")
      ("maxac,f", boost::program_options::value<int32_t>(&c.maxac)->default_value(10000), "max. allele count")
      ("bpdiff,b", boost::program_options::value<int32_t>(&c.bpdiff)->default_value(1000), "max. SV breakpoint offset")
      ("sizeratio,s", boost::program_options::value<float>(&c.sizeratio)->default_value(0.5), "min. SV size ratio")
      ("divergence,d", boost::program_options::value<float>(&c.divergence)->default_value(0.3), "max. SV allele divergence")
      ("outprefix,o", boost::program_options::value<std::string>(&c.outprefix)->default_value("out"), "output prefix")
      ("pass,p", "Filter sites for PASS")
      ("ignore,i", "Ignore duplicate IDs")
      ("ct,c", "Require matching CT value in addition to SV type")
      ;
    
    // Define hidden options
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "comparison VCF/BCF file")
      ;
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    // Set the visibility
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("base"))) {
      std::cerr << std::endl;
      std::cerr << "Usage: delly " << argv[0] << " [OPTIONS] -a <base.bcf> <input.bcf>" << std::endl;
      std::cerr << visible_options << "\n";
      return 0;
    }
    
    // Filter for PASS
    if (vm.count("pass")) c.filterForPass = true;
    else c.filterForPass = false;

    // Check duplicate IDs
    if (vm.count("ignore")) c.checkID = false;
    else c.checkID = true;

    // Check CT
    if (vm.count("ct")) c.checkCT = true;
    else c.checkCT = false;

    // Check base VCF file
    std::set<std::string> baseSamples;
    if (vm.count("base")) {
      if (!(boost::filesystem::exists(c.base) && boost::filesystem::is_regular_file(c.base) && boost::filesystem::file_size(c.base))) {
	std::cerr << "Input VCF/BCF file is missing: " << c.base.string() << std::endl;
	return 1;
      }
      htsFile* ifile = bcf_open(c.base.string().c_str(), "r");
      if (ifile == NULL) {
	std::cerr << "Fail to open file " << c.base.string() << std::endl;
      return 1;
      }
      hts_idx_t* bcfidx = NULL;
      tbx_t* tbx = NULL;
      if (hts_get_format(ifile)->format==vcf) tbx = tbx_index_load(c.base.string().c_str());
      else bcfidx = bcf_index_load(c.base.string().c_str());
      if ((bcfidx == NULL) && (tbx == NULL)) {
	std::cerr << "Fail to open index file for " << c.base.string() << std::endl;
	return 1;
      }
      bcf_hdr_t* hdr = bcf_hdr_read(ifile);
      if (hdr == NULL) {
	std::cerr << "Fail to header for " << c.base.string() << std::endl;
	return 1;
      }
      if (!(bcf_hdr_nsamples(hdr)>0)) {
	std::cerr << "Warning: BCF/VCF file has no sample genotypes " << c.base.string() << std::endl;
      }
      for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) baseSamples.insert(hdr->samples[i]);
      bcf_hdr_destroy(hdr);
      if (bcfidx) hts_idx_destroy(bcfidx);
      if (tbx) tbx_destroy(tbx);
      bcf_close(ifile);
    }
    
    // Check comparison VCF file
    std::set<std::string> compSamples;
    if (vm.count("input-file")) {
      if (!(boost::filesystem::exists(c.vcffile) && boost::filesystem::is_regular_file(c.vcffile) && boost::filesystem::file_size(c.vcffile))) {
	std::cerr << "Input VCF/BCF file is missing: " << c.vcffile.string() << std::endl;
	return 1;
      }
      htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
      if (ifile == NULL) {
	std::cerr << "Fail to open file " << c.vcffile.string() << std::endl;
      return 1;
      }
      hts_idx_t* bcfidx = NULL;
      tbx_t* tbx = NULL;
      if (hts_get_format(ifile)->format==vcf) tbx = tbx_index_load(c.vcffile.string().c_str());
      else bcfidx = bcf_index_load(c.vcffile.string().c_str());
      if ((bcfidx == NULL) && (tbx == NULL)) {
	std::cerr << "Fail to open index file for " << c.vcffile.string() << std::endl;
	return 1;
      }
      bcf_hdr_t* hdr = bcf_hdr_read(ifile);
      if (hdr == NULL) {
	std::cerr << "Fail to header for " << c.vcffile.string() << std::endl;
	return 1;
      }
      if (!(bcf_hdr_nsamples(hdr)>0)) {
	std::cerr << "Warning: BCF/VCF file has no sample genotypes " << c.vcffile.string() << std::endl;
      }
      for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) compSamples.insert(hdr->samples[i]);
      bcf_hdr_destroy(hdr);
      if (bcfidx) hts_idx_destroy(bcfidx);
      if (tbx) tbx_destroy(tbx);
      bcf_close(ifile);
    }
    
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cerr << "delly ";
    for(int i=0; i<argc; ++i) { std::cerr << argv[i] << ' '; }
    std::cerr << std::endl;

    // Common samples
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << baseSamples.size() << " base samples and " << compSamples.size() << " comparison samples" << std::endl;
    std::set_intersection(baseSamples.begin(), baseSamples.end(), compSamples.begin(), compSamples.end(), std::inserter(c.samples, c.samples.begin()));
    //for(uint32_t i = 0; i < c.samples.size(); ++i) std::cerr << c.samples[i] << std::endl;
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << c.samples.size() << " common samples" << std::endl;
    if (c.samples.size() < 1) std::cerr << "Warning: No common samples detected!" << std::endl;

    // Run comparison
    return compvcfRun(c);
  }
  
}

#endif
