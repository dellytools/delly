#ifndef MODVCF_H
#define MODVCF_H

#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "bolog.h"



namespace torali
{


void _remove_info_tag(bcf_hdr_t* hdr, bcf1_t* rec, std::string const& tag) {
  bcf_update_info(hdr, rec, tag.c_str(), NULL, 0, BCF_HT_INT);  // Type does not matter for n = 0
}

void _remove_format_tag(bcf_hdr_t* hdr, bcf1_t* rec, std::string const& tag) {
  bcf_update_format(hdr, rec, tag.c_str(), NULL, 0, BCF_HT_INT);  // Type does not matter for n = 0
}

void _remove_info(bcf_hdr_t* hdr, bcf1_t* rec) {
  std::string tmp[] = {"CT", "PRECISE", "IMPRECISE", "SVTYPE", "SVMETHOD", "CIEND", "CIPOS", "CHR2", "POS2", "END", "PE", "MAPQ", "SRMAPQ", "SR", "SRQ", "CONSENSUS"};
  std::set<std::string> keepInfo(tmp, tmp + sizeof(tmp)/sizeof(tmp[0]));

  if (!(rec->unpacked & BCF_UN_INFO)) bcf_unpack(rec, BCF_UN_INFO);
  
  for (uint32_t i = 0; i < rec->n_info; ++i){
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
  
  for(uint32_t i = 0; i<rec->n_fmt; ++i) {
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

inline int
_getInfoType(bcf_hdr_t const* hdr, std::string const& key) {
  return bcf_hdr_id2type(hdr, BCF_HL_INFO, bcf_hdr_id2int(hdr, BCF_DT_ID, key.c_str()));
}

inline int
_getFormatType(bcf_hdr_t const* hdr, std::string const& key) {
  return bcf_hdr_id2type(hdr, BCF_HL_FMT, bcf_hdr_id2int(hdr, BCF_DT_ID, key.c_str()));
}

inline bool _missing(bool const value) {
  return !value;
}

inline bool _missing(float const value) {
  return bcf_float_is_missing(value);
}

inline bool _missing(int8_t const value) {
  return (value == bcf_int8_missing);
}

inline bool _missing(int16_t const value) {
  return (value == bcf_int16_missing);
}

inline bool _missing(int32_t const value) {
  return (value == bcf_int32_missing);
}

inline bool _missing(std::string const& value) {
  return ((value.empty()) || (value == "."));
}

inline bool
_isKeyPresent(bcf_hdr_t const* hdr, std::string const& key) {
  return (bcf_hdr_id2int(hdr, BCF_DT_ID, key.c_str())>=0);
}

inline bool
_isDNA(std::string const& allele) {
  for(uint32_t i = 0; i<allele.size(); ++i) {
    if ((allele[i] != 'A') && (allele[i] != 'C') && (allele[i] != 'G') && (allele[i] != 'T') && (allele[i] != 'a') && (allele[i] != 'c') && (allele[i] != 'g') && (allele[i] != 't')) return false;
  }
  return true;
}
 
inline std::string
_replaceIUPAC(std::string const& alleles) {
  std::vector<char> out(alleles.size());
  int32_t inTag = 0;
  for(uint32_t i = 0; i<alleles.size(); ++i) {
    if ((inTag) || (alleles[i] == 'A') || (alleles[i] == 'C') || (alleles[i] == 'G') || (alleles[i] == 'T') || (alleles[i] == 'N') || (alleles[i] == 'a') || (alleles[i] == 'c') || (alleles[i] == 'g') || (alleles[i] == 't') || (alleles[i] == 'n') || (alleles[i] == '<') || (alleles[i] == '>') || (alleles[i] == ']') || (alleles[i] == '[') || (alleles[i] == ',')) {
      out[i] = alleles[i];
      if (alleles[i] == '<') inTag = 1;
      else if (alleles[i] == ']') inTag = 2;
      else if (alleles[i] == '[') inTag = 3;
      else if ((alleles[i] == '>') && (inTag == 1)) inTag = 0;
      else if ((alleles[i] == ']') && (inTag == 2)) inTag = 0;
      else if ((alleles[i] == '[') && (inTag == 3)) inTag = 0;
    } else {
      // Replace IUPAC
      if ((alleles[i] == 'U') || (alleles[i] == 'u')) out[i] = 'T';
      else if ((alleles[i] == 'R') || (alleles[i] == 'r')) out[i] = 'A';
      else if ((alleles[i] == 'Y') || (alleles[i] == 'y')) out[i] = 'C';
      else if ((alleles[i] == 'S') || (alleles[i] == 's')) out[i] = 'C';
      else if ((alleles[i] == 'W') || (alleles[i] == 'w')) out[i] = 'A';
      else if ((alleles[i] == 'K') || (alleles[i] == 'k')) out[i] = 'G';
      else if ((alleles[i] == 'M') || (alleles[i] == 'm')) out[i] = 'A';
      else if ((alleles[i] == 'B') || (alleles[i] == 'b')) out[i] = 'C';
      else if ((alleles[i] == 'D') || (alleles[i] == 'd')) out[i] = 'A';
      else if ((alleles[i] == 'H') || (alleles[i] == 'h')) out[i] = 'A';
      else if ((alleles[i] == 'V') || (alleles[i] == 'v')) out[i] = 'A';
      else out[i] = 'N';
    }
  }
  return std::string(out.begin(), out.end());
}     
  
 
 // Convert string to char*
struct cstyle_str {
  const char* operator ()(const std::string& s) {
    return s.c_str();
  }
};
 
// Parse Delly vcf file
template<typename TConfig, typename TStructuralVariantRecord>
inline void
vcfParse(TConfig const& c, bam_hdr_t* hd, std::vector<TStructuralVariantRecord>& svs) {
  // Load bcf file
  htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
  bcf_hdr_t* hdr = bcf_hdr_read(ifile);
  bcf1_t* rec = bcf_init();

  // Parse bcf
  int32_t nsvend = 0;
  int32_t* svend = NULL;
  int32_t npos2 = 0;
  int32_t* pos2 = NULL;
  int32_t nconsbp = 0;
  int32_t* consbp = NULL;
  int32_t nsvlen = 0;
  int32_t* svlen = NULL;
  int32_t npe = 0;
  int32_t* pe = NULL;
  int32_t ninslen = 0;
  int32_t* inslen = NULL;
  int32_t nhomlen = 0;
  int32_t* homlen = NULL;
  int32_t nsr = 0;
  int32_t* sr = NULL;
  int32_t ncipos = 0;
  int32_t* cipos = NULL;
  int32_t nmapq = 0;
  int32_t* mapq = NULL;
  int32_t nct = 0;
  char* ct = NULL;
  int32_t nsrq = 0;
  float* srq = NULL;
  int32_t nsvt = 0;
  char* svt = NULL;
  int32_t nmethod = 0;
  char* method = NULL;
  int32_t ncons = 0;
  char* cons = NULL;
  int32_t nchr2 = 0;
  char* chr2 = NULL;
  bool dellyVCF = false;
  while (bcf_read(ifile, hdr, rec) == 0) {
    bcf_unpack(rec, BCF_UN_INFO);

    // Delly BCF file?
    if (!dellyVCF) {
      if (bcf_get_info_string(hdr, rec, "SVMETHOD", &method, &nmethod) > 0) {
	std::string mstr = std::string(method);
	if ((mstr.size() >= 10) && (mstr.substr(0, 10)  == "EMBL.DELLY")) {
	  if (_isKeyPresent(hdr, "CONSBP")) dellyVCF = true;
	}
      }
    }

    // Delly
    if (dellyVCF) {
      // Fill SV record
      StructuralVariantRecord svRec;
      std::string chrName = bcf_hdr_id2name(hdr, rec->rid);
      int32_t tid = bam_name2id(hd, chrName.c_str());
      svRec.chr = tid;
      svRec.svStart = rec->pos + 1;
      svRec.id = svs.size();
      svRec.mapq = rec->qual;
      std::string refAllele = rec->d.allele[0];
      if ((!refAllele.empty()) && (refAllele != ".")) {
	std::string altAllele = rec->d.allele[1];
	if ((!altAllele.empty()) && (altAllele != ".")) {
	  svRec.alleles = refAllele + "," + altAllele;
	}
      }

      // Parse SV type
      if ((bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt) > 0) && (bcf_get_info_string(hdr, rec, "CT", &ct, &nct) > 0)) svRec.svt = _decodeOrientation(std::string(ct), std::string(svt));
      else continue;
      
      // Parse INFO
      if (bcf_get_info_flag(hdr, rec, "PRECISE", 0, 0) > 0) svRec.precise=true;
      else svRec.precise = false;
      if (bcf_get_info_int32(hdr, rec, "PE", &pe, &npe) > 0) svRec.peSupport = *pe;
      else {
	if (svRec.precise) svRec.peSupport = 0;
	else svRec.peSupport = 2;
      }
      if (svRec.svt != 4) {
	if (bcf_get_info_int32(hdr, rec, "INSLEN", &inslen, &ninslen) > 0) svRec.insLen = *inslen;
	else svRec.insLen = 0;
      } else {
	// Insertions must have INFO/SVLEN
	if (bcf_get_info_int32(hdr, rec, "SVLEN", &svlen, &nsvlen) > 0) svRec.insLen = *svlen;
	else continue;
      }
      if (bcf_get_info_int32(hdr, rec, "HOMLEN", &homlen, &nhomlen) > 0) svRec.homLen = *homlen;
      else svRec.homLen = 0;
      if (bcf_get_info_int32(hdr, rec, "SR", &sr, &nsr) > 0) svRec.srSupport = *sr;
      else svRec.srSupport = 0;

      // SV end assignment
      svRec.chr2 = tid;
      svRec.svEnd = rec->pos + 2;
      if (svRec.svt < DELLY_SVT_TRANS) {
	// Intra-chromosomal SV
	if (bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend) > 0) svRec.svEnd = *svend;
	svRec.svEnd += 1;
      } else {
	// Inter-chromosomal SV
	if (bcf_get_info_string(hdr, rec, "CHR2", &chr2, &nchr2) > 0) {
	  std::string chr2Name = std::string(chr2);
	  svRec.chr2 = bam_name2id(hd, chr2Name.c_str());
	}
	if (bcf_get_info_int32(hdr, rec, "POS2", &pos2, &npos2) > 0) svRec.svEnd = *pos2;
      }
      if (bcf_get_info_string(hdr, rec, "CONSENSUS", &cons, &ncons) > 0) {
	svRec.consensus = std::string(cons);
	if (bcf_get_info_int32(hdr, rec, "CONSBP", &consbp, &nconsbp) > 0) svRec.consBp = *consbp;
      }
      else svRec.precise = false;
      if (bcf_get_info_int32(hdr, rec, "CIPOS", &cipos, &ncipos) > 0) {
	svRec.ciposlow = cipos[0];
	svRec.ciposhigh = cipos[1];
      } else {
	svRec.ciposlow = -50;
	svRec.ciposhigh = 50;
      }
      if (bcf_get_info_int32(hdr, rec, "CIEND", &cipos, &ncipos) > 0) {
	svRec.ciendlow = cipos[0];
	svRec.ciendhigh = cipos[1];
      } else {
	svRec.ciendlow = -50;
	svRec.ciendhigh = 50;
      }
      if (bcf_get_info_int32(hdr, rec, "MAPQ", &mapq, &nmapq) > 0) svRec.peMapQuality = (uint8_t) *mapq;
      else svRec.peMapQuality = 0;
      if (bcf_get_info_int32(hdr, rec, "SRMAPQ", &mapq, &nmapq) > 0) svRec.srMapQuality = (uint8_t) *mapq;
      else svRec.srMapQuality = 0;
      if (bcf_get_info_float(hdr, rec, "SRQ", &srq, &nsrq) > 0) svRec.srAlignQuality = (double) *srq;
      else svRec.srAlignQuality = 0;

      svs.push_back(svRec);
    } else {
      std::cerr << "Error: Delly genotyping requires local SV assembly (INFO/CONSENSUS) and breakpoint (INFO/CONSBP) introduced in delly v1.1.7!" << std::endl;
      break;
    }
  }
  // Clean-up
  free(svend);
  free(pos2);
  free(consbp);
  free(svlen);
  free(svt);
  free(method);
  free(pe);
  free(inslen);
  free(homlen);
  free(sr);
  free(cons);
  free(cipos);
  free(mapq);
  free(srq);
  free(ct);
  free(chr2);

  // Close VCF
  bcf_hdr_destroy(hdr);
  bcf_close(ifile);
  bcf_destroy(rec);
}


template<typename TConfig, typename TStructuralVariantRecord, typename TJunctionCountMap, typename TReadCountMap, typename TCountMap>
inline void
vcfOutput(TConfig const& c, std::vector<TStructuralVariantRecord> const& svs, TJunctionCountMap const& jctCountMap, TReadCountMap const& readCountMap, TCountMap const& spanCountMap)
{
  // BoLog class
  BoLog<double> bl;

  // Open one bam file header
  samFile* samfile = sam_open(c.files[0].string().c_str(), "r");
  hts_set_fai_filename(samfile, c.genome.string().c_str());
  bam_hdr_t* bamhd = sam_hdr_read(samfile);

  // Output all structural variants
  std::string fmtout = "wb";
  if (c.outfile.string() == "-") fmtout = "w";
  htsFile *fp = hts_open(c.outfile.string().c_str(), fmtout.c_str());
  bcf_hdr_t *hdr = bcf_hdr_init("w");

  // Print vcf header
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  boost::gregorian::date today = now.date();
  std::string datestr("##fileDate=");
  datestr += boost::gregorian::to_iso_string(today);
  bcf_hdr_append(hdr, datestr.c_str());
  bcf_hdr_append(hdr, "##ALT=<ID=DEL,Description=\"Deletion\">");
  bcf_hdr_append(hdr, "##ALT=<ID=DUP,Description=\"Duplication\">");
  bcf_hdr_append(hdr, "##ALT=<ID=INV,Description=\"Inversion\">");
  bcf_hdr_append(hdr, "##ALT=<ID=BND,Description=\"Translocation\">");
  bcf_hdr_append(hdr, "##ALT=<ID=INS,Description=\"Insertion\">");
  bcf_hdr_append(hdr, "##FILTER=<ID=LowQual,Description=\"Poor quality and insufficient number of PEs and SRs.\">");
  bcf_hdr_append(hdr, "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"PE confidence interval around END\">");
  bcf_hdr_append(hdr, "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"PE confidence interval around POS\">");
  bcf_hdr_append(hdr, "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for POS2 coordinate in case of an inter-chromosomal translocation\">");
  bcf_hdr_append(hdr, "##INFO=<ID=POS2,Number=1,Type=Integer,Description=\"Genomic position for CHR2 in case of an inter-chromosomal translocation\">");
  bcf_hdr_append(hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">");
  bcf_hdr_append(hdr, "##INFO=<ID=PE,Number=1,Type=Integer,Description=\"Paired-end support of the structural variant\">");
  bcf_hdr_append(hdr, "##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of paired-ends\">");
  bcf_hdr_append(hdr, "##INFO=<ID=SRMAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of split-reads\">");
  bcf_hdr_append(hdr, "##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Split-read support\">");
  bcf_hdr_append(hdr, "##INFO=<ID=SRQ,Number=1,Type=Float,Description=\"Split-read consensus alignment quality\">");
  bcf_hdr_append(hdr, "##INFO=<ID=CONSENSUS,Number=1,Type=String,Description=\"Split-read consensus sequence\">");
  bcf_hdr_append(hdr, "##INFO=<ID=CONSBP,Number=1,Type=Integer,Description=\"Consensus SV breakpoint position\">");
  bcf_hdr_append(hdr, "##INFO=<ID=CE,Number=1,Type=Float,Description=\"Consensus sequence entropy\">");
  bcf_hdr_append(hdr, "##INFO=<ID=CT,Number=1,Type=String,Description=\"Paired-end signature induced connection type\">");
  bcf_hdr_append(hdr, "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Insertion length for SVTYPE=INS.\">");
  bcf_hdr_append(hdr, "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">");
  bcf_hdr_append(hdr, "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">");
  bcf_hdr_append(hdr, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">");
  bcf_hdr_append(hdr, "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">");
  bcf_hdr_append(hdr, "##INFO=<ID=INSLEN,Number=1,Type=Integer,Description=\"Predicted length of the insertion\">");
  bcf_hdr_append(hdr, "##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description=\"Predicted microhomology length using a max. edit distance of 2\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Log10-scaled genotype likelihoods for RR,RA,AA genotypes\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Per-sample genotype filter\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=RC,Number=1,Type=Integer,Description=\"Raw high-quality read counts or base counts for the SV\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=RCL,Number=1,Type=Integer,Description=\"Raw high-quality read counts or base counts for the left control region\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=RCR,Number=1,Type=Integer,Description=\"Raw high-quality read counts or base counts for the right control region\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=RDCN,Number=1,Type=Integer,Description=\"Read-depth based copy-number estimate for autosomal sites\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# high-quality reference pairs\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# high-quality variant pairs\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"# high-quality reference junction reads\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=RV,Number=1,Type=Integer,Description=\"# high-quality variant junction reads\">");

  // Add reference
  std::string refloc("##reference=");
  refloc += c.genome.string();
  bcf_hdr_append(hdr, refloc.c_str());
  for (int i = 0; i<bamhd->n_targets; ++i) {
    std::string refname("##contig=<ID=");
    refname += std::string(bamhd->target_name[i]) + ",length=" + boost::lexical_cast<std::string>(bamhd->target_len[i]) + ">";
    bcf_hdr_append(hdr, refname.c_str());
  }
  // Add samples
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) bcf_hdr_add_sample(hdr, c.sampleName[file_c].c_str());
  bcf_hdr_add_sample(hdr, NULL);
  if (bcf_hdr_write(fp, hdr) != 0) std::cerr << "Error: Failed to write BCF header!" << std::endl;

  if (!svs.empty()) {
    // Genotype arrays
    int32_t *gts = (int*) malloc(bcf_hdr_nsamples(hdr) * 2 * sizeof(int));
    float *gls = (float*) malloc(bcf_hdr_nsamples(hdr) * 3 * sizeof(float));
    int32_t *rcl = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    int32_t *rc = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    int32_t *rcr = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    int32_t *cnest = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    int32_t *drcount = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    int32_t *dvcount = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    int32_t *rrcount = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    int32_t *rvcount = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    int32_t *gqval = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    std::vector<std::string> ftarr;
    ftarr.resize(bcf_hdr_nsamples(hdr));
    
    // Iterate all structural variants
    typedef std::vector<TStructuralVariantRecord> TSVs;
    now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Genotyping" << std::endl;
    bcf1_t *rec = bcf_init();
    for(typename TSVs::const_iterator svIter = svs.begin(); svIter!=svs.end(); ++svIter) {
      if ((svIter->srSupport == 0) && (svIter->peSupport == 0)) continue;
      // In discovery mode, skip SVs that have less than 2 reads support after genotyping
      if (!c.hasVcfFile) {
	uint32_t totalGtSup = 0;
	for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	  totalGtSup += spanCountMap[file_c][svIter->id].alt.size() + jctCountMap[file_c][svIter->id].alt.size();
	}
	if (totalGtSup < 2) continue;
      }
      
      // Output main vcf fields
      int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
      if (svIter->chr == svIter->chr2) {
	// Intra-chromosomal
	if (((svIter->peSupport < 3) || (svIter->peMapQuality < 20)) && ((svIter->srSupport < 3) || (svIter->srMapQuality < 20))) tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "LowQual");
      } else {
	// Inter-chromosomal
	if (((svIter->peSupport < 5) || (svIter->peMapQuality < 20)) && ((svIter->srSupport < 5) || (svIter->srMapQuality < 20))) tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "LowQual");
      }
      rec->rid = bcf_hdr_name2id(hdr, bamhd->target_name[svIter->chr]);
      int32_t svStartPos = svIter->svStart - 1;
      if (svStartPos < 1) svStartPos = 1;
      int32_t svEndPos = svIter->svEnd;
      if (svIter->chr == svIter->chr2) --svEndPos; // To match VCF spec
      if (svEndPos < 1) svEndPos = 1;
      if (svEndPos >= (int32_t) bamhd->target_len[svIter->chr2]) svEndPos = bamhd->target_len[svIter->chr2] - 1;
      rec->pos = svStartPos;
      std::string id(_addID(svIter->svt));
      std::string padNumber = boost::lexical_cast<std::string>(svIter->id);
      padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
      id += padNumber;
      bcf_update_id(hdr, rec, id.c_str());
      std::string alleles = _replaceIUPAC(svIter->alleles);
      bcf_update_alleles_str(hdr, rec, alleles.c_str());
      bcf_update_filter(hdr, rec, &tmpi, 1);
      
      // Add INFO fields
      if (svIter->precise) bcf_update_info_flag(hdr, rec, "PRECISE", NULL, 1);
      else bcf_update_info_flag(hdr, rec, "IMPRECISE", NULL, 1);
      bcf_update_info_string(hdr, rec, "SVTYPE", _addID(svIter->svt).c_str());
      std::string dellyVersion("EMBL.DELLYv");
      dellyVersion += dellyVersionNumber;
      bcf_update_info_string(hdr,rec, "SVMETHOD", dellyVersion.c_str());
      if (svIter->svt < DELLY_SVT_TRANS) {
	tmpi = svEndPos;
	bcf_update_info_int32(hdr, rec, "END", &tmpi, 1);
      } else {
	tmpi = svStartPos + 2;
	bcf_update_info_int32(hdr, rec, "END", &tmpi, 1);
	bcf_update_info_string(hdr,rec, "CHR2", bamhd->target_name[svIter->chr2]);
	tmpi = svEndPos;
	bcf_update_info_int32(hdr, rec, "POS2", &tmpi, 1);
      }
      if (svIter->svt == 4) {
	tmpi = svIter->insLen;
	bcf_update_info_int32(hdr, rec, "SVLEN", &tmpi, 1);
      }
      tmpi = svIter->peSupport;
      bcf_update_info_int32(hdr, rec, "PE", &tmpi, 1);
      tmpi = svIter->peMapQuality;
      bcf_update_info_int32(hdr, rec, "MAPQ", &tmpi, 1);
      bcf_update_info_string(hdr, rec, "CT", _addOrientation(svIter->svt).c_str());
      int32_t ciend[2];
      ciend[0] = svIter->ciendlow;
      ciend[1] = svIter->ciendhigh;
      int32_t cipos[2];
      cipos[0] = svIter->ciposlow;
      cipos[1] = svIter->ciposhigh;
      bcf_update_info_int32(hdr, rec, "CIPOS", cipos, 2);
      bcf_update_info_int32(hdr, rec, "CIEND", ciend, 2);
      
      if (svIter->precise)  {
	tmpi = svIter->srMapQuality;
	bcf_update_info_int32(hdr, rec, "SRMAPQ", &tmpi, 1);
	tmpi = svIter->insLen;
	bcf_update_info_int32(hdr, rec, "INSLEN", &tmpi, 1);
	tmpi = svIter->homLen;
	bcf_update_info_int32(hdr, rec, "HOMLEN", &tmpi, 1);
	tmpi = svIter->srSupport;
	bcf_update_info_int32(hdr, rec, "SR", &tmpi, 1);
	float tmpf = svIter->srAlignQuality;
	bcf_update_info_float(hdr, rec, "SRQ", &tmpf, 1);
	if (svIter->consensus.size()) {
	  bcf_update_info_string(hdr, rec, "CONSENSUS", svIter->consensus.c_str());
	  tmpf = entropy(svIter->consensus);
	  bcf_update_info_float(hdr, rec, "CE", &tmpf, 1);
	  tmpi = svIter->consBp;
	  bcf_update_info_int32(hdr, rec, "CONSBP", &tmpi, 1);
	}
      }
      
      // Add genotype columns
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	// Counters
	rcl[file_c] = 0;
	rc[file_c] = 0;
	rcr[file_c] = 0;
	cnest[file_c] = 0;
	drcount[file_c] = 0;
	dvcount[file_c] = 0;
	rrcount[file_c] = 0;
	rvcount[file_c] = 0;
	drcount[file_c] = spanCountMap[file_c][svIter->id].ref.size();
	dvcount[file_c] = spanCountMap[file_c][svIter->id].alt.size();
	rrcount[file_c] = jctCountMap[file_c][svIter->id].ref.size();
	rvcount[file_c] = jctCountMap[file_c][svIter->id].alt.size();
	
	// Compute GLs
	if (svIter->precise) _computeGLs(bl, jctCountMap[file_c][svIter->id].ref, jctCountMap[file_c][svIter->id].alt, gls, gqval, gts, file_c);
	else _computeGLs(bl, spanCountMap[file_c][svIter->id].ref, spanCountMap[file_c][svIter->id].alt, gls, gqval, gts, file_c);
	
	// Compute RCs
	rcl[file_c] = readCountMap[file_c][svIter->id].leftRC;
	rc[file_c] = readCountMap[file_c][svIter->id].rc;
	rcr[file_c] = readCountMap[file_c][svIter->id].rightRC;
	cnest[file_c] = -1;
	if ((rcl[file_c] + rcr[file_c]) > 0) cnest[file_c] = boost::math::iround( 2.0 * (double) rc[file_c] / (double) (rcl[file_c] + rcr[file_c]) );
      
	// Genotype filter
	if (gqval[file_c] < 15) ftarr[file_c] = "LowQual";
	else ftarr[file_c] = "PASS";
      }
      int32_t qvalout = svIter->mapq;
      if (qvalout < 0) qvalout = 0;
      if (qvalout > 10000) qvalout = 10000;
      rec->qual = qvalout;
      
      bcf_update_genotypes(hdr, rec, gts, bcf_hdr_nsamples(hdr) * 2);
      bcf_update_format_float(hdr, rec, "GL",  gls, bcf_hdr_nsamples(hdr) * 3);
      bcf_update_format_int32(hdr, rec, "GQ", gqval, bcf_hdr_nsamples(hdr));
      std::vector<const char*> strp(bcf_hdr_nsamples(hdr));
      std::transform(ftarr.begin(), ftarr.end(), strp.begin(), cstyle_str());
      bcf_update_format_string(hdr, rec, "FT", &strp[0], bcf_hdr_nsamples(hdr));
      bcf_update_format_int32(hdr, rec, "RCL", rcl, bcf_hdr_nsamples(hdr));
      bcf_update_format_int32(hdr, rec, "RC", rc, bcf_hdr_nsamples(hdr));
      bcf_update_format_int32(hdr, rec, "RCR", rcr, bcf_hdr_nsamples(hdr));
      bcf_update_format_int32(hdr, rec, "RDCN", cnest, bcf_hdr_nsamples(hdr));
      bcf_update_format_int32(hdr, rec, "DR", drcount, bcf_hdr_nsamples(hdr));
      bcf_update_format_int32(hdr, rec, "DV", dvcount, bcf_hdr_nsamples(hdr));
      bcf_update_format_int32(hdr, rec, "RR", rrcount, bcf_hdr_nsamples(hdr));
      bcf_update_format_int32(hdr, rec, "RV", rvcount, bcf_hdr_nsamples(hdr));
      bcf_write1(fp, hdr, rec);
      bcf_clear1(rec);
    }
    bcf_destroy1(rec);
    
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
  }

  // Close BAM file
  bam_hdr_destroy(bamhd);
  sam_close(samfile);

  // Close VCF file
  bcf_hdr_destroy(hdr);
  hts_close(fp);

  // Build index
  if (c.outfile.string() != "-") bcf_index_build(c.outfile.string().c_str(), 14);
}


}

#endif
