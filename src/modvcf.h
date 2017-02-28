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
  std::string tmp[] = {"CT", "PRECISE", "IMPRECISE", "SVTYPE", "SVMETHOD", "CIEND", "CIPOS", "CHR2", "END", "PE", "MAPQ", "SR", "SRQ", "CONSENSUS"};
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


 // Convert string to char*
struct cstyle_str {
  const char* operator ()(const std::string& s) {
    return s.c_str();
  }
};


// Parse Delly vcf file
template<typename TConfig, typename TStructuralVariantRecord, typename TTag>
inline void
vcfParse(TConfig const& c, bam_hdr_t* hd, std::vector<TStructuralVariantRecord>& svs, SVType<TTag> svType)
{
  // Load bcf file
  htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
  bcf_hdr_t* hdr = bcf_hdr_read(ifile);
  bcf1_t* rec = bcf_init();

  // Parse bcf
  int32_t nsvend = 0;
  int32_t* svend = NULL;
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
  int32_t ncons = 0;
  char* cons = NULL;
  int32_t nchr2 = 0;
  char* chr2 = NULL;
  while (bcf_read(ifile, hdr, rec) == 0) {
    bcf_unpack(rec, BCF_UN_INFO);

    // Correct SV type
    if (bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt) > 0) {
      if (std::string(svt) != _addID(svType)) continue;
    } else continue;

    // Fill SV record
    StructuralVariantRecord svRec;
    std::string chrName = bcf_hdr_id2name(hdr, rec->rid);
    int32_t tid = bam_name2id(hd, chrName.c_str());
    svRec.chr = tid;
    svRec.svStart = rec->pos + 1;
    svRec.id = svs.size();
    std::string refAllele = rec->d.allele[0];
    std::string altAllele = rec->d.allele[1];
    svRec.alleles = refAllele + "," + altAllele;

    // Parse INFO
    if (bcf_get_info_flag(hdr, rec, "PRECISE", 0, 0) > 0) svRec.precise=true;
    else svRec.precise = false;
    if (bcf_get_info_int32(hdr, rec, "PE", &pe, &npe) > 0) svRec.peSupport = *pe;
    else {
      if (svRec.precise) svRec.peSupport = 0;
      else svRec.peSupport = 2;
    }
    if (bcf_get_info_int32(hdr, rec, "INSLEN", &inslen, &ninslen) > 0) svRec.insLen = *inslen;
    else svRec.insLen = 0;
    if (bcf_get_info_int32(hdr, rec, "HOMLEN", &homlen, &nhomlen) > 0) svRec.homLen = *homlen;
    else svRec.homLen = 0;
    if (bcf_get_info_int32(hdr, rec, "SR", &sr, &nsr) > 0) svRec.srSupport = *sr;
    else svRec.srSupport = 0;
    if (bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend) > 0) svRec.svEnd = *svend;
    else svRec.svEnd = rec->pos + 2;
    if (bcf_get_info_string(hdr, rec, "CONSENSUS", &cons, &ncons) > 0) svRec.consensus = std::string(cons);
    else svRec.precise = false;
    if (bcf_get_info_int32(hdr, rec, "CIPOS", &cipos, &ncipos) > 0) svRec.wiggle = cipos[1];
    else svRec.wiggle = 0;
    if (bcf_get_info_int32(hdr, rec, "MAPQ", &mapq, &nmapq) > 0) svRec.peMapQuality = (uint8_t) *mapq;
    else svRec.peMapQuality = 0;
    if (bcf_get_info_float(hdr, rec, "SRQ", &srq, &nsrq) > 0) svRec.srAlignQuality = (double) *srq;
    else svRec.srAlignQuality = 0;
    if (bcf_get_info_string(hdr, rec, "CT", &ct, &nct) > 0) svRec.ct = _decodeOrientation(std::string(ct));
    else {
      if (_addID(svType) == "DEL") svRec.ct = _decodeOrientation("3to5");
      else if (_addID(svType) == "DUP") svRec.ct = _decodeOrientation("5to3");
      else if (_addID(svType) == "INS") svRec.ct = _decodeOrientation("NtoN");
      else if (_addID(svType) == "INV") svRec.ct = _decodeOrientation("3to3"); // Insufficient
      else if (_addID(svType) == "TRA") svRec.ct = _decodeOrientation("3to5"); // Insufficient
    }
    if (bcf_get_info_string(hdr, rec, "CHR2", &chr2, &nchr2) > 0) {
      std::string chr2Name = std::string(chr2);
      svRec.chr2 = bam_name2id(hd, chr2Name.c_str());
    } else svRec.chr2 = tid;
    svs.push_back(svRec);
  }
  // Clean-up
  free(svend);
  free(svt);
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


template<typename TConfig, typename TStructuralVariantRecord, typename TJunctionCountMap, typename TReadCountMap, typename TCountMap, typename TTag>
inline void
vcfOutput(TConfig const& c, std::vector<TStructuralVariantRecord> const& svs, TJunctionCountMap const& jctCountMap, TReadCountMap const& readCountMap, TCountMap const& spanCountMap, SVType<TTag> svType) 
{
  // BoLog class
  BoLog<double> bl;

  // Open one bam file header
  samFile* samfile = sam_open(c.files[0].string().c_str(), "r");
  bam_hdr_t* bamhd = sam_hdr_read(samfile);

  // Output all structural variants
  htsFile *fp = hts_open(c.outfile.string().c_str(), "wb");
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
  bcf_hdr_append(hdr, "##ALT=<ID=TRA,Description=\"Translocation\">");
  bcf_hdr_append(hdr, "##ALT=<ID=INS,Description=\"Insertion\">");
  bcf_hdr_append(hdr, "##FILTER=<ID=LowQual,Description=\"PE/SR support below 3 or mapping quality below 20.\">");
  bcf_hdr_append(hdr, "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"PE confidence interval around END\">");
  bcf_hdr_append(hdr, "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"PE confidence interval around POS\">");
  bcf_hdr_append(hdr, "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">");
  bcf_hdr_append(hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">");
  bcf_hdr_append(hdr, "##INFO=<ID=PE,Number=1,Type=Integer,Description=\"Paired-end support of the structural variant\">");
  bcf_hdr_append(hdr, "##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of paired-ends\">");
  bcf_hdr_append(hdr, "##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Split-read support\">");
  bcf_hdr_append(hdr, "##INFO=<ID=SRQ,Number=1,Type=Float,Description=\"Split-read consensus alignment quality\">");
  bcf_hdr_append(hdr, "##INFO=<ID=CONSENSUS,Number=1,Type=String,Description=\"Split-read consensus sequence\">");
  bcf_hdr_append(hdr, "##INFO=<ID=CE,Number=1,Type=Float,Description=\"Consensus sequence entropy\">");
  bcf_hdr_append(hdr, "##INFO=<ID=CT,Number=1,Type=String,Description=\"Paired-end signature induced connection type\">");
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
  bcf_hdr_append(hdr, "##FORMAT=<ID=RC,Number=1,Type=Integer,Description=\"Raw high-quality read counts for the SV\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=RCL,Number=1,Type=Integer,Description=\"Raw high-quality read counts for the left control region\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=RCR,Number=1,Type=Integer,Description=\"Raw high-quality read counts for the right control region\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Read-depth based copy-number estimate for autosomal sites\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# high-quality reference pairs\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# high-quality variant pairs\">");
  if (c.isHaplotagged) {
    bcf_hdr_append(hdr, "##FORMAT=<ID=HP1DR,Number=1,Type=Integer,Description=\"# high-quality reference pairs on haplotype 1\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=HP2DR,Number=1,Type=Integer,Description=\"# high-quality reference pairs on haplotype 2\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=HP1DV,Number=1,Type=Integer,Description=\"# high-quality variant pairs on haplotype 1\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=HP2DV,Number=1,Type=Integer,Description=\"# high-quality variant pairs on haplotype 2\">");
  }
  bcf_hdr_append(hdr, "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"# high-quality reference junction reads\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=RV,Number=1,Type=Integer,Description=\"# high-quality variant junction reads\">");
  if (c.isHaplotagged) {
    bcf_hdr_append(hdr, "##FORMAT=<ID=HP1RR,Number=1,Type=Integer,Description=\"# high-quality reference junction reads on haplotype 1\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=HP2RR,Number=1,Type=Integer,Description=\"# high-quality reference junction reads on haplotype 2\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=HP1RV,Number=1,Type=Integer,Description=\"# high-quality variant junction reads on haplotype 1\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=HP2RV,Number=1,Type=Integer,Description=\"# high-quality variant junction reads on haplotype 2\">");
  }
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
  bcf_hdr_write(fp, hdr);

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
    int32_t *hp1drcount = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    int32_t *hp2drcount = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    int32_t *hp1dvcount = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    int32_t *hp2dvcount = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    int32_t *rrcount = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    int32_t *rvcount = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    int32_t *hp1rrcount = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    int32_t *hp2rrcount = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    int32_t *hp1rvcount = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    int32_t *hp2rvcount = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    int32_t *gqval = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    std::vector<std::string> ftarr;
    ftarr.resize(bcf_hdr_nsamples(hdr));
    
    // Iterate all structural variants
    typedef std::vector<TStructuralVariantRecord> TSVs;
    uint32_t lastId = svs.size();
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Genotyping" << std::endl;
    boost::progress_display show_progress( svs.size() );
    bcf1_t *rec = bcf_init();
    for(typename TSVs::const_iterator svIter = svs.begin(); svIter!=svs.end(); ++svIter) {
      ++show_progress;

      // Output main vcf fields
      int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
      if ((svIter->precise) && (svIter->chr == svIter->chr2)) {
	if ((svIter->srSupport < 3) || (svIter->peMapQuality < 20)) tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "LowQual");
      } else {
	if ((svIter->peSupport < 3) || (svIter->peMapQuality < 20) || ( (svIter->chr != svIter->chr2) && (svIter->peSupport < 5) ) ) tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "LowQual");
      }
      rec->rid = bcf_hdr_name2id(hdr, bamhd->target_name[svIter->chr]);
      rec->pos = svIter->svStart - 1;
      std::string id(_addID(svType));
      std::string padNumber = boost::lexical_cast<std::string>(svIter->id);
      padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
      id += padNumber;
      bcf_update_id(hdr, rec, id.c_str());
      bcf_update_alleles_str(hdr, rec, svIter->alleles.c_str());
      bcf_update_filter(hdr, rec, &tmpi, 1);
      
      // Add INFO fields
      if (svIter->precise) bcf_update_info_flag(hdr, rec, "PRECISE", NULL, 1);
      else bcf_update_info_flag(hdr, rec, "IMPRECISE", NULL, 1);
      bcf_update_info_string(hdr, rec, "SVTYPE", _addID(svType).c_str());
      std::string dellyVersion("EMBL.DELLYv");
      dellyVersion += dellyVersionNumber;
      bcf_update_info_string(hdr,rec, "SVMETHOD", dellyVersion.c_str());
      bcf_update_info_string(hdr,rec, "CHR2", bamhd->target_name[svIter->chr2]);
      tmpi = svIter->svEnd;
      bcf_update_info_int32(hdr, rec, "END", &tmpi, 1);
      tmpi = svIter->peSupport;
      bcf_update_info_int32(hdr, rec, "PE", &tmpi, 1);
      tmpi = svIter->peMapQuality;
      bcf_update_info_int32(hdr, rec, "MAPQ", &tmpi, 1);
      bcf_update_info_string(hdr, rec, "CT", _addOrientation(svIter->ct).c_str());
      int32_t ciend[2];
      ciend[0] = -svIter->wiggle;
      ciend[1] = svIter->wiggle;
      int32_t cipos[2];
      cipos[0] = -svIter->wiggle;
      cipos[1] = svIter->wiggle;
      bcf_update_info_int32(hdr, rec, "CIPOS", cipos, 2);
      bcf_update_info_int32(hdr, rec, "CIEND", ciend, 2);
      
      if (svIter->precise)  {
	tmpi = svIter->insLen;
	bcf_update_info_int32(hdr, rec, "INSLEN", &tmpi, 1);
	tmpi = svIter->homLen;
	bcf_update_info_int32(hdr, rec, "HOMLEN", &tmpi, 1);
	tmpi = svIter->srSupport;
	bcf_update_info_int32(hdr, rec, "SR", &tmpi, 1);
	float tmpf = svIter->srAlignQuality;
	bcf_update_info_float(hdr, rec, "SRQ", &tmpf, 1);
	bcf_update_info_string(hdr, rec, "CONSENSUS", svIter->consensus.c_str());
	tmpf = entropy(svIter->consensus);
	bcf_update_info_float(hdr, rec, "CE", &tmpf, 1);
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
	if (c.isHaplotagged) {
	  hp1drcount[file_c] = 0;
	  hp2drcount[file_c] = 0;
	  hp1dvcount[file_c] = 0;
	  hp2dvcount[file_c] = 0;
	}
	rrcount[file_c] = 0;
	rvcount[file_c] = 0;
	if (c.isHaplotagged) {
	  hp1rrcount[file_c] = 0;
	  hp2rrcount[file_c] = 0;
	  hp1rvcount[file_c] = 0;
	  hp2rvcount[file_c] = 0;
	}
	if (spanCountMap[file_c][svIter->id].ref.size() < spanCountMap[file_c][lastId + svIter->id].ref.size()) {
	  drcount[file_c] = spanCountMap[file_c][svIter->id].ref.size();
	  dvcount[file_c] = spanCountMap[file_c][svIter->id].alt.size();
	  if (c.isHaplotagged) {
	    hp1drcount[file_c] = spanCountMap[file_c][svIter->id].refh1;
	    hp2drcount[file_c] = spanCountMap[file_c][svIter->id].refh2;
	    hp1dvcount[file_c] = spanCountMap[file_c][svIter->id].alth1;
	    hp2dvcount[file_c] = spanCountMap[file_c][svIter->id].alth2;
	  }
	} else {
	  drcount[file_c] = spanCountMap[file_c][lastId + svIter->id].ref.size();
	  dvcount[file_c] = spanCountMap[file_c][lastId + svIter->id].alt.size();
	  if (c.isHaplotagged) {
	    hp1drcount[file_c] = spanCountMap[file_c][lastId + svIter->id].refh1;
	    hp2drcount[file_c] = spanCountMap[file_c][lastId + svIter->id].refh2;
	    hp1dvcount[file_c] = spanCountMap[file_c][lastId + svIter->id].alth1;
	    hp2dvcount[file_c] = spanCountMap[file_c][lastId + svIter->id].alth2;
	  }
	}
	rrcount[file_c] = jctCountMap[file_c][svIter->id].ref.size();
	rvcount[file_c] = jctCountMap[file_c][svIter->id].alt.size();
	if (c.isHaplotagged) {
	  hp1rrcount[file_c] = jctCountMap[file_c][svIter->id].refh1;
	  hp2rrcount[file_c] = jctCountMap[file_c][svIter->id].refh2;
	  hp1rvcount[file_c] = jctCountMap[file_c][svIter->id].alth1;
	  hp2rvcount[file_c] = jctCountMap[file_c][svIter->id].alth2;
	}
	
	// Compute GLs
	if (svIter->precise) _computeGLs(bl, jctCountMap[file_c][svIter->id].ref, jctCountMap[file_c][svIter->id].alt, gls, gqval, gts, file_c);
	else {  // Imprecise SVs
	  if (spanCountMap[file_c][svIter->id].ref.size() < spanCountMap[file_c][lastId + svIter->id].ref.size()) 
	    _computeGLs(bl, spanCountMap[file_c][svIter->id].ref, spanCountMap[file_c][svIter->id].alt, gls, gqval, gts, file_c);
	  else
	    _computeGLs(bl, spanCountMap[file_c][lastId + svIter->id].ref, spanCountMap[file_c][lastId + svIter->id].alt, gls, gqval, gts, file_c);
	}
	
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
      // ToDo
      //rec->qual = 0;
      
      
      bcf_update_genotypes(hdr, rec, gts, bcf_hdr_nsamples(hdr) * 2);
      bcf_update_format_float(hdr, rec, "GL",  gls, bcf_hdr_nsamples(hdr) * 3);
      bcf_update_format_int32(hdr, rec, "GQ", gqval, bcf_hdr_nsamples(hdr));
      std::vector<const char*> strp(bcf_hdr_nsamples(hdr));
      std::transform(ftarr.begin(), ftarr.end(), strp.begin(), cstyle_str());
      bcf_update_format_string(hdr, rec, "FT", &strp[0], bcf_hdr_nsamples(hdr));
      bcf_update_format_int32(hdr, rec, "RCL", rcl, bcf_hdr_nsamples(hdr));
      bcf_update_format_int32(hdr, rec, "RC", rc, bcf_hdr_nsamples(hdr));
      bcf_update_format_int32(hdr, rec, "RCR", rcr, bcf_hdr_nsamples(hdr));
      bcf_update_format_int32(hdr, rec, "CN", cnest, bcf_hdr_nsamples(hdr));
      bcf_update_format_int32(hdr, rec, "DR", drcount, bcf_hdr_nsamples(hdr));
      bcf_update_format_int32(hdr, rec, "DV", dvcount, bcf_hdr_nsamples(hdr));
      if (c.isHaplotagged) {
	bcf_update_format_int32(hdr, rec, "HP1DR", hp1drcount, bcf_hdr_nsamples(hdr));
	bcf_update_format_int32(hdr, rec, "HP2DR", hp2drcount, bcf_hdr_nsamples(hdr));
	bcf_update_format_int32(hdr, rec, "HP1DV", hp1dvcount, bcf_hdr_nsamples(hdr));
	bcf_update_format_int32(hdr, rec, "HP2DV", hp2dvcount, bcf_hdr_nsamples(hdr));
      }
      bcf_update_format_int32(hdr, rec, "RR", rrcount, bcf_hdr_nsamples(hdr));
      bcf_update_format_int32(hdr, rec, "RV", rvcount, bcf_hdr_nsamples(hdr));
      if (c.isHaplotagged) {
	bcf_update_format_int32(hdr, rec, "HP1RR", hp1rrcount, bcf_hdr_nsamples(hdr));
	bcf_update_format_int32(hdr, rec, "HP2RR", hp2rrcount, bcf_hdr_nsamples(hdr));
	bcf_update_format_int32(hdr, rec, "HP1RV", hp1rvcount, bcf_hdr_nsamples(hdr));
	bcf_update_format_int32(hdr, rec, "HP2RV", hp2rvcount, bcf_hdr_nsamples(hdr));
      }
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
    free(hp1drcount);
    free(hp2drcount);
    free(hp1dvcount);
    free(hp2dvcount);
    free(rrcount);
    free(rvcount);
    free(hp1rrcount);
    free(hp2rrcount);
    free(hp1rvcount);
    free(hp2rvcount);
    free(gqval);
  }

  // Close BAM file
  bam_hdr_destroy(bamhd);
  sam_close(samfile);

  // Close VCF file
  bcf_hdr_destroy(hdr);
  hts_close(fp);

  // Build index
  bcf_index_build(c.outfile.string().c_str(), 14);
}
 

}

#endif
