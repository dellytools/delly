#ifndef CORAL_H
#define CORAL_H

#include <limits>

#include <boost/icl/split_interval_map.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/math/special_functions/round.hpp>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "bed.h"
#include "scan.h"
#include "gcbias.h"
#include "cnv.h"
#include "version.h"

namespace torali
{

  struct CountDNAConfig {
    bool basecov;
    bool adaptive;
    bool hasMapFile;
    bool hasStatsFile;
    bool hasBedFile;
    bool hasScanFile;
    bool noScanWindowSelection;
    bool segmentation;
    bool hasGenoFile;
    bool hasVcfFile;
    uint32_t nchr;
    uint32_t meanisize;
    uint32_t window_size;
    uint32_t window_offset;
    uint32_t scanWindow;
    uint32_t minChrLen;
    uint32_t minCnvSize;
    uint32_t targetReads;
    double targetExpCov;
    uint16_t minQual;
    uint16_t mapqUniq;
    uint16_t mad;
    float ploidy;
    float ctrlPloidy;
    float expectedCN;
    float purity;
    float exclgc;
    float uniqueToTotalCovRatio;
    float fracWindow;
    float fragmentUnique;
    float controlMaf;
    float stringency;
    float cn_offset;
    std::string sampleName;
    boost::filesystem::path vcffile;
    boost::filesystem::path genofile;
    boost::filesystem::path outfile;
    boost::filesystem::path covfile;
    boost::filesystem::path genome;
    boost::filesystem::path statsFile;
    boost::filesystem::path mapFile;
    boost::filesystem::path bamFile;
    boost::filesystem::path bedFile;
    boost::filesystem::path scanFile;
    std::set<int32_t> refIdx;
  };

  struct CovWin {
    uint32_t start;
    uint32_t end;
    uint32_t winlen;
    double covsum;
    double expcov;
    double ucov;
    double tcov;
    bool valid;

    CovWin(uint32_t const s, uint32_t const e, uint32_t const w, double const cs, double const ec, double const uc, double const tc, bool vld) : start(s), end(e), winlen(w), covsum(cs), expcov(ec), ucov(uc), tcov(tc), valid(vld) {}
  };

  struct CountDNAConfigLib {
    uint16_t madCutoff;
    uint16_t madNormalCutoff;
    boost::filesystem::path genome;
    std::vector<boost::filesystem::path> files;
  };
  
  template<typename TConfig>
  inline int32_t
  bamCount(TConfig const& c, LibraryInfo const& li, std::vector<GcBias> const& gcbias, std::pair<uint32_t, uint32_t> const& gcbound) {
    // Load bam file
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // BED regions
    typedef std::set<std::pair<uint32_t, uint32_t> > TChrIntervals;
    typedef std::vector<TChrIntervals> TRegionsGenome;
    TRegionsGenome bedRegions;
    if (c.hasBedFile) {
      if (!_parsePotOverlappingIntervals(c.bedFile.string(), c.hasBedFile, hdr, bedRegions)) {
	std::cerr << "Couldn't parse BED intervals. Do the chromosome names match?" << std::endl;
	return 1;
      }
    }

    // Parse BAM file
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Count fragments" << std::endl;

    // Open output files
    boost::iostreams::filtering_ostream dataOut;
    if (!c.covfile.empty()) {
      dataOut.push(boost::iostreams::gzip_compressor());
      dataOut.push(boost::iostreams::file_sink(c.covfile.string(), std::ios_base::out | std::ios_base::binary));
      dataOut << "chr\tstart\tend\t" << c.sampleName << "_uniqfrac\t" << c.sampleName << "_logR\t" << c.sampleName << "_CN" << std::endl;
    }

    // CNVs
    std::vector<CNV> cnvs;
    if (c.hasGenoFile) parseVcfCNV(c, hdr, cnvs);

    // SVs for breakpoint refinement
    typedef std::vector<SVBreakpoint> TChrBreakpoints;
    typedef std::vector<TChrBreakpoints> TGenomicBreakpoints;
    TGenomicBreakpoints svbp(c.nchr, TChrBreakpoints());
    if (c.hasVcfFile) {
      std::vector<StructuralVariantRecord> svs;
      vcfParse(c, hdr, svs);
      for(uint32_t i = 0; i < svs.size(); ++i) {
	svbp[svs[i].chr].push_back(SVBreakpoint(svs[i].svStart, svs[i].ciposlow, svs[i].ciposhigh, svs[i].mapq));
	svbp[svs[i].chr2].push_back(SVBreakpoint(svs[i].svEnd, svs[i].ciendlow, svs[i].ciendhigh, svs[i].mapq));
      }
      for (uint32_t i = 0; i < svbp.size(); ++i) sort(svbp[i].begin(), svbp[i].end());
    }
    
    // Iterate chromosomes
    faidx_t* faiMap = NULL;
    if (c.hasMapFile) faiMap = fai_load(c.mapFile.string().c_str());
    faidx_t* faiRef = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      if ((!c.hasGenoFile) && (chrNoData(c, refIndex, idx))) continue;

      // Haploid chromosome?
      float chrCtrlPloidy = c.ctrlPloidy;
      float chrPloidy = c.ploidy;
      if (c.refIdx.find(refIndex) != c.refIdx.end()) {
	chrCtrlPloidy = c.ctrlPloidy - 1;
	chrPloidy = c.ploidy - 1;
      }
      
      // Reference sequence
      std::string tname(hdr->target_name[refIndex]);
      int32_t seqlen = faidx_seq_len(faiRef, tname.c_str());
      if (seqlen == - 1) continue;
      else seqlen = -1;
      char* ref = faidx_fetch_seq(faiRef, tname.c_str(), 0, faidx_seq_len(faiRef, tname.c_str()), &seqlen);
      if (ref == NULL) continue;

      // Mappability map
      char* seq = NULL;
      if (c.hasMapFile) {
	seqlen = faidx_seq_len(faiMap, tname.c_str());
	if (seqlen == - 1) {
	  if (ref != NULL) free(ref);
	  continue;
	} else seqlen = -1;
	seq = faidx_fetch_seq(faiMap, tname.c_str(), 0, faidx_seq_len(faiMap, tname.c_str()), &seqlen);
	if (seq == NULL) { free(ref); continue; }
      }

      // Get GC and Mappability
      std::vector<uint16_t> uniqContent(hdr->target_len[refIndex], 0);
      std::vector<uint16_t> gcContent(hdr->target_len[refIndex], 0);
      {
	// GC map
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet gcref(hdr->target_len[refIndex], false);
	for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	  if ((ref[i] == 'c') || (ref[i] == 'C') || (ref[i] == 'g') || (ref[i] == 'G')) gcref[i] = 1;
	}

	// Mappability map
	TBitSet uniq(hdr->target_len[refIndex], false);
	if (c.hasMapFile) {
	  for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	    if (seq[i] == 'C') uniq[i] = 1;
	  }
	}

	// Sum across fragment
	int32_t halfwin = (int32_t) (c.meanisize / 2);
	int32_t usum = 0;
	int32_t gcsum = 0;
	for(int32_t pos = halfwin; pos < (int32_t) hdr->target_len[refIndex] - halfwin; ++pos) {
	  if (pos == halfwin) {
	    for(int32_t i = pos - halfwin; i<=pos+halfwin; ++i) {
	      usum += uniq[i];
	      gcsum += gcref[i];
	    }
	  } else {
	    usum -= uniq[pos - halfwin - 1];
	    gcsum -= gcref[pos - halfwin - 1];
	    usum += uniq[pos + halfwin];
	    gcsum += gcref[pos + halfwin];
	  }
	  gcContent[pos] = gcsum;
	  if (c.hasMapFile) uniqContent[pos] = usum;
	}
      }

      // Coverage track
      typedef uint16_t TCount;
      uint32_t maxCoverage = std::numeric_limits<TCount>::max();
      typedef std::vector<TCount> TCoverage;
      TCoverage cov(hdr->target_len[refIndex], 0);
      TCoverage covUniq;
      if (!c.hasMapFile) covUniq.resize(hdr->target_len[refIndex], 0);

      {
	// Mate map
	typedef boost::unordered_map<std::size_t, bool> TMateMap;
	TMateMap mateMap;
	
	// Count reads
	hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	int32_t lastAlignedPos = 0;
	std::set<std::size_t> lastAlignedPosReads;
	while (sam_itr_next(samfile, iter, rec) >= 0) {
	  if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
	  if (rec->core.qual < c.minQual) continue;
	  if ((rec->core.flag & BAM_FPAIRED) && ((rec->core.flag & BAM_FMUNMAP) || (rec->core.tid != rec->core.mtid))) continue;

	  // Base coverage
	  if (c.basecov) {
	    if (c.hasMapFile) addBaseCoverage(rec, cov, hdr->target_len[refIndex], maxCoverage);
	    else addBaseCoverage(rec, cov, covUniq, c.mapqUniq, hdr->target_len[refIndex], maxCoverage);
	    continue;
	  }

	  // Fragment coverage
	  int32_t midPoint = rec->core.pos + halfAlignmentLength(rec);
	  if (rec->core.flag & BAM_FPAIRED) {
	    std::size_t seed = hash_sr(rec);
	    
	    // Clean-up the read store for identical alignment positions
	    if (rec->core.pos > lastAlignedPos) {
	      lastAlignedPosReads.clear();
	      lastAlignedPos = rec->core.pos;
	    }
	    
	    if (_firstPairObs(rec, lastAlignedPosReads)) {
	      // First read
	      lastAlignedPosReads.insert(seed);
	      std::size_t hv = hash_pair(rec);
	      mateMap[hv] = true;
	      continue;
	    } else {
	      // Second read
	      std::size_t hv = hash_pair_mate(rec);
	      auto itMM = mateMap.find(hv);
	      if ((itMM == mateMap.end()) || (!itMM->second)) continue; // Mate discarded
	      mateMap.erase(itMM);
	    }
	    
	    // update midpoint
	    int32_t isize = (rec->core.pos + alignmentLength(rec)) - rec->core.mpos;
	    if ((li.minNormalISize < isize) && (isize < li.maxNormalISize)) midPoint = rec->core.mpos + (int32_t) (isize/2);
	  }
	  
	  // Count fragment
	  if ((midPoint >= 0) && (midPoint < (int32_t) hdr->target_len[refIndex]) && (cov[midPoint] < maxCoverage - 1)) ++cov[midPoint];
	}
	// Clean-up
	if (seq != NULL) free(seq);
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }

      // No mappability: unique if mostly high-MAPQ reads
      if (!c.hasMapFile) {
	for(uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
	  bool u;
	  if (cov[pos] == 0) u = ((ref[pos] != 'N') && (ref[pos] != 'n'));
	  else u = (2 * (uint32_t) covUniq[pos] >= (uint32_t) cov[pos]);
	  uniqContent[pos] = (u ? (uint16_t) c.meanisize : 0);
	}
      }
      if (ref != NULL) free(ref);

      // CNV discovery
      if (!c.hasGenoFile) {
	// Call CNVs
	std::vector<CNV> chrcnv;
	callCNVs(c, gcbound, gcContent, uniqContent, gcbias, cov, hdr, refIndex, chrcnv);

	// Merge adjacent CNVs lacking read-depth shift
	mergeCNVs(c, chrcnv, cnvs);

	// Refine breakpoints
	if (c.hasVcfFile) breakpointRefinement(c, gcbound, gcContent, uniqContent, gcbias, cov, hdr, refIndex, svbp, cnvs);
      }
      
      // CNV genotyping
      genotypeCNVs(c, gcbound, gcContent, uniqContent, gcbias, cov, hdr, refIndex, cnvs);

      // BED File (target intervals)
      if (c.hasBedFile) {
	for(typename TChrIntervals::iterator it = bedRegions[refIndex].begin(); it != bedRegions[refIndex].end(); ++it) {
	  if ((it->first < it->second) && (it->second <= hdr->target_len[refIndex])) {
	    double covsum = 0;
	    double expcov = 0;
	    //double obsexp = 0;
	    uint32_t winlen = 0;
	    for(uint32_t pos = it->first; pos < it->second; ++pos) {
	      if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
		covsum += cov[pos];
		//obsexp += gcbias[gcContent[pos]].obsexp;
		expcov += gcbias[gcContent[pos]].coverage;
		++winlen;
	      }
	    }
	    if (winlen >= c.fracWindow * (it->second - it->first)) {
	      //obsexp /= (double) winlen;
	      // Normalized counts: double count = ((double) covsum / obsexp ) * (double) (it->second - it->first) / (double) winlen;
	      double cn = chrPloidy;
	      double logR = 0;
	      if (expcov > 0) {
		cn = (c.expectedCN * covsum / expcov - chrCtrlPloidy * (1 - c.purity)) / c.purity;
		logR = std::log2((covsum + 1.0) / (expcov + 1.0));
	      }
	      if (!c.covfile.empty()) dataOut << std::string(hdr->target_name[refIndex]) << "\t" << it->first << "\t" << it->second << "\t" << winlen << "\t" << logR << "\t" << cn << std::endl;
	    } else {
	      if (!c.covfile.empty()) dataOut << std::string(hdr->target_name[refIndex]) << "\t" << it->first << "\t" << it->second << "\tNA\tNA\tNA" << std::endl;
	    }
	  }
	}
      } else {
	// Genome-wide
	std::vector<CovWin> wins;
	if (c.adaptive) {
	  double covsum = 0;
	  double expcov = 0;
	  double ucov = 0;
	  double tcov = 0;
	  uint32_t winlen = 0;
	  uint32_t start = 0;
	  for(uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
	    if (!c.hasMapFile) {
	      ucov += covUniq[pos];
	      tcov += cov[pos];
	    }
	    if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	      covsum += cov[pos];
	      expcov += gcbias[gcContent[pos]].coverage;
	      ++winlen;
	      if (expcov >= c.targetExpCov) {
		wins.push_back(CovWin(start, pos + 1, winlen, covsum, expcov, ucov, tcov, true));
		covsum = 0;
		expcov = 0;
		ucov = 0;
		tcov = 0;
		winlen = 0;
		start = pos + 1;
	      }
	    }
	  }
	} else {
	  // Fixed windows (non-overlapping)
	  for(uint32_t start = 0; start < hdr->target_len[refIndex]; start = start + c.window_offset) {
	    if (start + c.window_size < hdr->target_len[refIndex]) {
	      double covsum = 0;
	      double expcov = 0;
	      double ucov = 0;
	      double tcov = 0;
	      uint32_t winlen = 0;
	      for(uint32_t pos = start; pos < start + c.window_size; ++pos) {
		if (!c.hasMapFile) {
		  ucov += covUniq[pos];
		  tcov += cov[pos];
		}
		if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
		  covsum += cov[pos];
		  expcov += gcbias[gcContent[pos]].coverage;
		  ++winlen;
		}
	      }
	      bool valid = (winlen >= c.fracWindow * c.window_size);
	      wins.push_back(CovWin(start, start + c.window_size, winlen, covsum, expcov, ucov, tcov, valid));
	    }
	  }
	}

	// Separate true hom. dels from unmappable
	uint32_t nw = wins.size();
	std::vector<bool> naFlag(nw, false);
	std::vector<bool> suspect(nw, false);
	std::vector<bool> strong(nw, false);
	double lowFrac = 0.1;         // close to CN0
	double flankFrac = 0.5;       // clear no CN0
	uint32_t maxHomDel = 1000000; // max size for hom. DEL
	for(uint32_t i = 0; i < nw; ++i) {
	  if ((!wins[i].valid) || (wins[i].expcov <= 0)) {
	    naFlag[i] = true;
	    continue;
	  }
	  double r = wins[i].covsum / wins[i].expcov;
	  suspect[i] = (r < lowFrac);
	  strong[i] = (r >= flankFrac);
	}
	for(uint32_t i = 0; i < nw; ) {
	  if (naFlag[i] || (!suspect[i])) {
	    ++i;
	    continue;
	  }
	  uint32_t a = i;
	  uint32_t b = i;
	  while ((b + 1 < nw) && (!naFlag[b+1]) && suspect[b+1]) ++b;
	  uint32_t runBp = wins[b].end - wins[a].start;
	  bool leftStrong = (a > 0) && (!naFlag[a-1]) && strong[a-1];
	  bool rightStrong = (b + 1 < nw) && (!naFlag[b+1]) && strong[b+1];
	  bool keepDel = leftStrong && rightStrong && (runBp <= maxHomDel);
	  if (!keepDel) {
	    for(uint32_t k = a; k <= b; ++k) naFlag[k] = true;
	  }
	  i = b + 1;
	}

	// Flag non-unique windows
	bool uniqGate = (c.basecov && (!c.hasMapFile));
	if (uniqGate) {
	  for(uint32_t i = 0; i < nw; ++i) {
	    if (naFlag[i]) continue;
	    if ((wins[i].tcov > 0) && (wins[i].ucov <= c.uniqueToTotalCovRatio * wins[i].tcov)) naFlag[i] = true;
	  }
	}

	// Write windows
	if (!c.covfile.empty()) {
	  std::string chrn(hdr->target_name[refIndex]);
	  for(uint32_t i = 0; i < nw; ++i) {
	    double uniqFrac;
	    if (uniqGate) uniqFrac = (wins[i].tcov > 0) ? (wins[i].ucov / wins[i].tcov) : -1.0;
	    else uniqFrac = (wins[i].end > wins[i].start) ? ((double) wins[i].winlen / (double) (wins[i].end - wins[i].start)) : -1.0;
	    if (naFlag[i]) {
	      dataOut << chrn << "\t" << wins[i].start << "\t" << wins[i].end << "\t" << uniqFrac << "\tNA\tNA" << std::endl;
	    } else {
	      double cn = chrPloidy;
	      double logR = 0;
	      if (wins[i].expcov > 0) {
		cn = (c.expectedCN * wins[i].covsum / wins[i].expcov - chrCtrlPloidy * (1 - c.purity)) / c.purity;
		logR = std::log2((wins[i].covsum + 1.0) / (wins[i].expcov + 1.0));
	      }
	      dataOut << chrn << "\t" << wins[i].start << "\t" << wins[i].end << "\t" << uniqFrac << "\t" << logR << "\t" << cn << std::endl;
	    }
	  }
	}
      }
    }

    // Sort CNVs
    sort(cnvs.begin(), cnvs.end());

    // Genotype CNVs
    cnvVCF(c, cnvs);

    // clean-up
    fai_destroy(faiRef);
    if (faiMap != NULL) fai_destroy(faiMap);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
    if (!c.covfile.empty()) {
      dataOut.pop();
      dataOut.pop();
    }
    
    return 0;
  }

  
  int coral(int argc, char **argv) {
    CountDNAConfig c;
    std::string haploidChr;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome file")
      ("quality,q", boost::program_options::value<uint16_t>(&c.minQual)->default_value(10), "min. mapping quality")
      ("mappability,m", boost::program_options::value<boost::filesystem::path>(&c.mapFile), "input mappability map (optional)")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "BCF output file")
      ("covfile,c", boost::program_options::value<boost::filesystem::path>(&c.covfile), "gzipped coverage file")
      ;

    boost::program_options::options_description cnv("CNV calling");
    cnv.add_options()
      ("sdrd,x", boost::program_options::value<float>(&c.stringency)->default_value(2), "min. SD read-depth shift")
      ("cn-offset,t", boost::program_options::value<float>(&c.cn_offset)->default_value(0.1), "min. CN offset")
      ("cnv-size,z", boost::program_options::value<uint32_t>(&c.minCnvSize)->default_value(1000), "min. CNV size")
      ("svfile,l", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "delly SV file for breakpoint refinement")
      ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.genofile), "input VCF/BCF file for re-genotyping")
      ("segmentation,u", "copy-number segmentation")
      ;

    boost::program_options::options_description cancer("Ploidy/purity correction");
    cancer.add_options()
      ("ploidy,y", boost::program_options::value<float>(&c.ploidy)->default_value(2), "sample ploidy")
      ("purity,p", boost::program_options::value<float>(&c.purity)->default_value(1), "sample purity [0.1, 1]")
      ("ctrl-ploidy", boost::program_options::value<float>(&c.ctrlPloidy)->default_value(2), "control ploidy")
      ("haploid-chr", boost::program_options::value<std::string>(&haploidChr), "haploid chromosomes, e.g. chrX,chrY")
      ;
    
    boost::program_options::options_description window("Read-depth windows");
    window.add_options()
      ("window,w", boost::program_options::value<uint32_t>(&c.window_size)->default_value(0), "window size in bp (0: automatic, coverage-adaptive)")
      ("bed-intervals,b", boost::program_options::value<boost::filesystem::path>(&c.bedFile), "input BED file (targeted intervals)")
      ("mapq-uniq", boost::program_options::value<uint16_t>(&c.mapqUniq)->default_value(20), "min. MAPQ for a uniquely-placed read")
      ("target-reads", boost::program_options::value<uint32_t>(&c.targetReads)->default_value(150), "target reads/window")

      ("fraction-unique", boost::program_options::value<float>(&c.uniqueToTotalCovRatio)->default_value(0.8), "uniqueness filter [0,1]")
      ("basecov", "force base-level read-depth counting")
      ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamFile), "input BAM/CRAM file")
      ("fragment", boost::program_options::value<float>(&c.fragmentUnique)->default_value(0.97), "min. fragment uniqueness [0,1]")
      ("statsfile", boost::program_options::value<boost::filesystem::path>(&c.statsFile), "gzipped stats output file (optional)")
      ("window-offset", boost::program_options::value<uint32_t>(&c.window_offset)->default_value(0), "window offset (0: window size)")
      ("fraction-window", boost::program_options::value<float>(&c.fracWindow)->default_value(0.25), "min. callable window fraction [0,1]")
      ("scan-window", boost::program_options::value<uint32_t>(&c.scanWindow)->default_value(10000), "GC scanning window size")
      ("scan-regions", boost::program_options::value<boost::filesystem::path>(&c.scanFile), "GC scanning regions in BED format")
      ("mad-cutoff", boost::program_options::value<uint16_t>(&c.mad)->default_value(3), "median + 3 * mad count cutoff")
      ("percentile", boost::program_options::value<float>(&c.exclgc)->default_value(0.0005), "excl. extreme GC fraction")
      ("no-window-selection", "no scan window selection")
      ;

    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);

    // Set the visibility
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(cnv).add(cancer).add(window).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(cnv).add(cancer).add(window);

    // Parse command-line
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
      std::cerr << std::endl;
      std::cerr << "Usage: delly " << argv[0] << " [OPTIONS] -g <genome.fa> <aligned.bam>" << std::endl;
      std::cerr << visible_options << "\n";
      return 1;
    }

    // Mappability map
    if (vm.count("mappability")) c.hasMapFile = true;
    else c.hasMapFile = false;

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cerr << "delly ";
    for(int i=0; i<argc; ++i) { std::cerr << argv[i] << ' '; }
    std::cerr << std::endl;

    // Stats file
    if (vm.count("statsfile")) c.hasStatsFile = true;
    else c.hasStatsFile = false;

    // BED intervals
    if (vm.count("bed-intervals")) c.hasBedFile = true;
    else c.hasBedFile = false;

    // Scan regions
    if (vm.count("scan-regions")) c.hasScanFile = true;
    else c.hasScanFile = false;

    // Scan window selection
    if (vm.count("no-window-selection")) c.noScanWindowSelection = true;
    else c.noScanWindowSelection = false;

    // Adaptive windows
    c.adaptive = false;
    if (c.window_size == 0) {
      if (c.hasBedFile) c.window_size = 10000;
      else c.adaptive = true;
    }
    if (c.targetReads == 0) c.targetReads = 150;

    // Purity
    if (c.purity > 1) c.purity = 1;
    if (c.purity < 0.1) c.purity = 0.1;

    // Expected ploidy
    c.expectedCN = c.purity * c.ploidy + (1.0 - c.purity) * c.ctrlPloidy;
    
    // Segmentation
    if (vm.count("segmentation")) c.segmentation = true;
    else c.segmentation = false;

    // Window offset
    if ((c.window_offset == 0) || (c.window_offset > c.window_size)) c.window_offset = c.window_size;

    // Check input VCF file (CNV genotyping)
    if (vm.count("vcffile")) {
      if (!(boost::filesystem::exists(c.genofile) && boost::filesystem::is_regular_file(c.genofile) && boost::filesystem::file_size(c.genofile))) {
	std::cerr << "Input VCF/BCF file is missing: " << c.genofile.string() << std::endl;
	return 1;
      }
      htsFile* ifile = bcf_open(c.genofile.string().c_str(), "r");
      if (ifile == NULL) {
	std::cerr << "Fail to open file " << c.genofile.string() << std::endl;
	return 1;
      }
      bcf_hdr_t* hdr = bcf_hdr_read(ifile);
      if (hdr == NULL) {
	std::cerr << "Fail to open index file " << c.genofile.string() << std::endl;
	return 1;
      }
      bcf_hdr_destroy(hdr);
      bcf_close(ifile);
      c.hasGenoFile = true;
    } else c.hasGenoFile = false;

    // Check input VCF file (delly SV file)
    if (vm.count("svfile")) {
      if (!(boost::filesystem::exists(c.vcffile) && boost::filesystem::is_regular_file(c.vcffile) && boost::filesystem::file_size(c.vcffile))) {
	std::cerr << "Input VCF/BCF file is missing: " << c.vcffile.string() << std::endl;
	return 1;
      }
      htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
      if (ifile == NULL) {
	std::cerr << "Fail to open file " << c.vcffile.string() << std::endl;
	return 1;
      }
      bcf_hdr_t* hdr = bcf_hdr_read(ifile);
      if (hdr == NULL) {
	std::cerr << "Fail to open index file " << c.vcffile.string() << std::endl;
	return 1;
      }
      bcf_hdr_destroy(hdr);
      bcf_close(ifile);
      c.hasVcfFile = true;
    } else c.hasVcfFile = false;

    // Check outfile
    if (!vm.count("outfile")) c.outfile = "-";
    else {
      if (c.outfile.string() != "-") {
	if (!_outfileValid(c.outfile)) return 1;
      }
    }
    
    // Check bam file
    LibraryInfo li;
    bool pairedLib = false;
    if (!(boost::filesystem::exists(c.bamFile) && boost::filesystem::is_regular_file(c.bamFile) && boost::filesystem::file_size(c.bamFile))) {
      std::cerr << "Alignment file is missing: " << c.bamFile.string() << std::endl;
      return 1;
    } else {
      // Get scan regions
      typedef boost::icl::interval_set<uint32_t> TChrIntervals;
      typedef typename TChrIntervals::interval_type TIVal;
      typedef std::vector<TChrIntervals> TRegionsGenome;
      TRegionsGenome scanRegions;

      // Open BAM file
      samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
      if (samfile == NULL) {
	std::cerr << "Fail to open file " << c.bamFile.string() << std::endl;
	return 1;
      }
      hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
      if (idx == NULL) {
	if (bam_index_build(c.bamFile.string().c_str(), 0) != 0) {
	  std::cerr << "Fail to open index for " << c.bamFile.string() << std::endl;
	  return 1;
	}
      }
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.bamFile.string() << std::endl;
	return 1;
      }
      c.nchr = hdr->n_targets;
      c.minChrLen = setMinChrLen(hdr, 0.95);
      std::string sampleName = "unknown";
      getSMTag(std::string(hdr->text), c.bamFile.stem().string(), sampleName);
      c.sampleName = sampleName;

      // Check special chromosomes
      if (haploidChr.size()) _selectedRefIndices(c, hdr, haploidChr);

      // Check matching chromosome names
      faidx_t* faiRef = fai_load(c.genome.string().c_str());
      faidx_t* faiMap = NULL;
      if (c.hasMapFile) faiMap = fai_load(c.mapFile.string().c_str());
      uint32_t mapFound = 0;
      uint32_t refFound = 0;
      for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
	std::string tname(hdr->target_name[refIndex]);
	if (c.hasMapFile && faidx_has_seq(faiMap, tname.c_str())) ++mapFound;
	if (faidx_has_seq(faiRef, tname.c_str())) ++refFound;
	else {
	  std::cerr << "Warning: BAM chromosome " << tname << " not present in reference genome!" << std::endl;
	}
      }
      fai_destroy(faiRef);
      if (faiMap != NULL) fai_destroy(faiMap);
      if (c.hasMapFile && !mapFound) {
	std::cerr << "Mappability map chromosome naming disagrees with BAM file!" << std::endl;
	return 1;
      }
      if (!refFound) {
	std::cerr << "Reference genome chromosome naming disagrees with BAM file!" << std::endl;
	return 1;
      }

      // Estimate library params
      if (c.hasScanFile) {
	if (!_parseBedIntervals(c.scanFile.string(), c.hasScanFile, hdr, scanRegions)) {
	  std::cerr << "Warning: Couldn't parse BED intervals. Do the chromosome names match?" << std::endl;
	  return 1;
	}
      } else {
	scanRegions.resize(hdr->n_targets);
	for (int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
	  scanRegions[refIndex].insert(TIVal::right_open(0, hdr->target_len[refIndex]));
	}
      }
      typedef std::vector<LibraryInfo> TSampleLibrary;
      TSampleLibrary sampleLib(1, LibraryInfo());
      CountDNAConfigLib dellyConf;
      dellyConf.genome = c.genome;
      dellyConf.files.push_back(c.bamFile);
      dellyConf.madCutoff = 9;
      dellyConf.madNormalCutoff = c.mad;
      getLibraryParams(dellyConf, scanRegions, sampleLib);
      li = sampleLib[0];
      pairedLib = (li.median > 0);
      if (!li.median) {
	li.median = 250;
	li.mad = 15;
	li.minNormalISize = 0;
	li.maxNormalISize = 400;
      }
      c.meanisize = ((int32_t) (li.median / 2)) * 2 + 1;

      // Clean-up
      bam_hdr_destroy(hdr);
      hts_idx_destroy(idx);
      sam_close(samfile);
    }

    // Counting model
    c.basecov = false;
    if (!c.hasMapFile) {
      // No mappability map req. base-level coverage arrays
      c.basecov = true;
    } else {
      // Base-level coverage
      if (vm.count("basecov")) c.basecov = true;
      else c.basecov = ((!pairedLib) && (li.rs >= 500));  // Single-end reads > 500bp 
    }

    // GC bias estimation
    typedef std::pair<uint32_t, uint32_t> TGCBound;
    TGCBound gcbound;
    std::vector<GcBias> gcbias(c.meanisize + 1, GcBias());
    {
      // Scan genomic windows
      typedef std::vector<ScanWindow> TWindowCounts;
      typedef std::vector<TWindowCounts> TGenomicWindowCounts;
      TGenomicWindowCounts scanCounts(c.nchr, TWindowCounts());
      scan(c, li, scanCounts);

      // Check coverage
      {
	std::vector<uint32_t> sampleScanVec;
	for(uint32_t i = 0; i < scanCounts.size(); ++i) {
	  for(uint32_t j = 0; j < scanCounts[i].size(); ++j) {
	    sampleScanVec.push_back(scanCounts[i][j].cov);
	  }
	  if (sampleScanVec.size() > 1000000) break;
	}
	std::sort(sampleScanVec.begin(), sampleScanVec.end());
	if (sampleScanVec.empty()) {
	  std::cerr << "Not enough windows!" << std::endl;
	  return 1;
	}
	if (sampleScanVec[sampleScanVec.size()/2] < 5) {
	  std::cerr << "Please increase the window size. Coverage is too low!" << std::endl;
	  return 1;
	}
      }
    
      // Select stable windows
      selectWindows(c, scanCounts);

      // Estimate GC bias
      gcBias(c, scanCounts, li, gcbias, gcbound);

      // Statistics output
      if (c.hasStatsFile) {
	// Open stats file
	boost::iostreams::filtering_ostream statsOut;
	statsOut.push(boost::iostreams::gzip_compressor());
	statsOut.push(boost::iostreams::file_sink(c.statsFile.string(), std::ios_base::out | std::ios_base::binary));
	
	// Library Info
	statsOut << "LP\t" << li.rs << ',' << li.median << ',' << li.mad << ',' << li.minNormalISize << ',' << li.maxNormalISize << std::endl;
	
	// Scan window summry
	samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
	bam_hdr_t* hdr = sam_hdr_read(samfile);
	statsOut << "SW\tchrom\tstart\tend\tselected\tcoverage\tuniqcov" <<  std::endl;
	for(uint32_t refIndex = 0; refIndex < (uint32_t) hdr->n_targets; ++refIndex) {
	  for(uint32_t i = 0; i < scanCounts[refIndex].size(); ++i) {
	    statsOut << "SW\t" <<  hdr->target_name[refIndex] << '\t' << scanCounts[refIndex][i].start << '\t' << scanCounts[refIndex][i].end << '\t' << scanCounts[refIndex][i].select << '\t' << scanCounts[refIndex][i].cov << '\t' << scanCounts[refIndex][i].uniqcov << std::endl;
	  }
	}
	bam_hdr_destroy(hdr);
	sam_close(samfile);
	
	// GC bias summary
	statsOut << "GC\tgcsum\tsample\treference\tpercentileSample\tpercentileReference\tfractionSample\tfractionReference\tobsexp\tmeancoverage" << std::endl;
	for(uint32_t i = 0; i < gcbias.size(); ++i) statsOut << "GC\t" << i << "\t" << gcbias[i].sample << "\t" << gcbias[i].reference << "\t" << gcbias[i].percentileSample << "\t" << gcbias[i].percentileReference << "\t" << gcbias[i].fractionSample << "\t" << gcbias[i].fractionReference << "\t" << gcbias[i].obsexp << "\t" << gcbias[i].coverage << std::endl;
	statsOut << "BoundsGC\t" << gcbound.first << "," << gcbound.second << std::endl;
	statsOut.pop();
	statsOut.pop();
      }
    }

    // Coverage-aware window size
    if (c.adaptive) {
      // Mean expected coverage at CN2 over the callable GC range
      double covMean = 0;
      uint64_t refCnt = 0;
      for(uint32_t i = gcbound.first + 1; i < gcbound.second; ++i) {
	covMean += gcbias[i].coverage * (double) gcbias[i].reference;
	refCnt += gcbias[i].reference;
      }
      if (refCnt) covMean /= (double) refCnt;
      if (covMean <= 0) {
	// Fixed windows
	c.adaptive = false;
	c.window_size = 10000;
	c.window_offset = c.window_size;
      } else {
	double readLen = (li.rs > 0) ? (double) li.rs : (double) c.meanisize;
	double molPerBp = c.basecov ? (covMean / readLen) : covMean;
	if (molPerBp <= 0) molPerBp = 1e-9;
	double winBp = (double) c.targetReads / molPerBp;
	double minWin = std::max(100.0, 4.0 * readLen);
	double maxWin = 2000000.0;
	if (winBp < minWin) winBp = minWin;
	if (winBp > maxWin) winBp = maxWin;
	c.targetExpCov = covMean * winBp;
	double effReads = molPerBp * winBp;
	double covDepth = c.basecov ? covMean : (covMean * readLen);
	now = boost::posix_time::second_clock::local_time();
	std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Auto window size: " << (uint32_t) winBp << " bp, " << (uint32_t) effReads << " reads/window (" << (c.basecov ? "base-level" : "fragment") << ", coverage " << (uint32_t) (covDepth + 0.5) << "x)" << std::endl;
      }
    }

    // Count reads
    if (bamCount(c, li, gcbias, gcbound)) {
      std::cerr << "Read counting error!" << std::endl;
      return 1;
    }

    // Done
    now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;

    return 0;
  }

  
}

#endif
