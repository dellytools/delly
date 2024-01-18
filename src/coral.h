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
    bool adaptive;
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
    uint16_t minQual;
    uint16_t mad;
    float ploidy;
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
      dataOut.push(boost::iostreams::file_sink(c.covfile.c_str(), std::ios_base::out | std::ios_base::binary));
      dataOut << "chr\tstart\tend\t" << c.sampleName << "_mappable\t" << c.sampleName << "_counts\t" << c.sampleName << "_CN" << std::endl;
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
      for (uint32_t i = 0; i < svbp.size(); ++i) sort(svbp[i].begin(), svbp[i].end(), SortSVBreakpoint<SVBreakpoint>());
    }
    
    // Iterate chromosomes
    faidx_t* faiMap = fai_load(c.mapFile.string().c_str());
    faidx_t* faiRef = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      if ((!c.hasGenoFile) && (chrNoData(c, refIndex, idx))) continue;
      
      // Check presence in mappability map
      std::string tname(hdr->target_name[refIndex]);
      int32_t seqlen = faidx_seq_len(faiMap, tname.c_str());
      if (seqlen == - 1) continue;
      else seqlen = -1;
      char* seq = faidx_fetch_seq(faiMap, tname.c_str(), 0, faidx_seq_len(faiMap, tname.c_str()), &seqlen);

      // Check presence in reference
      seqlen = faidx_seq_len(faiRef, tname.c_str());
      if (seqlen == - 1) continue;
      else seqlen = -1;
      char* ref = faidx_fetch_seq(faiRef, tname.c_str(), 0, faidx_seq_len(faiRef, tname.c_str()), &seqlen);

      // Get GC and Mappability
      std::vector<uint16_t> uniqContent(hdr->target_len[refIndex], 0);
      std::vector<uint16_t> gcContent(hdr->target_len[refIndex], 0);
      {
	// Mappability map
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet uniq(hdr->target_len[refIndex], false);
	for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	  if (seq[i] == 'C') uniq[i] = 1;
	}

	// GC map
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet gcref(hdr->target_len[refIndex], false);
	for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	  if ((ref[i] == 'c') || (ref[i] == 'C') || (ref[i] == 'g') || (ref[i] == 'G')) gcref[i] = 1;
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
	  uniqContent[pos] = usum;
	}
      }
      
      // Coverage track
      typedef uint16_t TCount;
      uint32_t maxCoverage = std::numeric_limits<TCount>::max();
      typedef std::vector<TCount> TCoverage;
      TCoverage cov(hdr->target_len[refIndex], 0);

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

	  int32_t midPoint = rec->core.pos + halfAlignmentLength(rec);
	  if (rec->core.flag & BAM_FPAIRED) {
	    // Clean-up the read store for identical alignment positions
	    if (rec->core.pos > lastAlignedPos) {
	      lastAlignedPosReads.clear();
	      lastAlignedPos = rec->core.pos;
	    }
	    
	    if ((rec->core.pos < rec->core.mpos) || ((rec->core.pos == rec->core.mpos) && (lastAlignedPosReads.find(hash_string(bam_get_qname(rec))) == lastAlignedPosReads.end()))) {
	      // First read
	      lastAlignedPosReads.insert(hash_string(bam_get_qname(rec)));
	      std::size_t hv = hash_pair(rec);
	      mateMap[hv] = true;
	      continue;
	    } else {
	      // Second read
	      std::size_t hv = hash_pair_mate(rec);
	      if ((mateMap.find(hv) == mateMap.end()) || (!mateMap[hv])) continue; // Mate discarded
	      mateMap[hv] = false;
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
	if (ref != NULL) free(ref);
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }

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
	if (c.adaptive) {
	  // Merge overlapping BED entries
	  TChrIntervals citv;
	  _mergeOverlappingBedEntries(bedRegions[refIndex], citv);

	  // Tile merged intervals
	  double covsum = 0;
	  double expcov = 0;
	  double obsexp = 0;
	  uint32_t winlen = 0;
	  uint32_t start = 0;
	  bool endOfWindow = true;
	  typename TChrIntervals::iterator it = citv.begin();
	  if (it != citv.end()) start = it->first;
	  while(endOfWindow) {
	    endOfWindow = false;
	    for(it = citv.begin(); ((it != citv.end()) && (!endOfWindow)); ++it) {
	      if ((it->first < it->second) && (it->second <= hdr->target_len[refIndex])) {
		if (start >= it->second) {
		  if (start == it->second) {
		    // Special case
		    typename TChrIntervals::iterator itNext = it;
		    ++itNext;
		    if (itNext != citv.end()) start = itNext->first;
		  }
		  continue;
		}
		for(uint32_t pos = it->first; ((pos < it->second) && (!endOfWindow)); ++pos) {
		  if (pos < start) continue;
		  if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
		    covsum += cov[pos];
		    obsexp += gcbias[gcContent[pos]].obsexp;
		    expcov += gcbias[gcContent[pos]].coverage;
		    ++winlen;
		    if (winlen == c.window_size) {
		      obsexp /= (double) winlen;
		      double count = ((double) covsum / obsexp ) * (double) c.window_size / (double) winlen;
		      double cn = c.ploidy;
		      if (expcov > 0) cn = c.ploidy * covsum / expcov;
		      if (!c.covfile.empty()) dataOut << std::string(hdr->target_name[refIndex]) << "\t" << start << "\t" << (pos + 1) << "\t" << winlen << "\t" << count << "\t" << cn << std::endl;
		      // reset
		      covsum = 0;
		      expcov = 0;
		      obsexp = 0;
		      winlen = 0;
		      if (c.window_offset == c.window_size) {
			// Move on
			start = pos + 1;
			endOfWindow = true;
		      } else {
			// Rewind
			for(typename TChrIntervals::iterator sit = citv.begin(); ((sit != citv.end()) && (!endOfWindow)); ++sit) {
			  if ((sit->first < sit->second) && (sit->second <= hdr->target_len[refIndex])) {
			    if (start >= sit->second) continue;
			    for(uint32_t k = sit->first; ((k < sit->second) && (!endOfWindow)); ++k) {
			      if (k < start) continue;
			      if ((gcContent[k] > gcbound.first) && (gcContent[k] < gcbound.second) && (uniqContent[k] >= c.fragmentUnique * c.meanisize)) {
				++winlen;
				if (winlen == c.window_offset) {
				  start = k + 1;
				  winlen = 0;
				  endOfWindow = true;
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	} else {
	  // Fixed Window Length
	  for(typename TChrIntervals::iterator it = bedRegions[refIndex].begin(); it != bedRegions[refIndex].end(); ++it) {
	    if ((it->first < it->second) && (it->second <= hdr->target_len[refIndex])) {
	      double covsum = 0;
	      double expcov = 0;
	      double obsexp = 0;
	      uint32_t winlen = 0;
	      for(uint32_t pos = it->first; pos < it->second; ++pos) {
		if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
		  covsum += cov[pos];
		  obsexp += gcbias[gcContent[pos]].obsexp;
		  expcov += gcbias[gcContent[pos]].coverage;
		  ++winlen;
		}
	      }
	      if (winlen >= c.fracWindow * (it->second - it->first)) {
		obsexp /= (double) winlen;
		double count = ((double) covsum / obsexp ) * (double) (it->second - it->first) / (double) winlen;
		double cn = c.ploidy;
		if (expcov > 0) cn = c.ploidy * covsum / expcov;
		if (!c.covfile.empty()) dataOut << std::string(hdr->target_name[refIndex]) << "\t" << it->first << "\t" << it->second << "\t" << winlen << "\t" << count << "\t" << cn << std::endl;
	      } else {
		if (!c.covfile.empty()) dataOut << std::string(hdr->target_name[refIndex]) << "\t" << it->first << "\t" << it->second << "\tNA\tNA\tNA" << std::endl;
	      }
	    }
	  }
	}
      } else {
	// Genome-wide
	if (c.adaptive) {
	  double covsum = 0;
	  double expcov = 0;
	  double obsexp = 0;
	  uint32_t winlen = 0;
	  uint32_t start = 0;
	  uint32_t pos = 0;
	  while(pos < hdr->target_len[refIndex]) {
	    if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	      covsum += cov[pos];
	      obsexp += gcbias[gcContent[pos]].obsexp;
	      expcov += gcbias[gcContent[pos]].coverage;
	      ++winlen;
	      if (winlen == c.window_size) {
		obsexp /= (double) winlen;
		double count = ((double) covsum / obsexp ) * (double) c.window_size / (double) winlen;
		double cn = c.ploidy;
		if (expcov > 0) cn = c.ploidy * covsum / expcov;
		if (!c.covfile.empty()) dataOut << std::string(hdr->target_name[refIndex]) << "\t" << start << "\t" << (pos + 1) << "\t" << winlen << "\t" << count << "\t" << cn << std::endl;
		// reset
		covsum = 0;
		expcov = 0;
		obsexp = 0;
		winlen = 0;
		if (c.window_offset == c.window_size) {
		  // Move on
		  start = pos + 1;
		} else {
		  // Rewind
		  for(uint32_t k = start; k < hdr->target_len[refIndex]; ++k) {
		    if ((gcContent[k] > gcbound.first) && (gcContent[k] < gcbound.second) && (uniqContent[k] >= c.fragmentUnique * c.meanisize)) {
		      ++winlen;
		      if (winlen == c.window_offset) {
			start = k + 1;
			pos = k;
			winlen = 0;
			break;
		      }
		    }
		  }
		}
	      }
	    }
	    ++pos;
	  }
	} else {
	  // Fixed windows (genomic tiling)
	  for(uint32_t start = 0; start < hdr->target_len[refIndex]; start = start + c.window_offset) {
	    if (start + c.window_size < hdr->target_len[refIndex]) {
	      double covsum = 0;
	      double expcov = 0;
	      double obsexp = 0;
	      uint32_t winlen = 0;
	      for(uint32_t pos = start; pos < start + c.window_size; ++pos) {
		if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
		  covsum += cov[pos];
		  obsexp += gcbias[gcContent[pos]].obsexp;
		  expcov += gcbias[gcContent[pos]].coverage;
		  ++winlen;
		}
	      }
	      if (winlen >= c.fracWindow * c.window_size) {
		obsexp /= (double) winlen;
		double count = ((double) covsum / obsexp ) * (double) c.window_size / (double) winlen;
		double cn = c.ploidy;
		if (expcov > 0) cn = c.ploidy * covsum / expcov;
		if (!c.covfile.empty()) dataOut << std::string(hdr->target_name[refIndex]) << "\t" << start << "\t" << (start + c.window_size) << "\t" << winlen << "\t" << count << "\t" << cn << std::endl;
	      }
	    }
	  }
	}
      }
    }

    // Sort CNVs
    sort(cnvs.begin(), cnvs.end(), SortCNVs<CNV>());

    // Genotype CNVs
    cnvVCF(c, cnvs);

    // clean-up
    fai_destroy(faiRef);
    fai_destroy(faiMap);
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

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome file")
      ("quality,q", boost::program_options::value<uint16_t>(&c.minQual)->default_value(10), "min. mapping quality")
      ("mappability,m", boost::program_options::value<boost::filesystem::path>(&c.mapFile), "input mappability map")
      ("ploidy,y", boost::program_options::value<float>(&c.ploidy)->default_value(2), "baseline ploidy")
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
    
    boost::program_options::options_description window("Read-depth windows");
    window.add_options()
      ("window-size,i", boost::program_options::value<uint32_t>(&c.window_size)->default_value(10000), "window size")
      ("window-offset,j", boost::program_options::value<uint32_t>(&c.window_offset)->default_value(10000), "window offset")
      ("bed-intervals,b", boost::program_options::value<boost::filesystem::path>(&c.bedFile), "input BED file")
      ("fraction-window,k", boost::program_options::value<float>(&c.fracWindow)->default_value(0.25), "min. callable window fraction [0,1]")
      ("adaptive-windowing,a", "use mappable bases for window size")
      ;

    boost::program_options::options_description gcopt("GC fragment normalization");
    gcopt.add_options()
      ("scan-window,w", boost::program_options::value<uint32_t>(&c.scanWindow)->default_value(10000), "scanning window size")
      ("fraction-unique,f", boost::program_options::value<float>(&c.uniqueToTotalCovRatio)->default_value(0.8), "uniqueness filter for scan windows [0,1]")
      ("scan-regions,r", boost::program_options::value<boost::filesystem::path>(&c.scanFile), "scanning regions in BED format")
      ("mad-cutoff,d", boost::program_options::value<uint16_t>(&c.mad)->default_value(3), "median + 3 * mad count cutoff")
      ("percentile,p", boost::program_options::value<float>(&c.exclgc)->default_value(0.0005), "excl. extreme GC fraction")
      ("no-window-selection,n", "no scan window selection")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamFile), "input bam file")
      ("fragment,e", boost::program_options::value<float>(&c.fragmentUnique)->default_value(0.97), "min. fragment uniqueness [0,1]")
      ("statsfile,s", boost::program_options::value<boost::filesystem::path>(&c.statsFile), "gzipped stats output file (optional)")
      ;

    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);

    // Set the visibility
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(cnv).add(window).add(gcopt).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(cnv).add(window).add(gcopt);

    // Parse command-line
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome")) || (!vm.count("mappability"))) {
      std::cerr << std::endl;
      std::cerr << "Usage: delly " << argv[0] << " [OPTIONS] -g <genome.fa> -m <genome.map> <aligned.bam>" << std::endl;
      std::cerr << visible_options << "\n";
      return 1;
    }

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

    // Adaptive windowing
    if (vm.count("adaptive-windowing")) c.adaptive = true;
    else c.adaptive = false;

    // Segmentation
    if (vm.count("segmentation")) c.segmentation = true;
    else c.segmentation = false;

    // Check window size
    if (c.window_offset > c.window_size) c.window_offset = c.window_size;
    if (c.window_size == 0) c.window_size = 1;
    if (c.window_offset == 0) c.window_offset = 1;

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

      // Check matching chromosome names
      faidx_t* faiRef = fai_load(c.genome.string().c_str());
      faidx_t* faiMap = fai_load(c.mapFile.string().c_str());
      uint32_t mapFound = 0;
      uint32_t refFound = 0;
      for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
	std::string tname(hdr->target_name[refIndex]);
	if (faidx_has_seq(faiMap, tname.c_str())) ++mapFound;
	if (faidx_has_seq(faiRef, tname.c_str())) ++refFound;
	else {
	  std::cerr << "Warning: BAM chromosome " << tname << " not present in reference genome!" << std::endl;
	}
      }
      fai_destroy(faiRef);
      fai_destroy(faiMap);
      if (!mapFound) {
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
	statsOut.push(boost::iostreams::file_sink(c.statsFile.string().c_str(), std::ios_base::out | std::ios_base::binary));
	
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
