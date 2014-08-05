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
#include <iomanip>
#include <fstream>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/progress.hpp>
#include "api/BamReader.h"
#include "api/BamIndex.h"

#include "version.h"
#include "intervaltree.h"
#include "util.h"

using namespace torali;

// Tags
struct TagMean {};
struct TagMin {};
struct TagMax {};

template<typename T>
struct Overlap {
  typedef T TTag;
};

template<typename T>
struct Overhang {
  typedef T TTag;
};


struct TOverlap {
  double overlapA;
  double overlapB;
  double overlapMean;
  std::string cargo;

  TOverlap(double a, double b, double mean, std::string& c) : overlapA(a), overlapB(b), overlapMean(mean), cargo(c) {}
  TOverlap(double a, double b, double mean, std::string const& c) : overlapA(a), overlapB(b), overlapMean(mean), cargo(c) {}
};



struct Config {
  unsigned int maxOverlaps;
  float cutoff;
  float threshold;
  std::string mode;
  std::string method;
  std::string fileSuffix;
  std::vector<std::string> anno_int;
  std::vector<unsigned int> fd;
  std::vector<bool> noIntervals;
  std::vector<boost::filesystem::path> files;
  std::vector<bool> queryNoIntervals;
};



template<typename TOverlapRecord>
struct SortRecordsMean : public std::binary_function<TOverlapRecord, TOverlapRecord, bool>
{
  inline bool operator()(TOverlapRecord const& rec1, TOverlapRecord const& rec2) {
    return ((rec1.overlapMean != rec2.overlapMean) ? (rec1.overlapMean > rec2.overlapMean) : (std::abs((int) (rec1.overlapA - rec1.overlapB) * 100) < std::abs((int) (rec2.overlapA - rec2.overlapB) * 100)));
  }
};

template<typename TOverlapRecord>
struct SortRecordsMin : public std::binary_function<TOverlapRecord, TOverlapRecord, bool>
{
  inline bool operator()(TOverlapRecord const& rec1, TOverlapRecord const& rec2) {
    return (std::min(rec1.overlapA, rec1.overlapB) != std::min(rec2.overlapA, rec2.overlapB)) ? (std::min(rec1.overlapA, rec1.overlapB) > std::min(rec2.overlapA, rec2.overlapB)) : (std::max(rec1.overlapA, rec1.overlapB) > std::max(rec2.overlapA, rec2.overlapB));
  }
};

template<typename TOverlapRecord>
struct SortRecordsMax : public std::binary_function<TOverlapRecord, TOverlapRecord, bool>
{
  inline bool operator()(TOverlapRecord const& rec1, TOverlapRecord const& rec2) {
    return (std::max(rec1.overlapA, rec1.overlapB) != std::max(rec2.overlapA, rec2.overlapB)) ? (std::max(rec1.overlapA, rec1.overlapB) > std::max(rec2.overlapA, rec2.overlapB)) : (std::min(rec1.overlapA, rec1.overlapB) > std::min(rec2.overlapA, rec2.overlapB));
  }
};


template<typename TQual, typename TThres>
inline bool
aboveThres(TQual const& qual, TThres th, TagMean const&) {
  return (qual.overlapMean >= th);
}

template<typename TQual, typename TThres>
inline bool
aboveThres(TQual const& qual, TThres th, TagMin const&) {
  return (std::min(qual.overlapA, qual.overlapB) >= th);
}

template<typename TQual, typename TThres>
inline bool
aboveThres(TQual const& qual, TThres th, TagMax const&) {
  return (std::max(qual.overlapA, qual.overlapB) >= th);
}


template<typename TQual>
inline void
sortMethod(TQual& qual, TagMean const&) {
  std::sort(qual.begin(), qual.end(), SortRecordsMean<TOverlap>());
}

template<typename TQual>
inline void
sortMethod(TQual& qual, TagMin const&) {
  std::sort(qual.begin(), qual.end(), SortRecordsMin<TOverlap>());
}

template<typename TQual>
inline void
sortMethod(TQual& qual, TagMax const&) {
  std::sort(qual.begin(), qual.end(), SortRecordsMax<TOverlap>());
}



// Overhang method
template<typename TFile, typename TConfig, typename TInterval, typename TResultVec, typename TTag>
inline void
sortResults(TFile& anno, TConfig const& c, TInterval const& searchInt, TResultVec const& results, Overhang<TTag> const&) {
  typename TResultVec::const_iterator vecBeg = results.begin();
  typename TResultVec::const_iterator vecEnd = results.end();
  typedef std::vector<TOverlap>  TOverlapQuality;
  TOverlapQuality qual;
  for(; vecBeg<vecEnd; ++vecBeg) {
    if (searchInt.low < vecBeg->low) {
      if (searchInt.high < vecBeg->high) qual.push_back(TOverlap(-1.0 * (vecBeg->low - searchInt.low), -1.0 * (vecBeg->high - searchInt.high), ((vecBeg->low - searchInt.low) + (vecBeg->high - searchInt.high)) / -2.0, vecBeg->cargo));
      else qual.push_back(TOverlap(-1.0 * (vecBeg->low - searchInt.low), -1.0 * (searchInt.high - vecBeg->high), ((vecBeg->low - searchInt.low) + (searchInt.high - vecBeg->high)) / -2.0, vecBeg->cargo));
    } else {
      if (searchInt.high < vecBeg->high) qual.push_back(TOverlap(-1.0 * (searchInt.low - vecBeg->low), -1.0 * (vecBeg->high - searchInt.high), ((searchInt.low - vecBeg->low) + (vecBeg->high - searchInt.high)) / -2.0, vecBeg->cargo));
      else qual.push_back(TOverlap(-1.0 * (searchInt.low - vecBeg->low), -1.0 * (searchInt.high - vecBeg->high), ((searchInt.low - vecBeg->low) + (searchInt.high - vecBeg->high)) / -2.0, vecBeg->cargo));
    }
  }
  sortMethod(qual, TTag());
  
  // Print n-best overlaps
  unsigned int numPrint = qual.size();
  if ((c.maxOverlaps>0) && (c.maxOverlaps<numPrint)) numPrint = c.maxOverlaps;
  TOverlapQuality::const_iterator qualBeg = qual.begin();
  TOverlapQuality::const_iterator qualEnd = qual.begin()+numPrint;
  if (c.cutoff == -1) {
    anno << "\t";
    for(;qualBeg < qualEnd; ++qualBeg) {
      if (aboveThres(*qualBeg, -c.threshold, TTag())) {
	anno << qualBeg->cargo << "[" << -qualBeg->overlapA << ',' << -qualBeg->overlapB << ',' << -qualBeg->overlapMean << "];";
      }
    }
  } else {
    if ((qualBeg != qualEnd) && (aboveThres(*qualBeg, -c.cutoff, TTag()))) anno << "1";
    else anno << "0";
  }
}


// Overlap method
template<typename TFile, typename TConfig, typename TInterval, typename TResultVec, typename TTag>
inline void
sortResults(TFile& anno, TConfig const& c, TInterval const& searchInt, TResultVec const& results, Overlap<TTag> const&) {
  typename TResultVec::const_iterator vecBeg = results.begin();
  typename TResultVec::const_iterator vecEnd = results.end();
  std::vector<int> mySort;
  typedef std::vector<TOverlap>  TOverlapQuality;
  TOverlapQuality qual;
  for(; vecBeg<vecEnd; ++vecBeg) {
    mySort.push_back(searchInt.low);
    mySort.push_back(searchInt.high);
    mySort.push_back(vecBeg->low);
    mySort.push_back(vecBeg->high);
    std::sort(mySort.begin(), mySort.end());
    qual.push_back(TOverlap(((double) (mySort[2] - mySort[1] + 1) / (double) (searchInt.high - searchInt.low + 1)),((double) (mySort[2] - mySort[1] + 1) / (double) (vecBeg->high - vecBeg->low + 1)),((((double) (mySort[2] - mySort[1] + 1) / (double) (searchInt.high - searchInt.low + 1)) + ((double) (mySort[2] - mySort[1] + 1) / (double) (vecBeg->high - vecBeg->low + 1))) / 2.0), vecBeg->cargo));
    mySort.clear();
  }
  sortMethod(qual, TTag());
  
  // Print n-best overlaps
  unsigned int numPrint = qual.size();
  if ((c.maxOverlaps>0) && (c.maxOverlaps<numPrint)) numPrint = c.maxOverlaps;
  TOverlapQuality::const_iterator qualBeg = qual.begin();
  TOverlapQuality::const_iterator qualEnd = qual.begin()+numPrint;
  if (c.cutoff == -1) {
    anno << "\t";
    for(;qualBeg < qualEnd; ++qualBeg) {
      if (aboveThres(*qualBeg, c.threshold, TTag())) {
	anno << qualBeg->cargo << std::setprecision(2) << "[" << qualBeg->overlapA << ',' << qualBeg->overlapB << ',' << qualBeg->overlapMean << "];";
      }
    }
  } else {
    if ((qualBeg != qualEnd) && (aboveThres(*qualBeg, c.cutoff, TTag()))) anno << "1";
    else anno << "0";
  }
}


template<typename TConfig, typename TMethod>
inline int
run(TConfig const& c, TMethod const&) {
  // Create for each file a hash of interval trees (one for each chromosome)
  typedef Interval<std::string> TInterval;
  typedef IntervalTree<TInterval> TIntervalTree;
  typedef std::map<std::string, TIntervalTree*> TChrIntervalTrees;
  typedef std::vector<TChrIntervalTrees> TIntervalTreeCollection;
  TIntervalTreeCollection iTreeColl;
  
  // Scan all annotation intervals
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Building annotation interval trees" << std::endl;
  boost::progress_display show_progress(c.anno_int.size());
  for(unsigned int anno_c = 0; anno_c<c.anno_int.size(); ++anno_c) {
    ++show_progress;
    TChrIntervalTrees chrIntervals;

    if (boost::filesystem::exists(c.anno_int[anno_c]) && boost::filesystem::is_regular_file(c.anno_int[anno_c]) && boost::filesystem::file_size(c.anno_int[anno_c])) {
      std::ifstream annoFile(c.anno_int[anno_c].c_str());
      unsigned int line_counter = 1;
      while (annoFile.good()) {
	std::string annoLine;
	getline(annoFile, annoLine);
	if (annoLine[0]=='#') continue;
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep("\t,;");
	Tokenizer tokens(annoLine, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter!=tokens.end()) {
	  std::string chrName=*tokIter++;
	  int32_t intStart = boost::lexical_cast<int32_t>(*tokIter++);
	  int32_t intEnd = intStart;
	  std::string key=*tokIter++;
	  if (!c.noIntervals[anno_c]) intEnd = boost::lexical_cast<int32_t>(key);
	  if (c.fd[anno_c] != 0) {
	    if (c.fd[anno_c] == 1) key = chrName;
	    else if (c.fd[anno_c] == 2) { 
	      std::stringstream out;
	      out << intStart;
	      key = out.str();
	    } else if (c.fd[anno_c] > 3) {
	      unsigned int k = 4;
	      for (;k < c.fd[anno_c];++k) tokIter++;
	      key=*tokIter++;
	    }
	  } else {
	    std::stringstream s; 
	    s << line_counter++;
	    key = s.str();
	  }
	  TChrIntervalTrees::iterator findTree = chrIntervals.find( chrName );
	  if (findTree == chrIntervals.end()) findTree = chrIntervals.insert(std::make_pair(chrName, new TIntervalTree())).first;
	  TIntervalTree* iTree = findTree->second;
	  TInterval newInt;
	  newInt.low = intStart;
	  newInt.high = intEnd;
	  newInt.cargo = key;
	  iTree->insertInterval(newInt);
	}
      }
      annoFile.close();
      iTreeColl.push_back(chrIntervals);
    }
  }

  // Scan all query intervals
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Searching query intervals" << std::endl;
  boost::progress_display show_progress2(c.files.size());
  for(unsigned int files_c = 0; files_c<c.files.size(); ++files_c) {
    ++show_progress2;

    if (!boost::filesystem::exists(c.files[files_c]) || !boost::filesystem::is_regular_file(c.files[files_c]) || !boost::filesystem::file_size(c.files[files_c])) continue;
    std::string outfile = c.files[files_c].string() + c.fileSuffix;
    std::ofstream anno(outfile.c_str());

    std::ifstream queryFile(c.files[files_c].string().c_str());
    while (queryFile.good()) {
      std::string queryLine;
      getline(queryFile, queryLine);
      if (queryLine[0]=='#') continue;
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep("\t,;");
      Tokenizer tokens(queryLine, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      if (tokIter!=tokens.end()) {
	std::string chrName=*tokIter++;
	int32_t intStart = boost::lexical_cast<int32_t>(*tokIter++);
	int32_t intEnd = intStart;
	if (!c.queryNoIntervals[files_c]) intEnd = boost::lexical_cast<int32_t>(*tokIter++);
	TInterval searchInt;
	searchInt.low = intStart;
	searchInt.high = intEnd;
	anno << queryLine;
	if (c.cutoff != -1) anno << "\t";

	// Find overlapping intervals
	for(unsigned int anno_c = 0; anno_c<c.anno_int.size(); ++anno_c) {
	  typedef std::vector<TInterval> TResultVec;
	  TResultVec results;
	  TChrIntervalTrees& chrIntervals = iTreeColl[anno_c];	
	  TChrIntervalTrees::iterator findTree = chrIntervals.find( chrName );
	  if (findTree != chrIntervals.end()) {
	    TIntervalTree* iTree = findTree->second;
	    iTree->enumOverlapInterval(searchInt, results);
	  }
	
	  // Sort results
	  sortResults(anno, c, searchInt, results, TMethod());
	}
	anno << std::endl;
      }
    }
    anno.close();
  }

  // Free all allocated interval trees
  for(unsigned int anno_c = 0; anno_c<c.anno_int.size(); ++anno_c) {
    TChrIntervalTrees& chrIntervals = iTreeColl[anno_c];	
    TChrIntervalTrees::iterator treeBeg = chrIntervals.begin();
    TChrIntervalTrees::iterator treeEnd = chrIntervals.end();
    for(;treeBeg != treeEnd; ++treeBeg) {
      delete treeBeg->second;
    }
  }

  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
  return 0;
}


int main(int argc, char **argv) {
  // Config object
  Config c;

  // Define required options
  boost::program_options::options_description required("Required options");
  required.add_options()
    ("annotation-intervals,a", boost::program_options::value<std::vector<std::string> >(&c.anno_int), "annotation intervals")
    ("field-identifier,f", boost::program_options::value<std::vector<unsigned int> >(&c.fd), "field identifier index")
    ("interval-or-point,i", boost::program_options::value<std::vector<bool> >(&c.noIntervals), "0: intervals, 1: point")
    ("query-interval-or-point,q", boost::program_options::value<std::vector<bool> >(&c.queryNoIntervals), "0: intervals, 1: point")
    ;

  // Define generic options
  boost::program_options::options_description thres("Threshold options");
  thres.add_options()
    ("method,e", boost::program_options::value<std::string>(&c.method)->default_value("overlap"), "threshold/cutoff is used for [overlap|overhang]")
    ("mode,d", boost::program_options::value<std::string>(&c.mode)->default_value("mean"), "threshold/cutoff apply mode [min|max|mean]")
    ("threshold,t", boost::program_options::value<float>(&c.threshold)->default_value(0), "threshold for overlap quality")
    ("cutoff,c", boost::program_options::value<float>(&c.cutoff)->default_value(-1), "presence/absence cutoff for c>= 0")
    ;
  
  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("max-overlaps,m", boost::program_options::value<unsigned int>(&c.maxOverlaps)->default_value(10), "max. overlaps to print, 0 means print all")
    ("output-suffix,o", boost::program_options::value<std::string>(&c.fileSuffix)->default_value(".aout"), "output file suffix")
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
  cmdline_options.add(required).add(thres).add(generic).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(required).add(thres).add(generic);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);


  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("annotation-intervals")) || (!vm.count("input-file"))) {
    printTitle("Interval overlap");
    if (vm.count("warranty")) {
      displayWarranty();
    } else if (vm.count("license")) {
      gplV3();
    } else {
      std::cout << "Usage: " << argv[0] << " [OPTIONS] <query_intervals1.txt> <query_intervals2.txt> ..." << std::endl;
      std::cout << visible_options << "\n"; 
    }
    return 1; 
  }

  // Resize field identifiers and noIntervals vectors if necessary
  if (c.anno_int.size() != c.fd.size()) c.fd.resize(c.anno_int.size());
  if (c.anno_int.size() != c.noIntervals.size()) c.noIntervals.resize(c.anno_int.size());
  if (c.files.size() != c.queryNoIntervals.size()) c.queryNoIntervals.resize(c.files.size());

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Which method
  if (!c.method.compare("overhang")) {
    if (!c.mode.compare("max")) return run(c, Overhang<TagMin>());
    else if (!c.mode.compare("min")) return run(c, Overhang<TagMax>());
    else return run(c, Overhang<TagMean>());
  } else {
    if (!c.mode.compare("max")) return run(c, Overlap<TagMax>());
    else if (!c.mode.compare("min")) return run(c, Overlap<TagMin>());
    else return run(c, Overlap<TagMean>());
  }
}
