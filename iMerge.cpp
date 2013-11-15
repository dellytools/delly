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
#include <iomanip>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/bron_kerbosch_all_cliques.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


#include "memory_mapped_file.h"
#include "version.h"
#include "util.h"
#include "intervaltree.h"

using namespace torali;

struct Config {
  float threshold;
  std::string fileSuffix;
  std::vector<std::string> files;
};

struct ConfidenceInterval {
  unsigned int id;
  unsigned int chr;
  unsigned int start;
  unsigned int end;
  int sLow;
  int sHigh;
  int eLow;
  int eHigh;
};


template<typename TOutputStream>
struct CliqueIterator {
  CliqueIterator(TOutputStream& stream, std::vector<unsigned int>& clIdVec, unsigned int& clCount) : os(stream), cliqueId(clIdVec), cliqCount(clCount) {
  }

  template<typename TClique, typename TGraph>
  void clique(const TClique& c, const TGraph& g) {
    typename TClique::const_iterator end = c.end();
    os << "Clique" << cliqCount << ": ";
    for(typename TClique::const_iterator i=c.begin(); i!=end; ++i) {
      if (!cliqueId[g[*i].name]) {
	cliqueId[g[*i].name]=cliqCount;
	os << g[*i].name << " ";
      } else {
	// Node has been used in a larger clique already
	os << g[*i].name << "(used) ";
      }
    }
    os << std::endl;
    ++cliqCount;
  }

  TOutputStream& os;
  std::vector<unsigned int>& cliqueId;
  unsigned int& cliqCount;
};

struct VertexName {
  unsigned int name;
};


template<typename TFileAccess>
inline int
run(Config const& c)
{
  // Define a chromosome map
  unsigned numChrom = 0;
  typedef std::map<std::string, unsigned int> TChrMap;
  TChrMap chrMap;

  // Define the interval tree
  unsigned int idInterval = 0;
  typedef Interval<unsigned int> TInterval;
  typedef IntervalTree<TInterval> TIntervalTree;
  typedef std::vector<TIntervalTree*> TChrIntervalTrees;
  TChrIntervalTrees chrIntervals;

  // All confidence intervals
  typedef std::vector<ConfidenceInterval> TConfVector;
  TConfVector confVec;
	
  // Read input confidence intervals
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    std::cout << "Confidence intervals: " << c.files[file_c] << std::endl;
    
    typedef Record<std::string, unsigned int, unsigned int, int, int, int, int, void, void, void, void, void> TRecord;
    Memory_mapped_file map_file(c.files[file_c]);
    char buffer[Memory_mapped_file::MAX_LINE_LENGTH];
    while (map_file.left_bytes() > 0) {
      map_file.read_line(buffer);
      Tokenizer token(buffer, Memory_mapped_file::MAX_LINE_LENGTH);
      TRecord line;
      addF0(token, line);
      addF1(token, line);
      addF2(token, line);
      addF3(token, line);
      addF4(token, line);
      addF5(token, line);
      addF6(token, line);

      // Add the maximum possible interval to the interval tree of that chromosome
      if (chrMap.find(line.f0) == chrMap.end()) chrMap.insert(std::make_pair(line.f0, numChrom++));
      unsigned int chr = chrMap.find(line.f0)->second;
      if (chr == chrIntervals.size()) chrIntervals.push_back(new TIntervalTree());
      TIntervalTree* iTree = chrIntervals[chr];
      TInterval newInt;
      newInt.low = ((int) line.f1 + (int) line.f3 < 0) ? 0 : ((int) line.f1 + (int) line.f3);
      newInt.high = (line.f2 + line.f6);
      newInt.cargo = idInterval;
      iTree->insertInterval(newInt);

      // Add the confidence interval to the vector
      ConfidenceInterval cI;
      cI.id = idInterval;
      cI.chr = chr;
      cI.start = line.f1;
      cI.end = line.f2;
      cI.sLow = ((int) line.f1 + (int) line.f3 < 0) ? 0 : ((int) line.f1 + (int) line.f3);
      cI.sHigh = (line.f1 + line.f4);
      cI.eLow = ((int) line.f2 + (int) line.f5 < 0) ? 0 : ((int) line.f2 + (int) line.f5);
      cI.eHigh = (line.f2 + line.f6);
      confVec.push_back(cI);

      // Next interval
      ++idInterval;
    }
  }
  

  // Create the merge graph

  // Define an undirected graph g
  typedef boost::undirected_graph<VertexName> Graph;
  Graph g;
      
  // Define the vertex property map
  typedef boost::property_map<Graph, unsigned int VertexName::*>::type NameMap;
  NameMap nm(get(&VertexName::name,g));
      
  // Define the reverse map
  typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
  typedef std::vector<Vertex> TId2Vertex;
  TId2Vertex idVertex;
      
  // Insert a vertex for each confidence interval
  typename TConfVector::const_iterator confBeg = confVec.begin();
  typename TConfVector::const_iterator confEnd = confVec.end();
  for(;confBeg != confEnd; ++confBeg) {
    // Add vertex 
    Vertex u = add_vertex(g);
    nm[u] = confBeg->id;
    idVertex.push_back(u);
  }

  // Iterate through all confidence intervals
  confBeg = confVec.begin();
  for(;confBeg != confEnd; ++confBeg) {
    //std::cout << "-----------" << std::endl;
    //std::cout << confBeg->chr << "," << confBeg->start << "," << confBeg->end << std::endl;

    typedef std::vector<TInterval> TResultVec;
    TResultVec results;
    if (confBeg->chr < chrIntervals.size()) {
      TIntervalTree* iTree = chrIntervals[confBeg->chr];
      TInterval searchInt;
      searchInt.low = confBeg->sLow;
      searchInt.high = confBeg->eHigh;
      iTree->enumOverlapInterval(searchInt, results);
    }
    // Iterate all overlapping intervals
    TResultVec::const_iterator vecBeg = results.begin();
    TResultVec::const_iterator vecEnd = results.end();
    for(; vecBeg<vecEnd; ++vecBeg) {
      if (confVec[vecBeg->cargo].id > confBeg->id) {
	// How good do both intervals agree
	unsigned int startOverlap = 0;
	if ((confBeg->sHigh > confVec[vecBeg->cargo].sLow) && (confBeg->sLow < confVec[vecBeg->cargo].sHigh)) {
	  std::vector<int> mySort;
	  mySort.push_back(confBeg->sLow);
	  mySort.push_back(confBeg->sHigh);
	  mySort.push_back(confVec[vecBeg->cargo].sLow);
	  mySort.push_back(confVec[vecBeg->cargo].sHigh);
	  std::sort(mySort.begin(), mySort.end());
	  startOverlap = (((mySort[2] - mySort[1]) * 100) / (mySort[3] - mySort[0]));
	}
	unsigned int endOverlap = 0;
	if ((confBeg->eHigh > confVec[vecBeg->cargo].eLow) && (confBeg->eLow < confVec[vecBeg->cargo].eHigh)) {
	  std::vector<int> mySort;
	  mySort.push_back(confBeg->eLow);
	  mySort.push_back(confBeg->eHigh);
	  mySort.push_back(confVec[vecBeg->cargo].eLow);
	  mySort.push_back(confVec[vecBeg->cargo].eHigh);
	  std::sort(mySort.begin(), mySort.end());
	  endOverlap = (((mySort[2] - mySort[1]) * 100) / (mySort[3] - mySort[0]));
	}
	//std::cout << confBeg->sLow << ',' << confBeg->sHigh << ':' << confVec[vecBeg->cargo].sLow << ',' << confVec[vecBeg->cargo].sHigh << ',' << startOverlap << std::endl;
	//std::cout << confBeg->eLow << ',' << confBeg->eHigh << ':' << confVec[vecBeg->cargo].eLow << ',' << confVec[vecBeg->cargo].eHigh << ',' << endOverlap << std::endl;
	// The weight of that overlap is the minimum percentage overlap in either the start or end confidence intervals
	int weight = std::min(startOverlap, endOverlap); 

	if ((weight > 0) && (weight >= c.threshold * 100)) {
	  Vertex u = idVertex[confBeg->id];
	  Vertex v = idVertex[vecBeg->cargo];
	  //std::cout << confBeg->id << ',' << vecBeg->cargo << std::endl;
	  //std::cout << vecBeg->cargo << ',' << confBeg->id << std::endl;

	  // Add the edge
	  add_edge(u,v, g);
	}
      }
    }
  }

  // Compute maximal cliques using bron-kerbosch algorithm, id 0 is used for all singletons
  unsigned int clCount = 1;
  std::vector<unsigned int> cliqueId(num_vertices(g));
  std::fill(cliqueId.begin(), cliqueId.end(), 0);
  CliqueIterator<std::ostream> clIt(std::cout, cliqueId, clCount);
  bron_kerbosch_all_cliques(g, clIt);

  // Parse input file a second time
  unsigned int localIntervalId = 0;
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    std::cout << "Components of confidence intervals: " << c.files[file_c] << std::endl;

    // Read files one-by-one
    std::string outfile = c.files[file_c] + c.fileSuffix;
    std::ofstream anno(outfile.c_str());
    Memory_mapped_file interval_file(c.files[file_c]);
    char interval_buffer[Memory_mapped_file::MAX_LINE_LENGTH];
    while (interval_file.left_bytes() > 0) {
      interval_file.read_line(interval_buffer);
      anno << interval_buffer << "\t";
      if (cliqueId[localIntervalId]>0) {
	anno << cliqueId[localIntervalId] << std::endl;
      } else {
	anno << clCount++ << std::endl;
      }
      ++localIntervalId;
    }
    interval_file.close();
    anno.close();
  }

  // Free all allocated interval trees
  typename TChrIntervalTrees::iterator chrIIt = chrIntervals.begin();
  typename TChrIntervalTrees::iterator chrEIt = chrIntervals.end();
  for(;chrIIt!=chrEIt; ++chrIIt) {
    delete *chrIIt;
  }
  return 0;
}


int main(int argc, char **argv) {
  Config c;

  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("output-file-suffix,o", boost::program_options::value<std::string>(&c.fileSuffix)->default_value(".comp.out"), "output file suffix")
    ("threshold,t", boost::program_options::value<float>(&c.threshold)->default_value(0), "threshold for confidence interval agreement")
    ;

  // Define hidden options
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value< std::vector<std::string> >(&c.files), "input file")
    ("warranty,w", "show warranty")
    ("license,l", "show license")
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
  if ((vm.count("help")) || (!vm.count("input-file"))) { 
    printTitle("Merge intervals");
    if (vm.count("warranty")) {
      displayWarranty();
    } else if (vm.count("license")) {
      gplV3();
    } else {
      std::cout << "Usage: " << argv[0] << " [OPTIONS] <confidence intervals1.txt> <confidence intervals2.txt> ..." << std::endl;
      std::cout << visible_options << "\n"; 
    }
    return 1; 
  }

  return run<Memory_mapped_file>(c);
}
