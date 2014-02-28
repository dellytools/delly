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
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "api/BamReader.h"
#include "api/BamIndex.h"


#include "version.h"
#include "extract.h"

#define MAX_CHROM_SIZE 250000000

using namespace torali;

int main(int argc, char **argv) {
  // Define command-line options
  ExtractConfig c;

  // Define required options
  boost::program_options::options_description required("Required options");
  required.add_options()
    ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome multi-fasta file")
    ;

  boost::program_options::options_description region("Region options");
  region.add_options()
    ("chr,c", boost::program_options::value<std::string>(&c.chr)->default_value("chrX"), "chromosome identifier")
    ("start,s", boost::program_options::value<unsigned int>(&c.start)->default_value(0), "region start (inclusive, 1-based)")
    ("end,e", boost::program_options::value<unsigned int>(&c.end)->default_value(1), "region end (exclusive, 1-based)")
    ;

  boost::program_options::options_description interval("Interval options");
  interval.add_options()
    ("interval,i", boost::program_options::value<boost::filesystem::path>(&c.intervals), "file with intervals")
    ("breaks,b", "parse list of breaks")
    ;

  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("linesize,z", boost::program_options::value<unsigned int>(&c.linesize)->default_value(60), "line size")
    ("closed-interval,d", "closed intervals [start, end]")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.fasta"), "FASTA output file")
    ;

  // Define hidden options
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("warranty,w", "show warranty")
    ("license,l", "show license")
    ;

  // Set the visibility
  boost::program_options::options_description cmdline_options;
  cmdline_options.add(required).add(generic).add(region).add(interval).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(required).add(generic).add(region).add(interval);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
  boost::program_options::notify(vm);

  // Check command line arguments
  if (vm.count("breaks")) c.breaks = true;
  else c.breaks = false;
  if (vm.count("closed-interval")) c.closed = true;
  else c.closed = false;

  // Help message
  if ((vm.count("help")) || (!vm.count("genome")) ) {
    printTitle("Region extraction");
    if (vm.count("warranty")) {
      displayWarranty();
    } else if (vm.count("license")) {
      gplV3();
    } else {
      std::cout << "Usage: " << argv[0] << " [OPTIONS]" << std::endl;
      std::cout << visible_options << "\n"; 
    }
    return 1; 
  }

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Run region extractor
  return runExtract(c);
}
