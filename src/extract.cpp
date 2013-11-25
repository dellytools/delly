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

#include "version.h"
#include "sam.h"
#include "fasta_reader.h"
#include "extract.h"

#define MAX_CHROM_SIZE 250000000

using namespace torali;

int main(int argc, char **argv) {
  // Define command-line options
  ExtractConfig c;
  c.breaks = false;
  c.closed = false;
  c.revComp = true;

  // Define required options
  boost::program_options::options_description required("Required options");
  required.add_options()
    ("genome,g", boost::program_options::value<std::string>(&c.genome), "genome multi-fasta file")
    ;

  boost::program_options::options_description region("Region options");
  region.add_options()
    ("chr,c", boost::program_options::value<std::string>(&c.chr)->default_value("chrX"), "chromosome identifier")
    ("start,s", boost::program_options::value<unsigned int>(&c.start)->default_value(0), "region start (inclusive, 0-based)")
    ("end,e", boost::program_options::value<unsigned int>(&c.end)->default_value(1), "region end (exclusive, 0-based)")
    ;

  boost::program_options::options_description interval("Interval options");
  interval.add_options()
    ("interval,i", boost::program_options::value<std::string>(&c.intervals), "file with intervals")
    ("field-identifier,f", boost::program_options::value<unsigned int>(&c.field_identifier)->default_value(0), "field identifier index")
    ("breaks,b", "parse list of breaks")
    ;

  // Define generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("linesize,z", boost::program_options::value<unsigned int>(&c.linesize)->default_value(60), "line size")
    ("closed-interval,o", "use closed interval [a,b], instead of [a,b[")
    ("no-reverse-complement,n", "do not reverse complement intervals [b,a] with b>a")
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
  if (vm.count("closed-interval")) c.closed = true;
  if (vm.count("no-reverse-complement")) c.revComp = false;
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

  // Show command line arguments
  std::cerr << "Region extraction" << std::endl;
  std::cerr << "Genome file: " << c.genome << std::endl;
  std::cerr << "Interval file: " << c.intervals << std::endl;
  std::cerr << "Chromosome identifier: " << c.chr << std::endl;
  std::cerr << "Start position: " << c.start << std::endl;
  std::cerr << "End position: " << c.end << std::endl;
  std::cerr << "Line size: " << c.linesize << std::endl;
  std::cerr << "Parse breaks: " << c.breaks << std::endl;
  std::cerr << "Closed-intervals: " << c.closed << std::endl;
  std::cerr << "Reverse-complement: " << c.revComp << std::endl;

  // Run region extractor
  if (c.breaks) runExtract<std::string>(c);
  else runExtract<void>(c);

  return 0;
}
