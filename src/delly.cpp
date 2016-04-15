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

#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>

#define BOOST_DISABLE_ASSERTS

#ifdef OPENMP
#include <omp.h>
#endif

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif

#include "version.h"
#include "delly.h"
#include "filter.h"
#include "merge.h"

using namespace torali;


inline void
displayUsage() {
  std::cout << "Usage: delly <command> <arguments>" << std::endl;
  std::cout << std::endl;
  std::cout << "Commands:" << std::endl;
  std::cout << std::endl;
  std::cout << "    call         discover and genotype structural variants" << std::endl;
  std::cout << "    merge        merge structural variants across VCF/BCF files and within a single VCF/BCF file" << std::endl;
  std::cout << "    filter       filter somatic or germline structural variants" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
}

int main(int argc, char **argv) {
    if (argc < 2) { 
      printTitle("Delly");
      displayUsage();
      return 0;
    }

    if ((std::string(argv[1]) == "version") || (std::string(argv[1]) == "--version") || (std::string(argv[1]) == "--version-only") || (std::string(argv[1]) == "-v")) {
      std::cout << "Delly version: v" << dellyVersionNumber << std::endl;
      return 0;
    }
    else if ((std::string(argv[1]) == "help") || (std::string(argv[1]) == "--help") || (std::string(argv[1]) == "-h") || (std::string(argv[1]) == "-?")) {
      printTitle("Delly");
      displayUsage();
      return 0;
    }
    else if ((std::string(argv[1]) == "warranty") || (std::string(argv[1]) == "--warranty") || (std::string(argv[1]) == "-w")) {
      displayWarranty();
      return 0;
    }
    else if ((std::string(argv[1]) == "license") || (std::string(argv[1]) == "--license") || (std::string(argv[1]) == "-l")) {
      gplV3();
      return 0;
    }
    else if ((std::string(argv[1]) == "call")) {
      return delly(argc-1,argv+1);
    }
    else if ((std::string(argv[1]) == "filter")) {
      return filter(argc-1,argv+1);
    }
    else if ((std::string(argv[1]) == "merge")) {
      return merge(argc-1,argv+1);
    }

    std::cerr << "Unrecognized command " << std::string(argv[1]) << std::endl;
    return 1;
}

