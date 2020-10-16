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
#include "tegua.h"
#include "coral.h"

using namespace torali;


inline void
displayUsage() {
  std::cout << "Usage: delly <command> <arguments>" << std::endl;
  std::cout << std::endl;
  std::cout << "Short-read commands:" << std::endl;
  std::cout << "    call         discover and genotype structural variants" << std::endl;
  std::cout << "    merge        merge structural variants across VCF/BCF files and within a single VCF/BCF file" << std::endl;
  std::cout << "    filter       filter somatic or germline structural variants" << std::endl;
  std::cout << std::endl;
  std::cout << "Long-read commands:" << std::endl;
  std::cout << "    lr           long-read SV discovery (currently, only INS and DEL are supported)" << std::endl;
  std::cout << std::endl;
  std::cout << "Read-depth commands:" << std::endl;
  std::cout << "    rd           read-depth normalization" << std::endl;
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
      std::cout << " using Boost: v" << BOOST_VERSION / 100000 << "." << BOOST_VERSION / 100 % 1000 << "." << BOOST_VERSION % 100 << std::endl;
      std::cout << " using HTSlib: v" << hts_version() << std::endl;
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
      bsd();
      return 0;
    }
    else if ((std::string(argv[1]) == "call")) {
      return delly(argc-1,argv+1);
    }
    else if ((std::string(argv[1]) == "lr")) {
      return tegua(argc-1,argv+1);
    }
    else if ((std::string(argv[1]) == "rd")) {
      return coral(argc-1,argv+1);
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

