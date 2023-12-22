#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>

#define BOOST_DISABLE_ASSERTS
#define BOOST_UUID_RANDOM_PROVIDER_FORCE_POSIX

#ifdef OPENMP
#include <omp.h>
#endif

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif

#include "version.h"
#include "delly.h"
#include "filter.h"
#include "classify.h"
#include "merge.h"
#include "tegua.h"
#include "coral.h"
#include "asmode.h"
#include "dpe.h"
#include "pangenome.h"
#include "markdup.h"
#include "compvcf.h"
#include "chimera.h"

using namespace torali;


inline void
displayUsage() {
  std::cerr << "Usage: delly <command> <arguments>" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Short-read SV calling:" << std::endl;
  std::cerr << "    call         discover and genotype structural variants" << std::endl;
  std::cerr << "    merge        merge structural variants across VCF/BCF files and within a single VCF/BCF file" << std::endl;
  std::cerr << "    filter       filter somatic or germline structural variants" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Long-read SV calling:" << std::endl;
  std::cerr << "    lr           long-read SV discovery" << std::endl;
  std::cerr << std::endl;
  //std::cerr << "Pan-genome based SV calling (work-in-progress):" << std::endl;
  //std::cerr << "    pg           pan-genome SV discovery" << std::endl;
  //std::cerr << std::endl;
  //std::cerr << "Assembly-based SV calling (work-in-progress):" << std::endl;
  //std::cerr << "    asm          assembly SV site discovery" << std::endl;
  //std::cerr << std::endl;
  std::cerr << "Copy-number variant calling:" << std::endl;
  std::cerr << "    cnv          discover and genotype copy-number variants" << std::endl;
  std::cerr << "    classify     classify somatic or germline copy-number variants" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Multi-sample VCF operations:" << std::endl;
  std::cerr << "    markdup      mark duplicate SV sites based on SV allele and GT concordance" << std::endl;
  std::cerr << "    compvcf      compare multi-sample VCF file to a ground truth VCF file" << std::endl;
  //std::cerr << "Deprecated:" << std::endl;
  //std::cerr << "    dpe          double paired-end signatures" << std::endl;
  //std::cerr << "    chimera      ONT chimera flagging" << std::endl;
  //std::cerr << std::endl;
  std::cerr << std::endl;
  std::cerr << std::endl;
}

int main(int argc, char **argv) {
    if (argc < 2) { 
      printTitle("Delly");
      displayUsage();
      return 0;
    }

    if ((std::string(argv[1]) == "version") || (std::string(argv[1]) == "--version") || (std::string(argv[1]) == "--version-only") || (std::string(argv[1]) == "-v")) {
      std::cerr << "Delly version: v" << dellyVersionNumber << std::endl;
      std::cerr << " using Boost: v" << BOOST_VERSION / 100000 << "." << BOOST_VERSION / 100 % 1000 << "." << BOOST_VERSION % 100 << std::endl;
      std::cerr << " using HTSlib: v" << hts_version() << std::endl;
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
    else if ((std::string(argv[1]) == "asm")) {
      return asmode(argc-1,argv+1);
    }
    else if ((std::string(argv[1]) == "markdup")) {
      return markdup(argc-1,argv+1);
    }
    else if ((std::string(argv[1]) == "compvcf")) {
      return compvcf(argc-1,argv+1);
    }
    else if ((std::string(argv[1]) == "pg")) {
      return pg(argc-1,argv+1);
    }
    else if ((std::string(argv[1]) == "dpe")) {
      return dpe(argc-1,argv+1);
    }
    else if ((std::string(argv[1]) == "chimera")) {
      return chimera(argc-1,argv+1);
    }
    else if ((std::string(argv[1]) == "cnv")) {
      return coral(argc-1,argv+1);
    }
    else if ((std::string(argv[1]) == "classify")) {
      return classify(argc-1,argv+1);
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

