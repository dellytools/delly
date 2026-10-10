#include "../src/delly.h"
#include <stdexcept>

namespace read_cap_test {
  struct Config : torali::Config {};

  // Exercise the real BAM collection and read-count bookkeeping independently
  // of whether a synthetic sequence can produce a biological consensus.
  inline int msa(Config const&, std::unordered_set<std::string> const& reads, std::string& consensus) {
    consensus = *reads.begin();
    return reads.size();
  }

  inline bool alignConsensus(Config const&, bam_hdr_t*, char const*, char const*, torali::StructuralVariantRecord&) {
    return true;
  }
}

void require(bool condition, char const* message) {
  if (!condition) throw std::runtime_error(message);
}

std::string sequence(unsigned value) {
  std::string result(150, 'A');
  for (unsigned i = 0; i < 12; ++i) {
    result[i] = "ACGT"[value % 4];
    value /= 4;
  }
  return result;
}

int main(int argc, char** argv) {
  require(argc == 2, "Provide a temporary directory");
  std::string directory(argv[1]);
  read_cap_test::Config config;
  config.genome = directory + "/ref.fa";
  config.maxThreads = 2;
  config.minMapQual = 1;
  {
    std::ofstream reference(config.genome.string());
    reference << ">chr1\n" << std::string(1000, 'A') << "\n>chr2\n" << std::string(1000, 'A') << '\n';
  }
  require(fai_build(config.genome.string().c_str()) == 0, "Cannot index reference");

  typedef boost::unordered_map<std::pair<int32_t, std::size_t>, int32_t> ReadSVs;
  std::vector<ReadSVs> readSVs(2);
  std::string headerText("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n");
  for (unsigned sample = 0; sample < 2; ++sample) {
    std::string file = directory + "/sample" + std::to_string(sample) + ".bam";
    config.files.push_back(file);
    samFile* output = sam_open(file.c_str(), "wb");
    bam_hdr_t* header = sam_hdr_parse(headerText.size(), headerText.c_str());
    require(output && header, "Cannot create BAM");
    require(sam_hdr_write(output, header) >= 0, "Cannot write BAM header");
    bam1_t* record = bam_init1();
    for (unsigned chromosome = 0; chromosome < 2; ++chromosome) {
      for (unsigned variant = chromosome == 0 ? 0 : 1; variant < 2; ++variant) {
        unsigned count = variant == 0 ? 30 : 25;
        for (unsigned read = 0; read < count; ++read) {
          unsigned id = sample * 1000 + chromosome * 100 + read + variant * 200;
          // Duplicate alignments with different names must not consume extra
          // slots when their sequence is already present.
          for (unsigned duplicate = 0; duplicate < 2; ++duplicate) {
            std::string row = "r" + std::to_string(id) + "_" + std::to_string(duplicate) +
              "\t0\tchr" + std::to_string(chromosome + 1) + "\t" + (variant == 0 ? "11" : "201") +
              "\t60\t150M\t*\t0\t0\t" + sequence(id) + "\t*";
            kstring_t line = {0, 0, NULL};
            kputs(row.c_str(), &line);
            require(sam_parse1(&line, header, record) >= 0, "Cannot parse synthetic read");
            readSVs[chromosome][std::make_pair(record->core.pos, torali::hash_sr(record))] = variant;
            require(sam_write1(output, header, record) >= 0, "Cannot write read");
            free(line.s);
          }
        }
      }
    }
    bam_destroy1(record);
    bam_hdr_destroy(header);
    require(sam_close(output) == 0, "Cannot close BAM");
    require(sam_index_build(file.c_str(), 0) == 0, "Cannot index BAM");
  }
  std::vector<boost::icl::interval_set<uint32_t>> regions(2);
  for (auto& region : regions)
    region.add(boost::icl::interval<uint32_t>::right_open(0, 1000));

  for (bool reverseSamples : {false, true}) {
    if (reverseSamples) std::reverse(config.files.begin(), config.files.end());
    for (unsigned limit : {20, 0, 1, 2, 21, 200}) {
      config.maxReadPerSV = limit;
      std::vector<torali::StructuralVariantRecord> variants(2);
      variants[0].id = 0;
      variants[0].svt = 2;
      variants[0].chr = variants[0].chr2 = 0;
      variants[0].svStart = 10;
      variants[0].svEnd = 100;
      variants[1].id = 1;
      variants[1].svt = 5;
      variants[1].chr = 1;
      variants[1].chr2 = 0;
      variants[1].svStart = variants[1].svEnd = 200;
      torali::assembleSplitReads(config, regions, readSVs, variants);
      int32_t deletionCount = std::min(limit, 60u);
      int32_t translocationCount = std::min(limit, 100u);
      require(variants[0].srSupport == (deletionCount > 1 ? deletionCount : 0), "Deletion read cap failed");
      require(variants[1].srSupport == (translocationCount > 1 ? translocationCount : 0), "Translocation read cap failed");
    }
  }
  std::cout << "Read caps passed for duplicates, both breakends and both sample orders\n";
}
