// Copyright (c) 2020 Robert Vaser

#include <getopt.h>

#include <cmath>
#include <iomanip>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/sequence.hpp"
#include "ram/minimizer_engine.hpp"

#include "spoa/spoa.hpp"

std::atomic<std::uint32_t> biosoup::Sequence::num_objects{0};
std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace {

static struct option options[] = {
  {"algorithm", required_argument, nullptr, 'l'},
  {"result", required_argument, nullptr, 'r'},
  {"dot", required_argument, nullptr, 'd'},
  {"strand-ambiguous", no_argument, nullptr, 's'},
  {"threads", required_argument, nullptr, 't'},
  {"version", no_argument, nullptr, 'v'},
  {"help", no_argument, nullptr, 'h'},
  {nullptr, 0, nullptr, 0}
};

std::unique_ptr<bioparser::Parser<biosoup::Sequence>> CreateParser(
    const std::string& path) {
  auto is_suffix = [] (const std::string& str, const std::string& suff) {
    return str.size() < suff.size() ? false :
        str.compare(str.size() - suff.size(), suff.size(), suff) == 0;
  };

  if (is_suffix(path, ".fasta") || is_suffix(path, ".fasta.gz") ||
      is_suffix(path, ".fna")   || is_suffix(path, ".fna.gz") ||
      is_suffix(path, ".faa")   || is_suffix(path, ".faa.gz") ||
      is_suffix(path, ".fa")    || is_suffix(path, ".fa.gz")) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastaParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }
  if (is_suffix(path, ".fastq") || is_suffix(path, ".fastq.gz") ||
      is_suffix(path, ".fq")    || is_suffix(path, ".fq.gz")) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastqParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  std::cerr << "[spoa::CreateParser] error: file " << path
            << " has unsupported format extension (valid extensions: .fasta, "
            << ".fasta.gz, .fna, .fna.gz, .faa, .faa.gz, .fa, .fa.gz, .fastq, "
            << ".fastq.gz, .fq, .fq.gz)"
            << std::endl;
  return nullptr;
}

void Help() {
  std::cout <<
      "usage: spoa [options ...] <backbone> <sequences> [<sequences>]\n"
      "\n"
      "  # default output is stdout\n"
      "  <backbone>/<sequences>\n"
      "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
      "\n"
      "  options:\n"
      "    -m <int>\n"
      "      default: 5\n"
      "      score for matching bases\n"
      "    -n <int>\n"
      "      default: -4\n"
      "      score for mismatching bases\n"
      "    -g <int>\n"
      "      default: -8\n"
      "      gap opening penalty (must be non-positive)\n"
      "    -e <int>\n"
      "      default: -8\n"
      "      gap extension penalty (must be non-positive)\n"
      "    -q <int>\n"
      "      default: -8\n"
      "      gap opening penalty of the second affine function\n"
      "      (must be non-positive)\n"
      "    -c <int>\n"
      "      default: -8\n"
      "      gap extension penalty of the second affine function\n"
      "      (must be non-positive)\n"
      "    -l, --algorithm <int>\n"
      "      default: 0\n"
      "      alignment mode:\n"
      "        0 - local (Smith-Waterman)\n"
      "        1 - global (Needleman-Wunsch)\n"
      "        2 - semi-global\n"
      "    -r, --result <int> (option can be used multiple times)\n"
      "      default: 1\n"
      "      result mode:\n"
      "        0 - consensus (FASTA)\n"
      "        1 - multiple sequence alignment (FASTA)\n"
      "        2 - 0 & 1 (FASTA)\n"
      "        3 - partial order graph (GFA)\n"
      "        4 - 0 & 3 (GFA)\n"
      "    -d, --dot <file>\n"
      "      output file for the partial order graph in DOT format\n"
      "    -s, --strand-ambiguous\n"
      "      for each sequence pick the strand with the better alignment\n"
      "    -t, --threads <int>\n"
      "      default: 1\n"
      "      number of threads\n"
      "    --version\n"
      "      prints the version number\n"
      "    -h, --help\n"
      "      prints the usage\n"
      "\n"
      "  gap mode:\n"
      "    linear if g >= e\n"
      "    affine if g <= q or e >= c\n"
      "    convex otherwise (default)\n";
}

void PrintGfa(
    const spoa::Graph& graph,
    const std::vector<std::string>& headers,
    const std::vector<bool>& is_reversed,
    bool include_consensus = false) {
  if (headers.size() < graph.sequences().size()) {
    std::cerr << "[spoa::PrintGfa] error: missing header(s)" << std::endl;
    return;
  }
  if (!is_reversed.empty() && is_reversed.size() < graph.sequences().size()) {
    std::cerr << "[spoa::PringGfa] error: missing reversion flag(s)" << std::endl;  // NOLINT
    return;
  }

  std::vector<bool> is_consensus_node(graph.nodes().size(), false);
  for (const auto& it : graph.consensus()) {
    is_consensus_node[it->id] = true;
  }

  std::cout << "H\tVN:Z:1.0" << std::endl;
  for (const auto& it : graph.nodes()) {
    std::cout << "S\t"
              << it->id + 1 << "\t"
              << static_cast<char>(graph.decoder(it->code));
    if (is_consensus_node[it->id]) {
      std::cout << "\tic:Z:true";
    }
    std::cout << std::endl;

    for (const auto& jt : it->outedges) {
      std::cout << "L\t"
                << it->id + 1 << "\t"
                << "+\t"
                << jt->head->id + 1 << "\t"
                << "+\t"
                << "OM\t"
                << "ew:f:" << jt->weight;
      if (is_consensus_node[it->id] &&
          is_consensus_node[jt->head->id]) {
        std::cout << "\tic:Z:true";
      }
      std::cout << std::endl;
    }
  }

  for (std::uint32_t i = 0; i < graph.sequences().size(); ++i) {
    std::cout << "P\t" << headers[i] << "\t";

    std::vector<std::uint32_t> path;
    auto curr = graph.sequences()[i];
    while (true) {
      path.emplace_back(curr->id + 1);
      if (!(curr = curr->Successor(i))) {
        break;
      }
    }

    bool ir = !is_reversed.empty() && is_reversed[i];
    if (ir) {
      std::reverse(path.begin(), path.end());
    }
    for (std::uint32_t j = 0; j < path.size(); ++j) {
      if (j != 0) {
        std::cout << ",";
      }
      std::cout << path[j] << (ir ? "-" : "+");
    }
    std::cout << "\t*" << std::endl;
  }

  if (include_consensus) {
    std::cout << "P\tConsensus\t";
    for (std::uint32_t i = 0; i < graph.consensus().size(); ++i) {
      if (i != 0) {
        std::cout << ",";
      }
      std::cout << graph.consensus()[i]->id + 1 << "+";
    }
    std::cout << "\t*" << std::endl;
  }
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> input_paths;

  std::int8_t m = 5;
  std::int8_t n = -4;
  std::int8_t g = -8;
  std::int8_t e = -8;
  std::int8_t q = -8;
  std::int8_t c = -8;

  std::uint8_t algorithm = 0;
  std::vector<std::uint8_t> results = { 1 };
  std::string dot_path{};
  bool is_strand_ambiguous = false;

  std::size_t num_threads = 1;

  std::string optstr = "m:n:g:e:q:c:l:r:d:st:h";
  int opt;
  while ((opt = getopt_long(argc, argv, optstr.c_str(), options, nullptr)) != -1) {  // NOLINT
    switch (opt) {
      case 'm': m = atoi(optarg); break;
      case 'n': n = atoi(optarg); break;
      case 'g': g = atoi(optarg); break;
      case 'e': e = atoi(optarg); break;
      case 'q': q = atoi(optarg); break;
      case 'c': c = atoi(optarg); break;
      case 'l': algorithm = atoi(optarg); break;
      case 'r': results.emplace_back(atoi(optarg)); break;
      case 'd': dot_path = optarg; break;
      case 's': is_strand_ambiguous = true; break;
      case 't': num_threads = atoi(optarg); break;
      case 'v': std::cout << VERSION << std::endl; return 0;
      case 'h': Help(); return 0;
      default: return 1;
    }
  }
  if (results.size() > 1) {
    results.erase(results.begin());
  }

  for (int32_t i = optind; i < argc; ++i) {
    input_paths.emplace_back(argv[i]);
  }

  if (input_paths.size() < 2) {
    std::cerr << "[spoa::] error: missing input file(s)!" << std::endl;
    Help();
    return 1;
  }

  std::vector<std::unique_ptr<biosoup::Sequence>> sequences;
  for (const auto& it : input_paths) {
    auto sparser = CreateParser(it);
    if (sparser == nullptr) {
      return 1;
    }
    auto chunk = sparser->Parse(-1);
    sequences.insert(
        sequences.end(),
        std::make_move_iterator(chunk.begin()),
        std::make_move_iterator(chunk.end()));
  }
  if (sequences.empty()) {
    std::cerr << "[spoa::] error: missing sequences!" << std::endl;
    return 1;
  }

  std::unique_ptr<spoa::AlignmentEngine> alignment_engine;
  try {
    alignment_engine = spoa::AlignmentEngine::Create(
        static_cast<spoa::AlignmentType>(algorithm), m, n, g, e, q, c);
  } catch(std::invalid_argument& exception) {
    std::cerr << exception.what() << std::endl;
    return 1;
  }

  {  // circular sequence preprocessing
    auto thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads);

    auto minimizer_engine = ram::MinimizerEngine(thread_pool);

    {
      std::vector<std::unique_ptr<biosoup::NucleicAcid>> compressed;
      compressed.emplace_back(std::unique_ptr<biosoup::NucleicAcid>(
          new biosoup::NucleicAcid("backbone", sequences.front()->data)));
      minimizer_engine.Minimize(compressed.begin(), compressed.end());
    }

    std::vector<std::future<std::vector<biosoup::Overlap>>> futures;
    for (std::uint32_t i = 1; i < sequences.size(); ++i) {
      futures.emplace_back(thread_pool->Submit(
          [&] (std::uint32_t j) -> std::vector<biosoup::Overlap> {
            auto na = std::unique_ptr<biosoup::NucleicAcid>(
                new biosoup::NucleicAcid("", sequences[j]->data));
            na->id = sequences[j]->id;
            return minimizer_engine.Map(na, false, false);
          },
          i));
    }

    std::size_t num_valid = 0;
    std::size_t num_longer = 0;
    std::size_t num_split = 0;
    std::size_t num_artefacts = 0;

    decltype(sequences) rotated;
    rotated.emplace_back(std::unique_ptr<biosoup::Sequence>(
        new biosoup::Sequence(*sequences.front())));

    for (auto& it : futures) {
      auto overlaps = it.get();
      if (overlaps.empty()) {  // junk
        continue;
      }

      if (overlaps.size() == 1) {  // one contiguous alignment
        rotated.emplace_back(std::unique_ptr<biosoup::Sequence>(
            new biosoup::Sequence(*sequences[overlaps.front().lhs_id])));
        rotated.back()->id = overlaps.front().lhs_id;
        rotated.back()->name += " "
            "[" +
            std::to_string(0) + ", " +
            std::to_string(rotated.back()->data.size()) +
            "]";
        if (overlaps.front().strand == 0) {
          rotated.back()->ReverseAndComplement();
        }
        ++num_valid;
        continue;
      }

      if (overlaps.size() == 3) {  // longer read
        std::sort(overlaps.begin(), overlaps.end(),
            [] (const biosoup::Overlap& lhs,
                const biosoup::Overlap& rhs) -> bool {
              return (lhs.rhs_end - lhs.rhs_begin) < (rhs.rhs_end - rhs.rhs_begin);  // NOLINT
            });

        auto o = overlaps.back();
        overlaps.pop_back();

        if (o.rhs_end - o.rhs_begin > 0.97 * sequences.front()->data.size()) {
          std::string name = sequences[o.lhs_id]->name + " ";
          name +=
              "[" +
              std::to_string(o.lhs_begin) + ", " +
              std::to_string(o.lhs_end) +
              "]";
          std::string data = sequences[o.lhs_id]->data.substr(
              o.lhs_begin,
              o.lhs_end - o.lhs_begin);

          rotated.emplace_back(std::unique_ptr<biosoup::Sequence>(
              new biosoup::Sequence(name, data)));
          rotated.back()->id = o.lhs_id;
          if (o.strand == 0) {
            rotated.back()->ReverseAndComplement();
          }
          ++num_longer;
        } else {
          // weird read
          ++num_artefacts;
          continue;
        }
      }

      if (overlaps.size() > 2 ||
          overlaps.front().strand != overlaps.back().strand) {
        // sequencing artefact
        ++num_artefacts;
        continue;
      }

      std::sort(overlaps.begin(), overlaps.end(),
          [] (const biosoup::Overlap& lhs,
              const biosoup::Overlap& rhs) -> bool {
            return lhs.rhs_begin < rhs.rhs_begin;
          });

      std::string first_name = sequences[overlaps.front().lhs_id]->name + " ";
      first_name +=
        "[" +
        std::to_string(overlaps.front().lhs_begin) + ", " +
        std::to_string(overlaps.front().lhs_end) +
        "]";
      std::string first_data = sequences[overlaps.front().lhs_id]->data.substr(
          overlaps.front().lhs_begin,
          overlaps.front().lhs_end);
      if (overlaps.front().strand == 0) {
        auto na = biosoup::NucleicAcid{"", first_data};
        na.ReverseAndComplement();
        first_data = na.InflateData();
      }

      std::string second_name = " ";
      second_name +=
          "[" +
          std::to_string(overlaps.back().lhs_begin) + ", " +
          std::to_string(overlaps.back().lhs_end) +
          "]";
      std::string second_data = sequences[overlaps.back().lhs_id]->data.substr(
          overlaps.back().lhs_begin,
          overlaps.back().lhs_end);
      if (overlaps.back().strand == 0) {
        auto na = biosoup::NucleicAcid{"", second_data};
        na.ReverseAndComplement();
        second_data = na.InflateData();
      }

      std::int32_t indel_size =
          overlaps.back().rhs_begin -
          overlaps.front().rhs_end;

      // add first
      rotated.emplace_back(std::unique_ptr<biosoup::Sequence>(
          new biosoup::Sequence(first_name, first_data)));
      rotated.back()->id = overlaps.front().lhs_id;

      if (indel_size > 0 && indel_size < 0.03 * sequences.front()->data.size()) {  // NOLINT
        // concatenate second
        rotated.back()->name += second_name;
        rotated.back()->data += second_data;
        ++num_valid;
      } else {
        // add second
        second_name = sequences[overlaps.back().lhs_id]->name + second_name;
        rotated.emplace_back(std::unique_ptr<biosoup::Sequence>(
            new biosoup::Sequence(second_name, second_data)));
        rotated.back()->id = overlaps.back().lhs_id;
        ++num_split;
      }
    }

    std::size_t md = log10(std::max(sequences.size(), rotated.size())) + 1;

    std::cerr << "Num sequences = " << std::setw(md) << sequences.size() << std::endl;  // NOLINT
    std::cerr << "Num valid     = " << std::setw(md) << num_valid << std::endl;
    std::cerr << "Num longer    = " << std::setw(md) << num_longer << std::endl;
    std::cerr << "Num split     = " << std::setw(md) << num_split << std::endl;
    std::cerr << "Num artefacts = " << std::setw(md) << num_artefacts << std::endl;  // NOLINT
    std::cerr << std::string(16 + md, '-') << std::endl;
    std::cerr << "Num modified  = " << std::setw(md) << rotated.size() << std::endl;  // NOLINT

    sequences.swap(rotated);

    std::sort(sequences.begin(), sequences.end(),
        [] (const std::unique_ptr<biosoup::Sequence>& lhs,
            const std::unique_ptr<biosoup::Sequence>& rhs) -> bool {
          return lhs->data.size() > rhs->data.size();
        });
  }

  if (sequences.empty()) {
    std::cerr << "[spoa::] error: filtered all sequences!" << std::endl;
    return 1;
  }

  std::size_t max_sequence_len = 0;
  for (const auto& it : sequences) {
    max_sequence_len = std::max(max_sequence_len, it->data.size());
  }
  try {
    alignment_engine->Prealloc(max_sequence_len, 4);
  } catch (std::invalid_argument& exception) {
    std::cerr << exception.what() << std::endl;
    return 1;
  }

  spoa::Graph graph{};
  std::vector<bool> is_reversed;
  for (const auto& it : sequences) {
    std::int32_t score = 0;
    spoa::Alignment alignment;
    try {
      alignment = alignment_engine->Align(it->data, graph, &score);
    } catch (std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return 1;
    }

    if (is_strand_ambiguous) {
      it->ReverseAndComplement();
      std::int32_t score_rev = 0;
      spoa::Alignment alignment_rev;
      try {
        alignment_rev = alignment_engine->Align(it->data, graph, &score_rev);
      } catch (std::invalid_argument& exception) {
        std::cerr << exception.what() << std::endl;
        return 1;
      }
      if (score >= score_rev) {
        it->ReverseAndComplement();
        is_reversed.push_back(false);
      } else {
        alignment = alignment_rev;
        is_reversed.push_back(true);
      }
    }

    try {
      if (it->quality.empty()) {
        graph.AddAlignment(alignment, it->data);
      } else {
        graph.AddAlignment(alignment, it->data, it->quality);
      }
    } catch(std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return 1;
    }
  }

  for (const auto& it : results) {
    switch (it) {
      case 0: {
        auto consensus = graph.GenerateConsensus();
        std::cout << ">Consensus LN:i:" << consensus.size() << std::endl
                  << consensus << std::endl;
        break;
      }
      case 1:
      case 2: {
        auto msa = graph.GenerateMultipleSequenceAlignment(it == 2);
        for (std::uint32_t i = 0; i < msa.size(); ++i) {
          std::string name = i < sequences.size() ? sequences[i]->name : "Consensus";  // NOLINT
          std::cout << ">" << name << std::endl
                    << msa[i] << std::endl;
        }
        break;
      }
      case 3:
      case 4: {
        std::vector<std::string> headers;
        for (const auto& it : sequences) {
          headers.emplace_back(it->name);
        }
        graph.GenerateConsensus();
        PrintGfa(graph, headers, is_reversed, it == 4);
        break;
      }
      default:
        break;
    }
  }

  graph.PrintDot(dot_path);

  return 0;
}
