#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "sequence.hpp"

#include "spoa/spoa.hpp"
#include "bioparser/bioparser.hpp"

static const char* version = "v2.0.1";

static struct option options[] = {
    {"match", required_argument, 0, 'm'},
    {"mismatch", required_argument, 0, 'x'},
    {"gap-open", required_argument, 0, 'g'},
    {"gap-extend", required_argument, 0, 'e'},
    {"algorithm", required_argument, 0, 'l'},
    {"result", required_argument, 0, 'r'},
    {"dot", required_argument, 0, 'd'},
    {"version", no_argument, 0, 'v'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

void help();

int main(int argc, char** argv) {

    int8_t match = 5;
    int8_t mismatch = -4;
    int8_t gap_open = -8;
    int8_t gap_extend = -6;

    uint8_t algorithm = 0;
    uint8_t result = 0;

    std::string dot_path = "";

    char opt;
    while ((opt = getopt_long(argc, argv, "m:x:g:e:l:r:d:h", options, nullptr)) != -1) {
        switch (opt) {
            case 'm':
                match = atoi(optarg);
                break;
            case 'x':
                mismatch = atoi(optarg);
                break;
            case 'g':
                gap_open = atoi(optarg);
                break;
            case 'e':
                gap_extend = atoi(optarg);
                break;
            case 'l':
                algorithm = atoi(optarg);
                break;
            case 'r':
                result = atoi(optarg);
                break;
            case 'd':
                dot_path = optarg;
                break;
            case 'v':
                printf("%s\n", version);
                exit(0);
            case 'h':
                help();
                exit(0);
            default:
                exit(1);
        }
    }

    if (optind >= argc) {
        fprintf(stderr, "[spoa::] error: missing input file!\n");
        help();
        exit(1);
    }

    std::string sequences_path = argv[optind];

    auto is_suffix = [](const std::string& src, const std::string& suffix) -> bool {
        if (src.size() < suffix.size()) {
            return false;
        }
        return src.compare(src.size() - suffix.size(), suffix.size(), suffix) == 0;
    };

    std::unique_ptr<bioparser::Parser<spoa::Sequence>> sparser = nullptr;

    if (is_suffix(sequences_path, ".fasta") || is_suffix(sequences_path, ".fa") ||
        is_suffix(sequences_path, ".fasta.gz") || is_suffix(sequences_path, ".fa.gz")) {
        sparser = bioparser::createParser<bioparser::FastaParser, spoa::Sequence>(
            sequences_path);
    } else if (is_suffix(sequences_path, ".fastq") || is_suffix(sequences_path, ".fq") ||
        is_suffix(sequences_path, ".fastq.gz") || is_suffix(sequences_path, ".fq.gz")) {
        sparser = bioparser::createParser<bioparser::FastqParser, spoa::Sequence>(
            sequences_path);
    } else {
        fprintf(stderr, "[spoa::] error: "
            "file %s has unsupported format extension (valid extensions: "
            ".fasta, .fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)!\n",
            sequences_path.c_str());
        exit(1);
    }

    auto alignment_engine = spoa::createAlignmentEngine(
        static_cast<spoa::AlignmentType>(algorithm), match, mismatch, gap_open,
        gap_extend);

    auto graph = spoa::createGraph();

    std::vector<std::unique_ptr<spoa::Sequence>> sequences;
    sparser->parse_objects(sequences, -1);

    size_t max_sequence_size = 0;
    for (const auto& it: sequences) {
        max_sequence_size = std::max(max_sequence_size, it->data().size());
    }
    alignment_engine->prealloc(max_sequence_size, 4);

    for (const auto& it: sequences) {
        auto alignment = (*alignment_engine)(it->data(), graph);
        graph->add_alignment(alignment, it->data(), it->quality());
    }

    if (result == 0 || result == 2) {
        std::string consensus = graph->generate_consensus();
        fprintf(stdout, "Consensus (%zu)\n", consensus.size());
        fprintf(stdout, "%s\n", consensus.c_str());
    }

    if (result == 1 || result == 2) {
        std::vector<std::string> msa;
        graph->generate_multiple_sequence_alignment(msa);
        fprintf(stdout, "Multiple sequence alignment\n");
        for (const auto& it: msa) {
            fprintf(stdout, "%s\n", it.c_str());
        }
    }

    graph->print_dot(dot_path);

    return 0;
}

void help() {
    printf(
        "usage: spoa [options ...] <sequences>\n"
        "\n"
        "    <sequences>\n"
        "        input file in FASTA/FASTQ format (can be compressed with gzip)\n"
        "        containing sequences\n"
        "\n"
        "    options:\n"
        "        -m, --match <int>\n"
        "            default: 5\n"
        "            score for matching bases\n"
        "        -x, --mismatch <int>\n"
        "            default: -4\n"
        "            score for mismatching bases\n"
        "        -g, --gap-open <int>\n"
        "            default: -8\n"
        "            gap opening penalty (must be negative)\n"
        "        -e, --gap-extend <int>\n"
        "            default: -6\n"
        "            gap extension penalty (must be negative)\n"
        "        -l, --algorithm <int>\n"
        "            default: 0\n"
        "            alignment mode:\n"
        "                0 - local (Smith-Waterman)\n"
        "                1 - global (Needleman-Wunsch)\n"
        "                2 - semi-global\n"
        "        -r, --result <int>\n"
        "            default: 0\n"
        "            result mode:\n"
        "                0 - consensus\n"
        "                1 - multiple sequence alignment\n"
        "                2 - 0 & 1\n"
        "        -d, --dot <file>\n"
        "            output file for the final POA graph in DOT format\n"
        "        --version\n"
        "            prints the version number\n"
        "        -h, --help\n"
        "            prints the usage\n");
}
