#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "chain.hpp"

#include "spoa/spoa.hpp"

static struct option options[] = {
    {"fasta-sequences", required_argument, 0, 'a'},
    {"fastq-sequences", required_argument, 0, 'q'},
    {"match", required_argument, 0, 'm'},
    {"mismatch", required_argument, 0, 'x'},
    {"gap-open", required_argument, 0, 'o'},
    {"gap-extend", required_argument, 0, 'e'},
    {"algorithm", required_argument, 0, 'l'},
    {"result", required_argument, 0, 'r'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

void help();

int main(int argc, char** argv) {

    std::string fasta_path;
    std::string fastq_path;

    int8_t match = 5;
    int8_t mismatch = -4;
    int8_t gap = -8;

    uint8_t algorithm = 0;
    uint8_t result = 0;

    while (true) {
        auto argument = getopt_long(argc, argv, "a:q:m:x:g:l:r:h", options,
            nullptr);
        if (argument == -1) {
            break;
        }

        switch (argument) {
            case 'a':
                fasta_path = optarg;
                break;
            case 'q':
                fastq_path = optarg;
                break;
            case 'm':
                match = atoi(optarg);
                break;
            case 'x':
                mismatch = atoi(optarg);
                break;
            case 'g':
                gap = atoi(optarg);
                break;
            case 'l':
                algorithm = atoi(optarg);
                break;
            case 'r':
                result = atoi(optarg);
                break;
            case 'h':
            default:
                help();
                return -1;
        }
    }

    if (fasta_path.empty() && fastq_path.empty()) {
        fprintf(stderr, "spoa:: error: "
            "missing option -a/-q (sequences file)\n");
        help();
        return -1;
    }

    auto alignment_engine = spoa::createAlignmentEngine(
        static_cast<spoa::AlignmentType>(algorithm), match, mismatch, gap);

    auto graph = spoa::createGraph();

    std::vector<std::unique_ptr<spoa::Chain>> chains;

    if (!fasta_path.empty()) {
        auto creader = bioparser::createReader<spoa::Chain,
            bioparser::FastaReader>(fasta_path);
        creader->read_objects(chains, -1);

        size_t max_sequence_size = 0;
        for (const auto& it: chains) {
            max_sequence_size = std::max(max_sequence_size, it->data().size());
        }
        alignment_engine->prealloc(max_sequence_size, 4);

        for (const auto& it: chains) {
            auto alignment = alignment_engine->align_sequence_with_graph(
                it->data(), graph);
            graph->add_alignment(alignment, it->data());
        }
    } else {
        auto creader = bioparser::createReader<spoa::Chain,
            bioparser::FastqReader>(fastq_path);
        creader->read_objects(chains, -1);

        size_t max_sequence_size = 0;
        for (const auto& it: chains) {
            max_sequence_size = std::max(max_sequence_size, it->data().size());
        }
        alignment_engine->prealloc(max_sequence_size, 4);

        for (const auto& it: chains) {
            auto alignment = alignment_engine->align_sequence_with_graph(
                it->data(), graph);
            graph->add_alignment(alignment, it->data(), it->quality());
        }
    }

    if (result == 0 || result == 2) {
        std::string consensus = graph->generate_consensus();
        fprintf(stderr, "Consensus (%zu)\n", consensus.size());
        fprintf(stderr, "%s\n", consensus.c_str());
    }

    if (result == 1 || result == 2) {
        std::vector<std::string> msa;
        graph->generate_multiple_sequence_alignment(msa);
        fprintf(stderr, "Multiple sequence alignment\n");
        for (const auto& it: msa) {
            fprintf(stderr, "%s\n", it.c_str());
        }
    }

    return 0;
}

void help() {
    printf(
        "usage: spoa -a <FASTA file> (or -q <FASTQ file>) [arguments ...]\n"
        "\n"
        "arguments:\n"
        "    -a, --fasta-sequences <file>\n"
        "        (required)\n"
        "        input FASTA file\n"
        "    -q, --fastq-sequences <file>\n"
        "        (required/alternative of -a)\n"
        "        input FASTQ file\n"
        "    -m, --match <int>\n"
        "        default: 5\n"
        "        score for matching bases\n"
        "    -x, --mismatch <int>\n"
        "        default: -4\n"
        "        score for mismatching bases\n"
        "    -g, --gap <int>\n"
        "        default: -8\n"
        "        gap penalty\n"
        "    -l, --algorithm <int>\n"
        "        default: 0\n"
        "        alignment mode: 0 - local (Smith-Waterman)\n"
        "                        1 - global (Needleman-Wunsch)\n"
        "                        2 - semi-global\n"
        "    -r, --result <int>\n"
        "        default: 0\n"
        "        result mode: 0 - consensus\n"
        "                     1 - multiple sequence alignment\n"
        "                     2 - 0 + 1\n"
        "    -h, --help\n"
        "        prints out the help\n");
}
