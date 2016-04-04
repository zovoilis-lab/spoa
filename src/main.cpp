#ifdef SPOA_MAIN_

#include <stdio.h>
#include <stdlib.h>

#include "spoa.hpp"
#include "chain.hpp"

using namespace SPOA;

int main(int argc, char** argv) {

    ChainSet reads;
    createChainSet(reads, argv[1]);

    std::vector<std::string> sequences;
    for (const auto& it: reads) {
        sequences.emplace_back(it->data());
    }

    auto params = AlignmentParams(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),
        atoi(argv[5]), (AlignmentType) atoi(argv[6]));

    std::string consensus = generate_consensus(sequences, params, true);

    fprintf(stderr, "Consensus (%zu)\n", consensus.size());
    fprintf(stderr, "%s\n", consensus.c_str());

    std::vector<std::string> msa;
    generate_msa(msa, sequences, params, true);

    fprintf(stderr, "Multiple sequence alignment\n");
    for (const auto& it: msa) {
        fprintf(stderr, "%s\n", it.c_str());
    }

    return 0;
}

#endif
