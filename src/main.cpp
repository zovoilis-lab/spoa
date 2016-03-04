#include <stdio.h>
#include <stdlib.h>

#include "alignment.hpp"
#include "graph.hpp"
#include "node.hpp"
#include "edge.hpp"

#include "chain.hpp"

int main(int argc, char** argv) {

    ChainSet reads;
    createChainSet(reads, argv[1]);

    GraphSharedPtr graph = createGraph(reads[0]->data());
    graph->topological_sort();
    //graph->print();

    for (uint32_t i = 1; i < reads.size(); ++i) {
        auto alignment = createAlignment(reads[i]->data(), graph,
            AlignmentParams(1, -1, -1, -1, atoi(argv[2]) == 0 ? AlignmentType::kNW : AlignmentType::kSW));

        alignment->align_sequence_to_graph();
        alignment->backtrack();

        const auto& node_ids = alignment->alignment_node_ids();
        const auto& seq_ids = alignment->alignment_seq_ids();

        graph->add_alignment(node_ids, seq_ids, reads[i]->data());
        //graph->print();
    }

    std::vector<std::string> msa;
    graph->generate_msa(msa);
    for (const auto& alignment_str: msa) {
        fprintf(stderr, "%s\n", alignment_str.c_str());
    }
    fprintf(stderr, "\n");
    std::string consensus = graph->generate_consensus();
    fprintf(stderr, "Consensus (%zu)\n", consensus.size());
    fprintf(stderr, "%s\n", consensus.c_str());

    return 0;
}
