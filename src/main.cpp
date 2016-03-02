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
    graph->print();

    for (uint32_t i = 1; i < reads.size(); ++i) {
        auto alignment = createAlignment(reads[i]->data(), graph,
            AlignmentParams(4, -2, -4, -2,  AlignmentType::kNW));

        alignment->align_sequence_to_graph();
        alignment->backtrack();

        const auto& node_ids = alignment->alignment_node_ids();
        const auto& seq_ids = alignment->alignment_seq_ids();

        for (const auto& id: node_ids) {
            if (id == -1) {
                fprintf(stderr, "-");
            } else {
                fprintf(stderr, "%c", graph->node(id)->letter());
            }
        }
        fprintf(stderr, "\n");
        for (const auto& id: seq_ids) {
            if (id == -1) {
                fprintf(stderr, "-");
            } else {
                fprintf(stderr, "%c", reads[i]->data()[id]);
            }
        }
        fprintf(stderr, "\n");
        graph->add_alignment(node_ids, seq_ids, reads[i]->data());
        graph->print();
    }

    return 0;
}
