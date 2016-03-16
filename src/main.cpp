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

    AlignmentType type = atoi(argv[2]) == 0 ? AlignmentType::kNW :
        (atoi(argv[2]) == 1 ? AlignmentType::kSW : AlignmentType::kOV);

    for (uint32_t i = 1; i < reads.size(); ++i) {
        auto alignment = createAlignment(reads[i]->data(), graph,
            AlignmentParams(1, -1, atoi(argv[3]), atoi(argv[4]), type));

        alignment->align_sequence_to_graph();
        alignment->backtrack();

        continue;

        const auto& node_ids = alignment->alignment_node_ids();
        const auto& seq_ids = alignment->alignment_seq_ids();

        if (seq_ids.empty() == false) {
            for (const auto& id: seq_ids) {
                if (id == -1) {
                    printf("-");
                } else {
                    printf("%c", reads[i]->data()[id] + 'A');
                }
            }
            printf("\n");
            for (const auto& id: node_ids) {
                if (id == -1) {
                    printf("-");
                } else {
                    printf("%c", graph->node(id)->letter() + 'A');
                }
            }
            printf("\n");
        }

        continue;

        graph->add_alignment(node_ids, seq_ids, reads[i]->data());
        //graph->print();
    }

    /*std::vector<std::string> msa;
    graph->generate_msa(msa);
    for (const auto& alignment_str: msa) {
        fprintf(stderr, "%s\n", alignment_str.c_str());
    }
    fprintf(stderr, "\n");
    std::string consensus = graph->generate_consensus();
    fprintf(stderr, "Consensus (%zu)\n", consensus.size());
    for (const auto& c: consensus) {
        fprintf(stderr, "%c", c + 'A');
    }
    fprintf(stderr, "\n");*/

    return 0;
}
