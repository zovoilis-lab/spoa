#include <stdio.h>
#include <stdlib.h>

#include "alignment.hpp"
#include "graph.hpp"
#include "node.hpp"
#include "edge.hpp"

int main(int argc, char** argv) {

    GraphSharedPtr graph = createGraph("huehueheh");
    graph->topological_sort();

    auto alignment = createAlignment("huehuehue", graph,
        AlignmentParams(1, -1, -1, -1,  AlignmentType::kNW));

    alignment->align_sequence_to_graph();

    NodeSharedPtr node1 = createNode(0, 'c');
    NodeSharedPtr node2 = createNode(1, 'c');
    EdgeSharedPtr edge = createEdge(node1, node2, 0);
    node1->add_out_edge(edge);
    node2->add_in_edge(edge);
    edge->add_sequence_label(1);
    edge->add_sequence_label(2);

    return 0;
}
