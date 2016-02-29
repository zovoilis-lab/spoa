#include <stdio.h>
#include <stdlib.h>

#include "graph.hpp"
#include "node.hpp"
#include "edge.hpp"

int main(int argc, char** argv) {

    auto graph = createGraph("huehuehehe");
    printf("%d %zu\n", graph->num_sequences(), graph->nodes().size());

    NodeSharedPtr node1 = createNode(0, 'c');
    NodeSharedPtr node2 = createNode(1, 'c');
    EdgeSharedPtr edge = createEdge(node1, node2, 0);
    node1->add_out_edge(edge);
    node2->add_in_edge(edge);
    edge->add_sequence_label(1);
    edge->add_sequence_label(2);

    //node1.reset();
    auto node3 = edge->begin_node();
    printf("%d\n", node3->id());

    edge.reset();
    const auto& edges = node1->out_edges();
    printf("%zu\n", edges[0]->sequence_labels().size());

    return 0;
}
