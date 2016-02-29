/*!
 * @file graph.cpp
 *
 * @brief Graph class source file
 */

#include <assert.h>

#include "node.hpp"
#include "edge.hpp"
#include "graph.hpp"

std::unique_ptr<Graph> createGraph(const std::string& sequence) {
    return std::unique_ptr<Graph>(new Graph(sequence));
}

Graph::Graph(const std::string& sequence) :
        num_sequences_(), nodes_() {

    assert(sequence.size() != 0);

    uint32_t node_id = 0;
    nodes_.emplace_back(createNode(node_id++, sequence[0]));

    for (uint32_t i = 1; i < sequence.size(); ++i) {
        nodes_.emplace_back(createNode(node_id++, sequence[i]));
        EdgeSharedPtr edge = createEdge(nodes_[i - 1], nodes_[i], num_sequences_);
        nodes_[i - 1]->add_out_edge(edge);
        nodes_[i]->add_in_edge(edge);
    }

    ++num_sequences_;
}

Graph::~Graph() {
}
