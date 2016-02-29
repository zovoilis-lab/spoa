/*!
 * @file graph.cpp
 *
 * @brief Graph class source file
 */

#include <assert.h>
#include <set>

#include "node.hpp"
#include "edge.hpp"
#include "graph.hpp"

std::unique_ptr<Graph> createGraph(const std::string& sequence) {
    return std::unique_ptr<Graph>(new Graph(sequence));
}

Graph::Graph(const std::string& sequence) :
        num_sequences_(), nodes_(), is_sorted_(false), sorted_nodes_ids_() {

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

void Graph::topological_sort() {

    sorted_nodes_ids_.clear();

    // 0 - unmarked, 1 - temporarily marked, 2 - permanently marked
    std::vector<uint8_t> marks(nodes_.size(), 0);

    while (true) {
        uint16_t i = 0;
        for (; i < nodes_.size(); ++i) {
            if (marks[nodes_[i]->id()] == 0) break;
        }
        if (i == nodes_.size()) break;

        visit_node(nodes_[i], marks);
    }

    assert(is_topologically_sorted() == true);
    is_sorted_ = true;
}

void Graph::visit_node(NodeSharedPtr node, std::vector<uint8_t>& marks) {
    assert(marks[node->id()] != 1 && "Graph is not a DAG!");

    if (marks[node->id()] == 0) {
        marks[node->id()] = 1;

        const auto& in_edges = node->in_edges();
        for (const auto& edge: in_edges) {
            visit_node(edge->begin_node(), marks);
        }

        marks[node->id()] = 2;
        sorted_nodes_ids_.emplace_back(node->id());
    }
}

bool Graph::is_topologically_sorted() const {
    assert(nodes_.size() == sorted_nodes_ids_.size());

    std::set<uint16_t> visited_nodes;
    for (const auto& id: sorted_nodes_ids_) {
        for (const auto& edge: nodes_[id]->in_edges()) {
            if (visited_nodes.count(edge->begin_node()->id()) == 0) {
                return false;
            }
        }
        visited_nodes.insert(id);
    }

    return true;
}
