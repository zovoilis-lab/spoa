/*!
 * @file node.cpp
 *
 * @brief Node class source file
 */

#include <unordered_set>

#include "edge.hpp"
#include "node.hpp"

namespace SPOA {

std::unique_ptr<Node> createNode(uint32_t id, char letter) {
    return std::unique_ptr<Node>(new Node(id, letter));
}

Node::Node(uint32_t id, char letter)
        : id_(id), letter_(letter), in_edges_(), out_edges_(), aligned_nodes_ids_() {
}

uint32_t Node::coverage() const {

    std::unordered_set<uint32_t> label_set;
    for (const auto& edge: in_edges_) {
        for (const auto& label: edge->sequence_labels()) {
            label_set.insert(label);
        }
    }
    for (const auto& edge: out_edges_) {
        for (const auto& label: edge->sequence_labels()) {
            label_set.insert(label);
        }
    }

    return label_set.size();
}

Node::~Node() {
}

}
