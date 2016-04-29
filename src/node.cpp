/*!
 * @file node.cpp
 *
 * @brief Node class source file
 */

#include "edge.hpp"
#include "node.hpp"

namespace SPOA {

std::unique_ptr<Node> createNode(uint32_t id, char letter, char type) {
    return std::unique_ptr<Node>(new Node(id, letter, type));
}

Node::Node(uint32_t id, char letter, char type) :
        id_(id), letter_(letter), type_(type), in_edges_(), out_edges_(), aligned_nodes_ids_() {
}

Node::~Node() {
}

}
