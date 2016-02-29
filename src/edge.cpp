/*!
 * @file edge.cpp
 *
 * @brief Edge class source file
 */

#include "edge.hpp"

std::unique_ptr<Edge> createEdge(NodeWeakPtr begin_node, NodeWeakPtr end_node,
    uint32_t label) {
    return std::unique_ptr<Edge>(new Edge(begin_node, end_node, label));
}

Edge::Edge(NodeWeakPtr begin_node, NodeWeakPtr end_node, uint32_t label) :
        begin_node_(begin_node), end_node_(end_node), sequence_labels_(1, label) {
}

Edge::~Edge() {
}
