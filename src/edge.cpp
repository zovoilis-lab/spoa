/*!
 * @file edge.cpp
 *
 * @brief Edge class source file
 */

#include "edge.hpp"

std::unique_ptr<Edge> createEdge(uint32_t begin_node_id, uint32_t end_node_id,
    uint32_t label) {
    return std::unique_ptr<Edge>(new Edge(begin_node_id, end_node_id, label));
}

Edge::Edge(uint32_t begin_node_id, uint32_t end_node_id, uint32_t label) :
        begin_node_id_(begin_node_id), end_node_id_(end_node_id),
        sequence_labels_(1, label) {
}

Edge::~Edge() {
}
