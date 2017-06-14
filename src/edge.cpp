/*!
 * @file edge.cpp
 *
 * @brief Edge class source file
 */

#include "spoa/graph.hpp"
#include "edge.hpp"

namespace spoa {

Edge::Edge(uint32_t begin_node_id, uint32_t end_node_id, uint32_t label,
    float weight)
        : begin_node_id_(begin_node_id), end_node_id_(end_node_id),
        sequence_labels_(1, label), sequence_weights_(1, weight),
        total_weight_(weight) {
}

Edge::~Edge() {
}

void Edge::add_sequence(uint32_t label, float weight) {
    sequence_labels_.emplace_back(label);
    sequence_weights_.emplace_back(weight);
    total_weight_ += weight;
}

}
