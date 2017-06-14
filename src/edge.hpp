/*!
 * @file edge.hpp
 *
 * @brief Edge class header file
 */

#pragma once

#include <stdint.h>
#include <vector>

namespace spoa {

class Graph;

class Graph::Edge {
public:
    ~Edge();

    void add_sequence(uint32_t label, float weight = 1.0);

    friend Graph;
private:
    Edge(uint32_t begin_node_id, uint32_t end_node_id, uint32_t label,
        float weight);
    Edge(const Edge&) = delete;
    const Edge& operator=(const Edge&) = delete;

    uint32_t begin_node_id_;
    uint32_t end_node_id_;
    std::vector<uint32_t> sequence_labels_;
    std::vector<float> sequence_weights_;
    float total_weight_;
};

}
