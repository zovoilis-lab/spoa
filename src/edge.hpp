/*!
 * @file edge.hpp
 *
 * @brief Edge class header file
 */

#pragma once

#include <assert.h>
#include <vector>
#include <memory>


class Edge;
std::unique_ptr<Edge> createEdge(uint32_t begin_node_id, uint32_t end_node_id,
    uint32_t label);

class Edge {
public:

    ~Edge();

    uint32_t begin_node_id() const {
        return begin_node_id_;
    }

    uint32_t end_node_id() const {
        return end_node_id_;
    }

    const std::vector<uint32_t>& sequence_labels() const {
        return sequence_labels_;
    }

    void add_sequence_label(uint32_t label) {
        sequence_labels_.emplace_back(label);
    }

    friend std::unique_ptr<Edge> createEdge(uint32_t begin_node_id,
        uint32_t end_node_id, uint32_t label);

private:

    Edge(uint32_t begin_node_id, uint32_t end_node_id, uint32_t label);
    Edge(const Edge&) = delete;
    const Edge& operator=(const Edge&) = delete;

    uint32_t begin_node_id_;
    uint32_t end_node_id_;
    std::vector<uint32_t> sequence_labels_;
};
