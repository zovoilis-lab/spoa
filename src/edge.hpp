/*!
 * @file edge.hpp
 *
 * @brief Edge class header file
 */

#pragma once

#include <assert.h>
#include <vector>
#include <memory>

class Node;
using NodeWeakPtr = std::weak_ptr<Node>;
using NodeSharedPtr = std::shared_ptr<Node>;

class Edge;
std::unique_ptr<Edge> createEdge(NodeWeakPtr begin_node, NodeWeakPtr end_node,
    uint32_t label);

class Edge {
public:

    ~Edge();

    NodeSharedPtr begin_node() const {
        assert(begin_node_.expired() == 0);
        return begin_node_.lock();
    }

    NodeSharedPtr end_node() const {
        assert(end_node_.expired() == 0);
        return end_node_.lock();
    }

    const std::vector<uint32_t>& sequence_labels() const {
        return sequence_labels_;
    }

    void add_sequence_label(uint32_t label) {
        sequence_labels_.emplace_back(label);
    }

    friend std::unique_ptr<Edge> createEdge(NodeWeakPtr begin_node,
        NodeWeakPtr end_node, uint32_t label);

private:

    Edge(NodeWeakPtr begin_node, NodeWeakPtr end_node, uint32_t label);
    Edge(const Edge&) = delete;
    const Edge& operator=(const Edge&) = delete;

    NodeWeakPtr begin_node_;
    NodeWeakPtr end_node_;
    std::vector<uint32_t> sequence_labels_;
};
