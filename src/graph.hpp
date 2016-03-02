/*!
 * @file graph.hpp
 *
 * @brief Graph class header file
 */

#pragma once

#include <assert.h>
#include <string>
#include <vector>
#include <memory>

class Node;
using NodeSharedPtr = std::shared_ptr<Node>;

class Graph;
std::unique_ptr<Graph> createGraph(const std::string& sequence);

class Graph {
public:

    ~Graph();

    uint32_t num_sequences() const {
        return num_sequences_;
    }

    const NodeSharedPtr node(uint32_t id) const {
        assert(id < num_nodes_);
        return nodes_[id];
    }

    const std::vector<NodeSharedPtr>& nodes() const {
        return nodes_;
    }

    const std::vector<uint32_t> sorted_nodes_ids() const {
        assert(is_sorted_ == true);
        return sorted_nodes_ids_;
    }

    void topological_sort();

    void add_alignment(const std::vector<int32_t>& node_ids,
        const std::vector<int32_t>& seq_ids, const std::string& sequence);

    void print() const;

    friend std::unique_ptr<Graph> createGraph(const std::string& sequence);

private:

    Graph(const std::string& sequence);
    Graph(const Graph&) = delete;
    const Graph& operator=(const Graph&) = delete;

    uint32_t add_node(char letter);

    void add_edge(uint32_t begin_node_id, uint32_t end_node_id);

    void visit_node(uint32_t node_id, std::vector<uint8_t>& marks);

    bool is_topologically_sorted() const;

    int32_t add_sequence(const std::string& sequence, uint32_t begin, uint32_t end);

    uint32_t num_sequences_;
    uint32_t num_nodes_;
    std::vector<NodeSharedPtr> nodes_;

    bool is_sorted_;
    std::vector<uint32_t> sorted_nodes_ids_;

    std::vector<uint32_t> sequences_start_nodes_ids_;
};
