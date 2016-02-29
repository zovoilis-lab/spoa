/*!
 * @file graph.hpp
 *
 * @brief Graph class header file
 */

#pragma once

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

    const std::vector<NodeSharedPtr>& nodes() const {
        return nodes_;
    }

    const std::vector<uint16_t> sorted_nodes_ids() const {
        return sorted_nodes_ids_;
    }

    void topological_sort();

    friend std::unique_ptr<Graph> createGraph(const std::string& sequence);

private:

    void visit_node(NodeSharedPtr node, std::vector<uint8_t>& marks);
    bool is_topologically_sorted() const;

    Graph(const std::string& sequence);
    Graph(const Graph&) = delete;
    const Graph& operator=(const Graph&) = delete;

    uint32_t num_sequences_;
    std::vector<NodeSharedPtr> nodes_;

    bool is_sorted_;
    std::vector<uint16_t> sorted_nodes_ids_;
};
