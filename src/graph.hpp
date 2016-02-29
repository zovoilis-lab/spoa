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

    friend std::unique_ptr<Graph> createGraph(const std::string& sequence);

private:

    Graph(const std::string& sequence);
    Graph(const Graph&) = delete;
    const Graph& operator=(const Graph&) = delete;

    uint32_t num_sequences_;
    std::vector<NodeSharedPtr> nodes_;
};
