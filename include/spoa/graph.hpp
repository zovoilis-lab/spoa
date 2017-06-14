/*!
 * @file graph.hpp
 *
 * @brief Graph class header file
 */

#pragma once

#include <memory>
#include <string>
#include <vector>
#include <unordered_set>

namespace spoa {

class Alignment;
class SisdAlignment;
class SimdAlignment;

class Graph;
std::unique_ptr<Graph> createGraph();

class Graph {
public:
    ~Graph();

    void add_alignment(const std::vector<std::pair<int32_t, int32_t>>& alignment,
        const std::string& sequence, float weight = 1.0);

    void add_alignment(const std::vector<std::pair<int32_t, int32_t>>& alignment,
        const std::string& sequence, const std::string& quality);

    void add_alignment(const std::vector<std::pair<int32_t, int32_t>>& alignment,
        const std::string& sequence, const std::vector<float>& weights);

    void generate_multiple_sequence_alignment(std::vector<std::string>& dst,
        bool include_consensus = false);

    std::string generate_consensus();
    // returns coverages
    std::string generate_consensus(std::vector<uint32_t>& dst);

    void print() const;

    friend Alignment;
    friend SisdAlignment;
    friend SimdAlignment;

    friend std::unique_ptr<Graph> createGraph();
private:
    Graph();
    Graph(const Graph&) = delete;
    const Graph& operator=(const Graph&) = delete;

    void topological_sort(bool rigorous = false);

    bool is_topologically_sorted() const;

    void traverse_heaviest_bundle();

    uint32_t branch_completion(std::vector<float>& scores,
        std::vector<int32_t>& predecessors,
        uint32_t rank);

    void check_msa(const std::vector<std::string>& msa,
        const std::vector<std::string>& sequences,
        const std::vector<uint32_t>& indices) const;

    std::unique_ptr<Graph> subgraph(uint32_t begin_node_id,
        uint32_t end_node_id, std::vector<int32_t>& map);

    void extract_subgraph_nodes(std::vector<bool>& dst, uint32_t current_node_id,
        uint32_t end_node_id) const;

    uint32_t add_node(char letter);

    void add_edge(uint32_t begin_node_id, uint32_t end_node_id, float weight);

    int32_t add_sequence(const std::string& sequence,
        const std::vector<float>& weights,
        uint32_t begin, uint32_t end);

    class Node;
    static std::unique_ptr<Node> createNode(uint32_t id, char letter);

    class Edge;
    static std::unique_ptr<Edge> createEdge(uint32_t begin_node_id,
        uint32_t end_node_id, uint32_t label, float weight);

    uint32_t num_sequences_;
    uint32_t num_nodes_;
    std::vector<std::shared_ptr<Node>> nodes_;
    std::unordered_set<uint8_t> alphabet_;
    bool is_sorted_;
    std::vector<uint32_t> sorted_nodes_ids_;
    std::vector<uint32_t> sequences_start_nodes_ids_;
    std::vector<uint32_t> consensus_;
};

}
