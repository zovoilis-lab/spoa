/*!
 * @file poa.cpp
 *
 * @brief Poa source file which encapsulates the implementation
 */

#include <stdlib.h>
#include <stdint.h>
#include <algorithm>

#include "graph.hpp"
#include "spoa.hpp"

void prepare_indices(std::vector<uint32_t>& dst, const std::vector<std::string>& sequences, bool sorted) {
    dst.resize(sequences.size());
    std::iota(dst.begin(), dst.end(), static_cast<uint32_t>(0));

    if (sorted) {
        std::sort(dst.begin(), dst.end(),
            [&](uint32_t lhs, uint32_t rhs) {
                return sequences[lhs].size() > sequences[rhs].size();
            }
        );
    }
}

namespace SPOA {

std::shared_ptr<Graph> construct_partial_order_graph(const std::vector<std::string>& sequences,
    AlignmentParams params, bool sorted) {

    std::vector<uint32_t> indices;
    prepare_indices(indices, sequences, sorted);

    std::shared_ptr<Graph> graph = createGraph(sequences[indices.front()]);

    for (uint32_t i = 1; i < sequences.size(); ++i) {
        auto alignment = createAlignment(sequences[indices[i]], graph, params);
        alignment->align_sequence_to_graph();
        alignment->backtrack();
        graph->add_alignment(std::move(alignment), sequences[indices[i]]);
    }

    return graph;
}

std::shared_ptr<Graph> construct_partial_order_graph(const std::vector<std::string>& sequences,
    const std::vector<std::string>& qualities, AlignmentParams params, bool sorted) {

    std::vector<uint32_t> indices;
    prepare_indices(indices, sequences, sorted);

    std::shared_ptr<Graph> graph = createGraph(sequences[indices.front()], qualities[indices.front()]);

    for (uint32_t i = 1; i < sequences.size(); ++i) {
        auto alignment = createAlignment(sequences[indices[i]], graph, params);
        alignment->align_sequence_to_graph();
        alignment->backtrack();
        graph->add_alignment(std::move(alignment), sequences[indices[i]], qualities[indices[i]]);
    }

    return graph;
}

std::shared_ptr<Graph> construct_partial_order_graph(const std::vector<std::string>& sequences,
    const std::vector<std::string>& qualities, const std::vector<uint32_t>& begin_positions,
    const std::vector<uint32_t>& end_positions, AlignmentParams params) {

    std::shared_ptr<Graph> graph = createGraph(sequences.front(), qualities.front());

    for (uint32_t i = 1; i < sequences.size(); ++i) {
        std::vector<int32_t> mapping;
        std::shared_ptr<Graph> subgraph = graph->subgraph(begin_positions[i], end_positions[i], mapping);

        auto alignment = createAlignment(sequences[i], subgraph, params);
        alignment->align_sequence_to_graph();
        alignment->backtrack();

        alignment->update_node_ids(mapping);

        graph->add_alignment(std::move(alignment), sequences[i], qualities[i]);
    }

    return graph;
}

std::string generate_consensus(const std::vector<std::string>& sequences,
    AlignmentParams params, bool sorted) {

    auto graph = construct_partial_order_graph(sequences, params, sorted);
    return graph->generate_consensus();
}

std::string generate_consensus(const std::vector<std::string>& sequences,
    const std::vector<std::string>& qualities, AlignmentParams params, bool sorted) {

    auto graph = construct_partial_order_graph(sequences, qualities, params, sorted);
    return graph->generate_consensus();
}

std::string generate_consensus(const std::vector<std::string>& sequences,
    const std::vector<std::string>& qualities, const std::vector<uint32_t>& begin_positions,
    const std::vector<uint32_t>& end_positions, AlignmentParams params) {

    auto graph = construct_partial_order_graph(sequences, qualities, begin_positions,
        end_positions, params);
    return graph->generate_consensus();
}

void generate_msa(std::vector<std::string>& dst, const std::vector<std::string>& sequences,
    AlignmentParams params, bool sorted) {

    std::vector<uint32_t> indices;
    prepare_indices(indices, sequences, sorted);

    auto graph = construct_partial_order_graph(sequences, params, sorted);
    graph->generate_msa(dst);
    graph->check_msa(dst, sequences, indices);
}

void generate_msa(std::vector<std::string>& dst, const std::vector<std::string>& sequences,
    const std::vector<std::string>& qualities, AlignmentParams params, bool sorted) {

    std::vector<uint32_t> indices;
    prepare_indices(indices, sequences, sorted);

    auto graph = construct_partial_order_graph(sequences, qualities, params, sorted);
    graph->generate_msa(dst);
    graph->check_msa(dst, sequences, indices);
}

}
