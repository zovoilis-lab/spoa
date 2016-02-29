/*!
 * @file alignment.cpp
 *
 * @brief Alignment class source file
 */

#include <limits>

#include "node.hpp"
#include "edge.hpp"
#include "graph.hpp"
#include "alignment.hpp"


AlignmentParams::AlignmentParams(int16_t m, int16_t mm, int16_t open,
    int16_t extend, AlignmentType t) :
        match(m), mismatch(mm), gap_open(open), gap_extend(extend), type(t) {
}

AlignmentParams::~AlignmentParams() {
}

Alignment::MatrixElement::MatrixElement(int32_t s, int32_t p_i, int32_t p_j,
    int16_t del, int16_t ins) :
        score(s), prev_i(p_i), prev_j(p_j), deletion_cost(del), insertion_cost(ins) {
}

Alignment::MatrixElement::~MatrixElement() {
}

std::unique_ptr<Alignment> createAlignment(const std::string& sequence,
    GraphSharedPtr graph, AlignmentParams&& params) {

    return std::unique_ptr<Alignment>(new Alignment(sequence, graph,
        std::move(params)));
}

Alignment::Alignment(const std::string& sequence, GraphSharedPtr graph,
    AlignmentParams&& params) :
        sequence_(sequence), graph_(graph), params_(std::move(params)) {

    assert(sequence_.size() != 0);

    int32_t min_score = std::numeric_limits<int32_t>::min();

    matrix_width_ = sequence_.size() + 1;
    matrix_height_ = graph_->nodes().size() + 1;

    matrix_.resize(matrix_width_ * matrix_height_, MatrixElement(min_score,
        -1, -1, params_.gap_open, params_.gap_open));
    matrix_[0].score = 0;

    graph_->topological_sort();
    const auto& sorted_nodes_ids = graph_->sorted_nodes_ids();
    node_id_to_graph_id_.resize(sorted_nodes_ids.size());

    for (uint32_t i = 0; i < sorted_nodes_ids.size(); ++i) {
        node_id_to_graph_id_[sorted_nodes_ids[i]] = i;
    }

    if (params_.type == AlignmentType::kNW) {

        for (uint32_t j = 1; j < matrix_width_; ++j) {
            matrix(0, j).score = j * params_.gap_extend;
            matrix(0, j).prev_i = 0;
            matrix(0, j).prev_j = j - 1;
        }

        for (uint32_t node_id: sorted_nodes_ids) {
            auto node = graph_->node(node_id);
            uint32_t graph_id = node_id_to_graph_id_[node_id];

            if (node->in_edges().size() == 0) {
                matrix(graph_id + 1, 0).score = params_.gap_extend;
                matrix(graph_id + 1, 0).prev_i = 0;
                matrix(graph_id + 1, 0).prev_j = 0;
            } else {
                for (const auto& edge: node->in_edges()) {
                    uint32_t pred_graph_id =
                        node_id_to_graph_id_[edge->begin_node()->id()];
                    if (matrix(graph_id + 1, 0).score < matrix(pred_graph_id + 1, 0).score) {
                        matrix(graph_id + 1, 0).score = matrix(pred_graph_id + 1, 0).score;
                        matrix(graph_id + 1, 0).prev_i = pred_graph_id + 1;
                    }
                }
                matrix(graph_id + 1, 0).prev_j = 0;
                matrix(graph_id + 1, 0).score += params_.gap_extend;
            }
        }
    }

    print_matrix();
}

Alignment::~Alignment() {
}

void Alignment::align_sequence_to_graph() {

    auto pair_score = [&](char lhs, char rhs) {
        return lhs == rhs ? params_.match : params_.mismatch;
    };

    graph_->topological_sort();
    const auto& sorted_nodes_ids = graph_->sorted_nodes_ids();

    for (uint32_t node_id: sorted_nodes_ids) {
        auto node = graph_->node(node_id);
        char graph_letter = node->letter();
        uint32_t i = node_id_to_graph_id_[node_id] + 1;

        for (uint32_t j = 1; j < matrix_width_; ++j) {
            char seq_letter = sequence_[j - 1];
            // match/mismatch
            int32_t score = pair_score(seq_letter, graph_letter);
            if (matrix(i - 1, j - 1).score + score > matrix(i, j).score) {
                matrix(i, j).score = matrix(i - 1, j - 1).score + score;
                matrix(i, j).prev_i = i - 1;
                matrix(i, j).prev_j = j - 1;
            }

            // insertion to sequence
            // deleteion from graph
        }
    }
    print_matrix();
}

void Alignment::print_matrix() {
    for (uint32_t i = 0; i < matrix_height_; ++i) {
        for (uint32_t j = 0; j < matrix_width_; ++j) {
            printf("(%3d, %3d, %3d) ", matrix(i, j).score, matrix(i, j).prev_i,
                matrix(i, j).prev_j);
            }
        printf("\n");
    }
}
