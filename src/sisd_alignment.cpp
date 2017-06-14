/*!
 * @file sisd_alignment.cpp
 *
 * @brief SisdAlignment class source file
 */

#include <limits>

#include "node.hpp"
#include "edge.hpp"
#include "spoa/graph.hpp"
#include "spoa/alignment.hpp"
#include "sisd_alignment.hpp"

namespace spoa {

constexpr uint32_t kMaxAlphabetSize = 256;
constexpr int32_t kBigNegativeValue = std::numeric_limits<int32_t>::min() + 1000;

std::unique_ptr<Alignment> createSisdAlignment(AlignmentType alignment_type,
    uint8_t match, uint8_t mismatch, uint8_t gap_open, uint8_t gap_extend) {

    return std::unique_ptr<Alignment>(new SisdAlignment(sequence, graph,
        std::move(params)));
}

struct SisdAlignment::Impl {
    std::vector<std::vector<int32_t>> sequence_profile;
    std::vector<uint32_t> node_id_to_node_rank;
    std::vector<int32_t> H;
    std::vector<int32_t> F;
    std::vector<int32_t> E;
};

SisdAlignment::SisdAlignment(AlignmentType alignment_type, uint8_t match,
    uint8_t mismatch, uint8_t gap_open, uint8_t gap_extend)
        : Alignment(alignment_type, match, mismatch, gap_open, gap_extend),
        pImpl() {
}

SisdAlignment::~SisdAlignment() {
}

std::vector<std::pair<int32_t, int32_t>> SisdAlignment::align_sequence_with_graph(
    const std::string& sequence, const std::unique_ptr<Graph>& graph) {

    uint32_t matrix_width = sequence.size() + 1;
    uint32_t matrix_height = graph->nodes().size() + 1;
    const auto& sorted_nodes_ids = graph->sorted_nodes_ids();

    auto print_matrix = [&]() -> void {
        for (uint32_t i = 0; i < matrix_height; ++i) {
            for (uint32_t j = 0; j < matrix_width; ++j) {
                printf("(%3d %3d %3d) ", pImpl_->H[i * matrix_width + j],
                    pImpl_->E[i * matrix_width + j],
                    pImpl_->F[i * matrix_width + j]);
            }
            printf("\n");
        }
        printf("\n");
    };

    /* realloc */
    if (pImpl_->sequence_profile.size() < kMaxAlphabetSize * matrix_width) {
        pImpl_->sequence_profile.resize(kMaxAlphabetSize * matrix_width, 0);
    }
    if (pImpl_->node_id_to_node_rank.size() < sorted_nodes_ids.size()) {
        pImpl_->node_id_to_node_rank.resize(sorted_nodes_ids.size(), 0);
    }
    if (pImpl_->H.size() < matrix_width * matrix_height) {
        pImpl_->H.resize(matrix_width * matrix_height, 0);
    }
    if (pImpl_->F.size() < matrix_width * matrix_height) {
        pImpl_->F.resize(matrix_width * matrix_height, kBigNegativeValue);
    }
    if (pImpl_->E.size() < matrix_width * matrix_height) {
        pImpl_->E.resize(matrix_width * matrix_height, kBigNegativeValue);
    }

    /* initialize */
    for (const auto& c: graph->alphabet()) {
        for (uint32_t i = 0; i < sequence.size(); ++i) {
            pImpl_->sequence_profile[(uint32_t) c * kMaxAlphabetSize + i] =
                c == sequence[i] ? match_ : mismatch_;
        }
    }
    for (uint32_t i = 0; i < sorted_nodes_ids.size(); ++i) {
        pImpl_->node_id_to_node_rank[sorted_nodes_ids[i]] = i;
    }

    int32_t max_score = 0;
    int32_t max_i = -1;
    int32_t max_j = -1;
    auto update_max_score = [&max_score, &max_i, &max_j](int32_t* H_row,
        uint32_t i, uint32_t j) -> void {

        if (max_score < H_row[j]) {
            max_score = H_row[j];
            max_i = i;
            max_j = j;
        }
        return;
    };

    if (alignment_type_ == AlignmentType::kNW) {
        max_score = kBigNegativeValue;
        for (const auto& node_id: sorted_nodes_ids) {
            uint32_t i = node_id_to_node_rank[node_id] + 1;
            const auto& node = graph->nodes_[node_id];
            if (node->in_edges_.empty()) {
                pImpl_->H[i * matrix_width] = gap_open_;
            } else {
                int32_t penalty = kBigNegativeValue;
                for (const auto& edge: node->in_edges_) {
                    uint32_t pred_i =
                        pImpl_->node_id_to_node_rank[edge->begin_node_id_] + 1;
                    penalty = std::max(penalty, pImpl_->H[pred_i * matrix_width]);
                }
                pImpl_->H[i * matrix_width] = penalty + gap_extend_;
            }
        }
        for (uint32_t j = 1; j < matrix_width_; ++j) {
            pImpl_->H[j] = gap_open_ + (j - 1) * gap_extend_;
        }

    } else if (alignment_type_ == AlignmentType::kOV) {
        max_score = kBigNegativeValue;
        for (uint32_t j = 1; j < matrix_width_; ++j) {
            pImpl_->H[j] = gap_open_ + (j - 1) * gap_extend_;
        }
    }

    // print_matrix();

    /* alignment */
    for (uint32_t node_id: sorted_nodes_ids) {
        const auto& node = pImpl_->graph->nodes_[node_id];
        const auto& char_profile = pImpl_->sequence_profile[node->letter_];
        uint32_t i = pImpl_->node_id_to_rank_id[node_id] + 1;

        int32_t* H_row = &(pImpl_->H[i * matrix_width]);
        int32_t* F_row = &(pImpl_->F[i * matrix_width]);

        uint32_t pred_i = node->in_edges_.empty() ? 0 :
            pImpl_->node_id_to_rank_id[node->in_edges_[0]->begin_node_id_] + 1;

        int32_t* H_pred_row = &(pImpl_->H[pred_i * matrix_width]);
        int32_t* F_pred_row = &(pImpl_->F[pred_i * matrix_width]);

        for (uint32_t j = 1; j < matrix_width; ++j) {
            // update F
            F_row[j] = std::max(H_pred_row[j] + gap_open_,
                F_pred_row[j] + gap_extend_);
            // update H
            H_row[j] = std::max(H_pred_row[j - 1] + char_profile[j - 1],
                F_row[j]);
        }

        // check other predeccessors
        for (uint32_t p = 1; p < node->in_edges_.size(); ++p) {
            pred_i = pImpl_->node_id_to_rank_id[node->in_edges_[p]->begin_node_id_] + 1;

            H_pred_row = &(pImpl_->H[pred_i * matrix_width]);
            F_pred_row = &(pImpl_->F[pred_i * matrix_width]);

            for (uint32_t j = 1; j < matrix_width; ++j) {
                // update F
                F_row[j] = std::max(F_row[j], std::max(H_pred_row[j] + gap_open_,
                    F_pred_row[j] + gap_extend_));
                // update H
                H_row[j] = std::max(H_row[j], std::max(H_pred_row[j - 1] +
                    char_profile[j - 1], F_row[j]));
            }
        }

        int32_t* E_row = &(pImpl_->E[i * matrix_width]);

        for (uint32_t j = 1; j < matrix_width; ++j) {
            // update E
            E_row[j] = std::max(H_row[j - 1] + gap_open_,
                E_row[j - 1] + gap_extend_);
            // update H
            H_row[j] = std::max(H_row[j], E_row[j]);

            if (alignment_type_ == AlignmentType::kSW) {
                H_row[j] = std::max(H_row[j], 0);
                update_max_score(H_row, i, j);

            } else if (alignment_type_ == AlignmentType::kNW &&
                (j == matrix_width_ - 1 && node->out_edges_.empty())) {
                update_max_score(H_row, i, j);

            } else if (alignment_type_ == AlignmentType::kOV &&
                node->out_edges_.empty())) {
                update_max_score(H_row, i, j);
            }
        }
    }

    /* backtrack */
    std::vector<std::pair<int32_t, int32_t>> alignment;

    uint32_t i = max_i_;
    uint32_t j = max_j_;

    auto sw_condition = [&]() {
        return (pImpl_->H[i * matrix_width + j] == 0) ? false : true;
    };
    auto nw_condition = [&]() {
        return (i == 0 && j == 0) ? false : true;
    };
    auto ov_condition = [&]() {
        return (i == 0 || j == 0) ? false : true;
    };

    uint32_t prev_i = 0;
    uint32_t prev_j = 0;

    while ((alignment_type_ == AlignmentType::kSW && sw_condition()) ||
        (alignment_type_ == AlignmentType::kNW && nw_condition()) ||
        (alignment_type_ == AlignmentType::kOV && ov_condition())) {

        // bloody backtrack
        auto H_ij = pImpl_->H[i * matrix_width_ + j];
        bool predecessor_found = false;

        if (i != 0) {
            const auto& node = graph_->nodes_[sorted_nodes_ids[i - 1]]);
            int32_t match_cost = j != 0 ?
                pImpl_->sequence_profile[node->letter_][j - 1] : 0;

            uint32_t pred_i = node->in_edges_.empty() ? 0 :
                pImpl_->node_id_to_rank_id[node->in_edges_[0]->begin_node_id_] + 1;

            if (j != 0 && H_ij == pImpl_->H[pred_i * matrix_width + (j - 1)] + match_cost) {
                prev_i = pred_i;
                prev_j = j - 1;
                predecessor_found = true;
            } else if ((H_ij == pImpl_->F[pred_i * matrix_width + j] + gap_extend_) ||
                (H_ij == pImpl_->H[pred_i * matrix_width_ + j] + gap_open_)) {
                prev_i = pred_i;
                prev_j = j;
                predecessor_found = true;
            }

            if (!predecessor_found) {
                const auto& edges = node->in_edges_;
                for (uint32_t p = 1; p < edges.size(); ++p) {
                    uint32_t pred_i =
                        pImpl_->node_id_to_node_rank[edges[p]->begin_node_id_] + 1;

                    if (j != 0 && H_ij == pImpl_->H[pred_i * matrix_width + (j - 1)] + match_cost) {
                        prev_i = pred_i;
                        prev_j = j - 1;
                        predecessor_found = true;
                        break;
                    }
                    if ((H_ij == pImpl_->F[pred_i * matrix_width + j] + gap_extend_) ||
                        (H_ij == pImpl_->H[pred_i * matrix_width + j] + gap_open_)){
                        prev_i = pred_i;
                        prev_j = j;
                        predecessor_found = true;
                        break;
                    }
                }
            }
        }

        if (!predecessor_found && H_ij == pImpl_->E[i * matrix_width + j]) {
            prev_i = i;
            prev_j = j - 1;
        }

        alignment.emplace_back(i == prev_i ? -1 : sorted_nodes_ids[i - 1],
            j == prev_j ? -1 : j - 1);

        i = prev_i;
        j = prev_j;
    }

    std::reverse(alignment.begin(), alignment.end());
    return alignment;
}

}
