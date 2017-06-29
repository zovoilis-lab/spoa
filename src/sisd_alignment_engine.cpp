/*!
 * @file sisd_alignment_engine.cpp
 *
 * @brief SisdAlignmentEngine class source file
 */

#include <limits>
#include <algorithm>

#include "spoa/graph.hpp"
#include "sisd_alignment_engine.hpp"

namespace spoa {

constexpr int32_t kNegativeInfinity = std::numeric_limits<int32_t>::min() + 1024;

std::unique_ptr<AlignmentEngine> createSisdAlignmentEngine(
    AlignmentType alignment_type, int8_t match, int8_t mismatch,
    int8_t gap_open, int8_t gap_extend) {

    return std::unique_ptr<AlignmentEngine>(new SisdAlignmentEngine(
        alignment_type, match, mismatch, gap_open, gap_extend));
}

struct SisdAlignmentEngine::Implementation {
    std::vector<uint32_t> node_id_to_rank;
    std::vector<int32_t> sequence_profile;
    std::vector<int32_t> X;
    int32_t* H;
    int32_t* F;
    int32_t* E;

    Implementation()
            : node_id_to_rank(), sequence_profile(), X(), H(nullptr),
            F(nullptr), E(nullptr) {
    }
};

SisdAlignmentEngine::SisdAlignmentEngine(AlignmentType alignment_type,
    int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend)
        : AlignmentEngine(alignment_type, match, mismatch, gap_open, gap_extend),
        pimpl_(new Implementation()) {
}

SisdAlignmentEngine::~SisdAlignmentEngine() {
}

void SisdAlignmentEngine::prealloc(uint32_t max_sequence_size,
    uint32_t alphabet_size) {

    this->realloc(alphabet_size * max_sequence_size, alphabet_size,
        max_sequence_size);
}

void SisdAlignmentEngine::realloc(uint32_t matrix_height, uint32_t num_codes,
    uint32_t matrix_width) {

    if (pimpl_->node_id_to_rank.size() < matrix_height - 1) {
        pimpl_->node_id_to_rank.resize(matrix_height - 1, 0);
    }
    if (pimpl_->sequence_profile.size() < num_codes * matrix_width) {
        pimpl_->sequence_profile.resize(num_codes * matrix_width, 0);
    }
    if (pimpl_->X.size() < 3 * matrix_height * matrix_width) {
        pimpl_->X.resize(3 * matrix_width * matrix_height, 0);
        pimpl_->H = &(pimpl_->X[0]);
        pimpl_->F = &(pimpl_->H[matrix_height * matrix_width]);
        pimpl_->E = &(pimpl_->F[matrix_height * matrix_width]);
    }
}

Alignment SisdAlignmentEngine::align_sequence_with_graph(
    const std::string& sequence, const std::unique_ptr<Graph>& graph) {

    if (graph->nodes().empty() || sequence.empty()) {
        return Alignment();
    }

    uint32_t matrix_width = sequence.size() + 1;
    uint32_t matrix_height = graph->nodes().size() + 1;
    const auto& sorted_nodes_ids = graph->sorted_nodes_ids();

    // realloc
    this->realloc(matrix_height, graph->num_codes(), matrix_width);

    // initialize
    for (uint32_t i = 0; i < graph->num_codes(); ++i) {
        char c = graph->decoder(i);
        pimpl_->sequence_profile[i * matrix_width] = 0;
        for (uint32_t j = 0; j < sequence.size(); ++j) {
            pimpl_->sequence_profile[i * matrix_width + (j + 1)] =
                (c == sequence[j] ? match_ : mismatch_);
        }
    }
    for (uint32_t i = 0; i < sorted_nodes_ids.size(); ++i) {
        pimpl_->node_id_to_rank[sorted_nodes_ids[i]] = i;
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

    for (uint32_t i = 0; i < matrix_height; ++i) {
        pimpl_->E[i * matrix_width] = kNegativeInfinity;
    }
    for (uint32_t j = 0; j < matrix_width; ++j) {
        pimpl_->F[j] = kNegativeInfinity;
    }

    if (alignment_type_ == AlignmentType::kSW) {
        for (uint32_t i = 0; i < matrix_height; ++i) {
            pimpl_->H[i * matrix_width] = 0;
        }
        for (uint32_t j = 0; j < matrix_width; ++j) {
            pimpl_->H[j] = 0;
        }
    }

    if (alignment_type_ == AlignmentType::kNW) {
        for (const auto& node_id: sorted_nodes_ids) {
            uint32_t i = pimpl_->node_id_to_rank[node_id] + 1;
            const auto& node = graph->nodes()[node_id];
            if (node->in_edges().empty()) {
                pimpl_->H[i * matrix_width] = gap_open_;
            } else {
                int32_t penalty = kNegativeInfinity;
                for (const auto& edge: node->in_edges()) {
                    uint32_t pred_i =
                        pimpl_->node_id_to_rank[edge->begin_node_id()] + 1;
                    penalty = std::max(penalty, pimpl_->H[pred_i * matrix_width]);
                }
                pimpl_->H[i * matrix_width] = penalty + gap_extend_;
            }
            // needed for easier NW backtrack
            pimpl_->F[i * matrix_width] = pimpl_->H[i * matrix_width];
        }
        for (uint32_t j = 1; j < matrix_width; ++j) {
            pimpl_->H[j] = gap_open_ + (j - 1) * gap_extend_;
            // needed for easier NW backtrack
            pimpl_->E[j] = pimpl_->H[j];
        }
        max_score = kNegativeInfinity;
    }

    if (alignment_type_ == AlignmentType::kOV) {
        for (uint32_t i = 0; i < matrix_height; ++i) {
            pimpl_->H[i * matrix_width] = 0;
        }
        for (uint32_t j = 1; j < matrix_width; ++j) {
            pimpl_->H[j] = gap_open_ + (j - 1) * gap_extend_;
        }
        max_score = kNegativeInfinity;
    }

    // alignment
    for (uint32_t node_id: sorted_nodes_ids) {
        const auto& node = graph->nodes()[node_id];
        const auto& char_profile =
            &(pimpl_->sequence_profile[node->code() * matrix_width]);

        uint32_t i = pimpl_->node_id_to_rank[node_id] + 1;

        int32_t* H_row = &(pimpl_->H[i * matrix_width]);
        int32_t* F_row = &(pimpl_->F[i * matrix_width]);

        uint32_t pred_i = node->in_edges().empty() ? 0 :
            pimpl_->node_id_to_rank[node->in_edges()[0]->begin_node_id()] + 1;

        int32_t* H_pred_row = &(pimpl_->H[pred_i * matrix_width]);
        int32_t* F_pred_row = &(pimpl_->F[pred_i * matrix_width]);

        for (uint32_t j = 1; j < matrix_width; ++j) {
            // update F
            F_row[j] = std::max(H_pred_row[j] + gap_open_,
                F_pred_row[j] + gap_extend_);
            // update H
            H_row[j] = H_pred_row[j - 1] + char_profile[j];
        }

        // check other predeccessors
        for (uint32_t p = 1; p < node->in_edges().size(); ++p) {
            pred_i = pimpl_->node_id_to_rank[node->in_edges()[p]->begin_node_id()] + 1;

            H_pred_row = &(pimpl_->H[pred_i * matrix_width]);
            F_pred_row = &(pimpl_->F[pred_i * matrix_width]);

            for (uint32_t j = 1; j < matrix_width; ++j) {
                // update F
                F_row[j] = std::max(F_row[j], std::max(H_pred_row[j] + gap_open_,
                    F_pred_row[j] + gap_extend_));
                // update H
                H_row[j] = std::max(H_row[j], H_pred_row[j - 1] + char_profile[j]);
            }
        }

        int32_t* E_row = &(pimpl_->E[i * matrix_width]);

        for (uint32_t j = 1; j < matrix_width; ++j) {
            // update E
            E_row[j] = std::max(H_row[j - 1] + gap_open_,
                E_row[j - 1] + gap_extend_);
            // update H
            H_row[j] = std::max(H_row[j], std::max(F_row[j], E_row[j]));

            if (alignment_type_ == AlignmentType::kSW) {
                H_row[j] = std::max(H_row[j], 0);
                update_max_score(H_row, i, j);

            } else if (alignment_type_ == AlignmentType::kNW &&
                (j == matrix_width - 1 && node->out_edges().empty())) {
                update_max_score(H_row, i, j);

            } else if (alignment_type_ == AlignmentType::kOV &&
                (node->out_edges().empty())) {
                update_max_score(H_row, i, j);
            }
        }
    }

    printf("%d| %d %d\n", max_score, max_i, max_j);

    // backtrack
    Alignment alignment;

    uint32_t i = max_i;
    uint32_t j = max_j;

    auto sw_condition = [this, &i, &j, &matrix_width]() {
        return (pimpl_->H[i * matrix_width + j] == 0) ? false : true;
    };
    auto nw_condition = [&i, &j]() {
        return (i == 0 && j == 0) ? false : true;
    };
    auto ov_condition = [&i, &j]() {
        return (i == 0 || j == 0) ? false : true;
    };

    uint32_t prev_i = 0;
    uint32_t prev_j = 0;

    while ((alignment_type_ == AlignmentType::kSW && sw_condition()) ||
        (alignment_type_ == AlignmentType::kNW && nw_condition()) ||
        (alignment_type_ == AlignmentType::kOV && ov_condition())) {

        auto H_ij = pimpl_->H[i * matrix_width + j];
        bool predecessor_found = false;

        if (i != 0) {
            const auto& node = graph->nodes()[sorted_nodes_ids[i - 1]];
            int32_t match_cost =
                pimpl_->sequence_profile[node->code() * matrix_width + j];

            uint32_t pred_i = node->in_edges().empty() ? 0 :
                pimpl_->node_id_to_rank[node->in_edges()[0]->begin_node_id()] + 1;

            if (j != 0 && H_ij == pimpl_->H[pred_i * matrix_width + (j - 1)] + match_cost) {
                prev_i = pred_i;
                prev_j = j - 1;
                predecessor_found = true;
            } else if ((H_ij == pimpl_->F[pred_i * matrix_width + j] + gap_extend_) ||
                (H_ij == pimpl_->H[pred_i * matrix_width + j] + gap_open_)) {
                prev_i = pred_i;
                prev_j = j;
                predecessor_found = true;
            }

            if (!predecessor_found) {
                const auto& edges = node->in_edges();
                for (uint32_t p = 1; p < edges.size(); ++p) {
                    uint32_t pred_i =
                        pimpl_->node_id_to_rank[edges[p]->begin_node_id()] + 1;

                    if (j != 0 && H_ij == pimpl_->H[pred_i * matrix_width + (j - 1)] + match_cost) {
                        prev_i = pred_i;
                        prev_j = j - 1;
                        predecessor_found = true;
                        break;
                    }
                    if ((H_ij == pimpl_->F[pred_i * matrix_width + j] + gap_extend_) ||
                        (H_ij == pimpl_->H[pred_i * matrix_width + j] + gap_open_)) {
                        prev_i = pred_i;
                        prev_j = j;
                        predecessor_found = true;
                        break;
                    }
                }
            }
        }

        if (!predecessor_found && H_ij == pimpl_->E[i * matrix_width + j]) {
            prev_i = i;
            prev_j = j - 1;
            predecessor_found = true;
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
