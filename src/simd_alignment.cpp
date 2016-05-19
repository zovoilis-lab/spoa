/*!
 * @file simd_alignment.cpp
 *
 * @brief SIMD alignment source file
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <limits>
#include <algorithm>

extern "C" {
    #include <immintrin.h> // AVX2 and lower
}

#include "node.hpp"
#include "edge.hpp"
#include "graph.hpp"
#include "alignment.hpp"
#include "simd_alignment.hpp"

namespace SPOA {

#ifdef __AVX2__

constexpr uint32_t kSimdRegisterSize = 256;
constexpr uint32_t kSimdNumVariables = 8;
constexpr uint32_t kLeftShiftSize = 4;
constexpr uint32_t kRightShiftSize = 28;
using __mxxxi = __m256i;

#define _mmxxx_load _mm256_load_si256
#define _mmxxx_store _mm256_store_si256
#define _mmxxx_and _mm256_and_si256
#define _mmxxx_or _mm256_or_si256
#define _mmxxx_testz _mm256_testz_si256
#define _mmxxx_lshift _mm256_slli_si256
#define _mmxxx_rshift _mm256_srli_si256

#define _mmxxx_add _mm256_add_epi32
#define _mmxxx_sub _mm256_sub_epi32
#define _mmxxx_min _mm256_min_epi32
#define _mmxxx_max _mm256_max_epi32
#define _mmxxx_set1 _mm256_set1_epi32
#define _mmxxx_cmpgt _mm256_cmpgt_epi32

#else // SSE4.1

constexpr uint32_t kSimdRegisterSize = 128;
constexpr uint32_t kSimdNumVariables = 8;
constexpr uint32_t kLeftShiftSize = 2;
constexpr uint32_t kRightShiftSize = 14;
using __mxxxi = __m128i;

#define _mmxxx_load _mm_load_si128
#define _mmxxx_store _mm_store_si128
#define _mmxxx_and _mm_and_si128
#define _mmxxx_or _mm_or_si128
#define _mmxxx_testz _mm_testz_si128
#define _mmxxx_lshift _mm_slli_si128
#define _mmxxx_rshift _mm_srli_si128

#define _mmxxx_add _mm_add_epi16
#define _mmxxx_sub _mm_sub_epi16
#define _mmxxx_min _mm_min_epi16
#define _mmxxx_max _mm_max_epi16
#define _mmxxx_set1 _mm_set1_epi16
#define _mmxxx_cmpgt _mm_cmpgt_epi16

#endif

inline bool check_if_gt(const __mxxxi& vec1, const __mxxxi& vec2) {
    __mxxxi vcmp = _mmxxx_cmpgt(vec1, vec2);
    return !(_mmxxx_testz(vcmp, vcmp));
}

void print_mxxxi(const __mxxxi& vec) {
    int16_t unpacked[kSimdNumVariables] __attribute__((aligned(kSimdRegisterSize / 8)));
    _mmxxx_store((__mxxxi*) unpacked, vec);
    for (uint32_t i = 0; i < kSimdNumVariables; i++) {
        printf("%d ", unpacked[i]);
    }
}

int16_t max_score_in_mxxxi(const __mxxxi& vec) {
    int16_t max_score = 0;
    int16_t unpacked[kSimdNumVariables] __attribute__((aligned(kSimdRegisterSize / 8)));
    _mmxxx_store((__mxxxi*) unpacked, vec);
    for (uint32_t i = 0; i < kSimdNumVariables; i++) {
        max_score = std::max(max_score, unpacked[i]);
    }
    return max_score;
}

int16_t score_of_element_in_mxxxi(const __mxxxi& vec, uint32_t pos) {
    int16_t unpacked[kSimdNumVariables] __attribute__((aligned(kSimdRegisterSize / 8)));
    _mmxxx_store((__mxxxi*) unpacked, vec);
    return unpacked[pos];
}

/* Taken from https://gcc.gnu.org/viewcvs/gcc?view=revision&revision=216149 */
inline void* align(size_t __align, size_t __size, void*& __ptr, size_t& __space) noexcept {

    const auto __intptr = reinterpret_cast<uintptr_t>(__ptr);
    const auto __aligned = (__intptr - 1u + __align) & -__align;
    const auto __diff = __aligned - __intptr;
    if ((__size + __diff) > __space)
        return nullptr;
    else {
        __space -= __diff;
        return __ptr = reinterpret_cast<void*>(__aligned);
    }
}

struct Cell {
    __mxxxi H;
    __mxxxi F;
};

void simd_align_sequence_to_graph(const std::string& sequence, std::shared_ptr<Graph> graph,
    AlignmentParams params) {

    // init stuff
    uint32_t actual_matrix_width = sequence.size();
    uint32_t matrix_width = actual_matrix_width + (actual_matrix_width % kSimdNumVariables == 0 ? 0 : kSimdNumVariables - actual_matrix_width % kSimdNumVariables);
    uint32_t matrix_height = graph->nodes().size() + 1;

    uint32_t alphabet_size = 256;
    uint32_t num_row_vectors = matrix_width / kSimdNumVariables;

    uint32_t alignment_size = kSimdRegisterSize / 8;

    __mxxxi* P_storage = new __mxxxi[num_row_vectors * alphabet_size + alignment_size - 1];
    void* storage_ptr = (void*) P_storage;
    size_t storage_size = (num_row_vectors * alphabet_size + alignment_size - 1) * sizeof(__mxxxi);
    __mxxxi* P = (__mxxxi*) align(alignment_size, sizeof(__mxxxi) * num_row_vectors * alphabet_size, storage_ptr, storage_size);

    for (const auto& c: graph->alphabet()) {
        std::vector<int16_t> profile(matrix_width, params.mismatch);
        for (uint32_t i = 0; i < sequence.size(); ++i) {
            char s = sequence[i];
            profile[i] = (c == s ? params.match : params.mismatch);
        }

        int16_t aligned_arr[kSimdNumVariables] __attribute__((aligned(kSimdRegisterSize / 8)));
        uint32_t j = 0;
        for (uint32_t i = 0; i < matrix_width; i += kSimdNumVariables) {
            std::copy(profile.begin() + i, profile.begin() + i + kSimdNumVariables, std::begin(aligned_arr));
            P[c * num_row_vectors + j++] = _mmxxx_load((__mxxxi const*) aligned_arr);
        }
    }

    const auto& sorted_nodes_ids = graph->sorted_nodes_ids();

    std::vector<uint32_t> node_id_to_graph_id(sorted_nodes_ids.size());
    for (uint32_t i = 0; i < sorted_nodes_ids.size(); ++i) {
        node_id_to_graph_id[sorted_nodes_ids[i]] = i;
    }

    Cell* matrix_storage = new Cell[num_row_vectors * matrix_height + alignment_size - 1];
    storage_ptr = (void*) matrix_storage;
    storage_size = (num_row_vectors * matrix_height + alignment_size - 1) * sizeof(Cell);
    Cell* matrix = (Cell*) align(alignment_size, sizeof(Cell) * num_row_vectors * matrix_height, storage_ptr, storage_size);

    int16_t big_negative_value = std::numeric_limits<int16_t>::min() + 1000;
    __mxxxi _big_negative = _mmxxx_set1(big_negative_value);
    __mxxxi _zeroes = _mmxxx_set1(0);

    for (uint32_t i = 0; i < num_row_vectors; ++i) {
        matrix[i].F = _big_negative;
    }

    std::vector<int32_t> first_column_values(matrix_height, 0);

    if (params.type == AlignmentType::kNW) {

        for (uint32_t i = 0; i < num_row_vectors; ++i) {
            matrix[i].H = _mmxxx_set1(params.insertion_open + i * kSimdNumVariables * params.insertion_extend);
            __mxxxi ext = _mmxxx_set1(params.insertion_extend);
            for (uint32_t j = 1; j < kSimdNumVariables; ++j) {
                ext = _mmxxx_lshift(ext, kLeftShiftSize);
                matrix[i].H = _mmxxx_add(matrix[i].H, ext);
            }
        }
        for (uint32_t node_id: sorted_nodes_ids) {
            const auto& node = graph->node(node_id);
            uint32_t i = node_id_to_graph_id[node_id] + 1;

            if (node->in_edges().size() == 0) {
                first_column_values[i] = params.deletion_open;
            } else {
                int32_t max_score = big_negative_value;
                for (const auto& edge: node->in_edges()) {
                    uint32_t pred_i = node_id_to_graph_id[edge->begin_node_id()] + 1;
                    max_score = std::max(max_score, first_column_values[pred_i]);
                }
                first_column_values[i] = max_score + params.deletion_extend;
            }
        }
    } else {
        for (uint32_t i = 0; i < num_row_vectors; ++i) {
            matrix[i].H = _zeroes;
        }
    }

    __mxxxi* hlp_storage = new __mxxxi[kSimdNumVariables + alignment_size - 1];
    storage_ptr = (void*) hlp_storage;
    storage_size = (kSimdNumVariables + alignment_size - 1) * sizeof(__mxxxi);
    __mxxxi* hlp = (__mxxxi*) align(alignment_size, sizeof(__mxxxi) * kSimdNumVariables, storage_ptr, storage_size);
    hlp[kSimdNumVariables - 1] = _big_negative;
    for (int32_t i = kSimdNumVariables - 2; i >= 0; --i) {
        hlp[i] = _mmxxx_rshift(hlp[i + 1], kLeftShiftSize);
    }

    // alignment
    params.insertion_open -= params.insertion_extend;
    params.deletion_open -= params.deletion_extend;

    __mxxxi _op = _mmxxx_set1(params.insertion_open);
    __mxxxi _ext = _mmxxx_set1(params.insertion_extend);

    int32_t max_score = params.type == AlignmentType::kNW ? big_negative_value : 0;
    int32_t max_i = -1;
    int32_t max_j = -1;

    uint32_t last_elem_pos = (actual_matrix_width - 1) % kSimdNumVariables;

    for (uint32_t node_id: sorted_nodes_ids) {
        const auto& node = graph->node(node_id);
        uint32_t i = node_id_to_graph_id[node_id] + 1;

        Cell* M_row = &matrix[i * num_row_vectors];
        __mxxxi* P_row = &P[node->letter() * num_row_vectors];

        uint32_t pred_i = node->in_edges().empty() ? 0 :
            node_id_to_graph_id[node->in_edges().front()->begin_node_id()] + 1;

        Cell* M_pred_row = &matrix[pred_i * num_row_vectors];

        __mxxxi _X = _mmxxx_rshift(_mmxxx_set1(first_column_values[pred_i]), kRightShiftSize);

        for (uint32_t j = 0; j < num_row_vectors; ++j) {
            // update F
            M_row[j].F = _mmxxx_add(_mmxxx_max(_mmxxx_add(M_pred_row[j].H, _op), M_pred_row[j].F), _ext);

            // get diagonal
            __mxxxi _T1 = _mmxxx_rshift(M_pred_row[j].H, kRightShiftSize);
            M_row[j].H = _mmxxx_or(_mmxxx_lshift(M_pred_row[j].H, kLeftShiftSize), _X);
            _X = _T1;

            // update H
            M_row[j].H = _mmxxx_max(_mmxxx_add(M_row[j].H, P_row[j]), M_row[j].F);
        }

        // check other predecessors
        for (uint32_t p = 1; p < node->in_edges().size(); ++p) {

            pred_i = node_id_to_graph_id[node->in_edges()[p]->begin_node_id()] + 1;
            Cell* M_pred_row = &matrix[pred_i * num_row_vectors];

            _X = _mmxxx_rshift(_mmxxx_set1(first_column_values[pred_i]), kRightShiftSize);

            for (uint32_t j = 0; j < num_row_vectors; ++j) {
                // update F
                M_row[j].F = _mmxxx_max(M_row[j].F, _mmxxx_add(_mmxxx_max(_mmxxx_add(M_pred_row[j].H, _op), M_pred_row[j].F), _ext));

                // get diagonal
                __mxxxi _T1 = _mmxxx_rshift(M_pred_row[j].H, kRightShiftSize);
                __mxxxi _H = _mmxxx_or(_mmxxx_lshift(M_pred_row[j].H, kLeftShiftSize), _X);
                _X = _T1;

                // updage H
                M_row[j].H = _mmxxx_max(M_row[j].H, _mmxxx_max(_mmxxx_add(_H, P_row[j]), M_row[j].F));
            }
        }

        __mxxxi E = _mmxxx_set1(first_column_values[i]);
        __mxxxi _score = _zeroes;

        for (uint32_t j = 0; j < num_row_vectors; ++j) {
            E = _mmxxx_add(_mmxxx_add(_mmxxx_or(_mmxxx_lshift(M_row[j].H, kLeftShiftSize), _mmxxx_rshift(E, kRightShiftSize)), _op), _ext);

            __mxxxi _T2 = E;
            __mxxxi ext = _ext;
            for (uint32_t k = 0; k < kSimdNumVariables - 1; ++k) {
                ext = _mmxxx_lshift(ext, kLeftShiftSize);
                _T2 = _mmxxx_add(_mmxxx_lshift(_T2, kLeftShiftSize), ext);
                E = _mmxxx_max(E, _mmxxx_or(hlp[k], _T2));
            }

            M_row[j].H = _mmxxx_max(M_row[j].H, E);
            E = _mmxxx_max(M_row[j].H, _mmxxx_sub(E, _op));

            if (params.type == AlignmentType::kSW) {
                M_row[j].H = _mmxxx_max(_zeroes, M_row[j].H);
            }
            _score = _mmxxx_max(_score, M_row[j].H);
        }


        if (params.type == AlignmentType::kSW) {
            int32_t max_row_score = max_score_in_mxxxi(_score);
            if (max_score < max_row_score) {
                max_score = max_row_score;
                max_i = i;
            }
        } else if (params.type == AlignmentType::kOV) {
            int32_t max_row_score = score_of_element_in_mxxxi(M_row[num_row_vectors-1].H, last_elem_pos);
            if (node->out_edges().empty()) {
                max_row_score = std::max(max_row_score, (int32_t) max_score_in_mxxxi(_score));
            }
            if (max_score < max_row_score) {
                max_score = max_row_score;
                max_i = i;
            }
        } else if (params.type == AlignmentType::kNW) {
            if (node->out_edges().empty()) {
                int32_t max_row_score = score_of_element_in_mxxxi(M_row[num_row_vectors-1].H, last_elem_pos);
                if (max_score < max_row_score) {
                    max_score = max_row_score;
                    max_i = i;
                }
            }
        }
    }

    if (max_i == -1 && max_j == -1) { // no alignment found
        printf("Simd  : %d| %d %d\n", max_score, max_i, max_j);
        delete[] hlp_storage;
        delete[] P_storage;
        delete[] matrix_storage;
        return;
    }

    if (params.type == AlignmentType::kSW) {
        Cell* best_row = &matrix[max_i * num_row_vectors];
        for (uint32_t j = 0; j < num_row_vectors; ++j) {
            int16_t unpacked[kSimdNumVariables] __attribute__((aligned(kSimdRegisterSize / 8)));
            _mmxxx_store((__mxxxi*) unpacked, best_row[j].H);
            for (uint32_t k = 0; k < kSimdNumVariables; k++) {
                if (unpacked[k] == max_score) {
                    max_j = j * kSimdNumVariables + k;
                    break;
                }
            }
        }
    } else if (params.type == AlignmentType::kOV) {
        if (graph->node(sorted_nodes_ids[max_i - 1])->out_edges().empty()) {
            Cell* best_row = &matrix[max_i * num_row_vectors];
            for (uint32_t j = 0; j < num_row_vectors; ++j) {
                int16_t unpacked[kSimdNumVariables] __attribute__((aligned(kSimdRegisterSize / 8)));
                _mmxxx_store((__mxxxi*) unpacked, best_row[j].H);
                for (uint32_t k = 0; k < kSimdNumVariables; k++) {
                    if (unpacked[k] == max_score) {
                        max_j = j * kSimdNumVariables + k;
                        break;
                    }
                }
            }
        } else {
            max_j = actual_matrix_width - 1;
        }
    } else if (params.type == AlignmentType::kNW) {
        max_j = actual_matrix_width - 1;
    }

    printf("Simd  : %d| %d %d\n", max_score, max_i, max_j + 1);

    params.deletion_open += params.deletion_extend;
    params.insertion_open += params.insertion_extend;

    // backtrack
    uint32_t max_num_predecessors = 0;
    for (uint32_t i = 0; i < (uint32_t) max_i; ++i) {
        max_num_predecessors = std::max(max_num_predecessors, (uint32_t) graph->node(sorted_nodes_ids[i])->in_edges().size());
    }

    int16_t Hn[kSimdNumVariables] __attribute__((aligned(kSimdRegisterSize / 8)));
    int16_t Hn_pred[kSimdNumVariables * max_num_predecessors] __attribute__((aligned(kSimdRegisterSize / 8)));
    int16_t Hn_diag_pred[kSimdNumVariables * max_num_predecessors] __attribute__((aligned(kSimdRegisterSize / 8)));
    int16_t Fn_pred[kSimdNumVariables * max_num_predecessors] __attribute__((aligned(kSimdRegisterSize / 8)));
    int16_t Pn[kSimdNumVariables] __attribute__((aligned(kSimdRegisterSize / 8)));

    std::vector<uint32_t> predecessors;

    std::vector<int32_t> alignment_node_ids;
    std::vector<int32_t> alignment_seq_ids;

    int32_t i = max_i;
    int32_t j = max_j;
    int32_t prev_i = 0, prev_j = 0;

    uint32_t j_div = j / kSimdNumVariables;
    uint32_t j_mod = j % kSimdNumVariables;

    bool load_next_segment = true;

    do {
        // check stop condition
        if (j == -1 || i == 0) {
            break;
        }

        const auto& node = graph->node(sorted_nodes_ids[i - 1]);
        // load everything
        if (load_next_segment) {
            predecessors.clear();

            // load current cells
            _mmxxx_store((__mxxxi*) Hn, matrix[i * num_row_vectors + j_div].H);
            // load predecessors cells
            if (node->in_edges().empty()) {
                predecessors.emplace_back(0);
                _mmxxx_store((__mxxxi*) Hn_pred, matrix[j_div].H);
                _mmxxx_store((__mxxxi*) Fn_pred, matrix[j_div].F);
            } else {
                uint32_t store_pos = 0;
                for (const auto& edge: node->in_edges()) {
                    predecessors.emplace_back(node_id_to_graph_id[edge->begin_node_id()] + 1);
                    _mmxxx_store((__mxxxi*) &Hn_pred[store_pos * kSimdNumVariables], matrix[predecessors.back() * num_row_vectors + j_div].H);
                    _mmxxx_store((__mxxxi*) &Fn_pred[store_pos * kSimdNumVariables], matrix[predecessors.back() * num_row_vectors + j_div].F);
                    ++store_pos;
                }
            }
            // load query profile cells
            _mmxxx_store((__mxxxi*) Pn, P[node->letter() * num_row_vectors + j_div]);
        }

        // check stop condition
        if (params.type == AlignmentType::kSW && Hn[j_mod] == 0) {
            break;
        }

        if (j_mod == 0) {
            // border case
            if (j_div > 0) {
                for (uint32_t i = 0; i < predecessors.size(); ++i) {
                    _mmxxx_store((__mxxxi*) &Hn_diag_pred[i * kSimdNumVariables], matrix[predecessors[i] * num_row_vectors + (j_div - 1)].H);
                }
            } else {
                for (uint32_t i = 0; i < predecessors.size(); ++i) {
                    Hn_diag_pred[(i + 1) * kSimdNumVariables - 1] = first_column_values[predecessors[i]];
                }
            }
        }

        // find best predecessor cell
        int32_t H_ij = Hn[j_mod];
        bool predecessor_found = false;

        if (i != 0) {
            int32_t match_cost = Pn[j_mod];

            for (uint32_t p = 0; p < predecessors.size(); ++p) {
                if ((j_mod == 0 && H_ij == Hn_diag_pred[(p + 1) * kSimdNumVariables - 1] + match_cost) ||
                    (j_mod != 0 && H_ij == Hn_pred[p * kSimdNumVariables + j_mod - 1] + match_cost)) {
                    prev_i = predecessors[p];
                    prev_j = j - 1;
                    predecessor_found = true;
                    break;
                }
                if ((H_ij == Fn_pred[p * kSimdNumVariables + j_mod] + params.insertion_extend) ||
                    (H_ij == Hn_pred[p * kSimdNumVariables + j_mod] + params.insertion_open)) {
                    prev_i = predecessors[p];
                    prev_j = j;
                    predecessor_found = true;
                    break;
                }
            }
        }

        if (!predecessor_found) {
            prev_i = i;
            prev_j = j - 1;
        }

        alignment_node_ids.emplace_back(i == prev_i ? -1 : sorted_nodes_ids[i - 1]);
        alignment_seq_ids.emplace_back(j == prev_j ? -1 : j);

        // update for next round
        load_next_segment = (i == prev_i ? false : true) || (j != prev_j && prev_j % kSimdNumVariables == kSimdNumVariables - 1 ? true : false );

        i = prev_i;
        j = prev_j;
        j_div = j / kSimdNumVariables;
        j_mod = j % kSimdNumVariables;

    } while (true);

    // update alignment for NW (backtrack stops on first row or column)
    if (params.type == AlignmentType::kNW) {
        while (i == 0 && j != -1) {
            alignment_node_ids.emplace_back(-1);
            alignment_seq_ids.emplace_back(j);
            --j;
        }
        while (i != 0 && j == -1) {
            alignment_node_ids.emplace_back(sorted_nodes_ids[i - 1]);
            alignment_seq_ids.emplace_back(-1);

            const auto& node = graph->node(sorted_nodes_ids[i - 1]);
            if (node->in_edges().empty()) {
                i = 0;
            } else {
                for (const auto& edge: node->in_edges()) {
                    uint32_t pred_i = node_id_to_graph_id[edge->begin_node_id()] + 1;
                    if (first_column_values[i] == first_column_values[pred_i] + params.insertion_extend) {
                        i = pred_i;
                        break;
                    }
                }
            }
        }
    }

    std::reverse(alignment_node_ids.begin(), alignment_node_ids.end());
    std::reverse(alignment_seq_ids.begin(), alignment_seq_ids.end());
    printf("%zu\n", alignment_node_ids.size());

    // free stuff
    delete[] hlp_storage;
    delete[] P_storage;
    delete[] matrix_storage;
}

}
