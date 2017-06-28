/*!
 * @file simd_alignment_engine.cpp
 *
 * @brief SimdAlignmentEngine class source file
 */

#include <limits>
#include <algorithm>
#include <cmath>

extern "C" {
    #include <immintrin.h> // AVX2 and lower
}

#include "spoa/graph.hpp"
#include "simd_alignment_engine.hpp"

namespace spoa {

// Taken from https://gcc.gnu.org/viewcvs/gcc?view=revision&revision=216149
inline void* align(size_t __align, size_t __size, void*& __ptr,
    size_t& __space) noexcept {

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

template<typename T>
T* allocateAlignedMemory(T** storage, uint32_t size, uint32_t alignment) {
    *storage = new T[size + alignment - 1];
    void* ptr = static_cast<void*>(*storage);
    size_t storage_size = (size + alignment - 1) * sizeof(T);
    return static_cast<T*>(align(alignment, size * sizeof(T), ptr, storage_size));
}

template<typename T>
struct InstructionSet;

#if defined(__AVX2__)

constexpr uint32_t kRegisterSize = 256;
constexpr uint32_t kAlignment = kRegisterSize / 8;
using __mxxxi = __m256i;

inline __mxxxi _mmxxx_load_si(__mxxxi const* mem_addr) {
    return _mm256_load_si256(mem_addr);
}

inline void _mmxxx_store_si(__mxxxi* mem_addr, const __mxxxi& a) {
    _mm256_store_si256(mem_addr, a);
}

inline __mxxxi _mmxxx_or_si(const __mxxxi& a, const __mxxxi& b) {
    return _mm256_or_si256(a, b);
}

template<>
struct InstructionSet<int16_t> {
    using type = int16_t;
    static constexpr uint32_t kNumVar = kRegisterSize / (8 * sizeof(type));
    inline __mxxxi _mmxxx_slli_si(const __mxxxi& a) {
        return _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a,
            _MM_SHUFFLE(0, 0, 2, 0)), 16 - sizeof(type));
    }
    inline __mxxxi _mmxxx_srli_si(const __mxxxi& a) {
        return _mm256_srli_si256(_mm256_permute2x128_si256(a, a,
            _MM_SHUFFLE(2, 0, 0, 1)), (kNumVar - 1) * sizeof(type) - 16);
    }
    static inline __mxxxi _mmxxx_add_epi(const __mxxxi& a, const __mxxxi& b) {
        return _mm256_add_epi16(a, b);
    };
    static inline __mxxxi _mmxxx_sub_epi(const __mxxxi& a, const __mxxxi& b) {
        return _mm256_sub_epi16(a, b);
    };
    static inline __mxxxi _mmxxx_min_epi(const __mxxxi& a, const __mxxxi& b) {
        return _mm256_min_epi16(a, b);
    };
    static inline __mxxxi _mmxxx_max_epi(const __mxxxi& a, const __mxxxi& b) {
        return _mm256_max_epi16(a, b);
    };
    static inline __mxxxi _mmxxx_set1_epi(int16_t a) {
        return _mm256_set1_epi16(a);
    }
};

template<>
struct InstructionSet<int32_t> {
    using type = int32_t;
    static constexpr uint32_t kNumVar = kRegisterSize / (8 * sizeof(type));
    inline __mxxxi _mmxxx_slli_si(const __mxxxi& a) {
        return _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a,
            _MM_SHUFFLE(0, 0, 2, 0)), 16 - sizeof(type));
    }
    inline __mxxxi _mmxxx_srli_si(const __mxxxi& a) {
        return _mm256_srli_si256(_mm256_permute2x128_si256(a, a,
            _MM_SHUFFLE(2, 0, 0, 1)), (kNumVar - 1) * sizeof(type) - 16);
    }
    static inline __mxxxi _mmxxx_add_epi(const __mxxxi& a, const __mxxxi& b) {
        return _mm256_add_epi32(a, b);
    };
    static inline __mxxxi _mmxxx_sub_epi(const __mxxxi& a, const __mxxxi& b) {
        return _mm256_sub_epi32(a, b);
    };
    static inline __mxxxi _mmxxx_min_epi(const __mxxxi& a, const __mxxxi& b) {
        return _mm256_min_epi32(a, b);
    };
    static inline __mxxxi _mmxxx_max_epi(const __mxxxi& a, const __mxxxi& b) {
        return _mm256_max_epi32(a, b);
    };
    static inline __mxxxi _mmxxx_set1_epi(int16_t a) {
        return _mm256_set1_epi32(a);
    }
};

#elif defined(__SSE4_1__)

constexpr uint32_t kRegisterSize = 128;
constexpr uint32_t kAlignment = kRegisterSize / 8;
using __mxxxi = __m128i;

inline __mxxxi _mmxxx_load_si(__mxxxi const* mem_addr) {
    return _mm_load_si128(mem_addr);
}

inline void _mmxxx_store_si(__mxxxi* mem_addr, const __mxxxi& a) {
    _mm_store_si128(mem_addr, a);
}

inline __mxxxi _mmxxx_or_si(const __mxxxi& a, const __mxxxi& b) {
    return _mm_or_si128(a, b);
}

template<>
struct InstructionSet<int16_t> {
    using type = int16_t;
    static constexpr uint32_t kNumVar = kRegisterSize / (8 * sizeof(type));
    inline __mxxxi _mmxxx_slli_si(const __mxxxi& a) {
        return _mm_slli_si128(a, sizeof(type));
    }
    inline __mxxxi _mmxxx_srli_si(const __mxxxi& a) {
        return _mm_srli_si128(a, (kNumVar - 1) * sizeof(type));
    }
    static inline __mxxxi _mmxxx_add_epi(const __mxxxi& a, const __mxxxi& b) {
        return _mm_add_epi16(a, b);
    };
    static inline __mxxxi _mmxxx_sub_epi(const __mxxxi& a, const __mxxxi& b) {
        return _mm_sub_epi16(a, b);
    };
    static inline __mxxxi _mmxxx_min_epi(const __mxxxi& a, const __mxxxi& b) {
        return _mm_min_epi16(a, b);
    };
    static inline __mxxxi _mmxxx_max_epi(const __mxxxi& a, const __mxxxi& b) {
        return _mm_max_epi16(a, b);
    };
    static inline __mxxxi _mmxxx_set1_epi(int16_t a) {
        return _mm_set1_epi16(a);
    }
};

template<>
struct InstructionSet<int32_t> {
    using type = int32_t;
    static constexpr uint32_t kNumVar = kRegisterSize / (8 * sizeof(type));
    inline __mxxxi _mmxxx_slli_si(const __mxxxi& a) {
        return _mm_slli_si128(a, sizeof(type));
    }
    inline __mxxxi _mmxxx_srli_si(const __mxxxi& a) {
        return _mm_srli_si128(a, (kNumVar - 1) * sizeof(type));
    }
    static inline __mxxxi _mmxxx_add_epi(const __mxxxi& a, const __mxxxi& b) {
        return _mm_add_epi32(a, b);
    };
    static inline __mxxxi _mmxxx_sub_epi(const __mxxxi& a, const __mxxxi& b) {
        return _mm_sub_epi32(a, b);
    };
    static inline __mxxxi _mmxxx_min_epi(const __mxxxi& a, const __mxxxi& b) {
        return _mm_min_epi32(a, b);
    };
    static inline __mxxxi _mmxxx_max_epi(const __mxxxi& a, const __mxxxi& b) {
        return _mm_max_epi32(a, b);
    };
    static inline __mxxxi _mmxxx_set1_epi(int16_t a) {
        return _mm_set1_epi32(a);
    }
};

#endif

#if defined(__AVX2__) || defined(__SSE4_1__)

template<typename T>
void _mmxxx_print(const __mxxxi& a) {

    typename T::type unpacked[T::kNumVar] __attribute__((aligned(kAlignment)));
    _mmxxx_store_si(reinterpret_cast<__mxxxi*>(unpacked), a);

    for (uint32_t i = 0; i < T::kNumVar; i++) {
        printf("%d ", unpacked[i]);
    }
}

template<typename T>
typename T::type _mmxxx_max_value(const __mxxxi& a) {

    typename T::type max_score = 0;
    typename T::type unpacked[T::kNumVar] __attribute__((aligned(kAlignment)));
    _mmxxx_store_si(reinterpret_cast<__mxxxi*>(unpacked), a);

    for (uint32_t i = 0; i < T::kNumVar; i++) {
        max_score = std::max(max_score, unpacked[i]);
    }

    return max_score;
}

template<typename T>
typename T::type _mmxxx_value_at(const __mxxxi& a, uint32_t i) {

    typename T::type unpacked[T::kNumVar] __attribute__((aligned(kAlignment)));
    _mmxxx_store_si(reinterpret_cast<__mxxxi*>(unpacked), a);

    return unpacked[i];
}

template<typename T>
int32_t _mmxxx_index_of(const __mxxxi* row, uint32_t row_width,
    typename T::type value) {

    for (uint32_t i = 0; i < row_width; ++i) {
        typename T::type unpacked[T::kNumVar] __attribute__((aligned(kAlignment)));
        _mmxxx_store_si(reinterpret_cast<__mxxxi*>(unpacked), row[i]);

        for (uint32_t j = 0; j < T::kNumVar; j++) {
            if (unpacked[j] == value) {
                return i * T::kNumVar + j;
            }
        }
    }

    return -1;
}

#endif

/*
template<typename Simd>
void alignSequenceToGraph(std::vector<std::vector<int32_t>>& sequence_profile, uint32_t sequence_size,
    std::shared_ptr<Graph> graph, AlignmentParams& params, std::vector<int32_t>& alignment_node_ids,
    std::vector<int32_t>& alignment_seq_ids) {

    // init stuff
    uint32_t actual_matrix_width = sequence_size;
    uint32_t matrix_width = actual_matrix_width + (actual_matrix_width % Simd::kNumVariables == 0 ?
        0 : Simd::kNumVariables - actual_matrix_width % Simd::kNumVariables);
    uint32_t matrix_height = graph->nodes().size() + 1;

    uint32_t alphabet_size = 256;
    uint32_t num_row_vectors = matrix_width / Simd::kNumVariables;

    __mxxxi* _P_storage = nullptr;
    __mxxxi* _P = allocateAlignedMemory(&_P_storage, num_row_vectors * alphabet_size, kRegisterSize / 8);

    uint32_t padding_penatly = std::max(params.deletion_open, params.insertion_open); // params.mismatch;
    for (const auto& c: graph->alphabet()) {
        while (sequence_profile[c].size() != matrix_width) {
            sequence_profile[c].emplace_back(padding_penatly);
        }

        typename Simd::type temp[Simd::kNumVariables] __attribute__((aligned(kRegisterSize / 8)));
        uint32_t j = 0;
        for (uint32_t i = 0; i < matrix_width; i += Simd::kNumVariables) {
            std::copy(sequence_profile[c].begin() + i,
                sequence_profile[c].begin() + i + Simd::kNumVariables,
                std::begin(temp));
            _P[c * num_row_vectors + j++] = _mmxxx_load((const __mxxxi*) temp);
        }
    }

    graph->topological_sort();
    const auto& sorted_nodes_ids = graph->sorted_nodes_ids();

    std::vector<uint32_t> node_id_to_graph_id(sorted_nodes_ids.size());
    for (uint32_t i = 0; i < sorted_nodes_ids.size(); ++i) {
        node_id_to_graph_id[sorted_nodes_ids[i]] = i;
    }

    Cell* _matrix_storage = nullptr;
    Cell* _matrix = allocateAlignedMemory(&_matrix_storage, num_row_vectors * matrix_height, kRegisterSize / 8);

    typename Simd::type big_negative_value = std::numeric_limits<typename Simd::type>::min() + 1000;
    __mxxxi _big_negative = Simd::set(big_negative_value);
    __mxxxi _zeroes = Simd::set(0);

    for (uint32_t i = 0; i < num_row_vectors; ++i) {
        _matrix[i].F = _big_negative;
    }

    std::vector<int32_t> first_column_values(matrix_height, 0);

    if (params.type == AlignmentType::kNW || params.type == AlignmentType::kOV) {

        for (uint32_t i = 0; i < num_row_vectors; ++i) {
            _matrix[i].H = Simd::set(params.insertion_open + i * Simd::kNumVariables * params.insertion_extend);
            __mxxxi _ext = Simd::set(params.insertion_extend);
            for (uint32_t j = 1; j < Simd::kNumVariables; ++j) {
                _ext = _mmxxx_lshift(_ext, Simd::kLeftShiftSize);
                _matrix[i].H = Simd::add(_matrix[i].H, _ext);
            }
        }
        if (params.type == AlignmentType::kNW) {
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
        }
    } else {
        for (uint32_t i = 0; i < num_row_vectors; ++i) {
            _matrix[i].H = _zeroes;
        }
    }

    __mxxxi* _masks_storage = nullptr;
    __mxxxi* _masks = allocateAlignedMemory(&_masks_storage, Simd::kNumVariables, kRegisterSize / 8);
    _masks[Simd::kNumVariables - 1] = _big_negative;
    for (int32_t i = Simd::kNumVariables - 2; i >= 0; --i) {
        _masks[i] = _mmxxx_rshift(_masks[i + 1], Simd::kLeftShiftSize);
    }

    // alignment
    params.insertion_open -= params.insertion_extend;
    params.deletion_open -= params.deletion_extend;

    __mxxxi _iop = Simd::set(params.insertion_open);
    __mxxxi _iext = Simd::set(params.insertion_extend);

    __mxxxi _dop = Simd::set(params.deletion_open);
    __mxxxi _dext = Simd::set(params.deletion_extend);

    typename Simd::type max_score = params.type == AlignmentType::kNW ? big_negative_value : 0;
    int32_t max_i = -1;
    int32_t max_j = -1;

    uint32_t last_elem_pos = (actual_matrix_width - 1) % Simd::kNumVariables;

    for (uint32_t node_id: sorted_nodes_ids) {
        const auto& node = graph->node(node_id);
        uint32_t i = node_id_to_graph_id[node_id] + 1;

        Cell* _M_row = &_matrix[i * num_row_vectors];
        __mxxxi* _P_row = &_P[node->letter() * num_row_vectors];

        uint32_t pred_i = node->in_edges().empty() ? 0 :
            node_id_to_graph_id[node->in_edges().front()->begin_node_id()] + 1;

        Cell* _M_pred_row = &_matrix[pred_i * num_row_vectors];

        __mxxxi _X = _mmxxx_rshift(Simd::set(first_column_values[pred_i]), Simd::kRightShiftSize);

        for (uint32_t j = 0; j < num_row_vectors; ++j) {
            // update F
            _M_row[j].F = Simd::add(Simd::max(Simd::add(_M_pred_row[j].H, _iop), _M_pred_row[j].F), _iext);

            // get diagonal
            __mxxxi _T1 = _mmxxx_rshift(_M_pred_row[j].H, Simd::kRightShiftSize);
            _M_row[j].H = _mmxxx_or(_mmxxx_lshift(_M_pred_row[j].H, Simd::kLeftShiftSize), _X);
            _X = _T1;

            // update H
            _M_row[j].H = Simd::max(Simd::add(_M_row[j].H, _P_row[j]), _M_row[j].F);
        }

        // check other predecessors
        for (uint32_t p = 1; p < node->in_edges().size(); ++p) {

            pred_i = node_id_to_graph_id[node->in_edges()[p]->begin_node_id()] + 1;
            Cell* _M_pred_row = &_matrix[pred_i * num_row_vectors];

            _X = _mmxxx_rshift(Simd::set(first_column_values[pred_i]), Simd::kRightShiftSize);

            for (uint32_t j = 0; j < num_row_vectors; ++j) {
                // update F
                _M_row[j].F = Simd::max(_M_row[j].F, Simd::add(Simd::max(Simd::add(_M_pred_row[j].H, _iop), _M_pred_row[j].F), _iext));

                // get diagonal
                __mxxxi _T1 = _mmxxx_rshift(_M_pred_row[j].H, Simd::kRightShiftSize);
                __mxxxi _H = _mmxxx_or(_mmxxx_lshift(_M_pred_row[j].H, Simd::kLeftShiftSize), _X);
                _X = _T1;

                // updage H
                _M_row[j].H = Simd::max(_M_row[j].H, Simd::max(Simd::add(_H, _P_row[j]), _M_row[j].F));
            }
        }

        __mxxxi _E = Simd::set(first_column_values[i]);
        __mxxxi _score = _zeroes;

        for (uint32_t j = 0; j < num_row_vectors; ++j) {
            _E = Simd::add(Simd::add(_mmxxx_or(_mmxxx_lshift(_M_row[j].H, Simd::kLeftShiftSize), _mmxxx_rshift(_E, Simd::kRightShiftSize)), _dop), _dext);

            __mxxxi _T2 = _E;
            __mxxxi _ext = _dext;
            for (uint32_t k = 0; k < Simd::kNumVariables - 1; ++k) {
                _ext = _mmxxx_lshift(_ext, Simd::kLeftShiftSize);
                _T2 = Simd::add(_mmxxx_lshift(_T2, Simd::kLeftShiftSize), _ext);
                _E = Simd::max(_E, _mmxxx_or(_masks[k], _T2));
            }

            _M_row[j].H = Simd::max(_M_row[j].H, _E);
            _E = Simd::max(_M_row[j].H, Simd::sub(_E, _dop));

            if (params.type == AlignmentType::kSW) {
                _M_row[j].H = Simd::max(_zeroes, _M_row[j].H);
            }
            _score = Simd::max(_score, _M_row[j].H);
        }

        if (params.type == AlignmentType::kSW) {
            int32_t max_row_score = maxValueInMxxxi<Simd>(_score);
            if (max_score < max_row_score) {
                max_score = max_row_score;
                max_i = i;
            }
        } else if (params.type == AlignmentType::kOV) {
            // int32_t max_row_score = valueOfElementInMxxxi<Simd>(_M_row[num_row_vectors - 1].H, last_elem_pos);
            if (node->out_edges().empty()) {
                // max_row_score = std::max(max_row_score, (int32_t) maxValueInMxxxi<Simd>(_score));
                int32_t max_row_score = maxValueInMxxxi<Simd>(_score);
                if (max_score < max_row_score) {
                    max_score = max_row_score;
                    max_i = i;
                }
            }
            //if (max_score < max_row_score) {
            //    max_score = max_row_score;
            //    max_i = i;
            //}
        } else if (params.type == AlignmentType::kNW) {
            if (node->out_edges().empty()) {
                int32_t max_row_score = valueOfElementInMxxxi<Simd>(_M_row[num_row_vectors - 1].H, last_elem_pos);
                if (max_score < max_row_score) {
                    max_score = max_row_score;
                    max_i = i;
                }
            }
        }
    }

    if (max_i == -1 && max_j == -1) { // no alignment found
        // printf("Simd  : %d| %d %d\n", max_score, max_i, max_j);
        delete[] _masks_storage;
        delete[] _matrix_storage;
        delete[] _P_storage;
        return;
    }

    if (params.type == AlignmentType::kSW) {
        max_j = findIndexOfValueInRow<Simd>(&_matrix[max_i * num_row_vectors], num_row_vectors, max_score);
    } else if (params.type == AlignmentType::kOV) {
        if (graph->node(sorted_nodes_ids[max_i - 1])->out_edges().empty()) {
            max_j = findIndexOfValueInRow<Simd>(&_matrix[max_i * num_row_vectors], num_row_vectors, max_score);
        } else {
            max_j = actual_matrix_width - 1;
        }
    } else if (params.type == AlignmentType::kNW) {
        max_j = actual_matrix_width - 1;
    }

    // printf("Simd  : %d| %d %d\n", max_score, max_i, max_j + 1);

    params.deletion_open += params.deletion_extend;
    params.insertion_open += params.insertion_extend;

    // backtrack
    uint32_t max_num_predecessors = 0;
    for (uint32_t i = 0; i < (uint32_t) max_i; ++i) {
        max_num_predecessors = std::max(max_num_predecessors, (uint32_t) graph->node(sorted_nodes_ids[i])->in_edges().size());
    }

    typename Simd::type H[Simd::kNumVariables] __attribute__((aligned(kRegisterSize / 8)));
    typename Simd::type H_pred[Simd::kNumVariables * max_num_predecessors] __attribute__((aligned(kRegisterSize / 8)));
    typename Simd::type H_diag_pred[Simd::kNumVariables * max_num_predecessors] __attribute__((aligned(kRegisterSize / 8)));
    typename Simd::type F_pred[Simd::kNumVariables * max_num_predecessors] __attribute__((aligned(kRegisterSize / 8)));
    typename Simd::type P[Simd::kNumVariables] __attribute__((aligned(kRegisterSize / 8)));

    std::vector<uint32_t> predecessors;

    int32_t i = max_i;
    int32_t j = max_j;
    int32_t prev_i = 0, prev_j = 0;

    uint32_t j_div = j / Simd::kNumVariables;
    uint32_t j_mod = j % Simd::kNumVariables;

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
            _mmxxx_store((__mxxxi*) H, _matrix[i * num_row_vectors + j_div].H);
            // load predecessors cells
            if (node->in_edges().empty()) {
                predecessors.emplace_back(0);
                _mmxxx_store((__mxxxi*) H_pred, _matrix[j_div].H);
                _mmxxx_store((__mxxxi*) F_pred, _matrix[j_div].F);
            } else {
                uint32_t store_pos = 0;
                for (const auto& edge: node->in_edges()) {
                    predecessors.emplace_back(node_id_to_graph_id[edge->begin_node_id()] + 1);
                    _mmxxx_store((__mxxxi*) &H_pred[store_pos * Simd::kNumVariables], _matrix[predecessors.back() * num_row_vectors + j_div].H);
                    _mmxxx_store((__mxxxi*) &F_pred[store_pos * Simd::kNumVariables], _matrix[predecessors.back() * num_row_vectors + j_div].F);
                    ++store_pos;
                }
            }
            // load query profile cells
            _mmxxx_store((__mxxxi*) P, _P[node->letter() * num_row_vectors + j_div]);
        }

        // check stop condition
        if (params.type == AlignmentType::kSW && H[j_mod] == 0) {
            break;
        }

        if (j_mod == 0) {
            // border case
            if (j_div > 0) {
                for (uint32_t i = 0; i < predecessors.size(); ++i) {
                    _mmxxx_store((__mxxxi*) &H_diag_pred[i * Simd::kNumVariables], _matrix[predecessors[i] * num_row_vectors + (j_div - 1)].H);
                }
            } else {
                for (uint32_t i = 0; i < predecessors.size(); ++i) {
                    H_diag_pred[(i + 1) * Simd::kNumVariables - 1] = first_column_values[predecessors[i]];
                }
            }
        }

        // find best predecessor cell
        int32_t H_ij = H[j_mod];
        bool predecessor_found = false;

        if (i != 0) {
            int32_t match_cost = P[j_mod];

            for (uint32_t p = 0; p < predecessors.size(); ++p) {
                if ((j_mod == 0 && H_ij == H_diag_pred[(p + 1) * Simd::kNumVariables - 1] + match_cost) ||
                    (j_mod != 0 && H_ij == H_pred[p * Simd::kNumVariables + j_mod - 1] + match_cost)) {
                    prev_i = predecessors[p];
                    prev_j = j - 1;
                    predecessor_found = true;
                    break;
                }
                if ((H_ij == F_pred[p * Simd::kNumVariables + j_mod] + params.insertion_extend) ||
                    (H_ij == H_pred[p * Simd::kNumVariables + j_mod] + params.insertion_open)) {
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
        load_next_segment = (i == prev_i ? false : true) || (j != prev_j && prev_j % Simd::kNumVariables == Simd::kNumVariables - 1 ? true : false );

        i = prev_i;
        j = prev_j;
        j_div = j / Simd::kNumVariables;
        j_mod = j % Simd::kNumVariables;

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

    // free stuff
    delete[] _masks_storage;
    delete[] _matrix_storage;
    delete[] _P_storage;
}

#endif

std::unique_ptr<Alignment> createSimdAlignment(const std::string& sequence,
    std::shared_ptr<Graph> graph, AlignmentParams params) {

#ifdef __SSE4_1__
    assert(sequence.size() != 0);
    return std::unique_ptr<Alignment>(new SimdAlignment(sequence, graph,
        std::move(params)));
#else
    return nullptr;
#endif
}
*/

std::unique_ptr<AlignmentEngine> createSimdAlignmentEngine(
    AlignmentType alignment_type, int8_t match, int8_t mismatch,
    int8_t gap_open, int8_t gap_extend) {

#if defined(__AVX2__) || defined(__SSE4_1__)

    return std::unique_ptr<AlignmentEngine>(new SimdAlignmentEngine(
        alignment_type, match, mismatch, gap_open, gap_extend));

#endif

    return nullptr;
}

struct SimdAlignmentEngine::Implementation {

#if defined(__AVX2__) || defined(__SSE4_1__)

    std::vector<uint32_t> node_id_to_rank;

    std::unique_ptr<__mxxxi[]> sequence_profile_storage;
    __mxxxi* sequence_profile;
    uint32_t sequence_profile_size;

    std::unique_ptr<__mxxxi[]> H_storage;
    __mxxxi* H;

    std::unique_ptr<__mxxxi[]> F_storage;
    __mxxxi* F;

    std::unique_ptr<__mxxxi[]> E_storage;
    __mxxxi* E;

    uint32_t X_size;

    Implementation()
            : node_id_to_rank(), sequence_profile_storage(nullptr),
            sequence_profile(nullptr), sequence_profile_size(0),
            H_storage(nullptr), H(nullptr),
            F_storage(nullptr), F(nullptr),
            E_storage(nullptr), E(nullptr),
            X_size(0) {
        }

#endif
};

SimdAlignmentEngine::SimdAlignmentEngine(AlignmentType alignment_type,
    int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend)
        : AlignmentEngine(alignment_type, match, mismatch, gap_open, gap_extend),
        pimpl_(new Implementation()) {
}

SimdAlignmentEngine::~SimdAlignmentEngine() {
}

Alignment SimdAlignmentEngine::align_sequence_with_graph(
    const std::string& sequence, const std::unique_ptr<Graph>& graph) {

    if (graph->nodes().empty() || sequence.empty()) {
        return Alignment();
    }

#if defined(__AVX2__) || defined(__SSE4_1__)

    uint32_t matrix_width = sequence.size() + InstructionSet<int16_t>::kNumVar;
    uint32_t matrix_height = graph->nodes().size() + 1;

    uint32_t max_short_value = std::numeric_limits<int16_t>::max();
    int8_t max_penalty = std::max(std::max(abs(match_), abs(mismatch_)),
        std::max(abs(gap_open_), abs(gap_extend_)));

    if (max_penalty * (matrix_width + matrix_height) < max_short_value) {
        return align<InstructionSet<int16_t>>(sequence, graph);
    } else {
        return align<InstructionSet<int32_t>>(sequence, graph);
    }

#else

    return Alignment();

#endif
}

template<typename T>
Alignment SimdAlignmentEngine::align(const std::string& sequence,
    const std::unique_ptr<Graph>& graph) {

    return Alignment();
}

/*

    for (const auto& c: graph->alphabet()) {
        sequence_profile_[c].reserve(sequence.size());
        for (const auto& s: sequence) {
            sequence_profile_[c].push_back(c == s ? params_.match : params_.mismatch);
        }
    }
}


// TODO: come up with a more elegant way for this function (its duplicate with SisdAlignment::adjust_node_ids)
void SimdAlignment::adjust_node_ids(const std::vector<int32_t>& mapping) {
    for (uint32_t i = 0; i < alignment_node_ids_.size(); ++i) {
        if (alignment_node_ids_[i] != -1) {
            alignment_node_ids_[i] = mapping[alignment_node_ids_[i]];
        }
    }
}*/

}
