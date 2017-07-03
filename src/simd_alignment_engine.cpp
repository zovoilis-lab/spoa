/*!
 * @file simd_alignment_engine.cpp
 *
 * @brief SimdAlignmentEngine class source file
 */

#include <stdlib.h>
#include <algorithm>
#include <limits>

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

#if 0 && defined(__AVX2__)

constexpr uint32_t kRegisterSize = 256;
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
    static inline __mxxxi _mmxxx_slli_si(const __mxxxi& a) {
        return _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a,
            _MM_SHUFFLE(0, 0, 2, 0)), 14);
    }
    static inline __mxxxi _mmxxx_srli_si(const __mxxxi& a) {
        return _mm256_srli_si256(_mm256_permute2x128_si256(a, a,
            _MM_SHUFFLE(2, 0, 0, 1)), 14);
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
    static inline __mxxxi _mmxxx_set1_epi(type a) {
        return _mm256_set1_epi16(a);
    }
};

template<>
struct InstructionSet<int32_t> {
    using type = int32_t;
    static constexpr uint32_t kNumVar = kRegisterSize / (8 * sizeof(type));
    static inline __mxxxi _mmxxx_slli_si(const __mxxxi& a) {
        return _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a,
            _MM_SHUFFLE(0, 0, 2, 0)), 12);
    }
    static inline __mxxxi _mmxxx_srli_si(const __mxxxi& a) {
        return _mm256_srli_si256(_mm256_permute2x128_si256(a, a,
            _MM_SHUFFLE(2, 0, 0, 1)), 12);
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
    static inline __mxxxi _mmxxx_set1_epi(type a) {
        return _mm256_set1_epi32(a);
    }
};

#elif defined(__SSE4_1__)

constexpr uint32_t kRegisterSize = 128;
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
    static inline __mxxxi _mmxxx_slli_si(const __mxxxi& a) {
        return _mm_slli_si128(a, 2);
    }
    static inline __mxxxi _mmxxx_srli_si(const __mxxxi& a) {
        return _mm_srli_si128(a, 14);
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
    static inline __mxxxi _mmxxx_set1_epi(type a) {
        return _mm_set1_epi16(a);
    }
};

template<>
struct InstructionSet<int32_t> {
    using type = int32_t;
    static constexpr uint32_t kNumVar = kRegisterSize / (8 * sizeof(type));
    static inline __mxxxi _mmxxx_slli_si(const __mxxxi& a) {
        return _mm_slli_si128(a, 4);
    }
    static inline __mxxxi _mmxxx_srli_si(const __mxxxi& a) {
        return _mm_srli_si128(a, 12);
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
    static inline __mxxxi _mmxxx_set1_epi(type a) {
        return _mm_set1_epi32(a);
    }
};

#endif

#if defined(__AVX2__) || defined(__SSE4_1__)

template<typename T>
void _mmxxx_print(const __mxxxi& a) {

    __attribute__((aligned(kRegisterSize / 8))) typename T::type
        unpacked[T::kNumVar];
    _mmxxx_store_si(reinterpret_cast<__mxxxi*>(unpacked), a);

    for (uint32_t i = 0; i < T::kNumVar; i++) {
        printf("%d ", unpacked[i]);
    }
}

template<typename T>
typename T::type _mmxxx_max_value(const __mxxxi& a) {

    typename T::type max_score = 0;
    __attribute__((aligned(kRegisterSize / 8))) typename T::type
        unpacked[T::kNumVar];
    _mmxxx_store_si(reinterpret_cast<__mxxxi*>(unpacked), a);

    for (uint32_t i = 0; i < T::kNumVar; i++) {
        max_score = std::max(max_score, unpacked[i]);
    }

    return max_score;
}

template<typename T>
typename T::type _mmxxx_value_at(const __mxxxi& a, uint32_t i) {

    __attribute__((aligned(kRegisterSize / 8))) typename T::type
        unpacked[T::kNumVar];
    _mmxxx_store_si(reinterpret_cast<__mxxxi*>(unpacked), a);

    return unpacked[i];
}

template<typename T>
int32_t _mmxxx_index_of(const __mxxxi* row, uint32_t row_width,
    typename T::type value) {

    for (uint32_t i = 0; i < row_width; ++i) {
        __attribute__((aligned(kRegisterSize / 8))) typename T::type
            unpacked[T::kNumVar];
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

std::unique_ptr<AlignmentEngine> createSimdAlignmentEngine(
    AlignmentType alignment_type, int8_t match, int8_t mismatch,
    int8_t gap_open, int8_t gap_extend) {

    return nullptr;

#if defined(__AVX2__) || defined(__SSE4_1__)

    return std::unique_ptr<AlignmentEngine>(new SimdAlignmentEngine(
        alignment_type, match, mismatch, gap_open, gap_extend));

#else

    return nullptr;

#endif
}

struct SimdAlignmentEngine::Implementation {

#if defined(__AVX2__) || defined(__SSE4_1__)

    std::vector<uint32_t> node_id_to_rank;

    std::unique_ptr<__mxxxi[]> sequence_profile_storage;
    uint32_t sequence_profile_size;
    __mxxxi* sequence_profile;

    std::vector<int32_t> first_column;
    std::unique_ptr<__mxxxi[]> X_storage;
    uint32_t X_size;
    __mxxxi* H;
    __mxxxi* F;
    __mxxxi* E;

    std::unique_ptr<__mxxxi[]> masks_storage;
    uint32_t masks_size;
    __mxxxi* masks;

    Implementation()
            : node_id_to_rank(), sequence_profile_storage(nullptr),
            sequence_profile_size(0), sequence_profile(nullptr), first_column(),
            X_storage(nullptr), X_size(0), H(nullptr), F(nullptr), E(nullptr),
            masks_storage(nullptr), masks_size(0), masks(nullptr) {
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

void SimdAlignmentEngine::prealloc(uint32_t max_sequence_size,
    uint32_t alphabet_size) {

#if defined(__AVX2__) || defined(__SSE4_1__)

    uint32_t longest_path = max_sequence_size * (alphabet_size + 1) + 1 +
        InstructionSet<int16_t>::kNumVar;

    uint32_t max_penalty = std::max(std::max(abs(match_), abs(mismatch_)),
        std::max(abs(gap_open_), abs(gap_extend_)));

    if (max_penalty * longest_path < std::numeric_limits<int16_t>::max()) {
        this->realloc((max_sequence_size / InstructionSet<int16_t>::kNumVar) + 1,
            alphabet_size * max_sequence_size, alphabet_size);
    } else {
        this->realloc((max_sequence_size / InstructionSet<int32_t>::kNumVar) + 1,
            alphabet_size * max_sequence_size, alphabet_size);
    }

#endif
}

void SimdAlignmentEngine::realloc(uint32_t matrix_width, uint32_t matrix_height,
    uint32_t num_codes) {

#if defined(__AVX2__) || defined(__SSE4_1__)

    if (pimpl_->node_id_to_rank.size() < matrix_height - 1) {
        pimpl_->node_id_to_rank.resize(matrix_height - 1, 0);
    }
    if (pimpl_->sequence_profile_size < num_codes * matrix_width) {
        __mxxxi* storage = nullptr;
        pimpl_->sequence_profile_size = num_codes * matrix_width;
        pimpl_->sequence_profile = allocateAlignedMemory(&storage,
            pimpl_->sequence_profile_size, kRegisterSize / 8);
        pimpl_->sequence_profile_storage.reset();
        pimpl_->sequence_profile_storage = std::unique_ptr<__mxxxi[]>(storage);
    }
    if (pimpl_->first_column.size() < matrix_height) {
        pimpl_->first_column.resize(matrix_height, 0);
    }
    if (gap_open_ == gap_extend_) {
        if (pimpl_->X_size < matrix_height * matrix_width) {
            __mxxxi* storage = nullptr;
            pimpl_->X_size = matrix_height * matrix_width;
            pimpl_->H = allocateAlignedMemory(&storage, pimpl_->X_size,
                kRegisterSize / 8);
            pimpl_->F = nullptr;
            pimpl_->E = nullptr;
            pimpl_->X_storage.reset();
            pimpl_->X_storage = std::unique_ptr<__mxxxi[]>(storage);
        }
    } else {
        if (pimpl_->X_size < 3 * matrix_height * matrix_width) {
            __mxxxi* storage = nullptr;
            pimpl_->X_size = 3 * matrix_height * matrix_width;
            pimpl_->H = allocateAlignedMemory(&storage, pimpl_->X_size,
                kRegisterSize / 8);
            pimpl_->F = &(pimpl_->H[matrix_height * matrix_width]);
            pimpl_->E = &(pimpl_->F[matrix_height * matrix_width]);
            pimpl_->X_storage.reset();
            pimpl_->X_storage = std::unique_ptr<__mxxxi[]>(storage);
        }
    }
    if (pimpl_->masks_size < InstructionSet<int16_t>::kNumVar) {
        __mxxxi* storage = nullptr;
        pimpl_->masks_size = InstructionSet<int16_t>::kNumVar;
        pimpl_->masks = allocateAlignedMemory(&storage,
            pimpl_->masks_size, kRegisterSize / 8);
        pimpl_->masks_storage.reset();
        pimpl_->masks_storage = std::unique_ptr<__mxxxi[]>(storage);
    }

#endif
}

template<typename T>
void SimdAlignmentEngine::initialize(const std::string& sequence,
    const std::unique_ptr<Graph>& graph, uint32_t normal_matrix_width,
    uint32_t matrix_width, uint32_t matrix_height) {

#if defined(__AVX2__) || defined(__SSE4_1__)

    int32_t padding_penatly = -1 * std::max(std::max(abs(match_), abs(mismatch_)),
        std::max(abs(gap_open_), abs(gap_extend_)));

    __attribute__((aligned(kRegisterSize / 8))) typename T::type
        unpacked[T::kNumVar] = {};

    for (uint32_t i = 0; i < graph->num_codes(); ++i) {
        char c = graph->decoder(i);
        for (uint32_t j = 0; j < matrix_width; ++j) {
            for (uint32_t k = 0; k < T::kNumVar; ++k) {
                unpacked[k] = (j * T::kNumVar + k) < normal_matrix_width ?
                    (c == sequence[j * T::kNumVar + k] ? match_ : mismatch_) :
                    padding_penatly;
            }
            pimpl_->sequence_profile[i * matrix_width + j] =
                _mmxxx_load_si(reinterpret_cast<const __mxxxi*>(unpacked));
        }
    }

    const auto& sorted_nodes_ids = graph->sorted_nodes_ids();

    for (uint32_t i = 0; i < sorted_nodes_ids.size(); ++i) {
        pimpl_->node_id_to_rank[sorted_nodes_ids[i]] = i;
    }

    typename T::type negative_infinity =
        std::numeric_limits<typename T::type>::min() + 1024;

    __mxxxi zeroes = T::_mmxxx_set1_epi(0);

    if (alignment_type_ == AlignmentType::kSW) {
        for (uint32_t i = 0; i < matrix_height; ++i) {
            pimpl_->first_column[i] = 0;
        }
        for (uint32_t j = 0; j < matrix_width; ++j) {
            pimpl_->H[j] = zeroes;
        }
    }

    if (alignment_type_ == AlignmentType::kNW) {
        pimpl_->first_column[0] = 0;
        for (const auto& node_id: sorted_nodes_ids) {
            uint32_t i = pimpl_->node_id_to_rank[node_id] + 1;
            const auto& node = graph->nodes()[node_id];
            if (node->in_edges().empty()) {
                pimpl_->first_column[i] = gap_open_;
            } else {
                int32_t penalty = negative_infinity;
                for (const auto& edge: node->in_edges()) {
                    uint32_t pred_i =
                        pimpl_->node_id_to_rank[edge->begin_node_id()] + 1;
                    penalty = std::max(penalty, pimpl_->first_column[pred_i]);
                }
                pimpl_->first_column[i] = penalty + gap_extend_;
            }
        }

        for (uint32_t j = 0; j < matrix_width; ++j) {
            pimpl_->H[j] = T::_mmxxx_set1_epi(gap_open_ +
                j * T::kNumVar * gap_extend_);
            __mxxxi ext = T::_mmxxx_set1_epi(gap_extend_);

            for (uint32_t k = 1; k < T::kNumVar; ++k) {
                ext = T::_mmxxx_slli_si(ext);
                pimpl_->H[j] = T::_mmxxx_add_epi(pimpl_->H[j], ext);
            }
        }
    }

    if (alignment_type_ == AlignmentType::kOV) {
        for (uint32_t i = 0; i < matrix_height; ++i) {
            pimpl_->first_column[i] = 0;
        }

        for (uint32_t j = 0; j < matrix_width; ++j) {
            pimpl_->H[j] = T::_mmxxx_set1_epi(gap_open_ +
                j * T::kNumVar * gap_extend_);
            __mxxxi ext = T::_mmxxx_set1_epi(gap_extend_);

            for (uint32_t k = 1; k < T::kNumVar; ++k) {
                ext = T::_mmxxx_slli_si(ext);
                pimpl_->H[j] = T::_mmxxx_add_epi(pimpl_->H[j], ext);
            }
        }
    }

#endif
}

Alignment SimdAlignmentEngine::align_sequence_with_graph(
    const std::string& sequence, const std::unique_ptr<Graph>& graph) {

    if (graph->nodes().empty() || sequence.empty()) {
        return Alignment();
    }

#if defined(__AVX2__) || defined(__SSE4_1__)

    uint32_t longest_path = graph->nodes().size() + 1 + sequence.size() +
        InstructionSet<int16_t>::kNumVar;

    uint32_t max_penalty = std::max(std::max(abs(match_), abs(mismatch_)),
        std::max(abs(gap_open_), abs(gap_extend_)));

    if (max_penalty * longest_path < std::numeric_limits<int16_t>::max()) {
        if (gap_open_ == gap_extend_) {
            return align_normal<InstructionSet<int16_t>>(sequence, graph);
        }
        return align_gotoh<InstructionSet<int16_t>>(sequence, graph);
    } else {
        if (gap_open_ == gap_extend_) {
            return align_normal<InstructionSet<int32_t>>(sequence, graph);
        }
        return align_gotoh<InstructionSet<int32_t>>(sequence, graph);
    }

#else

    return Alignment();

#endif
}

template<typename T>
Alignment SimdAlignmentEngine::align_normal(const std::string& sequence,
    const std::unique_ptr<Graph>& graph) {

#if defined(__AVX2__) || defined(__SSE4_1__)

    uint32_t normal_matrix_width = sequence.size();
    uint32_t matrix_width = (sequence.size() + (sequence.size() % T::kNumVar == 0 ?
        0 : T::kNumVar - sequence.size() % T::kNumVar)) / T::kNumVar;
    uint32_t matrix_height = graph->nodes().size() + 1;
    const auto& sorted_nodes_ids = graph->sorted_nodes_ids();

    // realloc
    this->realloc(matrix_width, matrix_height, graph->num_codes());

    // initialize
    this->initialize<T>(sequence, graph, normal_matrix_width, matrix_width,
        matrix_height);

    typename T::type negative_infinity =
        std::numeric_limits<typename T::type>::min() + 1024;

    __attribute__((aligned(kRegisterSize / 8))) typename T::type
        unpacked[T::kNumVar] = {};

    for (uint32_t i = 0; i < T::kNumVar; ++i) {
        unpacked[i] = 0;
    }
    for (uint32_t i = 0; i < T::kNumVar; ++i) {
        pimpl_->masks[i] =
            _mmxxx_load_si(reinterpret_cast<const __mxxxi*>(unpacked));
        unpacked[i] = negative_infinity;
    }

    typename T::type max_score = alignment_type_ == AlignmentType::kSW ? 0 :
        negative_infinity;
    int32_t max_i = -1;
    int32_t max_j = -1;
    uint32_t last_column_id = (normal_matrix_width - 1) % T::kNumVar;
    __mxxxi zeroes = T::_mmxxx_set1_epi(0);
    __mxxxi opn = T::_mmxxx_set1_epi(gap_open_);

    // alignment
    for (uint32_t node_id: sorted_nodes_ids) {
        const auto& node = graph->nodes()[node_id];
        __mxxxi* char_profile =
            &(pimpl_->sequence_profile[node->code() * matrix_width]);

        uint32_t i = pimpl_->node_id_to_rank[node_id] + 1;
        uint32_t pred_i = node->in_edges().empty() ? 0 :
            pimpl_->node_id_to_rank[node->in_edges()[0]->begin_node_id()] + 1;

        __mxxxi* H_row = &(pimpl_->H[i * matrix_width]);
        __mxxxi* H_pred_row = &(pimpl_->H[pred_i * matrix_width]);

        __mxxxi x = T::_mmxxx_srli_si(T::_mmxxx_set1_epi(
            pimpl_->first_column[pred_i]));

        for (uint32_t j = 0; j < matrix_width; ++j) {
            // get diagonal
            __mxxxi t1 = T::_mmxxx_srli_si(H_pred_row[j]);
            H_row[j] = _mmxxx_or_si(T::_mmxxx_slli_si(H_pred_row[j]), x);
            x = t1;

            // update H
            H_row[j] = T::_mmxxx_max_epi(T::_mmxxx_add_epi(H_row[j],
                char_profile[j]), T::_mmxxx_add_epi(H_pred_row[j], opn));
        }

        // check other predecessors
        for (uint32_t p = 1; p < node->in_edges().size(); ++p) {
            pred_i = pimpl_->node_id_to_rank[node->in_edges()[p]->begin_node_id()] + 1;

            H_pred_row = &(pimpl_->H[pred_i * matrix_width]);

            x = T::_mmxxx_srli_si(T::_mmxxx_set1_epi(
                pimpl_->first_column[pred_i]));

            for (uint32_t j = 0; j < matrix_width; ++j) {
                // get diagonal
                __mxxxi t1 = T::_mmxxx_srli_si(H_pred_row[j]);
                __mxxxi h = _mmxxx_or_si(T::_mmxxx_slli_si(H_pred_row[j]), x);
                x = t1;

                // updage H
                H_row[j] = T::_mmxxx_max_epi(H_row[j], T::_mmxxx_max_epi(
                    T::_mmxxx_add_epi(h, char_profile[j]),
                    T::_mmxxx_add_epi(H_pred_row[j], opn)));
            }
        }

        __mxxxi score = T::_mmxxx_set1_epi(negative_infinity);
        x = T::_mmxxx_set1_epi(pimpl_->first_column[i]);

        for (uint32_t j = 0; j < matrix_width; ++j) {
            __mxxxi t1 = T::_mmxxx_add_epi(_mmxxx_or_si(T::_mmxxx_slli_si(
                H_row[j]), T::_mmxxx_srli_si(x)), opn);
            __mxxxi t2 = opn;

            for (uint32_t k = 0; k < T::kNumVar; ++k) {
                H_row[j] = T::_mmxxx_max_epi(H_row[j], _mmxxx_or_si(t1,
                    pimpl_->masks[k]));
                t2 = T::_mmxxx_slli_si(t2);
                t1 = T::_mmxxx_add_epi(T::_mmxxx_slli_si(t1), t2);
            }

            x = H_row[j];

            if (alignment_type_ == AlignmentType::kSW) {
                H_row[j] = T::_mmxxx_max_epi(H_row[j], zeroes);
            }
            score = T::_mmxxx_max_epi(score, H_row[j]);
        }

        if (alignment_type_ == AlignmentType::kSW) {
            int32_t max_row_score = _mmxxx_max_value<T>(score);
            if (max_score < max_row_score) {
                max_score = max_row_score;
                max_i = i;
            }

        } else if (alignment_type_ == AlignmentType::kOV) {
            if (node->out_edges().empty()) {
                int32_t max_row_score = _mmxxx_max_value<T>(score);
                if (max_score < max_row_score) {
                    max_score = max_row_score;
                    max_i = i;
                }
            }

        } else if (alignment_type_ == AlignmentType::kNW) {
            if (node->out_edges().empty()) {
                int32_t max_row_score = _mmxxx_value_at<T>(
                    H_row[matrix_width - 1], last_column_id);
                if (max_score < max_row_score) {
                    max_score = max_row_score;
                    max_i = i;
                }
            }
        }
    }

    if (max_i == -1 && max_j == -1) { // no alignment found
        return Alignment();
    }

    if (alignment_type_ == AlignmentType::kSW) {
        max_j = _mmxxx_index_of<T>(&(pimpl_->H[max_i * matrix_width]),
            matrix_width, max_score);

    } else if (alignment_type_ == AlignmentType::kOV) {
        if (graph->nodes()[sorted_nodes_ids[max_i - 1]]->out_edges().empty()) {
            max_j = _mmxxx_index_of<T>(&(pimpl_->H[max_i * matrix_width]),
                matrix_width, max_score);
        } else {
            max_j = normal_matrix_width - 1;
        }

    } else if (alignment_type_ == AlignmentType::kNW) {
        max_j = normal_matrix_width - 1;
    }

    // backtrack
    uint32_t max_num_predecessors = 0;
    for (uint32_t i = 0; i < (uint32_t) max_i; ++i) {
        max_num_predecessors = std::max(max_num_predecessors,
            (uint32_t) graph->nodes()[sorted_nodes_ids[i]]->in_edges().size());
    }

    __attribute__((aligned(kRegisterSize / 8))) typename T::type H[T::kNumVar];
    __attribute__((aligned(kRegisterSize / 8))) typename T::type
        H_pred[T::kNumVar * max_num_predecessors];
    __attribute__((aligned(kRegisterSize / 8))) typename T::type
        H_diag_pred[T::kNumVar * max_num_predecessors];
    __attribute__((aligned(kRegisterSize / 8))) typename T::type
        H_left_pred[T::kNumVar];
    __attribute__((aligned(kRegisterSize / 8))) typename T::type
        profile[T::kNumVar];

    std::vector<uint32_t> predecessors;

    int32_t i = max_i;
    int32_t j = max_j;
    int32_t prev_i = 0, prev_j = 0;

    uint32_t j_div = j / T::kNumVar;
    uint32_t j_mod = j % T::kNumVar;

    bool load_next_segment = true;

    Alignment alignment;

    do {
        // check stop condition
        if (j == -1 || i == 0) {
            break;
        }

        const auto& node = graph->nodes()[sorted_nodes_ids[i - 1]];
        // load everything
        if (load_next_segment) {
            predecessors.clear();

            // load current cells
            _mmxxx_store_si(reinterpret_cast<__mxxxi*>(H),
                pimpl_->H[i * matrix_width + j_div]);

            // load predecessors cells
            if (node->in_edges().empty()) {
                predecessors.emplace_back(0);
                _mmxxx_store_si(reinterpret_cast<__mxxxi*>(H_pred),
                    pimpl_->H[j_div]);

            } else {
                uint32_t store_pos = 0;
                for (const auto& edge: node->in_edges()) {
                    predecessors.emplace_back(
                        pimpl_->node_id_to_rank[edge->begin_node_id()] + 1);
                    _mmxxx_store_si(
                        reinterpret_cast<__mxxxi*>(&H_pred[store_pos * T::kNumVar]),
                        pimpl_->H[predecessors.back() * matrix_width + j_div]);
                    ++store_pos;
                }
            }

            // load query profile cells
            _mmxxx_store_si(reinterpret_cast<__mxxxi*>(profile),
                pimpl_->sequence_profile[node->code() * matrix_width + j_div]);
        }

        // check stop condition
        if (alignment_type_ == AlignmentType::kSW && H[j_mod] == 0) {
            break;
        }

        if (j_mod == 0) {
            // border case
            if (j_div > 0) {
                _mmxxx_store_si(reinterpret_cast<__mxxxi*>(H_left_pred),
                    pimpl_->H[i * matrix_width + j_div - 1]);

                for (uint32_t p = 0; p < predecessors.size(); ++p) {
                    _mmxxx_store_si(
                        reinterpret_cast<__mxxxi*>(&H_diag_pred[p * T::kNumVar]),
                        pimpl_->H[predecessors[p] * matrix_width + (j_div - 1)]);
                }
            } else {
                H_left_pred[T::kNumVar - 1] = pimpl_->first_column[i];

                for (uint32_t p = 0; p < predecessors.size(); ++p) {
                    H_diag_pred[(p + 1) * T::kNumVar - 1] =
                        pimpl_->first_column[predecessors[p]];
                }
            }
        }

        // find best predecessor cell
        bool predecessor_found = false;

        if (i != 0) {
            for (uint32_t p = 0; p < predecessors.size(); ++p) {
                if (H[j_mod] == H_pred[p * T::kNumVar + j_mod] + gap_open_) {
                    prev_i = predecessors[p];
                    prev_j = j;
                    predecessor_found = true;
                    break;
                }
            }
        }

        if (!predecessor_found) {
            if ((j_mod == 0 && H[j_mod] == H_left_pred[T::kNumVar - 1] + gap_open_) ||
                (j_mod != 0 && H[j_mod] == H[j_mod - 1] + gap_open_)) {
                prev_i = i;
                prev_j = j - 1;
                predecessor_found = true;
            }
        }

        if (!predecessor_found && i != 0) {
            for (uint32_t p = 0; p < predecessors.size(); ++p) {
                if ((j_mod == 0 && H[j_mod] ==
                        H_diag_pred[(p + 1) * T::kNumVar - 1] + profile[j_mod]) ||
                    (j_mod != 0 && H[j_mod] ==
                        H_pred[p * T::kNumVar + j_mod - 1] + profile[j_mod])) {

                    prev_i = predecessors[p];
                    prev_j = j - 1;
                    predecessor_found = true;
                    break;
                }
            }
        }

        alignment.emplace_back(i == prev_i ? -1 : sorted_nodes_ids[i - 1],
            j == prev_j ? -1 : j);

        // update for next round
        load_next_segment = (i == prev_i ? false : true) ||
            (j != prev_j && prev_j % T::kNumVar == T::kNumVar - 1 ? true : false);

        i = prev_i;
        j = prev_j;
        j_div = j / T::kNumVar;
        j_mod = j % T::kNumVar;

    } while (true);

    // update alignment for NW (backtrack stops on first row or column)
    if (alignment_type_ == AlignmentType::kNW) {
        while (i == 0 && j != -1) {
            alignment.emplace_back(-1, j);
            --j;
        }
        while (i != 0 && j == -1) {
            alignment.emplace_back(sorted_nodes_ids[i - 1], -1);

            const auto& node = graph->nodes()[sorted_nodes_ids[i - 1]];
            if (node->in_edges().empty()) {
                i = 0;
            } else {
                for (const auto& edge: node->in_edges()) {
                    uint32_t pred_i =
                        pimpl_->node_id_to_rank[edge->begin_node_id()] + 1;
                    if (pimpl_->first_column[i] ==
                        pimpl_->first_column[pred_i] + gap_open_) {
                        i = pred_i;
                        break;
                    }
                }
            }
        }
    }

    std::reverse(alignment.begin(), alignment.end());
    return alignment;

#else

    return Alignment();

#endif
}

template<typename T>
Alignment SimdAlignmentEngine::align_gotoh(const std::string& sequence,
    const std::unique_ptr<Graph>& graph) {

#if defined(__AVX2__) || defined(__SSE4_1__)

    uint32_t normal_matrix_width = sequence.size();
    uint32_t matrix_width = (sequence.size() + (sequence.size() % T::kNumVar == 0 ?
        0 : T::kNumVar - sequence.size() % T::kNumVar)) / T::kNumVar;
    uint32_t matrix_height = graph->nodes().size() + 1;
    const auto& sorted_nodes_ids = graph->sorted_nodes_ids();

    // realloc
    this->realloc(matrix_width, matrix_height, graph->num_codes());

    // initialize
    this->initialize<T>(sequence, graph, normal_matrix_width, matrix_width,
        matrix_height);

    __attribute__((aligned(kRegisterSize / 8))) typename T::type
        unpacked[T::kNumVar] = {};

    typename T::type negative_infinity =
        std::numeric_limits<typename T::type>::min() + 1024;

    __mxxxi negative_infinities = T::_mmxxx_set1_epi(negative_infinity);

    for (uint32_t j = 0; j < matrix_width; ++j) {
        pimpl_->F[j] = negative_infinities;
    }

    for (uint32_t i = 0; i < T::kNumVar; ++i) {
        unpacked[i] = 0;
    }
    for (uint32_t i = 0; i < T::kNumVar; ++i) {
        unpacked[i] = negative_infinity;
        pimpl_->masks[i] =
            _mmxxx_load_si(reinterpret_cast<const __mxxxi*>(unpacked));
    }

    typename T::type max_score = alignment_type_ == AlignmentType::kSW ? 0 :
        negative_infinity;
    int32_t max_i = -1;
    int32_t max_j = -1;
    uint32_t last_column_id = (normal_matrix_width - 1) % T::kNumVar;
    __mxxxi zeroes = T::_mmxxx_set1_epi(0);
    __mxxxi opn = T::_mmxxx_set1_epi(gap_open_ - gap_extend_);
    __mxxxi ext = T::_mmxxx_set1_epi(gap_extend_);

    // alignment
    for (uint32_t node_id: sorted_nodes_ids) {
        const auto& node = graph->nodes()[node_id];
        __mxxxi* char_profile =
            &(pimpl_->sequence_profile[node->code() * matrix_width]);

        uint32_t i = pimpl_->node_id_to_rank[node_id] + 1;

        __mxxxi* H_row = &(pimpl_->H[i * matrix_width]);
        __mxxxi* F_row = &(pimpl_->F[i * matrix_width]);

        uint32_t pred_i = node->in_edges().empty() ? 0 :
            pimpl_->node_id_to_rank[node->in_edges()[0]->begin_node_id()] + 1;

        __mxxxi* H_pred_row = &(pimpl_->H[pred_i * matrix_width]);
        __mxxxi* F_pred_row = &(pimpl_->F[pred_i * matrix_width]);

        __mxxxi x = T::_mmxxx_srli_si(T::_mmxxx_set1_epi(
            pimpl_->first_column[pred_i]));

        for (uint32_t j = 0; j < matrix_width; ++j) {
            // update F
            F_row[j] = T::_mmxxx_add_epi(T::_mmxxx_max_epi(T::_mmxxx_add_epi(
                H_pred_row[j], opn), F_pred_row[j]), ext);

            // get diagonal
            __mxxxi t1 = T::_mmxxx_srli_si(H_pred_row[j]);
            H_row[j] = _mmxxx_or_si(T::_mmxxx_slli_si(H_pred_row[j]), x);
            x = t1;

            // update H
            H_row[j] = T::_mmxxx_max_epi(T::_mmxxx_add_epi(H_row[j],
                char_profile[j]), F_row[j]);
        }

        // check other predecessors
        for (uint32_t p = 1; p < node->in_edges().size(); ++p) {
            pred_i = pimpl_->node_id_to_rank[node->in_edges()[p]->begin_node_id()] + 1;

            H_pred_row = &(pimpl_->H[pred_i * matrix_width]);
            F_pred_row = &(pimpl_->F[pred_i * matrix_width]);

            x = T::_mmxxx_srli_si(T::_mmxxx_set1_epi(
                pimpl_->first_column[pred_i]));

            for (uint32_t j = 0; j < matrix_width; ++j) {
                // update F
                F_row[j] = T::_mmxxx_max_epi(F_row[j], T::_mmxxx_add_epi(
                    T::_mmxxx_max_epi(T::_mmxxx_add_epi(H_pred_row[j],
                    opn), F_pred_row[j]), ext));

                // get diagonal
                __mxxxi t1 = T::_mmxxx_srli_si(H_pred_row[j]);
                __mxxxi h = _mmxxx_or_si(T::_mmxxx_slli_si(H_pred_row[j]), x);
                x = t1;

                // updage H
                H_row[j] = T::_mmxxx_max_epi(H_row[j], T::_mmxxx_max_epi(
                    T::_mmxxx_add_epi(h, char_profile[j]), F_row[j]));
            }
        }

        __mxxxi* E_row = &(pimpl_->E[i * matrix_width]);
        __mxxxi score = negative_infinities;
        x = T::_mmxxx_set1_epi(pimpl_->first_column[i]);

        for (uint32_t j = 0; j < matrix_width; ++j) {
            E_row[j] = T::_mmxxx_add_epi(T::_mmxxx_add_epi(_mmxxx_or_si(
                T::_mmxxx_slli_si(H_row[j]), T::_mmxxx_srli_si(x)), opn), ext);

            __mxxxi t1 = E_row[j];
            __mxxxi t2 = ext;
            for (uint32_t k = 0; k < T::kNumVar - 1; ++k) {
                t2 = T::_mmxxx_slli_si(t2);
                t1 = T::_mmxxx_add_epi(T::_mmxxx_slli_si(t1), t2);
                E_row[j] = T::_mmxxx_max_epi(E_row[j], _mmxxx_or_si(t1,
                    pimpl_->masks[k]));
            }

            H_row[j] = T::_mmxxx_max_epi(H_row[j], E_row[j]);
            x = T::_mmxxx_max_epi(H_row[j], T::_mmxxx_sub_epi(E_row[j], opn));

            if (alignment_type_ == AlignmentType::kSW) {
                H_row[j] = T::_mmxxx_max_epi(H_row[j], zeroes);
            }
            score = T::_mmxxx_max_epi(score, H_row[j]);
        }

        if (alignment_type_ == AlignmentType::kSW) {
            int32_t max_row_score = _mmxxx_max_value<T>(score);
            if (max_score < max_row_score) {
                max_score = max_row_score;
                max_i = i;
            }

        } else if (alignment_type_ == AlignmentType::kOV) {
            if (node->out_edges().empty()) {
                int32_t max_row_score = _mmxxx_max_value<T>(score);
                if (max_score < max_row_score) {
                    max_score = max_row_score;
                    max_i = i;
                }
            }

        } else if (alignment_type_ == AlignmentType::kNW) {
            if (node->out_edges().empty()) {
                int32_t max_row_score = _mmxxx_value_at<T>(
                    H_row[matrix_width - 1], last_column_id);
                if (max_score < max_row_score) {
                    max_score = max_row_score;
                    max_i = i;
                }
            }
        }
    }

    if (max_i == -1 && max_j == -1) { // no alignment found
        return Alignment();
    }

    if (alignment_type_ == AlignmentType::kSW) {
        max_j = _mmxxx_index_of<T>(&(pimpl_->H[max_i * matrix_width]),
            matrix_width, max_score);

    } else if (alignment_type_ == AlignmentType::kOV) {
        if (graph->nodes()[sorted_nodes_ids[max_i - 1]]->out_edges().empty()) {
            max_j = _mmxxx_index_of<T>(&(pimpl_->H[max_i * matrix_width]),
                matrix_width, max_score);
        } else {
            max_j = normal_matrix_width - 1;
        }

    } else if (alignment_type_ == AlignmentType::kNW) {
        max_j = normal_matrix_width - 1;
    }

    // backtrack
    uint32_t max_num_predecessors = 0;
    for (uint32_t i = 0; i < (uint32_t) max_i; ++i) {
        max_num_predecessors = std::max(max_num_predecessors,
            (uint32_t) graph->nodes()[sorted_nodes_ids[i]]->in_edges().size());
    }

    __attribute__((aligned(kRegisterSize / 8))) typename T::type H[T::kNumVar];
    __attribute__((aligned(kRegisterSize / 8))) typename T::type
        H_pred[T::kNumVar * max_num_predecessors];
    __attribute__((aligned(kRegisterSize / 8))) typename T::type
        H_diag_pred[T::kNumVar * max_num_predecessors];
    __attribute__((aligned(kRegisterSize / 8))) typename T::type
        F_pred[T::kNumVar * max_num_predecessors];
    __attribute__((aligned(kRegisterSize / 8))) typename T::type E[T::kNumVar];
    __attribute__((aligned(kRegisterSize / 8))) typename T::type
        profile[T::kNumVar];

    std::vector<uint32_t> predecessors;

    int32_t i = max_i;
    int32_t j = max_j;
    int32_t prev_i = 0, prev_j = 0;

    uint32_t j_div = j / T::kNumVar;
    uint32_t j_mod = j % T::kNumVar;

    bool load_next_segment = true;

    Alignment alignment;

    do {
        // check stop condition
        if (j == -1 || i == 0) {
            break;
        }

        const auto& node = graph->nodes()[sorted_nodes_ids[i - 1]];
        // load everything
        if (load_next_segment) {
            predecessors.clear();

            // load current cells
            _mmxxx_store_si(reinterpret_cast<__mxxxi*>(H),
                pimpl_->H[i * matrix_width + j_div]);
            _mmxxx_store_si(reinterpret_cast<__mxxxi*>(E),
                pimpl_->E[i * matrix_width + j_div]);

            // load predecessors cells
            if (node->in_edges().empty()) {
                predecessors.emplace_back(0);
                _mmxxx_store_si(reinterpret_cast<__mxxxi*>(H_pred),
                    pimpl_->H[j_div]);
                _mmxxx_store_si(reinterpret_cast<__mxxxi*>(F_pred),
                    pimpl_->F[j_div]);

            } else {
                uint32_t store_pos = 0;
                for (const auto& edge: node->in_edges()) {
                    predecessors.emplace_back(
                        pimpl_->node_id_to_rank[edge->begin_node_id()] + 1);
                    _mmxxx_store_si(
                        reinterpret_cast<__mxxxi*>(&H_pred[store_pos * T::kNumVar]),
                        pimpl_->H[predecessors.back() * matrix_width + j_div]);
                    _mmxxx_store_si(
                        reinterpret_cast<__mxxxi*>(&F_pred[store_pos * T::kNumVar]),
                        pimpl_->F[predecessors.back() * matrix_width + j_div]);
                    ++store_pos;
                }
            }

            // load query profile cells
            _mmxxx_store_si(reinterpret_cast<__mxxxi*>(profile),
                pimpl_->sequence_profile[node->code() * matrix_width + j_div]);
        }

        // check stop condition
        if (alignment_type_ == AlignmentType::kSW && H[j_mod] == 0) {
            break;
        }

        if (j_mod == 0) {
            // border case
            if (j_div > 0) {
                for (uint32_t p = 0; p < predecessors.size(); ++p) {
                    _mmxxx_store_si(
                        reinterpret_cast<__mxxxi*>(&H_diag_pred[p * T::kNumVar]),
                        pimpl_->H[predecessors[p] * matrix_width + (j_div - 1)]);
                }
            } else {
                for (uint32_t p = 0; p < predecessors.size(); ++p) {
                    H_diag_pred[(p + 1) * T::kNumVar - 1] =
                        pimpl_->first_column[predecessors[p]];
                }
            }
        }

        // find best predecessor cell
        bool predecessor_found = false;

        if (i != 0) {
            for (uint32_t p = 0; p < predecessors.size(); ++p) {
                if ((H[j_mod] == F_pred[p * T::kNumVar + j_mod] + gap_extend_) ||
                    (H[j_mod] == H_pred[p * T::kNumVar + j_mod] + gap_open_)) {
                    prev_i = predecessors[p];
                    prev_j = j;
                    predecessor_found = true;
                    break;
                }
            }
        }

        if (!predecessor_found && H[j_mod] == E[j_mod]) {
            prev_i = i;
            prev_j = j - 1;
            predecessor_found = true;
        }

        if (!predecessor_found && i != 0) {
            for (uint32_t p = 0; p < predecessors.size(); ++p) {
                if ((j_mod == 0 && H[j_mod] ==
                        H_diag_pred[(p + 1) * T::kNumVar - 1] + profile[j_mod]) ||
                    (j_mod != 0 && H[j_mod] ==
                        H_pred[p * T::kNumVar + j_mod - 1] + profile[j_mod])) {

                    prev_i = predecessors[p];
                    prev_j = j - 1;
                    predecessor_found = true;
                    break;
                }
            }
        }

        alignment.emplace_back(i == prev_i ? -1 : sorted_nodes_ids[i - 1],
            j == prev_j ? -1 : j);

        // update for next round
        load_next_segment = (i == prev_i ? false : true) ||
            (j != prev_j && prev_j % T::kNumVar == T::kNumVar - 1 ? true : false);

        i = prev_i;
        j = prev_j;
        j_div = j / T::kNumVar;
        j_mod = j % T::kNumVar;

    } while (true);

    // update alignment for NW (backtrack stops on first row or column)
    if (alignment_type_ == AlignmentType::kNW) {
        while (i == 0 && j != -1) {
            alignment.emplace_back(-1, j);
            --j;
        }
        while (i != 0 && j == -1) {
            alignment.emplace_back(sorted_nodes_ids[i - 1], -1);

            const auto& node = graph->nodes()[sorted_nodes_ids[i - 1]];
            if (node->in_edges().empty()) {
                i = 0;
            } else {
                for (const auto& edge: node->in_edges()) {
                    uint32_t pred_i =
                        pimpl_->node_id_to_rank[edge->begin_node_id()] + 1;
                    if (pimpl_->first_column[i] ==
                        pimpl_->first_column[pred_i] + gap_extend_) {
                        i = pred_i;
                        break;
                    }
                }
            }
        }
    }

    std::reverse(alignment.begin(), alignment.end());
    return alignment;

#else

    return Alignment();

#endif
}

/*
// TODO: come up with a more elegant way for this function (its duplicate with SisdAlignment::adjust_node_ids)
void SimdAlignment::adjust_node_ids(const std::vector<int32_t>& mapping) {
    for (uint32_t i = 0; i < alignment_node_ids_.size(); ++i) {
        if (alignment_node_ids_[i] != -1) {
            alignment_node_ids_[i] = mapping[alignment_node_ids_[i]];
        }
    }
}*/

}
