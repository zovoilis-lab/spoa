/*!
 * @file simd_alignment_engine.hpp
 *
 * @brief SimdAlignmentEngine class header file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <string>
#include <vector>

#include "spoa/alignment_engine.hpp"

namespace spoa {

class Graph;

class SimdAlignmentEngine;
std::unique_ptr<AlignmentEngine> createSimdAlignmentEngine(AlignmentType type,
    AlignmentSubtype subtype, int8_t m, int8_t n, int8_t g, int8_t e);

class SimdAlignmentEngine: public AlignmentEngine {
public:
    ~SimdAlignmentEngine();

    void prealloc(uint32_t max_sequence_size, uint32_t alphabet_size) override;

    Alignment align(const char* sequence, uint32_t sequence_size,
        const std::unique_ptr<Graph>& graph) override;

    friend std::unique_ptr<AlignmentEngine> createSimdAlignmentEngine(
        AlignmentType type, AlignmentSubtype subtype, int8_t m, int8_t n,
        int8_t g, int8_t e);
private:
    SimdAlignmentEngine(AlignmentType type, AlignmentSubtype subtype,
        int8_t m, int8_t n, int8_t g, int8_t e);
    SimdAlignmentEngine(const SimdAlignmentEngine&) = delete;
    const SimdAlignmentEngine& operator=(const SimdAlignmentEngine&) = delete;

    template<typename T>
    Alignment linear(const char* sequence, uint32_t sequence_size,
        const std::unique_ptr<Graph>& graph) noexcept;

    template<typename T>
    Alignment affine(const char* sequence, uint32_t sequence_size,
        const std::unique_ptr<Graph>& graph) noexcept;

    void realloc(uint32_t matrix_width, uint32_t matrix_height,
        uint32_t num_codes);

    template<typename T>
    void initialize(const char* sequence, const std::unique_ptr<Graph>& graph,
        uint32_t normal_matrix_width, uint32_t matrix_width,
        uint32_t matrix_height) noexcept;

    struct Implementation;
    std::unique_ptr<Implementation> pimpl_;
};

}
