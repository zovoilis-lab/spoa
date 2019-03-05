/*!
 * @file sisd_alignment_engine.hpp
 *
 * @brief SisdAlignmentEngine class header file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <string>
#include <vector>

#include "spoa/alignment_engine.hpp"

namespace spoa {

class Graph;

class SisdAlignmentEngine;
std::unique_ptr<AlignmentEngine> createSisdAlignmentEngine(AlignmentType type,
    AlignmentSubtype subtype, int8_t m, int8_t n, int8_t g, int8_t e, int8_t q,
    int8_t c);

class SisdAlignmentEngine: public AlignmentEngine {
public:
    ~SisdAlignmentEngine();

    void prealloc(uint32_t max_sequence_size, uint32_t alphabet_size) override;

    Alignment align(const char* sequence, uint32_t sequence_size,
        const std::unique_ptr<Graph>& graph) override;

    friend std::unique_ptr<AlignmentEngine> createSisdAlignmentEngine(
        AlignmentType type, AlignmentSubtype subtype, int8_t m, int8_t n,
        int8_t g, int8_t e, int8_t c, int8_t q);
private:
    SisdAlignmentEngine(AlignmentType type, AlignmentSubtype subtype,
        int8_t m, int8_t n, int8_t g, int8_t e, int8_t q, int8_t c);
    SisdAlignmentEngine(const SisdAlignmentEngine&) = delete;
    const SisdAlignmentEngine& operator=(const SisdAlignmentEngine&) = delete;

    Alignment linear(const char* sequence, uint32_t sequence_size,
        const std::unique_ptr<Graph>& graph) noexcept;

    Alignment affine(const char* sequence, uint32_t sequence_size,
        const std::unique_ptr<Graph>& graph) noexcept;

    Alignment convex(const char* sequence, uint32_t sequence_size,
        const std::unique_ptr<Graph>& graph) noexcept;

    void realloc(uint32_t matrix_width, uint32_t matrix_height,
        uint32_t num_codes);

    void initialize(const char* sequence, uint32_t sequence_size,
        const std::unique_ptr<Graph>& graph) noexcept;

    struct Implementation;
    std::unique_ptr<Implementation> pimpl_;
};

}
