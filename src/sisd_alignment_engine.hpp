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
    AlignmentSubtype subtype, int8_t match, int8_t mismatch, int8_t gap_open,
    int8_t gap_extend);

class SisdAlignmentEngine: public AlignmentEngine {
public:
    ~SisdAlignmentEngine();

    void prealloc(uint32_t max_sequence_size, uint32_t alphabet_size) override;

    Alignment operator()(const char* sequence, uint32_t sequence_size,
        const std::unique_ptr<Graph>& graph) override;

    friend std::unique_ptr<AlignmentEngine> createSisdAlignmentEngine(
        AlignmentType type, AlignmentSubtype subtype, int8_t match,
        int8_t mismatch, int8_t gap_open, int8_t gap_extend);
private:
    SisdAlignmentEngine(AlignmentType type, AlignmentSubtype subtype,
        int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend);
    SisdAlignmentEngine(const SisdAlignmentEngine&) = delete;
    const SisdAlignmentEngine& operator=(const SisdAlignmentEngine&) = delete;

    Alignment linear(const char* sequence, uint32_t sequence_size,
        const std::unique_ptr<Graph>& graph) noexcept;

    Alignment affine(const char* sequence, uint32_t sequence_size,
        const std::unique_ptr<Graph>& graph) noexcept;

    void realloc(uint32_t matrix_width, uint32_t matrix_height,
        uint32_t num_codes);

    void initialize(const char* sequence, uint32_t sequence_size,
        const std::unique_ptr<Graph>& graph) noexcept;

    struct Implementation;
    std::unique_ptr<Implementation> pimpl_;
};

}
