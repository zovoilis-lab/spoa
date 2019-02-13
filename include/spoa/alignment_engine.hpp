/*!
 * @file alignment_engine.hpp
 *
 * @brief AlignmentEngine class header file
 */

#pragma once

#include <stdint.h>
#include <string>
#include <memory>
#include <vector>
#include <utility>

namespace spoa {

enum class AlignmentType {
    kSW, // Smith Waterman
    kNW, // Needleman Wunsch
    kOV // Overlap
};

enum class AlignmentSubtype {
    kLinear,
    kAffine,
    kConvex
};

class Graph;
using Alignment = std::vector<std::pair<int32_t, int32_t>>;

class AlignmentEngine;
std::unique_ptr<AlignmentEngine> createAlignmentEngine(AlignmentType type,
    int8_t match, int8_t mismatch, int8_t gap_open);

std::unique_ptr<AlignmentEngine> createAlignmentEngine(AlignmentType type,
    int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend);

class AlignmentEngine {
public:
    virtual ~AlignmentEngine() {}

    virtual void prealloc(uint32_t max_sequence_size, uint32_t alphabet_size) = 0;

    Alignment operator()(const std::string& sequence,
        const std::unique_ptr<Graph>& graph);

    virtual Alignment operator()(const char* sequence, uint32_t sequence_size,
        const std::unique_ptr<Graph>& graph) = 0;
protected:
    AlignmentEngine(AlignmentType type, AlignmentSubtype subtype, int8_t match,
        int8_t mismatch, int8_t gap_open, int8_t gap_extend);
    AlignmentEngine(const AlignmentEngine&) = delete;
    const AlignmentEngine& operator=(const AlignmentEngine&) = delete;

    AlignmentType type_;
    AlignmentSubtype subtype_;
    int8_t match_;
    int8_t mismatch_;
    int8_t gap_open_;
    int8_t gap_extend_;
};

}
