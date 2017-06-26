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

class Graph;
using Alignment = std::vector<std::pair<int32_t, int32_t>>;

class AlignmentEngine;
std::unique_ptr<AlignmentEngine> createAlignmentEngine(
    AlignmentType alignment_type, int8_t match, int8_t mismatch,
    int8_t gap_open, int8_t gap_extend);

class AlignmentEngine {
public:
    virtual ~AlignmentEngine() {}
    virtual Alignment align_sequence_with_graph(const std::string& sequence,
        const std::unique_ptr<Graph>& graph) = 0;
protected:
    AlignmentEngine(AlignmentType alignment_type, int8_t match, int8_t mismatch,
        int8_t gap_open, int8_t gap_extend);
    AlignmentEngine(const AlignmentEngine&) = delete;
    const AlignmentEngine& operator=(const AlignmentEngine&) = delete;

    AlignmentType alignment_type_;
    int8_t match_;
    int8_t mismatch_;
    int8_t gap_open_;
    int8_t gap_extend_;
};

}
