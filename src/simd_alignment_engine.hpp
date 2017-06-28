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
#include <utility>

#include "spoa/alignment_engine.hpp"

namespace spoa {

class Graph;

class SimdAlignmentEngine;
std::unique_ptr<AlignmentEngine> createSimdAlignmentEngine(
    AlignmentType alignment_type, int8_t match, int8_t mismatch,
    int8_t gap_open, int8_t gap_extend);

class SimdAlignmentEngine: public AlignmentEngine {
public:
    ~SimdAlignmentEngine();

    Alignment align_sequence_with_graph(const std::string& sequence,
        const std::unique_ptr<Graph>& graph) override;

    friend std::unique_ptr<AlignmentEngine> createSimdAlignmentEngine(
        AlignmentType alignment_type, int8_t match, int8_t mismatch,
        int8_t gap_open, int8_t gap_extend);
private:
    SimdAlignmentEngine(AlignmentType alignment_type, int8_t match,
        int8_t mismatch, int8_t gap_open, int8_t gap_extend);
    SimdAlignmentEngine(const SimdAlignmentEngine&) = delete;
    const SimdAlignmentEngine& operator=(const SimdAlignmentEngine&) = delete;

    template<typename T>
    Alignment align(const std::string& sequence,
        const std::unique_ptr<Graph>& graph);

    struct Implementation;
    std::unique_ptr<Implementation> pimpl_;
};

}
