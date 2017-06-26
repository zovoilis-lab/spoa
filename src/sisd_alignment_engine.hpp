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
#include <utility>

#include "spoa/alignment_engine.hpp"

namespace spoa {

class Graph;

class SisdAlignmentEngine;
std::unique_ptr<AlignmentEngine> createSisdAlignmentEngine(
    AlignmentType alignment_type, int8_t match, int8_t mismatch,
    int8_t gap_open, int8_t gap_extend);

class SisdAlignmentEngine: public AlignmentEngine {
public:
    ~SisdAlignmentEngine();

    Alignment align_sequence_with_graph(const std::string& sequence,
        const std::unique_ptr<Graph>& graph) override;

    friend std::unique_ptr<AlignmentEngine> createSisdAlignmentEngine(
        AlignmentType alignment_type, int8_t match, int8_t mismatch,
        int8_t gap_open, int8_t gap_extend);
private:
    SisdAlignmentEngine(AlignmentType alignment_type, int8_t match,
        int8_t mismatch, int8_t gap_open, int8_t gap_extend);
    SisdAlignmentEngine(const SisdAlignmentEngine&) = delete;
    const SisdAlignmentEngine& operator=(const SisdAlignmentEngine&) = delete;

    struct Implementation;
    std::unique_ptr<Implementation> pimpl_;
};

}
