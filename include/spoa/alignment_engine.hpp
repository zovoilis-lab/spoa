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
    int8_t m, int8_t n, int8_t g);

std::unique_ptr<AlignmentEngine> createAlignmentEngine(AlignmentType type,
    int8_t m, int8_t n, int8_t g, int8_t e);

std::unique_ptr<AlignmentEngine> createAlignmentEngine(AlignmentType type,
    int8_t m, int8_t n, int8_t g, int8_t e, int8_t q, int8_t c);

class AlignmentEngine {
public:
    virtual ~AlignmentEngine() {}

    virtual void prealloc(uint32_t max_sequence_size, uint32_t alphabet_size) = 0;

    Alignment align(const std::string& sequence,
        const std::unique_ptr<Graph>& graph);

    virtual Alignment align(const char* sequence, uint32_t sequence_size,
        const std::unique_ptr<Graph>& graph) = 0;
protected:
    AlignmentEngine(AlignmentType type, AlignmentSubtype subtype, int8_t m,
        int8_t n, int8_t g, int8_t e, int8_t q, int8_t c);
    AlignmentEngine(const AlignmentEngine&) = delete;
    const AlignmentEngine& operator=(const AlignmentEngine&) = delete;

    AlignmentType type_;
    AlignmentSubtype subtype_;
    int8_t m_;
    int8_t n_;
    int8_t g_;
    int8_t e_;
    int8_t q_;
    int8_t c_;
};

}
