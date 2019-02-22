/*!
 * @file alignment_engine.cpp
 *
 * @brief AlignmentEngine class source file
 */

#include <stdlib.h>
#include <limits>
#include <algorithm>

#include "sisd_alignment_engine.hpp"
#include "simd_alignment_engine.hpp"
#include "spoa/alignment_engine.hpp"

namespace spoa {

std::unique_ptr<AlignmentEngine> createAlignmentEngine(AlignmentType type,
    int8_t match, int8_t mismatch, int8_t gap_open) {

    return createAlignmentEngine(type, match, mismatch, gap_open, gap_open);
}

std::unique_ptr<AlignmentEngine> createAlignmentEngine(AlignmentType type,
    int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend) {

    if (type != AlignmentType::kSW && type != AlignmentType::kNW &&
        type != AlignmentType::kOV) {

        fprintf(stderr, "[spoa::createAlignmentEngine] error: "
            "invalid alignment type!\n");
        exit(1);
    }

    if (gap_open > 0) {
        fprintf(stderr, "[spoa::createAlignmentEngine] error: "
            "gap opening penalty must be non-positive!\n");
        exit(1);
    }
    if (gap_extend > 0) {
        fprintf(stderr, "[spoa::createAlignmentEngine] error: "
            "gap extension penalty must be non-positive!\n");
        exit(1);
    }

    AlignmentSubtype subtype = gap_open >= gap_extend ?
        AlignmentSubtype::kLinear : AlignmentSubtype::kAffine;

    if (subtype == AlignmentSubtype::kLinear) {
        gap_extend = gap_open;
    }

    auto alignment_engine = createSimdAlignmentEngine(type, subtype,
        match, mismatch, gap_open, gap_extend);

    if (alignment_engine == nullptr) {
        return createSisdAlignmentEngine(type, subtype, match, mismatch,
            gap_open, gap_extend);
    }

    return alignment_engine;
}

AlignmentEngine::AlignmentEngine(AlignmentType type, AlignmentSubtype subtype,
    int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend)
        : type_(type), subtype_(subtype), match_(match), mismatch_(mismatch),
        gap_open_(gap_open), gap_extend_(gap_extend) {
}

Alignment AlignmentEngine::operator()(const std::string& sequence,
    const std::unique_ptr<Graph>& graph) {

    return this->operator()(sequence.c_str(), sequence.size(), graph);
}

}
