/*!
 * @file alignment.hpp
 *
 * @brief Alignment class header file
 */

#pragma once

#include <string>
#include <memory>
#include <vector>
#include <pair>

namespace spoa {

enum class AlignmentType {
    kSW, // Smith Waterman
    kNW, // Needleman Wunsch
    kOV // Overlap
};

class Graph;

class Alignment;
std::unique_ptr<Alignment> createAlignment(AlignmentType alignment_type,
    uint8_t match, uint8_t mismatch, uint8_t gap_open, uint8_t gap_extend);

class Alignment {
public:
    virtual ~Alignment() {}
    virtual std::vector<std::pair<int32_t, int32_t>> align_sequence_with_graph(
        const std::string& sequence, const std::unique_ptr<Graph>& graph) = 0;
protected:
    Alignment(AlignmentType alignment_type, uint8_t match, uint8_t mismatch,
        uint8_t gap_open, uint8_t gap_extend));
    Alignment(const Alignment&) = delete;
    const Alignment& operator=(const Alignment&) = delete;

    AlignmentType alignment_type_;
    uint8_t match_;
    uint8_t mismatch_;
    uint8_t gap_open_;
    uint8_t gap_extend_;
};

}
