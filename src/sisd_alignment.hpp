/*!
 * @file sisd_alignment.hpp
 *
 * @brief SisdAlignment class header file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <string>
#include <vector>
#include <pair>

namespace spoa {

class Graph;
class Alignment;

class SisdAlignment;
std::unique_ptr<Alignment> createSisdAlignment(AlignmentType alignment_type,
    uint8_t match, uint8_t mismatch, uint8_t gap_open, uint8_t gap_extend);

class SisdAlignment: public Alignment {
public:

    ~SisdAlignment();

    std::vector<std::pair<int32_t, int32_t>> align_sequence_with_graph(
        const std::string& sequence, const std::unique_ptr<Graph>& graph) overide;

    friend std::unique_ptr<Alignment> createSisdAlignment(
        AlignmentType alignment_type, uint8_t match, uint8_t mismatch,
        uint8_t gap_open, uint8_t gap_extend);

private:

    SisdAlignment(AlignmentType alignment_type, uint8_t match, uint8_t mismatch,
        uint8_t gap_open, uint8_t gap_extend);
    SisdAlignment(const SisdAlignment&) = delete;
    const SisdAlignment& operator=(const SisdAlignment&) = delete;

    struct Impl;
    std::unique_ptr<Impl> pImpl_;
};

}
