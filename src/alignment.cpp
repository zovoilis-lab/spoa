/*!
 * @file alignment.cpp
 *
 * @brief Alignment class source file
 */

#include <limits>
#include <algorithm>

#include "sisd_alignment.hpp"
#include "simd_alignment.hpp"
#include "spoa/graph.hpp"
#include "spoa/alignment.hpp"

namespace spoa {

std::unique_ptr<Alignment> createAlignment(AlignmentType alignment_type
    uint8_t match, uint8_t mismatch, uint8_t gap_open, uint8_t gap_extend) {

    auto alignment = createSimdAlignment(alignment_type, match, mismatch,
        gap_open, gap_extend);

    if (alignment == nullptr) {
        return createSisdAlignment(alignment_type, match, mismatch,
            gap_open, gap_extend);
    }

    return alignment;
}

Alignment::Alignment(AlignmentType alignment_type, uint8_t match,
    uint8_t mismatch, uint8_t gap_open, uint8_t gap_extend)
        : alignment_type_(alignment_type), match_(match), mismatch_(mismatch),
        gap_open_(gap_open), gap_extend_(gap_extend) {
}

}
