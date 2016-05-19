/*!
 * @file simd_alignment.hpp
 *
 * @brief SIMD alignment header file
 */

#pragma once

#include <assert.h>
#include <string>
#include <memory>
#include <vector>

namespace SPOA {

class AlignmentParams;

void simd_align_sequence_to_graph(const std::string& sequence, std::shared_ptr<Graph> graph,
    AlignmentParams params);

}
