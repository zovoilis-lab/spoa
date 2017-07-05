/*!
 * @file chain.cpp
 *
 * @brief Chain class source file
 */

#include "chain.hpp"

namespace spoa {

Chain::Chain(uint64_t id, const char* name, uint32_t name_length,
    const char* data, uint32_t data_length)
        : id_(id), name_(name, name_length), data_(data, data_length),
        quality_() {
}

Chain::Chain(uint64_t id, const char* name, uint32_t name_length,
    const char* data, uint32_t data_length, const char* quality,
    uint32_t quality_length)
        : id_(id), name_(name, name_length), data_(data, data_length),
        quality_(quality, quality_length) {
}

}
