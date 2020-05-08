/*!
 * @file dispatcher.cpp
 *
 * @brief CPU dispatching mechanism that also covers non-dispatching case
 */

#ifdef GEN_DISPATCH

// FeatureDetector - can be substituted for google's cpu_features
#include "cpu_x86.h"

#endif
#include "simd_alignment_engine_impl.hpp"


namespace spoa{

#ifndef GEN_DISPATCH
template class SimdAlignmentEngine<Arch::automatic>;

template
std::unique_ptr<AlignmentEngine> createSimdAlignmentEngine<Arch::automatic>(AlignmentType type,
    AlignmentSubtype subtype, std::int8_t m, std::int8_t n, std::int8_t g,
    std::int8_t e, std::int8_t q, std::int8_t c);
#endif


std::unique_ptr<AlignmentEngine> createSimdAlignmentEngine(AlignmentType type,
    AlignmentSubtype subtype, std::int8_t m, std::int8_t n, std::int8_t g,
    std::int8_t e, std::int8_t q, std::int8_t c) {

#ifdef GEN_DISPATCH
        FeatureDetector::cpu_x86 cpu_features;
        cpu_features.detect_host();

    if (cpu_features.HW_AVX2)
    {
        //std::cout<<"AVX2"<<std::endl;
        return createSimdAlignmentEngine<Arch::avx2>(type,
            subtype, m, n, g, e, q, c);
    }
    else if (cpu_features.HW_SSE41){

        //std::cout<<"SSE4"<<std::endl;
        return createSimdAlignmentEngine<Arch::sse4_1>(type,
            subtype, m, n, g, e, q, c);
    }
    else {
        //std::cout<<"SSE2"<<std::endl;
        return createSimdAlignmentEngine<Arch::sse2>(type,
            subtype, m, n, g, e, q, c);
    }
#else
    return createSimdAlignmentEngine<Arch::automatic>(type,
            subtype, m, n, g, e, q, c);
#endif


}

}

