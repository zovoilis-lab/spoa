#include <iostream>

#ifdef GEN_DISPATCH

//#include "cpuinfo_x86.h"
#include "cpu_x86.h"

#endif
#include "simd_alignment_engine_impl.hpp"

//static const cpu_features::X86Features features = cpu_features::GetX86Info().features;

//static const struct cpu_x86 features(cpu_x86.detect_host());

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
        //std::cout<<features.sse4_1<<std::endl;
#ifdef GEN_DISPATCH
        FeatureDetector::cpu_x86 cpu_features;
        cpu_features.detect_host();

    if (cpu_features.HW_AVX2)
    {
        std::cout<<"AVX2"<<std::endl;
        return createSimdAlignmentEngine<Arch::avx2>(type,
            subtype, m, n, g, e, q, c);
    }
    else if (cpu_features.HW_SSE41){

        std::cout<<"SSE4"<<std::endl;
        return createSimdAlignmentEngine<Arch::sse4_1>(type,
            subtype, m, n, g, e, q, c);
    }
    else {
        std::cout<<"SSE2"<<std::endl;
        return createSimdAlignmentEngine<Arch::sse2>(type,
            subtype, m, n, g, e, q, c);
    }
#else
    return createSimdAlignmentEngine<Arch::automatic>(type,
            subtype, m, n, g, e, q, c);
#endif


}

}

