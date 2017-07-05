# spoa

![image](https://travis-ci.org/rvaser/spoa.svg?branch=master)

Spoa (SIMD POA) is a c++ implementation of the partial order alignment (POA) algorithm (as described in 10.1093/bioinformatics/18.3.452) which is used to generate consensus sequences (as described in 10.1093/bioinformatics/btg109). It supports three alignment modes: local (Smith-Waterman), global (Needleman-Wunsh) and semi-global alignment (overlap). It supports Intel SSE4.1+ vectorization for all alignment modes. AVX2 is implemented as well and can be enabled manually in `simd_alignment_engine.cpp` file. It is disabled by default as it does not provide adequate speedup for shorter sequences due to high latency shifts.

## Dependencies

### Linux

Application uses following software:

1. gcc 4.8+ or clang 3.4+
2. cmake 3.2+

## Installation

CmakeLists is provided in the project root folder. By running the following commands:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
an executable named `spoa` will appear in `build/bin` directory and a library named `libspoa.a` in the `build/lib` directory.

Optionally, you can run `sudo make install` to install spoa executable and library to your machine.

## Usage

Usage of spoa is as following:

    spoa -a <FASTA file> (or -q <FASTQ file>) [arguments ...]

    arguments:
    -m, --match <int>
        default: 5
        score for matching bases (must be positive integer)
    -x, --mismatch <int>
        default: -4
        score for mismatching bases (must be negative integer)
    -o, --gap-open <int>
        default: -8
        gap opening penalty (must be negative integer)
    -e, --gap-extend <int>
        default: -6
        gap extension penalty (must be negative integer)
    -l, --algorithm <int>
        default: 0
        alignment mode: 0 - local, 1 - global, 2 - semi-global
    -r, --result <int>
        default: 0
        result mode: 0 - consensus, 1 - multiple sequence alignment, 2 - 0 + 1
    -h, -help
        prints out the help

###Library

Simple library usage can be seen in the following `example.cpp` file. This code shows how to get consensus and multiple sequence alignment for a set of sequences without quality values.

```cpp
#include "spoa/spoa.hpp"

int main(int argc, char** argv) {

    std::vector<std::string> sequences = {
        "CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT",
        "ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT",
        "ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG",
        "CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA",
        "GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT",
        "CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT"
    };

    auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(atoi(argv[1])),
        atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));

    auto graph = spoa::createGraph();

    for (const auto& it: sequences) {
        auto alignment = alignment_engine->align_sequence_with_graph(it, graph);
        graph->add_alignment(alignment, it);
    }

    std::string consensus = graph->generate_consensus();

    fprintf(stderr, "Consensus (%zu)\n", consensus.size());
    fprintf(stderr, "%s\n", consensus.c_str());

    std::vector<std::string> msa;
    graph->generate_multiple_sequence_alignment(msa);

    fprintf(stderr, "Multiple sequence alignment\n");
    for (const auto& it: msa) {
        fprintf(stderr, "%s\n", it.c_str());
    }

    return 0;
}
```

This code can be compiled from spoa root with:
```bash
g++ example.cpp -std=c++11 -Iinclude/ -Lbuild/lib/ -lspoa -o example
```
If spoa was installed, the following modified command will compile as well:
```bash
g++ example.cpp -std=c++11 -lspoa -o example
```

The executable can be run with:
```bash
./example 0 5 -4 -8 -6
```

On the other hand, if you are using `cmake` you can add spoa to your project by adding commands `add_subdirectory(vendor/spoa EXCLUDE_FROM_ALL)` and `target_link_libraries(your_exe spoa)` to your main CMakeLists file.


## Contact information

For additional information, help and bug reports please send an email to: robert.vaser@fer.hr.

## Acknowledgement

This work has been supported in part by Croatian Science Foundation under the project UIP-11-2013-7353.
