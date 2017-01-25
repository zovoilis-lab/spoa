# spoa

Spoa (SIMD POA) is a c++ implementation of the partial order alignment (POA) algorithm (as described in 10.1093/bioinformatics/18.3.452) which is used to generate consensus sequences (as described in 10.1093/bioinformatics/btg109). It supports three alignment modes: local (Smith-Waterman), global (Needleman-Wunsh) and semi-global alignment (overlap). It supports Intel SSE4.1+ vectorization for all alignment modes.

## Dependencies

### Linux

Application uses following software:

1. gcc 4.8+

## Installation

### Linux

Makefile is provided in the project root folder. Inside spoa root, run:

    make

After running make, an executable named spoa will appear in the current directory.

If you would like to add spoa to your project as a static library, run:

    make install

and add -Iinclude/ -Llib/ -lspoa while compiling your project. Please look usage for a
detailed example.

If you would like to add spoa source files to you project, just include the src/ directory.
Unnecessary files (chain.* and main.cpp) have a macro which excludes them from compiling in other projects).

## Usage

## Executables

All examples assume that make has been run and that spoa was successfully compiled.

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

To remove spoa executable, run:

    make clean

###Library

Simple library usage can be seen in the following example.c file. This code shows
how to get consensus and multiple sequence alignment for a set of sequences.

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

        auto params = spoa::AlignmentParams(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]),
            atoi(argv[4]), (spoa::AlignmentType) atoi(argv[5]));

        std::string consensus = spoa::generate_consensus(sequences, params, true);

        fprintf(stderr, "Consensus (%zu)\n", consensus.size());
        fprintf(stderr, "%s\n", consensus.c_str());

        std::vector<std::string> msa;
        spoa::generate_msa(msa, sequences, params, true);

        fprintf(stderr, "Multiple sequence alignment\n");
        for (const auto& it: msa) {
            fprintf(stderr, "%s\n", it.c_str());
        }

        return 0;
    }

This code can be compiled with:

    g++ example.cpp -std=c++11 -Iinclude/ -Llib/ -lspoa -o example

And the executable can be run with:

    ./example 5 -4 -8 -6 0

To remove include/ and lib/ directories created for the library, run:

    make remove

## Contact information

For additional information, help and bug reports please send an email to: robert.vaser@fer.hr.

## Acknowledgement

This work has been supported in part by Croatian Science Foundation under the project UIP-11-2013-7353.
