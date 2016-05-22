# spoa

Spoa (SIMD POA) is a c++ implementation of the partial order alignment (POA) algorithm (as described in 10.1093/bioinformatics/18.3.452) which is used to generate consensus sequences (as described in 10.1093/bioinformatics/btg109). It supports three alignment modes: local (Smith-Waterman), global (Needleman-Wunsh) and semi-global alignment (overlap). It supports Intel SSE4.1+ vectorization for all alignment modes.

## DEPENDENCIES

### LINUX

Application uses following software:

1. gcc 4.*

## INSTALLATION

### LINUX

Makefile is provided in the project root folder. Inside spoa root, run:

    make

After running make, an executable named spoa will appear in the current directory.

## USAGE

All examples assume that make has been run and that spoa was successfully compiled.

Spoa must be run as following:

    ./spoa <sequences.fasta> <match> <missmatch> <gap_open> <gap_extend> <algorithm>

Missmatch, gap open and gap extend must be negative values. Supported algorithms are: 0 (local), 1 (global) and 2 (semi-global).

To remove spoa executable, run:

    make clean

## Contact information

For additional information, help and bug reports please send an email to: robert.vaser@fer.hr.

## Acknowledgment

This work has been supported in part by Croatian Science Foundation under the project UIP-11-2013-7353.
