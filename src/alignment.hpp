/*!
 * @file alignment.hpp
 *
 * @brief Alignment class header file
 */

#pragma once

#include <assert.h>
#include <string>
#include <memory>
#include <vector>

enum class AlignmentType {
    kSW,
    kNW
};

class AlignmentParams {
public:

    AlignmentParams(int16_t match, int16_t mismatch, int16_t gap_open,
         int16_t gap_extend, AlignmentType type);
    ~AlignmentParams();

    int16_t match;
    int16_t mismatch;
    int16_t gap_open;
    int16_t gap_extend;
    AlignmentType type;
};

class Graph;
using GraphSharedPtr = std::shared_ptr<Graph>;

class Alignment;
std::unique_ptr<Alignment> createAlignment(const std::string& sequence,
    GraphSharedPtr graph, AlignmentParams&& params);

class Alignment {
public:

    ~Alignment();

    void align_sequence_to_graph();

    friend std::unique_ptr<Alignment> createAlignment(const std::string& sequence,
        GraphSharedPtr graph, AlignmentParams&& params);

private:

    Alignment(const std::string& sequence, GraphSharedPtr graph,
        AlignmentParams&& params);
    Alignment(const Alignment&) = delete;
    const Alignment& operator=(const Alignment&) = delete;

    class MatrixElement {
    public:
        MatrixElement(int32_t score, int32_t prev_i, int32_t prev_j,
            int16_t deletion_cost, int16_t insertion_cost);
        ~MatrixElement();

        int32_t score;
        int32_t prev_i;
        int32_t prev_j;
        int16_t deletion_cost;
        int16_t insertion_cost;
    };

    inline MatrixElement& matrix(uint32_t i, uint32_t j) {
        assert(i < matrix_height_);
        assert(j < matrix_width_);
        return matrix_[i * matrix_width_ + j];
    }

    void print_matrix();

    std::string sequence_;
    GraphSharedPtr graph_;
    AlignmentParams params_;

    uint32_t matrix_width_;
    uint32_t matrix_height_;
    std::vector<MatrixElement> matrix_;

    std::vector<uint32_t> node_id_to_graph_id_;
};
