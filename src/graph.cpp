/*!
 * @file graph.cpp
 *
 * @brief Graph class source file
 */

#include <set>

#include "node.hpp"
#include "edge.hpp"
#include "graph.hpp"

std::unique_ptr<Graph> createGraph(const std::string& sequence) {
    return std::unique_ptr<Graph>(new Graph(sequence));
}

Graph::Graph(const std::string& sequence) :
        num_sequences_(), num_nodes_(), nodes_(), is_sorted_(false),
        sorted_nodes_ids_(), sequences_start_nodes_ids_() {

    assert(sequence.size() != 0);

    int32_t start_node_id = this->add_sequence(sequence, 0, sequence.size());

    sequences_start_nodes_ids_.emplace_back(start_node_id);
    ++num_sequences_;
}

Graph::~Graph() {
}

uint32_t Graph::add_node(char letter) {
    nodes_.emplace_back(createNode(num_nodes_, letter));
    return num_nodes_++;
}

void Graph::add_edge(uint32_t begin_node_id, uint32_t end_node_id) {
    if (begin_node_id >= num_nodes_ || end_node_id >= num_nodes_) {
        fprintf(stderr, "%u %u\n", begin_node_id, end_node_id);
    }
    assert(begin_node_id < num_nodes_ && end_node_id < num_nodes_);

    for (const auto& edge: nodes_[begin_node_id]->out_edges()) {
        if (edge->end_node_id() == end_node_id) {
            edge->add_sequence_label(num_sequences_);
            return;
        }
    }

    EdgeSharedPtr edge = createEdge(begin_node_id, end_node_id, num_sequences_);
    nodes_[begin_node_id]->add_out_edge(edge);
    nodes_[end_node_id]->add_in_edge(edge);
}

void Graph::topological_sort() {

    if (is_sorted_) {
        return;
    }
    sorted_nodes_ids_.clear();

    // 0 - unmarked, 1 - temporarily marked, 2 - permanently marked
    std::vector<uint8_t> marks(num_nodes_, 0);

    uint32_t i = 0;
    while (true) {
        for (; i < num_nodes_; ++i) {
            if (marks[nodes_[i]->id()] == 0) {
                break;
            }
        }
        if (i == nodes_.size()) {
            break;
        }

        this->visit_node(i, marks);
    }

    assert(this->is_topologically_sorted() == true);
    is_sorted_ = true;
}

void Graph::visit_node(uint32_t node_id, std::vector<uint8_t>& marks) {
    assert(marks[node_id] != 1 && "Graph is not a DAG!");

    if (marks[node_id] == 0) {
        marks[node_id] = 1;

        const auto& in_edges = nodes_[node_id]->in_edges();
        for (const auto& edge: in_edges) {
            this->visit_node(edge->begin_node_id(), marks);
        }

        marks[node_id] = 2;
        sorted_nodes_ids_.emplace_back(node_id);
    }
}

bool Graph::is_topologically_sorted() const {
    assert(nodes_.size() == sorted_nodes_ids_.size());

    std::set<uint32_t> visited_nodes;
    for (uint32_t node_id: sorted_nodes_ids_) {
        for (const auto& edge: nodes_[node_id]->in_edges()) {
            if (visited_nodes.count(edge->begin_node_id()) == 0) {
                return false;
            }
        }
        visited_nodes.insert(node_id);
    }

    return true;
}

void Graph::add_alignment(const std::vector<int32_t>& node_ids,
    const std::vector<int32_t>& seq_ids, const std::string& sequence) {

    assert(sequence.size() != 0);
    assert(node_ids.size() == seq_ids.size());

    std::vector<uint32_t> valid_seq_ids;
    for (const auto& id: seq_ids) {
        if (id != -1) {
            valid_seq_ids.emplace_back(id);
        }
    }
    assert(valid_seq_ids.size() != 0);

    uint32_t tmp = num_nodes_;
    int32_t start_node_id = this->add_sequence(sequence, 0, valid_seq_ids.front());
    int32_t head_node_id = tmp == num_nodes_ ? -1 : num_nodes_ - 1;

    fprintf(stderr, "%d %d\n", start_node_id, head_node_id);

    int32_t tail_node_id = this->add_sequence(sequence, valid_seq_ids.back() + 1, sequence.size());

    fprintf(stderr, "%d\n", tail_node_id);

    int32_t new_node_id = -1;

    for (uint32_t i = 0; i < seq_ids.size(); ++i) {
        //fprintf(stderr, "%d| %d %d %d\n", seq_ids[i], start_node_id, head_node_id, new_node_id);
        if (seq_ids[i] == -1) {
            continue;
        }

        char letter = sequence[seq_ids[i]];
        //fprintf(stderr, "%c\n", letter);
        if (node_ids[i] == -1) {
            new_node_id = this->add_node(letter);

        } else {
            auto node = nodes_[node_ids[i]];
            if (node->letter() == letter) {
                new_node_id = node_ids[i];

            } else {
                int32_t aligned_to_node_id = -1;
                for (const auto& aid: node->aligned_nodes_ids()) {
                    if (nodes_[aid]->letter() == letter) {
                        aligned_to_node_id = aid;
                        break;
                    }
                }

                if (aligned_to_node_id == -1) {
                    new_node_id = this->add_node(letter);

                    for (const auto& aid: node->aligned_nodes_ids()) {
                        nodes_[new_node_id]->add_aligned_node_id(aid);
                        nodes_[aid]->add_aligned_node_id(new_node_id);
                    }

                    nodes_[new_node_id]->add_aligned_node_id(node_ids[i]);
                    node->add_aligned_node_id(new_node_id);

                } else {
                    new_node_id = aligned_to_node_id;
                }
            }
        }

        if (start_node_id == -1) {
            start_node_id = new_node_id;
        }
        if (head_node_id != -1) {
            this->add_edge(head_node_id, new_node_id);
        }
        head_node_id = new_node_id;
    }
    //fprintf(stderr, "%d %d %d\n", start_node_id, head_node_id, new_node_id);

    if (tail_node_id != -1) {
        this->add_edge(head_node_id, tail_node_id);
    }

    ++num_sequences_;
    sequences_start_nodes_ids_.emplace_back(start_node_id);

    is_sorted_ = false;
    this->topological_sort();
}

int32_t Graph::add_sequence(const std::string& sequence, uint32_t begin, uint32_t end) {

    if (begin == end) {
        return -1;
    }

    assert(begin < sequence.size() && end <= sequence.size());

    int32_t first_node_id = this->add_node(sequence[begin]);

    uint32_t node_id;
    for (uint32_t i = begin + 1; i < end; ++i) {
        node_id = this->add_node(sequence[i]);
        this->add_edge(node_id - 1, node_id);
    }

    return first_node_id;
}

void Graph::print() const {
    printf("digraph %d {\n", num_sequences_);
    printf("    graph [rankdir=LR]\n");
    for (uint32_t i = 0; i < num_nodes_; ++i) {
        printf("    %d [label = \"%d|%c\"]\n", i, i, nodes_[i]->letter());
        for (const auto& edge: nodes_[i]->out_edges()) {
            printf("    %d -> %d [label = \"%zu\"]\n", i, edge->end_node_id(),
                edge->sequence_labels().size());
        }
        for (const auto& aid: nodes_[i]->aligned_nodes_ids()) {
            if (aid > i) {
                printf("    %d -> %d [style = dotted, arrowhead = none]\n", i, aid);
            }
        }
    }
    printf("}\n");
}
