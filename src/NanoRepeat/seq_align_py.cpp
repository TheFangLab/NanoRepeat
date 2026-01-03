#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cstdint>

struct Sequence {
    std::string id;
    std::string sequence;
    int length;
};

struct AlignmentResult {
    std::string read_id;
    std::string template_id;
    int score;
    int read_start_pos;
    int read_end_pos;
    int template_start_pos;
    int template_end_pos;
    int num_match;
    int num_mismatch;
    int num_gap;
    int read_seq_len;
    int template_seq_len;
};

std::string trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\n\r");
    if (start == std::string::npos) return "";
    size_t end = str.find_last_not_of(" \t\n\r");
    return str.substr(start, end - start + 1);
}

void read_fasta(const std::string& file_path, std::vector<Sequence>& sequences) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + file_path);
    }

    std::string line;
    std::string current_id;
    std::string current_seq;

    while (std::getline(file, line)) {
        std::string trimmed_line = trim(line);
        
        if (!trimmed_line.empty() && trimmed_line[0] == '>') {
            if (!current_id.empty() && !current_seq.empty()) {
                Sequence seq;
                seq.id = current_id;
                seq.sequence = current_seq;
                seq.length = current_seq.length();
                sequences.push_back(seq);
            }

            current_id = trimmed_line.substr(1);
            current_seq.clear();
        } else if (!trimmed_line.empty()) {
            current_seq += trimmed_line;
        }
    }

    if (!current_id.empty() && !current_seq.empty()) {
        Sequence seq;
        seq.id = current_id;
        seq.sequence = current_seq;
        seq.length = current_seq.length();
        sequences.push_back(seq);
    }

    file.close();
}

AlignmentResult smith_waterman(const std::string& read_id, const std::string& template_id,
                               const std::string& seq1, const std::string& seq2,
                               int match_score, int mismatch_penalty, int gap_penalty) {
    int m = seq1.length();
    int n = seq2.length();
    size_t row_stride = n + 1;
    std::vector<int> H;
    try {
        H.resize((size_t)(m + 1) * (n + 1), 0);
    } catch (const std::bad_alloc&) {
        throw std::runtime_error("Memory allocation failed. Sequences are too long.");
    }
    int max_score = 0;
    int max_i = 0;
    int max_j = 0;

    for (int i = 1; i <= m; ++i) {
        size_t current_row_offset = (size_t)i * row_stride;
        size_t prev_row_offset    = (size_t)(i - 1) * row_stride;

        for (int j = 1; j <= n; ++j) { 
            int score_diag = H[prev_row_offset + j - 1] + (seq1[i-1] == seq2[j-1] ? match_score : mismatch_penalty);
            int score_up   = H[prev_row_offset + j] + gap_penalty;   
            int score_left = H[current_row_offset + j - 1] + gap_penalty; 
            
            int current_val = std::max({0, score_diag, score_up, score_left});
            H[current_row_offset + j] = current_val;

            if (current_val > max_score) {
                max_score = current_val;
                max_i = i;
                max_j = j;
            }
        }
    }

    AlignmentResult result;
    result.read_id = read_id;
    result.template_id = template_id;
    result.score = max_score;
    result.read_seq_len = m;
    result.template_seq_len = n;
    
    result.read_end_pos = max_i;
    result.template_end_pos = max_j;
    
    result.num_match = 0;
    result.num_mismatch = 0;
    result.num_gap = 0;

    int i = max_i;
    int j = max_j;

    while (i > 0 && j > 0) {
        size_t current_idx = (size_t)i * row_stride + j;
        int current_score = H[current_idx];
        
        if (current_score == 0) break; 
        size_t prev_row_idx = (size_t)(i - 1) * row_stride;
        int score_from_diag = H[prev_row_idx + j - 1] + (seq1[i-1] == seq2[j-1] ? match_score : mismatch_penalty);
        int score_from_up = H[prev_row_idx + j] + gap_penalty;
        int score_from_left = H[current_idx - 1] + gap_penalty;

        if (current_score == score_from_diag) {
            if (seq1[i-1] == seq2[j-1]) {
                result.num_match++;
            } else {
                result.num_mismatch++;
            }
            i--;
            j--;
        }
        else if (current_score == score_from_up) {
            result.num_gap++;
            i--;
        }
        else {
            result.num_gap++;
            j--; 
        }
    }
    result.read_start_pos = i;
    result.template_start_pos = j;
    return result;
}

std::vector<AlignmentResult> fasta_sw_align(
    const std::string& reads_file, 
    const std::string& template_file,
    int match_score, 
    int mismatch_penalty, 
    int gap_penalty,
    int min_align_score) {
    
    std::vector<Sequence> reads_seq;
    std::vector<Sequence> template_seq;
    
    try {
        read_fasta(reads_file, reads_seq);
        read_fasta(template_file, template_seq);
    } catch (const std::runtime_error& e) {
        throw std::runtime_error("Error reading FASTA files: " + std::string(e.what()));
    }

    if (reads_seq.empty()) {
        throw std::runtime_error("No sequences found in reads file: " + reads_file);
    }

    if (template_seq.empty()) {
        throw std::runtime_error("No sequences found in template file: " + template_file);
    }

    std::vector<AlignmentResult> results;
    results.reserve(reads_seq.size());

    for (const auto& read : reads_seq) {
        for (const auto& tmpl : template_seq) {             
            AlignmentResult result = smith_waterman(
                read.id, tmpl.id,
                read.sequence, tmpl.sequence,
                match_score, mismatch_penalty, gap_penalty
            );
            
            if (result.score >= min_align_score) {
                results.push_back(result);
            }
        }
    }
    
    return results;
}

PYBIND11_MODULE(seq_align, m) {
    m.doc() = "Smith-Waterman alignment for FASTA sequences (1D Array Optimized)";
    
    pybind11::class_<AlignmentResult>(m, "AlignmentResult")
        .def_readwrite("read_id", &AlignmentResult::read_id)
        .def_readwrite("template_id", &AlignmentResult::template_id)
        .def_readwrite("score", &AlignmentResult::score)
        .def_readwrite("read_start_pos", &AlignmentResult::read_start_pos)
        .def_readwrite("read_end_pos", &AlignmentResult::read_end_pos)
        .def_readwrite("template_start_pos", &AlignmentResult::template_start_pos)
        .def_readwrite("template_end_pos", &AlignmentResult::template_end_pos)
        .def_readwrite("num_match", &AlignmentResult::num_match)
        .def_readwrite("num_mismatch", &AlignmentResult::num_mismatch)
        .def_readwrite("num_gap", &AlignmentResult::num_gap)
        .def_readwrite("read_seq_len", &AlignmentResult::read_seq_len)
        .def_readwrite("template_seq_len", &AlignmentResult::template_seq_len);

    m.def("fasta_sw_align", &fasta_sw_align, 
          "Perform Smith-Waterman alignment",
          pybind11::arg("reads_file"),
          pybind11::arg("template_file"),
          pybind11::arg("match_score") = 2,
          pybind11::arg("mismatch_penalty") = -4,
          pybind11::arg("gap_penalty") = -4,
          pybind11::arg("min_align_score") = 10);
}