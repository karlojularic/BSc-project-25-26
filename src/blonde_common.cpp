#include "blonde_common.hpp"
#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <climits>

namespace blonde {

size_t calculate_N50(const std::vector<size_t> &lengths) {
    std::vector<size_t> sorted = lengths;
    std::sort(sorted.begin(), sorted.end(), std::greater<size_t>());
    size_t total = 0;
    for (auto l : sorted) {
        total += l;
    }

    size_t half_total = total / 2;
    size_t running_total = 0;
    for (auto l : sorted) {
        running_total += l;
        if (running_total >= half_total) {
            return l;
        }
    }

    return 0;
}

void PrintStats(const std::string& file1, const std::string& file2) {
    using namespace std;

    vector<Sequence> references;
    auto ref_parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(file1);
    auto ref_seqs = ref_parser->Parse(-1);  // -1 znaci parsiranje cijele datoteke

    for (auto &seq_ptr : ref_seqs) {
        references.push_back(*seq_ptr);
    }

    vector<Sequence> fragments;
    if (file2.find(".fastq") != string::npos) {
        auto frag_parser = bioparser::Parser<Sequence>::Create<bioparser::FastqParser>(file2);
        auto frag_seqs = frag_parser->Parse(-1);
        for (auto &seq_ptr : frag_seqs) {
            fragments.push_back(*seq_ptr);
        }
    } else {
        auto frag_parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(file2);
        auto frag_seqs = frag_parser->Parse(-1);
        for (auto &seq_ptr : frag_seqs) {
            fragments.push_back(*seq_ptr);
        }
    }

    cerr << "Reference sequences:\n";
    for (auto &seq : references) {
        cerr << seq.name << " : " << seq.seq.length() << "\n";
    }

    size_t total_len = 0;
    size_t min_len = SIZE_MAX;
    size_t max_len = 0;
    vector<size_t> lengths;
    for (auto &f : fragments) {
        size_t len = f.seq.length();
        lengths.push_back(len);
        total_len += len;
        if (len < min_len) min_len = len;
        if (len > max_len) max_len = len;
    }

    double avg_len = fragments.empty() ? 0 : static_cast<double>(total_len)/fragments.size();
    size_t n50 = calculate_N50(lengths);

    cerr << "\nFragment statistics:\n";
    cerr << "Number of fragments: " << fragments.size() << "\n";
    cerr << "Average length: " << avg_len << "\n";
    cerr << "Min length: " << min_len << "\n";
    cerr << "Max length: " << max_len << "\n";
    cerr << "N50 length: " << n50 << "\n";
}

} // namespace blonde