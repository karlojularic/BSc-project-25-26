#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <getopt.h>
#include "version.h"

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"

using namespace std;

struct Sequence {
    std::string name;
    std::string seq;
    std::string qual;  // za FASTQ

    // Fasta konstruktor
    Sequence(const char* n, std::uint32_t n_len,
             const char* s, std::uint32_t s_len)
        : name(n, n_len), seq(s, s_len) {}

    // Fastq konstruktor
    Sequence(const char* n, std::uint32_t n_len,
             const char* s, std::uint32_t s_len,
             const char* q, std::uint32_t q_len)
        : name(n, n_len), seq(s, s_len), qual(q, q_len) {}

    friend bioparser::FastaParser<Sequence>;
    friend bioparser::FastqParser<Sequence>;
};


void print_help() {
    cout << "Usage: blonde_mapper [options] <file1> <file2>\n";
    cout << "Options:\n";
    cout << "  -h, --help       Prikaz pomoći\n";
    cout << "  -v, --version    Prikaz verzije programa\n";
}

size_t calculate_N50(const vector<size_t> &lengths) {
    vector<size_t> sorted = lengths;
    sort(sorted.begin(), sorted.end(), greater<size_t>());
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

int main(int argc, char *argv[]) {

    const struct option long_options[] = {
        {"help",    no_argument, 0, 'h'},
        {"version", no_argument, 0, 'v'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;

    while((opt = getopt_long(argc, argv, "hv", long_options, &option_index)) != -1) {
        switch(opt) {
            case 'h':
                print_help();
                return 0;
            case 'v':
                cout << "v" << PROJECT_VERSION << "\n";
                return 0;
            default:
                print_help();
                return 1;
        }
    }

    if (optind + 2 != argc) {
        cerr << "Error: two arguments are expected.\n";
        print_help();
        return 1;
    }

    string file1 = argv[optind];
    string file2 = argv[optind + 1];
    cout << "Input file 1: " << file1 << "\n";
    cout << "Input file 2: " << file2 << "\n";

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

    return 0;
    
}