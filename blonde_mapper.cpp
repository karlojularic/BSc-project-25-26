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

    return 0;
    
}