#include <iostream>
#include <getopt.h>
#include <vector>
#include <random>
#include <algorithm>
#include "version.h"
#include "blonde_common.hpp"
#include "blonde_alignment.hpp"
#include "blonde_mapper.hpp" 

using namespace std;

void print_help() {
    cout << "Usage: blonde_mapper [options] <file1> <file2>\n";
    cout << "Options:\n";
    cout << "  -h, --help       Prikaz pomoći\n";
    cout << "  -v, --version    Prikaz verzije programa\n";
    cout << "  -k <int>         K-mer size (default: 15)\n";
    cout << "  -w <int>         Window size (default: 5)\n";
    cout << "  -f <double>      Filter fraction (default: 0.001)\n";
    cout << "  -t <int>         Number of threads (default: 1)\n";
    cout << "  -c               Print CIGAR strings\n";
    cout << "  -a <type>        Alignment type: global|local|semiglobal (default: semiglobal)\n";
    cout << "  -m <int>         Match score (default: 2)\n";
    cout << "  -n <int>         Mismatch cost (default: 1)\n";
    cout << "  -g <int>         Gap cost (default: 2)\n";
}

static blonde::AlignmentType parse_aln_type(const string& s) {
    string t = s;
    transform(t.begin(), t.end(), t.begin(), ::tolower);
    if (t == "global") return blonde::AlignmentType::GLOBAL;
    if (t == "local") return blonde::AlignmentType::LOCAL;
    if (t == "semiglobal" || t == "semi-global" || t == "semi") return blonde::AlignmentType::SEMIGLOBAL;

    cerr << "Error: unknown alignment type '" << s << "'. Use global|local|semiglobal.\n";
    return blonde::AlignmentType::SEMIGLOBAL;
}

static void report_random_alignment(
    const vector<blonde::Sequence>& fragments,
    blonde::AlignmentType aln_type,
    int match,
    int mismatch,
    int gap
) {
    vector<size_t> idx;
    idx.reserve(fragments.size());
    for (size_t i = 0; i < fragments.size(); ++i) {
        if (fragments[i].seq.size() <= 5000) idx.push_back(i);
    }

    if (idx.size() < 2) {
        cerr << "[random-align] Not enough sequences (len<=5000) in second input file.\n";
        return;
    }

    std::mt19937 rng(std::random_device{}());
    std::shuffle(idx.begin(), idx.end(), rng);

    const auto& A = fragments[idx[0]];
    const auto& B = fragments[idx[1]];

    string cigar;
    unsigned int target_begin = 0;

    int score = blonde::Align(
        A.seq.c_str(), (unsigned int)A.seq.size(),
        B.seq.c_str(), (unsigned int)B.seq.size(),
        aln_type,
        match, mismatch, gap,
        &cigar, &target_begin
    );

    cerr << "-----------------------------\n";
    cerr << "RANDOM ALIGNMENT (from 2nd input file)\n";
    cerr << "seq1: " << A.name << " len=" << A.seq.size() << "\n";
    cerr << "seq2: " << B.name << " len=" << B.seq.size() << "\n";
    cerr << "score: " << score << "\n";
    cerr << "target_begin: " << target_begin << "\n";
    cerr << "cigar: " << cigar << "\n";
    cerr << "-----------------------------\n";
}

int main(int argc, char *argv[]) {
    unsigned int k = 15;
    unsigned int w = 5;
    double f = 0.001;
    int num_threads = 3;
    bool print_cigar = false;

    blonde::AlignmentType aln_type = blonde::AlignmentType::SEMIGLOBAL;
    int match = 2;
    int mismatch = -1;
    int gap = -2;

    const struct option long_options[] = {
        {"help",    no_argument, 0, 'h'},
        {"version", no_argument, 0, 'v'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;
    while((opt = getopt_long(argc, argv, "hvk:w:f:t:ca:m:n:g:", long_options, &option_index)) != -1) {
        switch(opt) {
            case 'h':
                print_help();
                return 0;
            case 'v':
                cout << "v" << PROJECT_VERSION << "\n";
                return 0;
            case 'k': k = stoi(optarg); break;
            case 'w': w = stoi(optarg); break;
            case 'f': f = stod(optarg); break;
            case 't': num_threads = stoi(optarg); break;
            case 'c': print_cigar = true; break;
            case 'a': aln_type = parse_aln_type(optarg); break;
            case 'm': match = stoi(optarg); break;
            case 'n': mismatch = stoi(optarg); break;
            case 'g': gap = stoi(optarg); break;

            default: print_help(); return 1;
        }
    }

    if (mismatch > 0) mismatch = -mismatch;
    if (gap > 0) gap = -gap;

    if (optind + 2 != argc) {
        cerr << "Error: reference and fragments files are expected.\n";
        print_help();
        return 1;
    }

    string file1 = argv[optind];
    string file2 = argv[optind + 1];

    blonde::PrintStats(file1, file2);
    
    vector<blonde::Sequence> references = blonde::LoadSequences(file1);
    vector<blonde::Sequence> fragments = blonde::LoadSequences(file2);

    if (references.empty() || fragments.empty()) {
        cerr << "Error: Could not load sequences from provided files." << endl;
        return 1;
    }

    // Align two random sequences from second input file (len <= 5000) and report
    // zadano u projektu pod Alignment
    report_random_alignment(fragments, aln_type, match, mismatch, gap);

    blonde::RunMapper(
        references, 
        fragments, 
        k, w, f, 
        aln_type,
        match, mismatch, gap,
        print_cigar,
        num_threads
    );

    return 0;
}
