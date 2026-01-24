#include <iostream>
#include <getopt.h>
#include <vector>
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
}

int main(int argc, char *argv[]) {
    unsigned int k = 15;
    unsigned int w = 5;
    double f = 0.001;
    int num_threads = 3;
    bool print_cigar = false;

    const struct option long_options[] = {
        {"help",    no_argument, 0, 'h'},
        {"version", no_argument, 0, 'v'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;
    while((opt = getopt_long(argc, argv, "hvk:w:f:t:c", long_options, &option_index)) != -1) {
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
            default: print_help(); return 1;
        }
    }

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

    blonde::RunMapper(
        references, 
        fragments, 
        k, w, f, 
        blonde::AlignmentType::LOCAL,
        2, -1, -2,
        print_cigar,
        num_threads
    );

    return 0;
}