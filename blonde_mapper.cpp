#include <iostream>
#include <unistd.h>
#include "version.h"

using namespace std;

void print_help() {
    cout << "Usage: blonde_mapper [options] <file1> <file2>\n";
    cout << "Options:\n";
    cout << "  -h, --help       Prikaz pomoći\n";
    cout << "  -v, --version    Prikaz verzije programa\n";
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
                cout << "blonde_mapper version " << PROJECT_VERSION << "\n";
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