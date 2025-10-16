#include <iostream>
#include <unistd.h>
#include "version.h"

using namespace std;

int main(int argc, char *argv[]) {

    int opt;
    
    while ((opt = getopt(argc, argv, ":hv")) != -1) {
        switch(opt) {
            case 'h':
                cout << "Help: This program does something useful." << endl;
                break;
            case 'v':
                cout << "Version: " << PROJECT_VERSION << endl;
                break;
        }
    }

    return 0;
    
}