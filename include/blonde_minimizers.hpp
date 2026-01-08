#ifndef BLONDE_MINIMIZERS_HPP
#define BLONDE_MINIMIZERS_HPP

#include <vector>
#include <tuple>
#include <string>

namespace blonde {

std::vector<std::tuple<unsigned int, unsigned int, bool>> Minimize(
    const char* sequence, 
    unsigned int sequence_len,
    unsigned int kmer_len,
    unsigned int window_len);

inline unsigned char charTo2Bit(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 0;
    }
}

inline unsigned int getReverseComplement(unsigned int kmer, unsigned int k) {
    unsigned int rc = 0;
    for (unsigned int i = 0; i < k; ++i) {
        unsigned int baza = kmer & 0x03; 
    
        unsigned int komplement = 3 - baza;
        
        rc = (rc << 2) | komplement;

        kmer >>= 2;
    }
    return rc;
}

}

#endif