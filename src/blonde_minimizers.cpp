#include "blonde_minimizers.hpp"
#include <limits> // Za std::numeric_limits

namespace blonde {

std::vector<std::tuple<unsigned int, unsigned int, bool>> Minimize(
    const char* sequence, 
    unsigned int sequence_len,
    unsigned int kmer_len,
    unsigned int window_len) {

    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers;

    if (sequence_len < kmer_len) return minimizers;

    int num_kmers = sequence_len - kmer_len + 1;
    
    struct KmerInfo {
        unsigned int value;
        unsigned int pos;
        bool is_original;
    };
    std::vector<KmerInfo> canonical_kmers;

    for (int i = 0; i < num_kmers; ++i) {
        unsigned int current_kmer = 0;
        for (unsigned int j = 0; j < kmer_len; ++j) {
            current_kmer = (current_kmer << 2) | charTo2Bit(sequence[i + j]);
        }

        unsigned int rc_kmer = getReverseComplement(current_kmer, kmer_len);

        if (current_kmer <= rc_kmer) {
            canonical_kmers.push_back({current_kmer, (unsigned int)i, true});
        } else {
            canonical_kmers.push_back({rc_kmer, (unsigned int)i, false});
        }
    }

    if (canonical_kmers.size() < window_len) {
        window_len = canonical_kmers.size();
    }

    int last_added_pos = -1; // izbjegavanje dodavanja istog minimizera vise puta

    for (int i = 0; i <= (int)canonical_kmers.size() - (int)window_len; ++i) {
        KmerInfo min_in_window = canonical_kmers[i];

        for (unsigned int j = 1; j < window_len; ++j) {
            if (canonical_kmers[i + j].value < min_in_window.value) {
                min_in_window = canonical_kmers[i + j];
            }
        }

        // dodaj minimizer samo ako vec nije dodan
        if (last_added_pos != (int)min_in_window.pos) {
            minimizers.push_back(std::make_tuple(
                min_in_window.value, 
                min_in_window.pos, 
                min_in_window.is_original
            ));
            last_added_pos = (int)min_in_window.pos;
        }
    }

    return minimizers;
}

}