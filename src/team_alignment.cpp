#include "team_alignment.hpp"
#include <iostream>

namespace blonde {

int Align(
    const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    AlignmentType type,
    int match,
    int mismatch,
    int gap,
    std::string* cigar,
    unsigned int* target_begin)
{
    std::cerr << "[DEBUG] Called Align() with sequences of length "
              << query_len << " and " << target_len << "\n";
    std::cerr << "  Type: "
              << (type == AlignmentType::GLOBAL ? "GLOBAL" :
                  type == AlignmentType::LOCAL ? "LOCAL" : "SEMIGLOBAL")
              << "\n";
    std::cerr << "  Parameters: match=" << match
              << ", mismatch=" << mismatch
              << ", gap=" << gap << "\n";

    if (cigar) *cigar = "10M"; // mock vrijednost
    if (target_begin) *target_begin = 0;

    // privremeno
    return 42;
}

}