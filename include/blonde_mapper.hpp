#ifndef BLONDE_MAPPER_HPP
#define BLONDE_MAPPER_HPP

#include "blonde_common.hpp"
#include "blonde_alignment.hpp"

#include <vector>

namespace blonde {

void RunMapper(
    const std::vector<Sequence>& references,
    const std::vector<Sequence>& fragments,
    unsigned int k,
    unsigned int w,
    double f,
    AlignmentType aln_type,
    int match,
    int mismatch,
    int gap,
    bool print_cigar
);

}

#endif