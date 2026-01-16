#pragma once
#include <string>
#include <vector>
#include <cstdint>

namespace blonde {

struct Sequence {
    std::string name;
    std::string seq;
    std::string qual;

    Sequence(const char* n, std::uint32_t n_len,
             const char* s, std::uint32_t s_len)
        : name(n, n_len), seq(s, s_len) {}

    Sequence(const char* n, std::uint32_t n_len,
             const char* s, std::uint32_t s_len,
             const char* q, std::uint32_t q_len)
        : name(n, n_len), seq(s, s_len), qual(q, q_len) {}
};

void PrintStats(const std::string& file1, const std::string& file2);

std::vector<Sequence> LoadSequences(const std::string& path);

}