#include "blonde_mapper.hpp"
#include "blonde_common.hpp"
#include "blonde_minimizers.hpp"
#include "blonde_alignment.hpp"
#include <unordered_map>
#include <vector>
#include <tuple>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif


namespace blonde {

// jedan seed hit
struct Seed {
    uint32_t frag_pos;
    uint32_t ref_id;
    uint32_t ref_pos;
    bool is_reverse;
};

using MinimizerIndex =
    std::unordered_map<uint32_t, std::vector<std::pair<uint32_t, uint32_t>>>;

// ---------- 1. izgradnja indeksa ---------- 

MinimizerIndex BuildReferenceIndex(
    const std::vector<Sequence>& references,
    unsigned int k,
    unsigned int w
) {
    MinimizerIndex index;

    for (uint32_t rid = 0; rid < references.size(); ++rid) {
        const auto& ref = references[rid];

        auto mins = Minimize(
            ref.seq.c_str(),
            ref.seq.size(),
            k,
            w
        );

        for (auto& m : mins) {
            uint32_t hash = std::get<0>(m);
            uint32_t pos = std::get<1>(m);

            index[hash].push_back({rid, pos});
        }
    }

    return index;
}

//pomocna fija1
std::string ReverseComplement(const std::string& s) {
    std::string rc;
    rc.reserve(s.size());
    for (auto it = s.rbegin(); it != s.rend(); ++it) {
        switch (*it) {
            case 'A': case 'a': rc.push_back('T'); break;
            case 'C': case 'c': rc.push_back('G'); break;
            case 'G': case 'g': rc.push_back('C'); break;
            case 'T': case 't': rc.push_back('A'); break;
            default: rc.push_back('N');
        }
    }
    return rc;
}

//pomocna fija2
void FilterFrequentMinimizers(MinimizerIndex& index, double f) { //pomocna fija koja uklanja minimizere iz indeksa
    if (index.empty() || f <= 0.0) return;
    std::vector<std::pair<uint32_t, size_t>> freqs;
    freqs.reserve(index.size());

    for (const auto& kv : index) {
        freqs.push_back({kv.first, kv.second.size()});
    }

    std::sort(freqs.begin(), freqs.end(),
              [](const auto& a, const auto& b) {
                  return a.second < b.second;
              });

    size_t keep = static_cast<size_t>((1.0 - f) * freqs.size());

    size_t min_to_keep = std::min<size_t>(index.size(), 100);
    if (keep < min_to_keep) {
        keep = min_to_keep;
    }
  

    for (size_t i = keep; i < freqs.size(); ++i) {
        index.erase(freqs[i].first);
    }
}

//pomocna fija3 LIS
std::vector<Seed> ComputeLIS(const std::vector<Seed>& hits, unsigned int k) {
    int n = hits.size();
    if (n == 0) return {};

    std::vector<int> dp(n, 1);
    std::vector<int> parent(n, -1);

    int best_len = 1;
    int best_end = 0;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            if (hits[j].frag_pos < hits[i].frag_pos &&
                hits[j].ref_pos < hits[i].ref_pos) {

                int d1 = (int)hits[j].ref_pos - (int)hits[j].frag_pos;
                int d2 = (int)hits[i].ref_pos - (int)hits[i].frag_pos;

                if (std::abs(d1 - d2) > (int)k) continue;

                if (dp[j] + 1 > dp[i]) {
                    dp[i] = dp[j] + 1;
                    parent[i] = j;
                }
            }
        }

        if (dp[i] > best_len) {
            best_len = dp[i];
            best_end = i;
        }
    }

    std::vector<Seed> chain;
        for (int i = best_end; i != -1; i = parent[i]) {
            chain.push_back(hits[i]);
        }
    std::reverse(chain.begin(), chain.end());
    return chain;
}

//pomocna fija4 Chain Score
static long long ChainScore(const std::vector<Seed>& chain, unsigned int k) {
    if (chain.empty()) return std::numeric_limits<long long>::min();

    // Tunables
    const long long SEED_BONUS  = 100; // bodovi
    const long long GAP_PENALTY = 1;   // kazna
    const long long DIAG_PENALTY= 2;   // kazna
dr-dq
    long long score = 0;
    score += SEED_BONUS * (long long)chain.size();

    for (size_t i = 1; i < chain.size(); ++i) {
        int dq = (int)chain[i].frag_pos - (int)chain[i-1].frag_pos;
        int dr = (int)chain[i].ref_pos  - (int)chain[i-1].ref_pos;

        // (mali skokovi)
        int step = std::max(dq, dr);
        int step_bins = (k ? (step + (int)k - 1) / (int)k : step);
        score -= GAP_PENALTY * (long long)step_bins;

        // konzistentna dijagonala (dr ~ dq)
        int diag = std::abs(dr - dq);
        int diag_bins = (k ? (diag + (int)k - 1) / (int)k : diag);
        score -= DIAG_PENALTY * (long long)diag_bins;
    }

    return score;
}

// ---------- 2. mapiranje jednog fragmenta ---------- 

void MapFragment(
    const Sequence& fragment,
    const std::vector<Sequence>& references,
    const MinimizerIndex& index,
    unsigned int k,
    unsigned int w,
    AlignmentType aln_type,
    int match,
    int mismatch,
    int gap,
    bool print_cigar
) {
    // ---------- 2.1. minimizeri fragmenta ----------

    auto frag_mins_fwd = Minimize(
        fragment.seq.c_str(), 
        fragment.seq.size(), 
        k, 
        w
    );

    std::string frag_rc = ReverseComplement(fragment.seq);
    auto frag_mins_rev = Minimize(
        frag_rc.c_str(), 
        frag_rc.size(), 
        k, 
        w
    );

    // ---------- 2.2. seed hitovi ----------
    std::vector<Seed> seeds;

    for (auto& m : frag_mins_fwd) {
        uint32_t hash = std::get<0>(m);
        uint32_t frag_pos = std::get<1>(m);

        auto it = index.find(hash);
        if (it != index.end()) {
            for (auto& hit : it->second) 
                seeds.push_back({frag_pos, hit.first, hit.second, false});
        }
    }

    for (auto& m : frag_mins_rev) {
        uint32_t hash = std::get<0>(m);
        uint32_t frag_pos_rc = std::get<1>(m);

        // frag_pos ostaje u RC koordinatama (indeks u frag_rc)
        uint32_t frag_pos = frag_pos_rc;

        auto it = index.find(hash);
        if (it != index.end()) {
            for (auto& hit : it->second)
                seeds.push_back({frag_pos, hit.first, hit.second, true});
        }
    }

    if (seeds.empty()) return;

    // ---------- 2.3. chaining (LIS) ----------
    const int DIAG_BAND = 2*k;

    std::sort(seeds.begin(), seeds.end(),
    [](const Seed& a, const Seed& b) {
        if (a.is_reverse != b.is_reverse) return a.is_reverse < b.is_reverse;
        if (a.ref_id != b.ref_id) return a.ref_id < b.ref_id;
        int da = (int)a.ref_pos - (int)a.frag_pos;
        int db = (int)b.ref_pos - (int)b.frag_pos;
        if (da != db) return da < db;
        return a.frag_pos < b.frag_pos;
    });

    std::vector<std::vector<Seed>> diag_groups;

    for (const auto& s : seeds) {
        int d = (int)s.ref_pos - (int)s.frag_pos;

        if (diag_groups.empty()) {
            diag_groups.push_back({s});
            continue;
        }

        auto& last_group = diag_groups.back();
        int last_d = (int)last_group.back().ref_pos
                - (int)last_group.back().frag_pos;

        if (std::abs(d - last_d) <= DIAG_BAND &&
            s.ref_id == last_group.back().ref_id &&
            s.is_reverse == last_group.back().is_reverse) {

            last_group.push_back(s);
        } else {
            diag_groups.push_back({s});
        }
    }

    std::vector<std::vector<Seed>> chains;

    for (const auto& group : diag_groups) {
        if (group.size() < 2) continue;

        auto c = ComputeLIS(group, k);
        if (c.size() >= 1) { //ovo mijenjati po potrebi za kratke primjere, originalno >=3
            chains.push_back(c);
        }
    }

    if (chains.empty()) return;

    // uzimam najbolji chain po chaining score
    const std::vector<Seed>* best_plus = nullptr;
    const std::vector<Seed>* best_minus = nullptr;
    long long best_plus_score = std::numeric_limits<long long>::min();
    long long best_minus_score = std::numeric_limits<long long>::min();

    for (const auto& c : chains) {
        if (c.empty()) continue;
        long long sc = ChainScore(c, k);

        if (!c[0].is_reverse) {
            if (!best_plus || sc > best_plus_score || (sc == best_plus_score && c.size() > best_plus->size())) {
                best_plus = &c;
                best_plus_score = sc;
            }
        } else {
            if (!best_minus || sc > best_minus_score || (sc == best_minus_score && c.size() > best_minus->size())) {
                best_minus = &c;
                best_minus_score = sc;
            }
        }
    }

    if (!best_plus && !best_minus) return;
    // std::cerr << fragment.name
    //       << " best_plus_score=" << best_plus_score
    //       << " best_minus_score=" << best_minus_score << "\n";

    // odabir strand-a po chain-scoreu
    const std::vector<Seed>* best_chain = nullptr;

    if (!best_plus) best_chain = best_minus;
    else if (!best_minus) best_chain = best_plus;
    else {
        // ako su scoreovi gotovo jednaki, uzmi onaj s vise seedova (stabilnije)
        if (best_plus_score != best_minus_score) {
            best_chain = (best_plus_score > best_minus_score) ? best_plus : best_minus;
        } else {
            best_chain = (best_plus->size() >= best_minus->size()) ? best_plus : best_minus;
        }
    }

    auto& chain = *best_chain;

    bool is_reverse = chain[0].is_reverse;
    const std::string& frag_seq_used = is_reverse ? frag_rc : fragment.seq;

    // ---------- 2.4. prozor za aligment ----------
    uint32_t ref_id = chain[0].ref_id;

    uint32_t frag_min = UINT32_MAX;
    uint32_t frag_max = 0;

    for (const auto& s : chain) {
        frag_min = std::min(frag_min, s.frag_pos);
        frag_max = std::max(frag_max, s.frag_pos + k);
    }


    uint32_t ref_min = UINT32_MAX;
    uint32_t ref_max = 0;

    for (const auto& s : chain) {
        ref_min = std::min(ref_min, s.ref_pos);
        ref_max = std::max(ref_max, s.ref_pos + k);
    }

    // padding
    const int PAD = 3 * k;

    int fs = std::max<int>(0, frag_min - PAD);
    int fe = std::min<int>(frag_seq_used.size(), frag_max + PAD);

    int rs = std::max<int>(0, ref_min - PAD);
    int re = std::min<int>(references[ref_id].seq.size(), ref_max + PAD);

    std::string frag_sub = frag_seq_used.substr(fs, fe - fs);
    std::string ref_sub = references[ref_id].seq.substr(rs, re - rs);

    // ---------- 2.5. alignment ----------
    std::string cigar;
    unsigned int target_begin = 0;

    int score = Align(
        frag_sub.data(), frag_sub.size(),
        ref_sub.data(), ref_sub.size(),
        aln_type,
        match,
        mismatch,
        gap,
        &cigar,
        &target_begin
    );

    if (cigar.empty()) return;

    // ---------- 2.6. PAF ----------  

    uint32_t query_aligned = 0;
    uint32_t target_aligned = 0;
    uint32_t aln_len = 0;
    uint32_t nmatch = 0;

    uint32_t num = 0;
    for (char c : cigar) {
        if (std::isdigit(c)) {
            num = num * 10 + (c - '0');
        } else {
            switch (c) {
                case '=':
                    query_aligned += num;
                    target_aligned += num;
                    aln_len += num;
                    nmatch += num;
                    break;
                case 'X':
                case 'M': //M zapravo nema, sve je =/X
                    query_aligned += num;
                    target_aligned += num;
                    aln_len += num;
                    break;
                case 'I':
                    query_aligned += num;
                    aln_len += num;
                    break;
                case 'D':
                    target_aligned += num;
                    aln_len += num;
                    break;
            }
            num = 0;
        }
    }

    uint32_t q_start = fs;
    uint32_t q_end = fs + query_aligned;

    uint32_t t_start = rs + target_begin;
    uint32_t t_end = t_start + target_aligned;

    if (is_reverse) {
        uint32_t new_q_start = fragment.seq.size() - q_end;
        uint32_t new_q_end = fragment.seq.size() - q_start;
        q_start = new_q_start;
        q_end = new_q_end;
    }

    uint32_t left_clip  = q_start;
    uint32_t right_clip = fragment.seq.size() - q_end;

    std::string final_cigar = cigar;
    // if (left_clip > 0)  final_cigar = std::to_string(left_clip) + "S" + final_cigar;
    // if (right_clip > 0) final_cigar += std::to_string(right_clip) + "S";

    int mapq = 255;

    #pragma omp critical

    std::cout
        << fragment.name << "\t"
        << fragment.seq.size() << "\t"
        << q_start << "\t"
        << q_end << "\t"
        << (is_reverse ? "-" : "+") << "\t"
        << references[ref_id].name << "\t"
        << references[ref_id].seq.size() << "\t"
        << t_start << "\t"
        << t_end << "\t"
        << nmatch << "\t"
        << aln_len << "\t"
        << mapq;
        if (print_cigar) {
            std::cout << "\tcg:Z:" << final_cigar;
        }
        std::cout << "\n";
}

// ---------- 3. glavni ulaz mappera ----------

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
    bool print_cigar,
    int num_threads
) {
    auto index = BuildReferenceIndex(references, k, w);
    FilterFrequentMinimizers(index, f);

    omp_set_num_threads(num_threads);

    #pragma omp parallel for schedule(dynamic)

    for (const auto& frag : fragments) {
        MapFragment(
            frag,
            references,
            index,
            k, w,
            aln_type,
            match, mismatch, gap,
            print_cigar
        );
    }
}

}