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
                hits[j].ref_pos < hits[i].ref_pos &&
                (hits[i].frag_pos - hits[j].frag_pos) >= 1 && // minimalni pomak
                (hits[i].ref_pos - hits[j].ref_pos) >= 1) {

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

// ---------- 2. mapiranje jednog fragmenta ---------- 

void MapFragment(
    const Sequence& fragment,
    const std::vector<Sequence>& references,
    const MinimizerIndex& index,
    unsigned int k,
    unsigned int w,
    double f,
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
        if (it == index.end()) continue;

        for (auto& hit : it->second) {
            uint32_t ref_id = hit.first;
            seeds.push_back({frag_pos, ref_id, hit.second, false});
        }
    }

    for (auto& m : frag_mins_rev) {
        uint32_t hash = std::get<0>(m);
        uint32_t frag_pos_rc = std::get<1>(m);
        uint32_t frag_pos = fragment.seq.size() - frag_pos_rc - k;


        auto it = index.find(hash);
        if (it == index.end()) continue;

        for (auto& hit : it->second) {
            uint32_t ref_id = hit.first;
            seeds.push_back({frag_pos, ref_id, hit.second, true});
        }
    }

    if (seeds.empty()) return;

    // ---------- 2.3. chaining (LIS) ----------
    const int DIAG_BAND = 2*k;

    std::sort(seeds.begin(), seeds.end(),
    [](const Seed& a, const Seed& b) {
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
        if (c.size() >= 3) {
            chains.push_back(c);
        }
    }

    if (chains.empty()) return;

    //mapq
    std::vector<int> chain_scores;
    for (const auto& c : chains) {
        chain_scores.push_back((int)c.size());
    }

    std::sort(chain_scores.begin(), chain_scores.end(), std::greater<int>());

    int best = chain_scores[0];
    int second = (chain_scores.size() > 1) ? chain_scores[1] : 0;

    int mapq = 0;

    if (best > 0) {
        if (second == 0) {
            mapq = 255;
        } else {
            double ratio = (double)second / (double)best;
            mapq = (int)(255.0 * (1.0 - ratio));
        }
    }

    if (mapq < 0) mapq = 0;
    if (mapq > 255) mapq = 255;

    // privremeno: uzima najduži chain
    auto& chain = *std::max_element(
        chains.begin(), chains.end(),
        [](const auto& a, const auto& b) {
            return a.size() < b.size();
        });


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
                case 'M': //M zapravo nema, sve je trenutno =/X
                    query_aligned += num;
                    target_aligned += num;
                    aln_len += num;
                    break;
                case 'I':
                    aln_len += num;
                    break;
                case 'S':
                    query_aligned += num;
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
        << mapq; //moze mozda jednostavnije 
        if (print_cigar) {
            std::cout << "\tcg:Z:" << cigar;
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
    bool print_cigar
) {
    auto index = BuildReferenceIndex(references, k, w);
    FilterFrequentMinimizers(index, f);

    for (const auto& frag : fragments) {
        MapFragment(
            frag,
            references,
            index,
            k, w, f,
            aln_type,
            match, mismatch, gap,
            print_cigar
        );
    }
}

}