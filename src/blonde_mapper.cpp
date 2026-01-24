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
#include <cctype>

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
        uint32_t frag_pos = fragment.seq.size() - frag_pos_rc - k;

        auto it = index.find(hash);
        if (it != index.end()) {
            for (auto& hit : it->second) 
                seeds.push_back({frag_pos, hit.first, hit.second, true});
        }
    }

    if (seeds.empty()) return; //prepraviti

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
        if (c.size() >= 1) { //ovo mijenjati po potrebi za kratke primjere, orignalno >=3
            chains.push_back(c);
        }
    }

    if (chains.empty()) return;

    //novo
    const std::vector<Seed>* best_plus = nullptr;
    const std::vector<Seed>* best_minus = nullptr;

    for (const auto& c : chains) {
        if (c.empty()) continue;
        if (!c[0].is_reverse) {
            if (!best_plus || c.size() > best_plus->size()) best_plus = &c;
        } else {
            if (!best_minus || c.size() > best_minus->size()) best_minus = &c;
        }
    }

    if (!best_plus && !best_minus) return;
    //novo kraj

    //na stari nacin se pregledavaju chainovi, ali ih je vise
    //na novi nacin se bira najbolji
    auto run_extend = [&](const std::vector<Seed>* chain_ptr, bool want_reverse,
                        // output
                        int& out_score,
                        uint32_t& out_ref_id,
                        uint32_t& out_q_start, uint32_t& out_q_end,
                        uint32_t& out_t_start, uint32_t& out_t_end,
                        uint32_t& out_nmatch, uint32_t& out_aln_len,
                        std::string& out_cigar) -> bool
    {
        if (!chain_ptr || chain_ptr->empty()) return false;

        const auto& chain = *chain_ptr;
        out_ref_id = chain[0].ref_id;

        const std::string& frag_seq_used = want_reverse ? frag_rc : fragment.seq;

        //dosl sve staro
        uint32_t frag_min = UINT32_MAX, frag_max = 0;
        uint32_t ref_min = UINT32_MAX, ref_max = 0;

        for (const auto& s : chain) {
            frag_min = std::min(frag_min, s.frag_pos);
            frag_max = std::max(frag_max, s.frag_pos + k);
            ref_min = std::min(ref_min, s.ref_pos);
            ref_max = std::max(ref_max, s.ref_pos + k);
        }

        const int PAD = 3 * (int)k;

        int fs = std::max<int>(0, (int)frag_min - PAD);
        int fe = std::min<int>((int)frag_seq_used.size(), (int)frag_max + PAD);

        int rs = std::max<int>(0, (int)ref_min - PAD);
        int re = std::min<int>((int)references[out_ref_id].seq.size(), (int)ref_max + PAD);

        std::string frag_sub = frag_seq_used.substr(fs, fe - fs);
        std::string ref_sub = references[out_ref_id].seq.substr(rs, re - rs);

        unsigned int target_begin = 0;
        std::string cigar;

        int score = Align(
            frag_sub.data(), frag_sub.size(),
            ref_sub.data(), ref_sub.size(),
            aln_type, match, mismatch, gap,
            &cigar, &target_begin
        );

        if (cigar.empty()) return false;

        uint32_t query_aligned = 0, target_aligned = 0, aln_len = 0, nmatch = 0;
        uint32_t num = 0;

        for (char c : cigar) {
            if (std::isdigit((unsigned char)c)) {
                num = num * 10 + (c - '0');
            } else {
                switch (c) {
                    case '=': query_aligned += num; target_aligned += num; aln_len += num; nmatch += num; break;
                    case 'X':
                    case 'M': query_aligned += num; target_aligned += num; aln_len += num; break;
                    case 'I': query_aligned += num; aln_len += num; break;
                    case 'D': target_aligned += num; aln_len += num; break;
                }
                num = 0;
            }
        }

        uint32_t q_start = fs;
        uint32_t q_end = fs + query_aligned;

        uint32_t t_start = rs + target_begin;
        uint32_t t_end = t_start + target_aligned;

        if (want_reverse) {
            uint32_t new_q_start = (uint32_t)fragment.seq.size() - q_end;
            uint32_t new_q_end = (uint32_t)fragment.seq.size() - q_start;
            q_start = new_q_start;
            q_end = new_q_end;
        }
        //dosl sve staro kraj

        out_score = score;
        out_q_start = q_start; out_q_end = q_end;
        out_t_start = t_start; out_t_end = t_end;
        out_nmatch = nmatch; out_aln_len = aln_len;
        out_cigar = cigar;

        return true;
    };

    // rezultati za plus i minus //novo
    bool ok_p=false, ok_m=false;

    int score_p=0, score_m=0;
    uint32_t ref_p=0, ref_m=0;
    uint32_t qs_p=0, qe_p=0, ts_p=0, te_p=0, nm_p=0, al_p=0;
    uint32_t qs_m=0, qe_m=0, ts_m=0, te_m=0, nm_m=0, al_m=0;
    std::string cigar_p, cigar_m;

    if (best_plus) ok_p = run_extend(best_plus, false, score_p, ref_p, qs_p, qe_p, ts_p, te_p, nm_p, al_p, cigar_p);
    if (best_minus) ok_m = run_extend(best_minus, true, score_m, ref_m, qs_m, qe_m, ts_m, te_m, nm_m, al_m, cigar_m);

    if (!ok_p && !ok_m) return;

    // odabir najboljeg (minimap2, nakon extend-a) //novo
    bool take_minus = false;
    if (!ok_p) take_minus = true;
    else if (!ok_m) take_minus = false;
    else {
        if (score_m != score_p) take_minus = (score_m > score_p);
        else if (nm_m != nm_p) take_minus = (nm_m > nm_p);
        else if (al_m != al_p) take_minus = (al_m < al_p);
        else take_minus = false; // total tie
    }

    int mapq = 255;
    // mapq je stari
    if (ok_p && ok_m) {
        int best = take_minus ? score_m : score_p;
        int second = take_minus ? score_p : score_m;
        double ratio = (double)second / (double)best;
        mapq = (int)(254.0 * (1.0 - ratio));
    } else {
        mapq = 254;
    }

    //novo - samo odabir
    bool is_reverse = take_minus;

    uint32_t ref_id = is_reverse ? ref_m : ref_p;
    uint32_t q_start = is_reverse ? qs_m : qs_p;
    uint32_t q_end = is_reverse ? qe_m : qe_p;
    uint32_t t_start = is_reverse ? ts_m : ts_p;
    uint32_t t_end = is_reverse ? te_m : te_p;
    uint32_t nmatch = is_reverse ? nm_m : nm_p;
    uint32_t aln_len = is_reverse ? al_m : al_p;
    std::string cigar = is_reverse ? cigar_m : cigar_p;

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

    if (print_cigar) std::cout << "\tcg:Z:" << cigar;
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