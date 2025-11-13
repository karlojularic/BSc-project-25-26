#include "blonde_alignment.hpp"
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>
#include <limits>

namespace blonde {

namespace {

enum class Parent : std::uint8_t {
    NONE  = 0,
    DIAG  = 1,  // match/mismatch (M)
    UP    = 2,  // gap u targetu  -> I 
    LEFT  = 3   // gap u queryju -> D 
};

struct Cell {
    int score;
    Parent parent;
};

inline int score_match(char q, char t, int match, int mismatch) {
    return (q == t) ? match : mismatch;
}

std::string BuildCigar(const std::string& ops) {
    if (ops.empty()) {
        return "*"; // nema poravnanja
    }

    std::string cigar;
    char current = ops[0];
    unsigned int count = 1;

    for (std::size_t i = 1; i < ops.size(); ++i) {
        if (ops[i] == current) {
            ++count;
        } else {
            cigar += std::to_string(count);
            cigar += current;
            current = ops[i];
            count = 1;
        }
    }
    cigar += std::to_string(count);
    cigar += current;

    return cigar;
}

} 

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
    const unsigned int n = query_len;
    const unsigned int m = target_len;

    if (n == 0 || m == 0) {
        if (cigar) *cigar = "*";
        if (target_begin) *target_begin = 0;
        return 0;
    }

    // --------------------- GLOBAL ALIGNMENT ---------------------------
    if (type == AlignmentType::GLOBAL) {

        std::vector<std::vector<Cell>> dp(n+1, std::vector<Cell>(m+1));

        dp[0][0] = {0, Parent::NONE};

        for (unsigned int i=1;i<=n;i++) {
            dp[i][0].score  = dp[i-1][0].score + gap;
            dp[i][0].parent = Parent::UP;
        }
        for (unsigned int j=1;j<=m;j++) {
            dp[0][j].score  = dp[0][j-1].score + gap;
            dp[0][j].parent = Parent::LEFT;
        }

        // Fill
        for (unsigned int i=1;i<=n;i++) {
            for (unsigned int j=1;j<=m;j++) {

                int diag = dp[i-1][j-1].score +
                           score_match(query[i-1],target[j-1],match,mismatch);
                int up   = dp[i-1][j].score + gap;
                int left = dp[i][j-1].score + gap;

                int best = diag;
                Parent p = Parent::DIAG;

                if (up > best) { best = up; p = Parent::UP; }
                if (left > best) { best = left; p = Parent::LEFT; }

                dp[i][j].score = best;
                dp[i][j].parent = p;
            }
        }

        // Traceback (0,0) → (n,m)
        unsigned int i = n, j = m;
        std::string ops_rev;

        while (i>0 || j>0) {
            Parent p = dp[i][j].parent;
            if (p == Parent::DIAG) { ops_rev.push_back('M'); --i; --j; }
            else if (p == Parent::UP) { ops_rev.push_back('I'); --i; }
            else if (p == Parent::LEFT) { ops_rev.push_back('D'); --j; }
            else break;
        }

        if (target_begin) *target_begin = j;

        if (cigar) {
            std::string ops(ops_rev.rbegin(), ops_rev.rend());
            *cigar = BuildCigar(ops);
        }

        return dp[n][m].score;
    }

    // -------------------- SEMI-GLOBAL ALIGNMENT -----------------------
    if (type == AlignmentType::SEMIGLOBAL) {

        std::vector<std::vector<Cell>> dp(n+1, std::vector<Cell>(m+1));

        for (unsigned int i=0;i<=n;i++) {
            dp[i][0] = {0, Parent::NONE};
        }
        for (unsigned int j=0;j<=m;j++) {
            dp[0][j] = {0, Parent::NONE};
        }

        // Fill
        for (unsigned int i=1;i<=n;i++) {
            for (unsigned int j=1;j<=m;j++) {
                int diag = dp[i-1][j-1].score +
                           score_match(query[i-1],target[j-1],match,mismatch);
                int up   = dp[i-1][j].score + gap;
                int left = dp[i][j-1].score + gap;

                int best = diag;
                Parent p = Parent::DIAG;

                if (up > best)  { best = up;  p = Parent::UP; }
                if (left > best){ best = left; p = Parent::LEFT; }

                dp[i][j].score = best;
                dp[i][j].parent = p;
            }
        }

        //  goal cell: max last row iili last column
        int best = std::numeric_limits<int>::min();
        unsigned int gi = n, gj = m;

        for (unsigned int j=0;j<=m;j++) {
            if (dp[n][j].score > best) {
                best = dp[n][j].score;
                gi = n;
                gj = j;
            }
        }
        for (unsigned int i=0;i<=n;i++) {
            if (dp[i][m].score > best) {
                best = dp[i][m].score;
                gi = i;
                gj = m;
            }
        }

        // Traceback until parent==NONE
        unsigned int i = gi, j = gj;
        std::string ops_rev;

        while ((i>0 || j>0) && dp[i][j].parent != Parent::NONE) {
            if (i == 0 || j == 0)
                break;

            Parent p = dp[i][j].parent;
            if (p == Parent::DIAG) { ops_rev.push_back('M'); --i; --j; }
            else if (p == Parent::UP) { ops_rev.push_back('I'); --i; }
            else if (p == Parent::LEFT) { ops_rev.push_back('D'); --j; }
        }

        if (target_begin) *target_begin = j;

        if (cigar) {
            std::string ops(ops_rev.rbegin(), ops_rev.rend());
            *cigar = BuildCigar(ops);
        }

        return best;
    }

    // ----------------------- LOCAL ALIGNMENT --------------------------
    if (type == AlignmentType::LOCAL) {

        std::vector<std::vector<Cell>> dp(n+1, std::vector<Cell>(m+1));

        for (unsigned int i=0;i<=n;i++) {
            dp[i][0] = {0, Parent::NONE};
        }
        for (unsigned int j=0;j<=m;j++) {
            dp[0][j] = {0, Parent::NONE};
        }

        int best_score = 0;
        unsigned int bi = 0, bj = 0;

        // Fill
        for (unsigned int i=1;i<=n;i++) {
            for (unsigned int j=1;j<=m;j++) {
                int diag = dp[i-1][j-1].score +
                           score_match(query[i-1],target[j-1],match,mismatch);
                int up   = dp[i-1][j].score + gap;
                int left = dp[i][j-1].score + gap;

                int cell = std::max({0, diag, up, left});
                Parent p = Parent::NONE;

                if (cell == diag) p = Parent::DIAG;
                else if (cell == up) p = Parent::UP;
                else if (cell == left) p = Parent::LEFT;

                dp[i][j].score = cell;
                dp[i][j].parent = p;

                if (cell > best_score) {
                    best_score = cell;
                    bi = i;
                    bj = j;
                }
            }
        }

        // v1 traceback
        if (cigar) {
            std::string raw;
            unsigned int i = bi, j = bj;
            unsigned int j_start = j;

            while (i > 0 && j > 0) {
                Parent p = dp[i][j].parent;

                if (dp[i][j].score == 0) {
                    j_start = j;
                    break;
                }

                if (p == Parent::DIAG) {
                    raw.push_back(query[i - 1] == target[j - 1] ? 'M' : 'X');
                    --i; 
                    --j;
                } 
                else if (p == Parent::UP) {
                    raw.push_back('I');
                    --i;
                } 
                else if (p == Parent::LEFT) {
                    raw.push_back('D');
                    --j;
                } 
                else {
                    break;
                }
            }

            std::reverse(raw.begin(), raw.end());
            *cigar = BuildCigar(raw);

            if (target_begin)
                *target_begin = j_start;
        }

        return best_score;

    }

    return 0;
}

} // namespace blonde