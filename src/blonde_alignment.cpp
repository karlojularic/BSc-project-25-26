#include "blonde_alignment.hpp"
#include <vector>
#include <algorithm>
#include <iostream>

namespace blonde {

static std::string compress_cigar(const std::string &s) {
    if (s.empty()) return "";
    std::string out;
    char cur = s[0];
    int cnt = 1;
    for (size_t i = 1; i < s.size(); i++) {
        if (s[i] == cur) cnt++;
        else {
            out += std::to_string(cnt) + cur;
            cur = s[i];
            cnt = 1;
        }
    }
    out += std::to_string(cnt) + cur;
    return out;
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
    const int n = query_len;
    const int m = target_len;

    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));
    std::vector<std::vector<char>> trace(n + 1, std::vector<char>(m + 1, 'X'));

    int max_i = n, max_j = m;
    int max_score = 0;

    if (type == AlignmentType::GLOBAL) {
        for (int i = 0; i <= n; i++) {
            dp[i][0] = i * gap;
            trace[i][0] = 'U';
        }
        for (int j = 0; j <= m; j++) {
            dp[0][j] = j * gap;
            trace[0][j] = 'L';
        }
    } else if (type == AlignmentType::SEMIGLOBAL) {
        //bez penalizacije rubova
        for (int i = 0; i <= n; i++) dp[i][0] = 0;
        for (int j = 0; j <= m; j++) dp[0][j] = 0;
    }

    //dinamicko prog
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int score_diag = dp[i - 1][j - 1] + (query[i - 1] == target[j - 1] ? match : mismatch);
            int score_up   = dp[i - 1][j] + gap;
            int score_left = dp[i][j - 1] + gap;

            int best;
            if (type == AlignmentType::LOCAL)
                best = std::max({0, score_diag, score_up, score_left});
            else
                best = std::max({score_diag, score_up, score_left});

            dp[i][j] = best;

            if (best == score_diag) trace[i][j] = 'D';
            else if (best == score_up) trace[i][j] = 'U';
            else if (best == score_left) trace[i][j] = 'L';
            else trace[i][j] = 'X';

            if (type == AlignmentType::LOCAL && best > max_score) {
                max_score = best;
                max_i = i;
                max_j = j;
            }
        }
    }

    if (type == AlignmentType::GLOBAL) {
        max_i = n; max_j = m;
        max_score = dp[n][m];
    } else if (type == AlignmentType::SEMIGLOBAL) {
        max_score = dp[n][0];
        max_i = n; max_j = 0;
        for (int j = 0; j <= m; j++) {
            if (dp[n][j] > max_score) { max_score = dp[n][j]; max_i = n; max_j = j; }
        }
        for (int i = 0; i <= n; i++) {
            if (dp[i][m] > max_score) { max_score = dp[i][m]; max_i = i; max_j = m; }
        }
    }

    //rekonstrukcija
    if (cigar) {
        std::string raw;
        int i = max_i, j = max_j;

        while (i > 0 && j > 0) {
            char move = trace[i][j];
            if (type == AlignmentType::LOCAL && dp[i][j] == 0)
                break;

            if (move == 'D') { 
                raw.push_back(query[i - 1] == target[j - 1] ? 'M' : 'X');
                i--; j--;
            } else if (move == 'U') { 
                raw.push_back('I'); //insertion u TARGETU (deletion u queryu)
                i--;
            } else if (move == 'L') { 
                raw.push_back('D'); //deletion u TARGETU (insertion u queryu)
                j--;
            } else break;
        }

        std::reverse(raw.begin(), raw.end());
        *cigar = compress_cigar(raw);

        if (target_begin) *target_begin = j;  //j pokazuje na početak segmenta u targetu !
    }

    return max_score;
}

} 