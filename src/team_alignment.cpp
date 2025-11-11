#include "team_alignment.hpp"
#include <vector>
#include <algorithm>
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
    if (type != AlignmentType::GLOBAL) {
        std::cerr << "[WARN] Only GLOBAL alignment implemented so far.\n";
    }

    const int n = query_len;
    const int m = target_len;

    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));
    std::vector<std::vector<char>> trace(n + 1, std::vector<char>(m + 1, 'X')); 
    //'D' = dijagonalno, 'U' = gore, 'L' = lijevo

    //rubovi
    for (int i = 0; i <= n; i++) {
        dp[i][0] = i * gap;
        trace[i][0] = 'U';
    }
    for (int j = 0; j <= m; j++) {
        dp[0][j] = j * gap;
        trace[0][j] = 'L';
    }
    trace[0][0] = 'X';

    //popunjavanje 
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int score_diag = dp[i - 1][j - 1] +
                (query[i - 1] == target[j - 1] ? match : mismatch);
            int score_up = dp[i - 1][j] + gap;
            int score_left = dp[i][j - 1] + gap;

            dp[i][j] = std::max({score_diag, score_up, score_left});

            if (dp[i][j] == score_diag)
                trace[i][j] = 'D';
            else if (dp[i][j] == score_up)
                trace[i][j] = 'U';
            else
                trace[i][j] = 'L';
        }
    }

    int score = dp[n][m]; 

    //dodatno
    //rekonstrukcija CIGAR stringa
    if (cigar) {
        std::string result;
        int i = n, j = m;
        while (i > 0 || j > 0) {
            if (trace[i][j] == 'D') {
                result.push_back(query[i - 1] == target[j - 1] ? 'M' : 'X');
                i--; j--;
            } else if (trace[i][j] == 'U') {
                result.push_back('D'); //deletion u targetu
                i--;
            } else if (trace[i][j] == 'L') {
                result.push_back('I'); //insertion u targetu
                j--;
            } else break;
        }
        std::reverse(result.begin(), result.end());
        *cigar = result;
    }

    if (target_begin) *target_begin = 0; 
    return score;
}

}