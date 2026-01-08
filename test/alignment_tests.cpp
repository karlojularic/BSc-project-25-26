#include "gtest/gtest.h"
#include "blonde_alignment.hpp" 
#include <string>

using namespace blonde; 

// test 1: savrseno poklapanje
TEST(GlobalAlignmentTest, PerfectMatch) {
    const char* query = "GATTACA";
    const char* target = "GATTACA";
    int match = 2, mismatch = -1, gap = -3;
    std::string cigar;
    unsigned int begin_pos = 999;
    
    // ocekivano: 7 * 2 = 14
    int expected_score = 14;
    // ocekivani neformatirani CIGAR: M M M M M M M
    std::string expected_raw_cigar = "MMMMMMM"; 

    int actual_score = Align(query, 7, target, 7, 
                             AlignmentType::GLOBAL, 
                             match, mismatch, gap, 
                             &cigar, &begin_pos);

    ASSERT_EQ(actual_score, expected_score);
    ASSERT_EQ(cigar, expected_raw_cigar);
    ASSERT_EQ(begin_pos, 0); // globalno poravnanje uvijek pocinje na 0
}

// test 2: mismatch
TEST(GlobalAlignmentTest, SimpleMismatch) {
    const char* query = "GATTACA";
    const char* target = "GACTACA"; // mismatch na poziciji 2 (T vs C)
    int match = 2, mismatch = -1, gap = -3;
    std::string cigar;
    
    // ocekivano: 6 match * 2 + 1 mismatch * (-1) = 12 - 1 = 11
    int expected_score = 11;
    // ocekivani neformatirani CIGAR (koristeci 'X' za mismatch): M M X M M M M
    std::string expected_raw_cigar = "MMXMMMM";

    int actual_score = Align(query, 7, target, 7, 
                             AlignmentType::GLOBAL, 
                             match, mismatch, gap, 
                             &cigar);

    ASSERT_EQ(actual_score, expected_score);
    ASSERT_EQ(cigar, expected_raw_cigar);
}

// test 3: umetanje/brisanje (gap)
TEST(GlobalAlignmentTest, SimpleGap) {
    const char* query = "ACGTA"; // N=5
    const char* target = "ACGATA"; // M=6, (G je dodano u target)
    int match = 1, mismatch = -2, gap = -3;
    std::string cigar;
    
    int expected_score = 2;
    std::string expected_raw_cigar = "MMMIM";

    int actual_score = Align(query, 5, target, 6, 
                             AlignmentType::GLOBAL, 
                             match, mismatch, gap, 
                             &cigar);

    ASSERT_EQ(actual_score, expected_score);
    ASSERT_EQ(cigar, expected_raw_cigar);
}

// test 4: ugradeni sub-string
TEST(LocalAlignmentTest, EmbeddedMatch) {
    const char* query = "GTCAG";
    const char* target = "GATTACAGA"; // T C A G je unutra (index 3)
    int match = 3, mismatch = -1, gap = -2;
    std::string cigar;
    unsigned int begin_pos = 999;
    
    // Query T C A G (4 bp) vs Target T A C A G (4 bp)
    // ocekivano: A C A G (4 bp)
    // Score: 4 Match * 3 = 12 (ako je GTCAG query, a TACA dio targeta)
    // ako se uzme najbolje poklapanje "TCA" unutar GATTACAGA vs "GTCAG":
    // Target: A C A
    // Query:  T C A
    // Score: 2M*3 + 1X*(-1) = 6 - 1 = 5.
    
    const char* q = "CAT";
    const char* t = "GGGCATTT"; // CAT je na indexu 3 (j=3)
    
    // ocekivano: CAT / CAT. Score = 3 * 3 = 9.
    int expected_score = 9;
    unsigned int expected_begin = 3; 
    std::string expected_raw_cigar = "MMM";

    int actual_score = Align(q, 3, t, 8, 
                             AlignmentType::LOCAL, 
                             match, mismatch, gap, 
                             &cigar, &begin_pos);

    ASSERT_EQ(actual_score, expected_score);
    ASSERT_EQ(cigar, expected_raw_cigar);
    ASSERT_EQ(begin_pos, expected_begin); 
}

// test 5: prefix alignment
TEST(SemiGlobalAlignmentTest, PrefixTarget) {
    const char* query = "AATT"; // N=4
    const char* target = "AATTGGCC"; // M=8
    int match = 2, mismatch = -1, gap = -3;
    std::string cigar;
    unsigned int begin_pos = 999;

    // ocekivano: AATT na pocetku targeta, ostatak targeta (GGCC) se ne kaznjava.
    // Score: 4 Match * 2 = 8
    int expected_score = 8;
    unsigned int expected_begin = 0; 
    std::string expected_raw_cigar = "MMMM"; 

    int actual_score = Align(query, 4, target, 8, 
                             AlignmentType::SEMIGLOBAL, 
                             match, mismatch, gap, 
                             &cigar, &begin_pos);

    ASSERT_EQ(actual_score, expected_score);
    ASSERT_EQ(cigar, expected_raw_cigar);
    ASSERT_EQ(begin_pos, expected_begin);
}

// test 6: suffix alignment
TEST(SemiGlobalAlignmentTest, SuffixTarget) {
    const char* query = "GGCC"; // N=4
    const char* target = "AATTGGCC"; // M=8
    int match = 2, mismatch = -1, gap = -3;
    std::string cigar;
    unsigned int begin_pos = 999;

    // ocekivano: GGCC na kraju targeta, pocetak targeta (AATT) se ne kaznjava
    // Score: 4 Match * 2 = 8
    int expected_score = 8;
    unsigned int expected_begin = 4; // target index 4 je G
    std::string expected_raw_cigar = "MMMM"; 

    int actual_score = Align(query, 4, target, 8, 
                             AlignmentType::SEMIGLOBAL, 
                             match, mismatch, gap, 
                             &cigar, &begin_pos);

    ASSERT_EQ(actual_score, expected_score);
    ASSERT_EQ(cigar, expected_raw_cigar);
    ASSERT_EQ(begin_pos, expected_begin);
}

// dodatni test za provjeru null pointera
TEST(AlignmentTest, NullPointers) {
    const char* query = "ACGT";
    const char* target = "ACGT";
    
    // Očekivano: 4 * 1 = 4
    int actual_score = Align(query, 4, target, 4, 
                             AlignmentType::GLOBAL, 
                             1, -1, -1, 
                             nullptr, nullptr); // nema CIGAR-a, nema pocetne pozicije

    ASSERT_EQ(actual_score, 4);
}