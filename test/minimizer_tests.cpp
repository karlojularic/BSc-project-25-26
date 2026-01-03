#include "gtest/gtest.h"
#include "team_minimizers.hpp"
#include <tuple>
#include <vector>

using namespace blonde;

// Test 1: Provjera na vrlo jednostavnoj sekvenci
TEST(MinimizerTest, SimpleSequence) {
    // Sekvenca: AAAAA GGGGG (dužina 10)
    // k = 5, w = 2
    // k-meri:
    // 0: AAAAA (vrijednost 0)
    // 1: AAAAG
    // 2: AAAGG
    // 3: AAGGG
    // 4: AGGGG
    // 5: GGGGG
    const char* seq = "AAAAAGGGGG";
    unsigned int k = 5;
    unsigned int w = 2;

    auto minimizers = Minimize(seq, 10, k, w);

    // Očekujemo da će AAAAA biti minimizer jer je leksikografski najmanji (sve nule)
    ASSERT_FALSE(minimizers.empty());
    
    // Prvi minimizer bi trebao biti AAAAA (vrijednost 0) na poziciji 0
    EXPECT_EQ(std::get<0>(minimizers[0]), 0); 
    EXPECT_EQ(std::get<1>(minimizers[0]), 0);
}

// Test 2: Provjera Reverse Complement logike
TEST(MinimizerTest, ReverseComplement) {
    // Sekvenca: TTTTT (k=5, w=1)
    // Originalni k-mer: TTTTT (vrijednost 1023 za k=5)
    // Reverse complement: AAAAA (vrijednost 0)
    // Program mora odabrati AAAAA jer je 0 < 1023
    const char* seq = "TTTTT";
    auto minimizers = Minimize(seq, 5, 5, 1);

    ASSERT_EQ(minimizers.size(), 1);
    EXPECT_EQ(std::get<0>(minimizers[0]), 0); // Vrijednost AAAAA
    EXPECT_FALSE(std::get<2>(minimizers[0])); // is_original mora biti false jer je uzeo RC
}