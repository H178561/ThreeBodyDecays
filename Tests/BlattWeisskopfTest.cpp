#include <gtest/gtest.h>
#include "../BlattWeisskopf.hh"
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace ThreeBodyDecays;

// Test fixture for BlattWeisskopf tests
class BlattWeisskopfTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Set precision for floating point comparisons
        tolerance = 1e-10;
    }

    double tolerance;

    // Helper function to compare doubles with tolerance
    void expectNear(double expected, double actual, double tol = 1e-10)
    {
        EXPECT_NEAR(expected, actual, tol);
    }
};

// Test S-wave (L=0) - should always return 1.0
TEST_F(BlattWeisskopfTest, SWave_AlwaysReturnsOne)
{
    double d = 1.5; // radius parameter

    // Test various momentum values
    std::vector<double> momenta = {0.0, 0.1, 0.5, 1.0, 2.0, 10.0};

    for (double q : momenta)
    {
        double result = BlattWeisskopf(q, 0, d);
        expectNear(1.0, result, tolerance);
    }
}

// Test P-wave (L=1) at specific known values
TEST_F(BlattWeisskopfTest, PWave_KnownValues)
{
    double d = 1.5;

    // Test q = 0 (should give 0)
    double result_zero = BlattWeisskopf(0.0, 1, d);
    expectNear(0.0, result_zero, tolerance);

    // Test q = d (z = 1, should give sqrt(1/(1+1)) = sqrt(1/2))
    // When q = d, then z² = (q*d)² = d² = (1.5)² = 2.25, not 1!
    double result_unit = BlattWeisskopf(d, 1, d);
    double z2_actual = (d * d) * (d * d); // This is wrong in the comment
    // Correct: when q = d, z² = (q*d)² = d²
    double z2_correct = d * d; // z² = 2.25
    double expected_correct = std::sqrt(z2_correct / (1.0 + z2_correct));
    expectNear(expected_correct, result_unit, tolerance);

    // Test a specific case: q = 0.5, d = 1.0
    double q = 0.5;
    d = 1.0;
    double z2 = (q * d) * (q * d); // z² = 0.25
    double expected = std::sqrt(z2 / (1.0 + z2)); // sqrt(0.25/1.25) = sqrt(0.2)
    double result = BlattWeisskopf(q, 1, d);
    expectNear(expected, result, tolerance);
}

// Test D-wave (L=2) at specific known values
TEST_F(BlattWeisskopfTest, DWave_KnownValues)
{
    double d = 1.5;

    // Test q = 0 (should give 0)
    double result_zero = BlattWeisskopf(0.0, 2, d);
    expectNear(0.0, result_zero, tolerance);

    // Test a specific case: q = 1.0, d = 1.0
    double q = 1.0;
    d = 1.0;
    double z2 = (q * d) * (q * d); // z² = 1.0
    double z4 = z2 * z2; // z⁴ = 1.0
    double expected = std::sqrt(z4 / (9.0 + 3.0 * z2 + z4)); // sqrt(1/(9+3+1)) = sqrt(1/13)
    double result = BlattWeisskopf(q, 2, d);
    expectNear(expected, result, tolerance);
}

// Test unsupported L values (should return 1.0)
TEST_F(BlattWeisskopfTest, UnsupportedL_ReturnsOne)
{
    double q = 1.0;
    double d = 1.5;

    // Test L values not implemented (3, 4, 5, etc.)
    std::vector<int> unsupported_L = {3, 4, 5, 10, -1};

    for (int L : unsupported_L)
    {
        double result = BlattWeisskopf(q, L, d);
        expectNear(1.0, result, tolerance);
    }
}

// Test monotonic behavior: for fixed d, BlattWeisskopf should increase with q for L > 0
TEST_F(BlattWeisskopfTest, MonotonicBehavior)
{
    double d = 1.0;
    std::vector<double> momenta = {0.1, 0.5, 1.0, 1.5, 2.0};

    // Test P-wave monotonicity
    for (size_t i = 1; i < momenta.size(); ++i)
    {
        double bf_prev = BlattWeisskopf(momenta[i-1], 1, d);
        double bf_curr = BlattWeisskopf(momenta[i], 1, d);
        EXPECT_GT(bf_curr, bf_prev) << "P-wave should be monotonically increasing";
    }

    // Test D-wave monotonicity
    for (size_t i = 1; i < momenta.size(); ++i)
    {
        double bf_prev = BlattWeisskopf(momenta[i-1], 2, d);
        double bf_curr = BlattWeisskopf(momenta[i], 2, d);
        EXPECT_GT(bf_curr, bf_prev) << "D-wave should be monotonically increasing";
    }
}

// Test asymptotic behavior: for large q, BlattWeisskopf should approach 1
TEST_F(BlattWeisskopfTest, AsymptoticBehavior)
{
    double d = 1.0;
    double large_q = 100.0;

    // For P-wave: sqrt(z²/(1+z²)) → 1 as z → ∞
    double result_p = BlattWeisskopf(large_q, 1, d);
    expectNear(1.0, result_p, 1e-4); // Relaxed tolerance for asymptotic behavior

    // For D-wave: sqrt(z⁴/(9+3z²+z⁴)) → 1 as z → ∞
    double result_d = BlattWeisskopf(large_q, 2, d);
    expectNear(1.0, result_d, 1e-4); // Relaxed tolerance for asymptotic behavior
}

// Test the breakup momentum function
TEST_F(BlattWeisskopfTest, BreakupMomentum_KnownValues)
{
    // Test case where M < m1 + m2 (should return 0)
    double result_below_threshold = breakup(1.0, 0.6, 0.5); // 1.0 < 0.6 + 0.5
    expectNear(0.0, result_below_threshold, tolerance);

    // Test case at threshold (should return 0)
    double result_at_threshold = breakup(1.1, 0.6, 0.5); // 1.1 = 0.6 + 0.5
    expectNear(0.0, result_at_threshold, tolerance);

    // Test symmetric case: M = 2, m1 = m2 = 0.5
    double M = 2.0;
    double m1 = 0.5, m2 = 0.5;
    double expected = 0.5 * std::sqrt((M*M - (m1+m2)*(m1+m2)) * (M*M - (m1-m2)*(m1-m2))) / M;
    // = 0.5 * sqrt((4 - 1) * (4 - 0)) / 2 = 0.5 * sqrt(12) / 2 = sqrt(3)/2
    double result_symmetric = breakup(M, m1, m2);
    expectNear(std::sqrt(3.0)/2.0, result_symmetric, tolerance);

    // Test asymmetric case: M = 2, m1 = 0.8, m2 = 0.2
    M = 2.0;
    m1 = 0.8; m2 = 0.2;
    double lambda = (M*M - (m1+m2)*(m1+m2)) * (M*M - (m1-m2)*(m1-m2));
    // = (4 - 1) * (4 - 0.36) = 3 * 3.64 = 10.92
    expected = 0.5 * std::sqrt(lambda) / M;
    double result_asymmetric = breakup(M, m1, m2);
    expectNear(expected, result_asymmetric, tolerance);
}

// Test edge cases
TEST_F(BlattWeisskopfTest, EdgeCases)
{
    double d = 1.0;

    // Test very small momentum
    double small_q = 1e-10;
    double result_p_small = BlattWeisskopf(small_q, 1, d);
    EXPECT_GT(result_p_small, 0.0);
    EXPECT_LT(result_p_small, 0.1); // Should be very small

    double result_d_small = BlattWeisskopf(small_q, 2, d);
    EXPECT_GT(result_d_small, 0.0);
    EXPECT_LT(result_d_small, 0.01); // Should be very small

    // Test very small d (should not cause division by zero issues)
    double small_d = 1e-10;
    double result_small_d = BlattWeisskopf(1.0, 1, small_d);
    EXPECT_GT(result_small_d, 0.0);
    EXPECT_LT(result_small_d, 1.0);
}

// Test specific physics values (typical hadron physics parameters)
TEST_F(BlattWeisskopfTest, TypicalPhysicsValues)
{
    // Typical values in hadron physics
    double d_typical = 1.5; // GeV^-1, typical hadron size

    // Test with typical momentum values in GeV
    std::vector<double> typical_momenta = {0.1, 0.3, 0.5, 1.0, 2.0};

    for (double q : typical_momenta)
    {
        // S-wave should always be 1
        double s_wave = BlattWeisskopf(q, 0, d_typical);
        expectNear(1.0, s_wave, tolerance);

        // P-wave should be between 0 and 1
        double p_wave = BlattWeisskopf(q, 1, d_typical);
        EXPECT_GE(p_wave, 0.0);
        EXPECT_LE(p_wave, 1.0);

        // D-wave should be between 0 and 1
        double d_wave = BlattWeisskopf(q, 2, d_typical);
        EXPECT_GE(d_wave, 0.0);
        EXPECT_LE(d_wave, 1.0);

        // Higher L should suppress more at low momentum
        if (q < 1.0)
        {
            EXPECT_GE(p_wave, d_wave) << "P-wave should be larger than D-wave at low momentum";
        }
    }
}

// Test consistency with analytical formulas
TEST_F(BlattWeisskopfTest, AnalyticalConsistency)
{
    double q = 0.7;
    double d = 1.2;
    double z2 = (q * d) * (q * d);

    // Test P-wave formula explicitly
    double expected_p = std::sqrt(z2 / (1.0 + z2));
    double result_p = BlattWeisskopf(q, 1, d);
    expectNear(expected_p, result_p, tolerance);

    // Test D-wave formula explicitly
    double expected_d = std::sqrt(z2 * z2 / (9.0 + 3.0 * z2 + z2 * z2));
    double result_d = BlattWeisskopf(q, 2, d);
    expectNear(expected_d, result_d, tolerance);
}

// Test that the function works with the typical usage pattern from your JSON model
TEST_F(BlattWeisskopfTest, TypicalUsagePattern)
{
    // From your JSON: "BlattWeisskopf_resonance_l1" with radius 1.5, l=1
    double radius = 1.5;
    int L = 1;

    // Typical momentum values that might occur in three-body decays
    std::vector<double> test_momenta = {0.0, 0.2, 0.5, 0.8, 1.0, 1.5, 2.0};

    for (double q : test_momenta)
    {
        double bf = BlattWeisskopf(q, L, radius);

        // Basic sanity checks
        EXPECT_GE(bf, 0.0);
        EXPECT_LE(bf, 1.0);

        // At q=0, should be 0 for L>0
        if (q == 0.0 && L > 0)
        {
            expectNear(0.0, bf, tolerance);
        }

        // Should approach 1 for large q
        if (q > 5.0)
        {
            EXPECT_GT(bf, 0.9);
        }
    }

    // Test L=2 case (from "BlattWeisskopf_resonance_l2")
    L = 2;
    for (double q : test_momenta)
    {
        double bf = BlattWeisskopf(q, L, radius);
        EXPECT_GE(bf, 0.0);
        EXPECT_LE(bf, 1.0);
    }
}

// Performance test (optional, for large-scale usage)
TEST_F(BlattWeisskopfTest, PerformanceTest)
{
    double d = 1.5;
    int num_iterations = 100000;

    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < num_iterations; ++i)
    {
        double q = 0.001 * i; // varying momentum
        BlattWeisskopf(q, 1, d);
        BlattWeisskopf(q, 2, d);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Performance test: " << num_iterations << " iterations took "
              << duration.count() << " microseconds" << std::endl;

    // Should complete in reasonable time (less than 1 second)
    EXPECT_LT(duration.count(), 1000000);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
