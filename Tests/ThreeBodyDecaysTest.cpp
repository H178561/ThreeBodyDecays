#include <gtest/gtest.h>
#include "include/ThreeBodyDecays.hh"
#include <iostream>
#include <iomanip>

// Test fixture for ThreeBodyDecays tests
class ThreeBodyDecaysTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Common setup for all tests
        p_p = {13.6499, 11.4971, 26.0864, 31.621};
        p_K = {36.2116, 30.5006, 68.9512, 83.643};
        p_pi = {19.252, 16.2157, 37.0456, 44.7882};
        p_L = {69.1135, 58.2135, 132.083, 160.052};

        // Setup masses
        ms = {1.0, 2, 3, 15.0};
        mssquared = {1.0, 4, 9, 15.0 * 15.0};

        // Initialize spins
        spins = {1, 0, 0, 1};
    }

    // Common variables
    std::vector<double> p_p, p_K, p_pi, p_L;
    ThreeBodyMasses ms, mssquared;
    ThreeBodySpins spins;
    ThreeBodyDecays decays;
};

/*
// Test for amplitude calculations
TEST_F(ThreeBodyDecaysTest, AmplitudeCalculation) {
    std::vector<double> p_Lc(p_p.size());
    std::vector<complex> amplitudes;

    // Sum of vectors
    for (size_t i = 0; i < p_p.size(); ++i) {
        p_Lc[i] = p_p[i] + p_K[i] + p_pi[i];
    }

    // Possible helicities (-1/2, 1/2)
    std::array<double, 2> helicities = { -0.5, 0.5 };

    for (double lambda_p : helicities) {
        for (double lambda_Lc : helicities) {
            complex amplitude = 0.3 * lambda_Lc * min_scalar_product(p_p, p_pi) +
                                complex(0.0, 0.7) * lambda_p * min_scalar_product(p_K, p_pi) +
                                complex(1.1, -1.2) * (lambda_p + lambda_Lc) * min_scalar_product(p_K, p_p);
            amplitudes.push_back(amplitude);
        }
    }

    // Expected values for comparison
    std::vector<complex> expected_amplitudes = {
        complex(-1.44803, 1.42474),
        complex(0.0959718, -0.0502259),
        complex(-0.0959718, 0.0502259),
        complex(1.44803, -1.42474)
    };

    // Verify each amplitude
    ASSERT_EQ(amplitudes.size(), 4);
    for (size_t i = 0; i < amplitudes.size(); ++i) {
        EXPECT_NEAR(amplitudes[i].real(), expected_amplitudes[i].real(), 1e-3);
        EXPECT_NEAR(amplitudes[i].imag(), expected_amplitudes[i].imag(), 1e-3);
    }
}*/

// Test for Mandelstam variables and cosine angles
TEST_F(ThreeBodyDecaysTest, MandelstamVariablesAndCosines)
{
    auto σs = decays.x2σs({0.3, 0.3}, ms, 1);

    // Test Mandelstam variables
    double cos23 = decays.cosθ23(σs, mssquared);
    double cos31 = decays.cosθ31(σs, mssquared);
    double cos12 = decays.cosθ12(σs, mssquared);

    // Expected values
    EXPECT_NEAR(cos23, -0.40, 1e-5);
    EXPECT_NEAR(cos31, 0.115476, 1e-5);
    EXPECT_NEAR(cos12, 0.342264, 1e-5);
}

// Test for Wigner rotations
TEST_F(ThreeBodyDecaysTest, WignerRotations)
{
    auto σs = decays.x2σs({0.3, 0.3}, ms, 1);

    // Test cosζ values for different rotations
    double cosζ110 = wr(1, 1, 0)->cos_zeta(σs, mssquared);
    double cosζ210 = wr(2, 1, 0)->cos_zeta(σs, mssquared);
    double cosζ310 = wr(3, 1, 0)->cos_zeta(σs, mssquared);

    // Test more cosζ values
    double cosζ120 = wr(1, 2, 0)->cos_zeta(σs, mssquared);
    double cosζ220 = wr(2, 2, 0)->cos_zeta(σs, mssquared);
    double cosζ320 = wr(3, 2, 0)->cos_zeta(σs, mssquared);

    // Expected values for first set
    EXPECT_NEAR(cosζ110, 1.0, 1e-5);
    EXPECT_NEAR(cosζ210, -0.19649, 1e-5);
    EXPECT_NEAR(cosζ310, -0.792381, 1e-5);

    // Expected values for second set
    EXPECT_NEAR(cosζ120, -0.19649, 1e-5);
    EXPECT_NEAR(cosζ220, 1.0, 1e-5);
    EXPECT_NEAR(cosζ320, -0.442439, 1e-5);

    // Test a few more rotations from different systems
    double cosζ130 = wr(1, 3, 0)->cos_zeta(σs, mssquared);
    double cosζ230 = wr(2, 3, 0)->cos_zeta(σs, mssquared);
    double cosζ330 = wr(3, 3, 0)->cos_zeta(σs, mssquared);

    EXPECT_NEAR(cosζ130, -0.792381, 1e-5);
    EXPECT_NEAR(cosζ230, -0.442439, 1e-5);
    EXPECT_NEAR(cosζ330, 1.0, 1e-5);
}

// Test for decay chain creation with NoRecoupling
TEST_F(ThreeBodyDecaysTest, DecayChainNoRecoupling)
{
    auto σs = decays.x2σs({0.3, 0.3}, ms, 1);

    // Create lineshape function
    auto Xlineshape = [](double sigma) -> complex
    {
        return complex(2.2 / sigma, 0.0);
    };

    // Create decay chain with NoRecoupling
    auto dc = createDecayChainCoupling(
        1,                                     // k-value
        Xlineshape,                            // Lineshape function
        "2+",                                  // jp (Spin-parity)
        //ThreeBodyParities{'+', '+', '+', '+'}, // Parities
        ThreeBodySystem(ms, spins),            // ThreeBodySystem
        RecouplingType::NoRecoupling, {0, 0}   // NoRecoupling function
    );

    // Reference values for decay chain
    std::vector<int> refζs = {1, 2, 3, 1};

    // Calculate amplitude tensor
    Tensor4D result4d = decays.amplitude4d(*dc, σs, refζs);

    // Verify tensor dimensions
    EXPECT_GT(result4d.size(), 0);
    if (result4d.size() > 0)
    {
        EXPECT_GT(result4d[0].size(), 0);
        if (result4d[0].size() > 0)
        {
            EXPECT_GT(result4d[0][0].size(), 0);
            if (result4d[0][0].size() > 0)
            {
                EXPECT_GT(result4d[0][0][0].size(), 0);
            }
        }
    }

    std::cout << "Tensor4D result4d : " << std::endl;
    for (int i = 0; i < result4d.size(); ++i)
    {
        for (int j = 0; j < result4d[0].size(); ++j)
        {
            for (int k = 0; k < result4d[0][0].size(); ++k)
            {
                for (int z = 0; z < result4d[0][0][0].size(); ++z)
                {
                    std::cout << result4d[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }
}

// Test for decay chain with ParityRecoupling
TEST_F(ThreeBodyDecaysTest, DecayChainParityRecoupling)
{
    auto σs = decays.x2σs({0.3, 0.3}, ms, 1);

    // Create lineshape function
    auto Xlineshape = [](double sigma) -> complex
    {
        return complex(2.2 / sigma, 0.0);
    };

    // Create decay chain with ParityRecoupling for HRk
    auto dc2 = createDecayChainCoupling(
        1,                                             // k-value
        Xlineshape,                                    // Lineshape function
        "2+",                                          // jp (Spin-parity)
        //ThreeBodyParities{'+', '+', '+', '+'},         // Parities
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {0, 0}, false,   // NoRecoupling for Hij
        RecouplingType::ParityRecoupling, {2, 0}, true // ParityRecoupling for HRk
    );

    // Reference values
    std::vector<int> refζs = {1, 2, 3, 1};

    // Calculate amplitude tensor
    Tensor4D resultdc2 = decays.amplitude4d(*dc2, σs, refζs);

    // Verify tensor dimensions
    EXPECT_GT(resultdc2.size(), 0);
    if (resultdc2.size() > 0)
    {
        EXPECT_GT(resultdc2[0].size(), 0);
    }
}

// Test for decay chain with double ParityRecoupling
TEST_F(ThreeBodyDecaysTest, DecayChainDoubleParityRecoupling)
{
    auto σs = decays.x2σs({0.3, 0.3}, ms, 1);

    // Create lineshape function
    auto Xlineshape = [](double sigma) -> complex
    {
        return complex(2.2 / sigma, 0.0);
    };

    // Create decay chain with ParityRecoupling for both HRk and Hij
    auto dc3 = createDecayChainCoupling(
        3,                                              // k-value
        Xlineshape,                                     // Lineshape function
        "2+",                                           // jp (Spin-parity)
        //ThreeBodyParities{'+', '+', '+', '+'},          // Parities
        ThreeBodySystem(ms, spins),                     // ThreeBodySystem
        RecouplingType::ParityRecoupling, {2, 0}, true, // ParityRecoupling for HRk
        RecouplingType::ParityRecoupling, {2, 0}, false // ParityRecoupling for Hij
    );

    // Reference values
    std::vector<int> refζs = {3, 3, 3, 3};

    // Calculate amplitude tensor
    Tensor4D resultdc3 = decays.amplitude4d(*dc3, σs, refζs);

    // Verify tensor dimensions
    EXPECT_GT(resultdc3.size(), 0);
    if (resultdc3.size() > 0)
    {
        EXPECT_GT(resultdc3[0].size(), 0);
    }
}

// Main function to run the tests
#ifndef RUNNING_COMBINED_TESTS
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
#endif

void ThreeBodyDecays::A_test()
{
    std::vector<double> p_p = {13.6499, 11.4971, 26.0864, 31.621};
    std::vector<double> p_K = {36.2116, 30.5006, 68.9512, 83.643};
    std::vector<double> p_pi = {19.252, 16.2157, 37.0456, 44.7882};
    std::vector<double> p_L{69.1135, 58.2135, 132.083, 160.052};

    std::vector<double> p_Lc(p_p.size());

    std::vector<complex> amplitudes;

    // Sum of vectors
    for (size_t i = 0; i < p_p.size(); ++i)
    {
        p_Lc[i] = p_p[i] + p_K[i] + p_pi[i];
    }
    // std::cout << "Momentum: " << p_Lc[0] << " " << p_Lc[1] << " " << p_Lc[2] << " " << p_Lc[3] << std::endl;

    // Possible helicities
    // std::vector<std::pair<double, double>> possible_helicities = { {-0.5, -0.5}, {-0.5, 0.5}, {0.5, -0.5}, {0.5, 0.5} };

    // Possible helicities (-1/2, 1/2)
    std::array<double, 2> helicities = {-0.5, 0.5};

    for (double lambda_p : helicities)
    {
        for (double lambda_Lc : helicities)
        {
            complex amplitude = 0.3 * lambda_Lc *
                                    min_scalar_product(p_p, p_pi) +
                                complex(0.0, 0.7) * lambda_p *
                                    min_scalar_product(p_K, p_pi) +
                                complex(1.1, -1.2) *
                                    (lambda_p + lambda_Lc) *
                                    min_scalar_product(p_K, p_p);
            amplitudes.push_back(amplitude);
        }
    }

    // return amplitudes;

    // Return or use amplitudes as needed

    std::cout
        << "Amplitudes: " << std::setprecision(8) << amplitudes[0] << " "
        << amplitudes[1] << " " << amplitudes[2] << " " << amplitudes[3]
        << std::endl;

    std::cout
        << "Should match to "
        << "-1.44803+1.42474im, 0.0959718-0.0502259im, -0.0959718+0.0502259im, 1.44803-1.42474im "
        << std::endl;

    // namespace {using output_t = std::tuple<LHCb::Particles, LHCb::Particles, LHCb::ProtoParticles, LHCb::Tracks>;}

    ThreeBodyMasses ms = {1.0, 2, 3, 15.0};
    ThreeBodyMasses mssquared = {1.0, 4, 9, 15.0 * 15.0};
    auto σs = x2σs({0.3, 0.3}, ms, 1); // just a mapping function

    std::cout
        << "Mandelstam variables: " << σs[0] << " " << σs[1] << " " << σs[2]
        << std::endl;

    double cos23 = cosθ23(σs, mssquared);
    double cos31 = cosθ31(σs, mssquared);
    double cos12 = cosθ12(σs, mssquared);

    std::cout
        << "Cosine values: " << "cos23: " << cos23 << " cos31: " << cos31
        << " cos12: " << cos12 << std::endl;
    std::cout
        << "Should match to "
        << "cosθ_23: -0.42 │ cosθ_31: 0.115476 │ cosθ_12: 0.342264" << std::endl;

    double cosζ110 = wr(1, 1, 0)->cos_zeta(σs, mssquared);
    double cosζ210 = wr(2, 1, 0)->cos_zeta(σs, mssquared);
    double cosζ310 = wr(3, 1, 0)->cos_zeta(σs, mssquared);
    double cosζ120 = wr(1, 2, 0)->cos_zeta(σs, mssquared);
    double cosζ220 = wr(2, 2, 0)->cos_zeta(σs, mssquared);
    double cosζ320 = wr(3, 2, 0)->cos_zeta(σs, mssquared);
    double cosζ130 = wr(1, 3, 0)->cos_zeta(σs, mssquared);
    double cosζ230 = wr(2, 3, 0)->cos_zeta(σs, mssquared);
    double cosζ330 = wr(3, 3, 0)->cos_zeta(σs, mssquared);
    double cosζ111 = wr(1, 1, 1)->cos_zeta(σs, mssquared);
    double cosζ211 = wr(2, 1, 1)->cos_zeta(σs, mssquared);
    double cosζ311 = wr(3, 1, 1)->cos_zeta(σs, mssquared);
    double cosζ121 = wr(1, 2, 1)->cos_zeta(σs, mssquared);
    double cosζ221 = wr(2, 2, 1)->cos_zeta(σs, mssquared);
    double cosζ321 = wr(3, 2, 1)->cos_zeta(σs, mssquared);
    double cosζ131 = wr(1, 3, 1)->cos_zeta(σs, mssquared);
    double cosζ231 = wr(2, 3, 1)->cos_zeta(σs, mssquared);
    double cosζ331 = wr(3, 3, 1)->cos_zeta(σs, mssquared);
    double cosζ112 = wr(1, 1, 2)->cos_zeta(σs, mssquared);
    double cosζ212 = wr(2, 1, 2)->cos_zeta(σs, mssquared);
    double cosζ312 = wr(3, 1, 2)->cos_zeta(σs, mssquared);
    double cosζ122 = wr(1, 2, 2)->cos_zeta(σs, mssquared);
    double cosζ222 = wr(2, 2, 2)->cos_zeta(σs, mssquared);
    double cosζ322 = wr(3, 2, 2)->cos_zeta(σs, mssquared);
    double cosζ132 = wr(1, 3, 2)->cos_zeta(σs, mssquared);
    double cosζ232 = wr(2, 3, 2)->cos_zeta(σs, mssquared);
    double cosζ332 = wr(3, 3, 2)->cos_zeta(σs, mssquared);
    double cosζ113 = wr(1, 1, 3)->cos_zeta(σs, mssquared);
    double cosζ213 = wr(2, 1, 3)->cos_zeta(σs, mssquared);
    double cosζ313 = wr(3, 1, 3)->cos_zeta(σs, mssquared);
    double cosζ123 = wr(1, 2, 3)->cos_zeta(σs, mssquared);
    double cosζ223 = wr(2, 2, 3)->cos_zeta(σs, mssquared);
    double cosζ323 = wr(3, 2, 3)->cos_zeta(σs, mssquared);
    double cosζ133 = wr(1, 3, 3)->cos_zeta(σs, mssquared);
    double cosζ233 = wr(2, 3, 3)->cos_zeta(σs, mssquared);
    double cosζ333 = wr(3, 3, 3)->cos_zeta(σs, mssquared);

    /*
    1 │ cosζ_1(1)_for_0   1.0
    2 │ cosζ_2(1)_for_0  -0.19649
    3 │ cosζ_3(1)_for_0  -0.792381
    4 │ cosζ_1(2)_for_0  -0.19649
    5 │ cosζ_2(2)_for_0   1.0
    6 │ cosζ_3(2)_for_0  -0.442439
    7 │ cosζ_1(3)_for_0  -0.792381
    8 │ cosζ_2(3)_for_0  -0.442439
    9 │ cosζ_3(3)_for_0   1.0
   10 │ cosζ_1(1)_for_1   1.0
   11 │ cosζ_2(1)_for_1   0.99797
   12 │ cosζ_3(1)_for_1   0.989415
   13 │ cosζ_1(2)_for_1   0.99797
   14 │ cosζ_2(2)_for_1   1.0
   15 │ cosζ_3(2)_for_1   0.978166
   16 │ cosζ_1(3)_for_1   0.989415
   17 │ cosζ_2(3)_for_1   0.978166
   18 │ cosζ_3(3)_for_1   1.0
   19 │ cosζ_1(1)_for_2   1.0
   20 │ cosζ_2(1)_for_2   0.951225
   21 │ cosζ_3(1)_for_2   0.728673
   22 │ cosζ_1(2)_for_2   0.951225
   23 │ cosζ_2(2)_for_2   1.0
   24 │ cosζ_3(2)_for_2   0.904411
   25 │ cosζ_1(3)_for_2   0.728673
   26 │ cosζ_2(3)_for_2   0.904411
   27 │ cosζ_3(3)_for_2   1.0
   28 │ cosζ_1(1)_for_3   1.0
   29 │ cosζ_2(1)_for_3   0.892623
   30 │ cosζ_3(1)_for_3   0.95766
   31 │ cosζ_1(2)_for_3   0.892623
   32 │ cosζ_2(2)_for_3   1.0
   33 │ cosζ_3(2)_for_3   0.984616
   34 │ cosζ_1(3)_for_3   0.95766
   35 │ cosζ_2(3)_for_3   0.984616
   36 │ cosζ_3(3)_for_3   1.0
   */

    std::cout
        << "cosζ_1(1)_for_0: " << cosζ110 << "                    " << " real: 1.0 "
        << approx_equal(cosζ110, 1.0) << "\n "
        << "cosζ_2(1)_for_0: " << cosζ210 << "                    "
        << " real: -0.19649 " << approx_equal(cosζ210, -0.19649)
        << "cosζ_3(1)_for_0: " << cosζ310 << "                    "
        << " real: -0.792381 " << approx_equal(cosζ310, -0.792381) << "\n "
        << "cosζ_1(2)_for_0: " << cosζ120 << "                    "
        << " real: -0.19649 " << approx_equal(cosζ120, -0.19649) << "\n "
        << "cosζ_2(2)_for_0: " << cosζ220 << "                    "
        << " real: 1.0 " << approx_equal(cosζ220, 1.0) << "\n "
        << "cosζ_3(2)_for_0: " << cosζ320 << "                    "
        << " real: -0.442439 " << approx_equal(cosζ320, -0.442439) << "\n "
        << "cosζ_1(3)_for_0: " << cosζ130 << "                    "
        << " real: -0.792381 " << approx_equal(cosζ130, -0.792381) << "\n "
        << "cosζ_2(3)_for_0: " << cosζ230 << "                    "
        << " real: -0.442439 " << approx_equal(cosζ230, -0.442439) << "\n "
        << "cosζ_3(3)_for_0: " << cosζ330 << "                    " << " real: 1.0 "
        << approx_equal(cosζ330, 1.0) << "\n "
        << "cosζ_1(1)_for_1: " << cosζ111 << "                    " << " real: 1.0 "
        << approx_equal(cosζ111, 1.0) << "\n "
        << "cosζ_2(1)_for_1: " << cosζ211 << "                    "
        << " real: 0.99797 " << approx_equal(cosζ211, 0.99797) << "\n "
        << "cosζ_3(1)_for_1: " << cosζ311 << "                    "
        << " real: 0.989415 " << approx_equal(cosζ311, 0.989415) << "\n "
        << "cosζ_1(2)_for_1: " << cosζ121 << "                    "
        << " real: 0.99797 " << approx_equal(cosζ121, 0.99797) << "\n "
        << "cosζ_2(2)_for_1: " << cosζ221 << "                    " << " real: 1.0 "
        << approx_equal(cosζ221, 1.0) << "\n "
        << "cosζ_3(2)_for_1: " << cosζ321 << "                    "
        << " real: 0.978166 " << approx_equal(cosζ321, 0.978166) << "\n "
        << "cosζ_1(3)_for_1: " << cosζ131 << "                    "
        << " real: 0.989415 " << approx_equal(cosζ131, 0.989415) << "\n "
        << "cosζ_2(3)_for_1: " << cosζ231 << "                    "
        << " real: 0.978166 " << approx_equal(cosζ231, 0.978166) << "\n "
        << "cosζ_3(3)_for_1: " << cosζ331 << "                    " << " real: 1.0 "
        << approx_equal(cosζ331, 1.0) << "\n "
        << "cosζ_1(1)_for_2: " << cosζ112 << "                    " << " real: 1.0 "
        << approx_equal(cosζ112, 1.0) << "\n "
        << "cosζ_2(1)_for_2: " << cosζ212 << "                    "
        << " real: 0.951225 " << approx_equal(cosζ212, 0.951225) << "\n "
        << "cosζ_3(1)_for_2: " << cosζ312 << "                    "
        << " real: 0.728673 " << approx_equal(cosζ312, 0.728673) << "\n "
        << "cosζ_1(2)_for_2: " << cosζ122 << "                    "
        << " real: 0.951225 " << approx_equal(cosζ122, 0.951225) << "\n "
        << "cosζ_2(2)_for_2: " << cosζ222 << "                    " << " real: 1.0 "
        << approx_equal(cosζ222, 1.0) << "\n "
        << "cosζ_3(2)_for_2: " << cosζ322 << "                    "
        << " real: 0.904411 " << approx_equal(cosζ322, 0.904411) << "\n "
        << "cosζ_1(3)_for_2: " << cosζ132 << "                    "
        << " real: 0.728673 " << approx_equal(cosζ132, 0.728673) << "\n "
        << "cosζ_2(3)_for_2: " << cosζ232 << "                    "
        << " real: 0.904411 " << approx_equal(cosζ232, 0.904411) << "\n "
        << "cosζ_3(3)_for_2: " << cosζ332 << "                    " << " real: 1.0 "
        << approx_equal(cosζ332, 1.0) << "\n "
        << "cosζ_1(1)_for_3: " << cosζ113 << "                    " << " real: 1.0 "
        << approx_equal(cosζ113, 1.0) << "\n "
        << "cosζ_2(1)_for_3: " << cosζ213 << "                    "
        << " real: 0.892623 " << approx_equal(cosζ213, 0.892623) << "\n "
        << "cosζ_3(1)_for_3: " << cosζ313 << "                    "
        << " real: 0.95766 " << approx_equal(cosζ313, 0.95766) << "\n "
        << "cosζ_1(2)_for_3: " << cosζ123 << "                    "
        << " real: 0.892623 " << approx_equal(cosζ123, 0.892623) << "\n "
        << "cosζ_2(2)_for_3: " << cosζ223 << "                    " << " real: 1.0 "
        << approx_equal(cosζ223, 1.0) << "\n "
        << "cosζ_3(2)_for_3: " << cosζ323 << "                    "
        << " real: 0.984616 " << approx_equal(cosζ323, 0.984616) << "\n "
        << "cosζ_1(3)_for_3: " << cosζ133 << "                    "
        << " real: 0.95766 " << approx_equal(cosζ133, 0.95766) << "\n "
        << "cosζ_2(3)_for_3: " << cosζ233 << "                    "
        << " real: 0.984616 " << approx_equal(cosζ233, 0.984616) << "\n "
        << "cosζ_3(3)_for_3: " << cosζ333 << "                    " << " real: 1.0 "
        << approx_equal(cosζ333, 1.0) << std::endl;

    ThreeBodySpins spins = {1, 0, 0, 1}; // h0=1 bezieht sich auf den Spin des Elternteilchens

    auto tbs = ThreeBodySystem(ms, spins);

    auto Xlineshape = [](double sigma) -> complex
    {
        return complex(2.2 / sigma, 0.0);
    };

    // NoRecoupling-Funktionen erstellen
    HelicityFunction noRecoupling = [](const std::array<int, 2> &twoMs, const std::array<int, 3> &twoJs) -> complex
    {
        return complex(1.0, 0.0); // Einfache Implementierung von NoRecoupling
    };

    // DecayChain erstellen
    // DecayChain erstellen

    auto dc = createDecayChainLS(
        1,                                     // k-Wert
        Xlineshape,                            // Lineshape-Funktion
        "2+",                                  // jp (Spin-Parität, hier 2+)
        ThreeBodyParities{'+', '+', '+', '+'}, // Paritäten, ändern Sie diese entsprechend
        tbs                                    // ThreeBodySystem as reference
    );

    auto dchain = DecayChain(
        1, // k-Wert
        // ms,
        // spins,
        2,                                            // two_j (for "2+")
        Xlineshape,                                   // Lineshape-Funktion
        noRecoupling,                                 // NoRecoupling-Funktion
        noRecoupling,                                 // NoRecoupling-Funktion
        *std::make_shared<ThreeBodySystem>(ms, spins) // ThreeBodySystem by value, not pointer
    );

    // MandelstamTuple erstellen (σs)
    // Schon erstellt
    // auto σs = x2σs({ 0.3, 0.3 }, ms, 1);

    // Berechnen Sie die Amplitude
    // Gehen wir davon aus, dass Sie für beide Teilchen Helizitätswerte von 1/2 verwenden möchten
    std::vector<int> two_λs = {1, 1, 1, 1}; // Für 1/2 Helizität verwenden wir 1 in doubled representation
    std::vector<int> refζs = {1, 2, 3, 1};

    // complex result = amplitude0(*dc, σs, two_λs);

    // Ausgabe
    // std::cout << "Amplitude result: " << result << std::endl;

    /*
    Tensor4Dcomp result4dold = amplitude4dold(*dc, σs, refζs);

    // Print a summary of the Tensor4D instead of trying to stream it directly
    std::cout << "Tensor4D dimensions: "
        << result4dold.size() << " x "
        << (result4dold.empty() ? 0 : result4dold[0].size()) << " x "
        << (result4dold.empty() || result4dold[0].empty() ? 0 : result4dold[0][0].size()) << " x "
        << (result4dold.empty() || result4dold[0].empty() || result4dold[0][0].empty() ? 0 : result4dold[0][0][0].size())
        << std::endl;

    std::cout << " Tensor4D result4dold 0000: " << result4dold[0][0][0][0] << std::endl;
    std::cout << "Tensor4D result4dold 0001: " << result4dold[0][0][0][1] << std::endl;
    std::cout << "Tensor4D result4dold 0002: " << result4dold[1][0][0][2] << std::endl;
    std::cout << "Tensor4D result4dold 1000: " << result4dold[1][0][0][0] << std::endl;
    std::cout << "Tensor4D result4dold 1001: " << result4dold[1][0][0][1] << std::endl;
    std::cout << "Tensor4D result4dold 1002: " << result4dold[1][0][0][2] << std::endl;
    std::cout << "Tensor4D result4dold 2000: " << result4dold[2][0][0][0] << std::endl;
    std::cout << "Tensor4D result4dold 2001: " << result4dold[2][0][0][1] << std::endl;
    std::cout << "Tensor4D result4dold 2002: " << result4dold[2][0][0][2] << std::endl;
    */
    // Tensor4D result4d = amplitude4d(dchain, σs, refζs);
    Tensor4D result4d = amplitude4d(*dc, σs, refζs);

    // Print a summary of the Tensor4D instead of trying to stream it directly
    std::cout << "Tensor4D dimensions: "
              << result4d.size() << " x "
              << (result4d.empty() ? 0 : result4d[0].size()) << " x "
              << (result4d.empty() || result4d[0].empty() ? 0 : result4d[0][0].size()) << " x "
              << (result4d.empty() || result4d[0].empty() || result4d[0][0].empty() ? 0 : result4d[0][0][0].size())
              << std::endl;

    std::cout << "Tensor4D result4ddc : " << std::endl;
    for (int i = 0; i < result4d.size(); ++i)
    {
        for (int j = 0; j < result4d[0].size(); ++j)
        {
            for (int k = 0; k < result4d[0][0].size(); ++k)
            {
                for (int z = 0; z < result4d[0][0][0].size(); ++z)
                {
                    std::cout << result4d[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }

    // dc = DecayChain(; k=1, two_j=2, Xlineshape=σ->2.2/σ,
    //     HRk=ParityRecoupling(2,0,true), Hij=NoRecoupling(0,0),
    //     tbs=ThreeBodySystem(ms, ThreeBodySpins(1,0,0; h0=1)))

    auto dc2 = createDecayChainLS(
        1,                                     // k-Wert
        Xlineshape,                            // Lineshape-Funktion
        "2+",                                  // jp (Spin-Parität, hier 2+)
        ThreeBodyParities{'+', '+', '+', '+'}, // Paritäten, ändern Sie diese entsprechend
        tbs                                    // ThreeBodySystem as reference
    );
    Tensor4D resultdc2 = amplitude4d(*dc2, σs, refζs);

    std::cout << "Tensor4D resultdc2 : " << std::endl;
    for (int i = 0; i < resultdc2.size(); ++i)
    {
        for (int j = 0; j < resultdc2[0].size(); ++j)
        {
            for (int k = 0; k < resultdc2[0][0].size(); ++k)
            {
                for (int z = 0; z < resultdc2[0][0][0].size(); ++z)
                {
                    std::cout << resultdc2[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }

    // dc = DecayChain(; k=3, two_j=2, Xlineshape=σ->2.2/σ,
    // HRk=ParityRecoupling(2,0,true), Hij=ParityRecoupling(2,0, false),
    // tbs=ThreeBodySystem(ms, ThreeBodySpins(1,0,0; h0=1)))

    auto dc3 = createDecayChainLS(
        3,                                     // k-Wert
        Xlineshape,                            // Lineshape-Funktion
        "2+",                                  // jp (Spin-Parität, hier 2+)
        ThreeBodyParities{'+', '+', '+', '+'}, // Paritäten, ändern Sie diese entsprechend
        tbs                                    // ThreeBodySystem as reference
    );
    Tensor4D resultdc3 = amplitude4d(*dc3, σs, {3, 3, 3, 3});

    std::cout << "Tensor4D resultdc3 : " << std::endl;
    for (int i = 0; i < resultdc3.size(); ++i)
    {
        for (int j = 0; j < resultdc3[0].size(); ++j)
        {
            for (int k = 0; k < resultdc3[0][0].size(); ++k)
            {
                for (int z = 0; z < resultdc3[0][0][0].size(); ++z)
                {
                    std::cout << resultdc3[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }
}
