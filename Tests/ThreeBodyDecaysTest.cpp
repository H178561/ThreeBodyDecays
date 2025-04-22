#include <gtest/gtest.h>
#include "ThreeBodyDecays.hh"
#include <iostream>
#include <iomanip>

// Test fixture for ThreeBodyDecays tests
class ThreeBodyDecaysTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup for all tests
        p_p = { 13.6499, 11.4971, 26.0864, 31.621 };
        p_K = { 36.2116, 30.5006, 68.9512, 83.643 };
        p_pi = { 19.252, 16.2157, 37.0456, 44.7882 };
        p_L = { 69.1135, 58.2135, 132.083, 160.052 };

        // Setup masses
        ms = { 1.0, 2, 3, 15.0 };
        mssquared = { 1.0, 4, 9, 15.0 * 15.0 };

        // Initialize spins
        spins = { 1, 0, 0, 1 };
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
TEST_F(ThreeBodyDecaysTest, MandelstamVariablesAndCosines) {
    auto σs = decays.x2σs({ 0.3, 0.3 }, ms, 1);

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
TEST_F(ThreeBodyDecaysTest, WignerRotations) {
    auto σs = decays.x2σs({ 0.3, 0.3 }, ms, 1);

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
TEST_F(ThreeBodyDecaysTest, DecayChainNoRecoupling) {
    auto σs = decays.x2σs({ 0.3, 0.3 }, ms, 1);

    // Create lineshape function
    auto Xlineshape = [](double sigma) -> complex {
        return complex(2.2 / sigma, 0.0);
    };

    // Create decay chain with NoRecoupling
    auto dc = createDecayChainLS(
        1,                // k-value
        Xlineshape,       // Lineshape function
        "2+",             // jp (Spin-parity)
        ThreeBodyParities{'+', '+', '+', '+'}, // Parities
        std::make_shared<ThreeBodySystem>(ms, spins),  // ThreeBodySystem
        RecouplingType::NoRecoupling, {0,0} // NoRecoupling function
    );

    // Reference values for decay chain
    std::vector<int> refζs = {1, 2, 3, 1};

    // Calculate amplitude tensor
    Tensor4D result4d = decays.amplitude4d(*dc, σs, refζs);

    // Verify tensor dimensions
    EXPECT_GT(result4d.size(), 0);
    if (result4d.size() > 0) {
        EXPECT_GT(result4d[0].size(), 0);
        if (result4d[0].size() > 0) {
            EXPECT_GT(result4d[0][0].size(), 0);
            if (result4d[0][0].size() > 0) {
                EXPECT_GT(result4d[0][0][0].size(), 0);
            }
        }
    }
}

// Test for decay chain with ParityRecoupling
TEST_F(ThreeBodyDecaysTest, DecayChainParityRecoupling) {
    auto σs = decays.x2σs({ 0.3, 0.3 }, ms, 1);

    // Create lineshape function
    auto Xlineshape = [](double sigma) -> complex {
        return complex(2.2 / sigma, 0.0);
    };

    // Create decay chain with ParityRecoupling for HRk
    auto dc2 = createDecayChainLS(
        1,                // k-value
        Xlineshape,       // Lineshape function
        "2+",             // jp (Spin-parity)
        ThreeBodyParities{'+', '+', '+', '+'}, // Parities
        std::make_shared<ThreeBodySystem>(ms, spins),  // ThreeBodySystem
        RecouplingType::NoRecoupling, {0,0}, false, // NoRecoupling for Hij
        RecouplingType::ParityRecoupling, {2,0}, true // ParityRecoupling for HRk
    );

    // Reference values
    std::vector<int> refζs = {1, 2, 3, 1};

    // Calculate amplitude tensor
    Tensor4D resultdc2 = decays.amplitude4d(*dc2, σs, refζs);

    // Verify tensor dimensions
    EXPECT_GT(resultdc2.size(), 0);
    if (resultdc2.size() > 0) {
        EXPECT_GT(resultdc2[0].size(), 0);
    }
}

// Test for decay chain with double ParityRecoupling
TEST_F(ThreeBodyDecaysTest, DecayChainDoubleParityRecoupling) {
    auto σs = decays.x2σs({ 0.3, 0.3 }, ms, 1);

    // Create lineshape function
    auto Xlineshape = [](double sigma) -> complex {
        return complex(2.2 / sigma, 0.0);
    };

    // Create decay chain with ParityRecoupling for both HRk and Hij
    auto dc3 = createDecayChainLS(
        3,                // k-value
        Xlineshape,       // Lineshape function
        "2+",             // jp (Spin-parity)
        ThreeBodyParities{'+', '+', '+', '+'}, // Parities
        std::make_shared<ThreeBodySystem>(ms, spins),  // ThreeBodySystem
        RecouplingType::ParityRecoupling, {2,0}, true, // ParityRecoupling for HRk
        RecouplingType::ParityRecoupling, {2,0}, false // ParityRecoupling for Hij
    );

    // Reference values
    std::vector<int> refζs = {3, 3, 3, 3};

    // Calculate amplitude tensor
    Tensor4D resultdc3 = decays.amplitude4d(*dc3, σs, refζs);

    // Verify tensor dimensions
    EXPECT_GT(resultdc3.size(), 0);
    if (resultdc3.size() > 0) {
        EXPECT_GT(resultdc3[0].size(), 0);
    }
}

// Main function to run the tests
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
