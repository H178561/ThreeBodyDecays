#include <gtest/gtest.h>
#include "ThreeBodyDecays.hh"
#include "ThreeBodyAmplitudeModel.hh"
#include <iostream>
#include <iomanip>

// Test fixture for ThreeBodyDecays tests
class JsonModelTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Common setup for all tests


        // Setup masses
        ms = {0.938272046, 0.13957018, 0.493677, 2.28646};
        // mssquared = {1.0, 4, 9, 15.0 * 15.0};

        // Initialize spins
        spins = {1, 0, 0, 1};
    }

    // Create Breit-Wigner lineshape function
    std::function<complex(double)> createBW(double mass, double width)
    {
        return [mass, width](double sigma) -> complex
        {
            return complex(1.0, 0.0) / complex(mass * mass - sigma, -mass * width);
        };
    }

    // Common variables
    std::vector<double> p_p, p_K, p_pi, p_L;
    ThreeBodyMasses ms, mssquared;
    ThreeBodySpins spins;
    ThreeBodyDecays decays;
};

/*
L1520XLineshape: (0.061,0.001)
L1520Xlineshapeold: (0.434,0.004)
L1520FF 0.9210.152
L1520 k 2 Lineshape (0.061,0.001)3/2- helicity -10 parity 010
(3.605,0.544)
Single amp tensor:
(-0.084,0.003)	(-0.084,0.003)
(0.128,-0.004)	(0.128,-0.004)
L1520XLineshape: (0.061,0.001)
L1520Xlineshapeold: (0.434,0.004)
L1520FF 0.9210.152
L1520 k 2 Lineshape (0.061,0.001)3/2- helicity 10 parity 010
(-1.970,18.380)
Single amp tensor:
(-0.128,0.004)	(0.128,-0.004)
(-0.084,0.003)	(0.084,-0.003)
*/

TEST_F(JsonModelTest, L1520_Test)
{
    MandelstamTuple σs = {1, 3, 2};

    // Create Breit-Wigner lineshape function for L1520
    auto originalBreitWigner = createBW(1.518, 0.015);
    complex formFactor1 = 0.921;
    complex formFactor2 = 0.152;
    auto Xlineshape = [originalBreitWigner, formFactor1, formFactor2](double s) -> complex
    {
        return originalBreitWigner(s) * formFactor1 * formFactor2;
    };

    // std::cout << "L1520XLineshape: " << Xlineshape(0.0) << std::endl;
    //  Create decay chain with ParityRecoupling for HRk
    auto L15201 = createDecayChainCoupling(
        2,                                              // k-value
        Xlineshape,                                     // Lineshape function
        "3/2-",                                         // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                     // ThreeBodySystem
        RecouplingType::NoRecoupling, {-1, 0}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, false // ParityRecoupling for Hij
    );
    // std::cout << "L1520Xlineshapeold: " << Xlineshape(0.0) << std::endl;

    // Calculate amplitude tensor
    Tensor4Dcomp result1 = decays.amplitude4dcomp(*L15201, σs, 1);
    double prec = 1e-3;

    EXPECT_NEAR(result1[0][0][0][0].real(), -0.084, prec);
    EXPECT_NEAR(result1[0][0][0][0].imag(), 0.003, prec);
    EXPECT_NEAR(result1[0][0][0][1].real(), -0.084, prec);
    EXPECT_NEAR(result1[0][0][0][1].imag(), 0.003, prec);
    EXPECT_NEAR(result1[1][0][0][0].real(), 0.128, prec);
    EXPECT_NEAR(result1[1][0][0][0].imag(), -0.004, prec);
    EXPECT_NEAR(result1[1][0][0][1].real(), 0.128, prec);
    EXPECT_NEAR(result1[1][0][0][1].imag(), -0.004, prec);

    /*
    std::cout << "Single amp tensor: " << std::endl;
    for (int i = 0; i < result1.size(); ++i)
    {
        for (int j = 0; j < result1[0].size(); ++j)
        {
            for (int k = 0; k < result1[0][0].size(); ++k)
            {
                for (int z = 0; z < result1[0][0][0].size(); ++z)
                {
                    std::cout << std::fixed << std::setprecision(3) << result1[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }*/

    // Create decay chain with ParityRecoupling for HRk
    auto L15202 = createDecayChainCoupling(
        2,                                              // k-value
        Xlineshape,                                     // Lineshape function
        "3/2-",                                         // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                     // ThreeBodySystem
        RecouplingType::NoRecoupling, {1, 0}, false,    // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, false // ParityRecoupling for Hij
    );

    // Calculate amplitude tensor
    Tensor4Dcomp result2 = decays.amplitude4dcomp(*L15202, σs, 1);
    EXPECT_NEAR(result2[0][0][0][0].real(), -0.128, prec);
    EXPECT_NEAR(result2[0][0][0][0].imag(), 0.004, prec);
    EXPECT_NEAR(result2[0][0][0][1].real(), 0.128, prec);
    EXPECT_NEAR(result2[0][0][0][1].imag(), -0.004, prec);
    EXPECT_NEAR(result2[1][0][0][0].real(), -0.084, prec);
    EXPECT_NEAR(result2[1][0][0][0].imag(), 0.003, prec);
    EXPECT_NEAR(result2[1][0][0][1].real(), 0.084, prec);
    EXPECT_NEAR(result2[1][0][0][1].imag(), -0.003, prec);
}

TEST_F(JsonModelTest, L1600_Test)
{

    /*
    L1600XLineshape: (0.212,0.033)
    L1600Xlineshapeold: (0.368,0.056) 1.630 0.250
    L1600FF 1.0000.577
    L1600 k 2 Lineshape (0.212,0.033)1/2- helicity 10 parity 010
    (10.063,-1.216)
    10 -10
    10 10
    Single amp tensor:
    (0.595,-0.707)	(-0.593,0.704)
    (-0.366,0.435)	(0.365,-0.434)
    L1600XLineshape: (0.212,0.033)
    L1600Xlineshapeold: (0.368,0.056) 1.630 0.250
    L1600FF 1.0000.577
    L1600 k 2 Lineshape (0.212,0.033)1/2- helicity -10 parity 010
    (-6.987,-4.450)
    -10 -10
    -10 10
    Single amp tensor:
    (-0.365,0.434)	(-0.366,0.435)
    (-0.593,0.704)	(-0.595,0.707)
    */

    MandelstamTuple σs = {1, 3, 2};
    // Create Breit-Wigner lineshape function for L1600
    auto originalBreitWigner = createBW(1.630, 0.250);
    complex formFactor1 = 1.000;
    complex formFactor2 = 0.577;
    auto Xlineshape = [originalBreitWigner, formFactor1, formFactor2](double s) -> complex
    {
        return originalBreitWigner(s) * formFactor1 * formFactor2;
    };

    auto L16001 = createDecayChainCoupling(
        2,                                              // k-value
        Xlineshape,                                     // Lineshape function
        "1/2-",                                         // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                     // ThreeBodySystem
        RecouplingType::NoRecoupling, {1, 0}, false,    // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, false // ParityRecoupling for Hij
    );

    Tensor4Dcomp result1 = decays.amplitude4dcomp(*L16001, σs, 1);
    double prec = 1e-3;
    EXPECT_NEAR(result1[0][0][0][0].real(), 0.595, prec);
    EXPECT_NEAR(result1[0][0][0][0].imag(), -0.707, prec);
    EXPECT_NEAR(result1[0][0][0][1].real(), -0.593, prec);
    EXPECT_NEAR(result1[0][0][0][1].imag(), 0.704, prec);
    EXPECT_NEAR(result1[1][0][0][0].real(), -0.366, prec);
    EXPECT_NEAR(result1[1][0][0][0].imag(), 0.435, prec);
    EXPECT_NEAR(result1[1][0][0][1].real(), 0.365, prec);
    EXPECT_NEAR(result1[1][0][0][1].imag(), -0.434, prec);

    auto L16002 = createDecayChainCoupling(
        2,                                              // k-value
        Xlineshape,                                     // Lineshape function
        "1/2-",                                         // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                     // ThreeBodySystem
        RecouplingType::NoRecoupling, {-1, 0}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, false // ParityRecoupling for Hij
    );
    Tensor4Dcomp result2 = decays.amplitude4dcomp(*L16002, σs, 1);
    EXPECT_NEAR(result2[0][0][0][0].real(), -0.365, prec);
    EXPECT_NEAR(result2[0][0][0][0].imag(), 0.434, prec);
    EXPECT_NEAR(result2[0][0][0][1].real(), -0.366, prec);
    EXPECT_NEAR(result2[0][0][0][1].imag(), 0.435, prec);
    EXPECT_NEAR(result2[1][0][0][0].real(), -0.593, prec);
    EXPECT_NEAR(result2[1][0][0][0].imag(), 0.704, prec);
    EXPECT_NEAR(result2[1][0][0][1].real(), -0.595, prec);
    EXPECT_NEAR(result2[1][0][0][1].imag(), 0.707, prec);
}

TEST_F(JsonModelTest, L1670)
{
    /*
    """
L1670XLineshape: (0.358,0.006)
L1670Xlineshapeold: (0.358,0.006)
BW(1.670, 0.030)
L1670FF 1.0001.000
Creating DecayChain with k=2, two_j=1
L1670 k 2 Lineshape (0.358,0.006)1/2+ helicity -10 parity 011
(-0.240,-0.102)
Single amp tensor:
(4.245,-1.007)	(4.261,-1.011)
(-1.420,0.337)	(-1.425,0.338)
L1670XLineshape: (0.358,0.006)
L1670Xlineshapeold: (0.358,0.006)
L1670FF 1.0001.000
Creating DecayChain with k=2, two_j=1
L1670 k 2 Lineshape (0.358,0.006)1/2+ helicity 10 parity 011
(-0.404,0.715)
Single amp tensor:
(-1.425,0.338)	(1.420,-0.337)
(-4.261,1.011)	(4.245,-1.007)
"""
*/

    MandelstamTuple σs = {1, 3, 2};

    // Create Breit-Wigner lineshape function for L1670
    auto originalBreitWigner = createBW(1.670, 0.030);
    complex formFactor1 = 1.000;
    complex formFactor2 = 1.000;
    auto Xlineshape = [originalBreitWigner, formFactor1, formFactor2](double s) -> complex
    {
        return originalBreitWigner(s) * formFactor1 * formFactor2;
    };

    auto L16701 = createDecayChainCoupling(
        2,                                             // k-value
        Xlineshape,                                    // Lineshape function
        "1/2+",                                        // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {-1, 0}, false,  // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, true // ParityRecoupling for Hij
    );

    Tensor4Dcomp result1 = decays.amplitude4dcomp(*L16701, σs, 1);
    double prec = 1e-3;
    EXPECT_NEAR(result1[0][0][0][0].real(), 4.245, prec);
    EXPECT_NEAR(result1[0][0][0][0].imag(), -1.007, prec);
    EXPECT_NEAR(result1[0][0][0][1].real(), 4.261, prec);
    EXPECT_NEAR(result1[0][0][0][1].imag(), -1.011, prec);
    EXPECT_NEAR(result1[1][0][0][0].real(), -1.420, prec);
    EXPECT_NEAR(result1[1][0][0][0].imag(), 0.337, prec);
    EXPECT_NEAR(result1[1][0][0][1].real(), -1.425, prec);
    EXPECT_NEAR(result1[1][0][0][1].imag(), 0.338, prec);

    auto L16702 = createDecayChainCoupling(
        2,                                             // k-value
        Xlineshape,                                    // Lineshape function
        "1/2+",                                        // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {1, 0}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, true // ParityRecoupling for Hij
    );
    Tensor4Dcomp result2 = decays.amplitude4dcomp(*L16702, σs, 1);
    EXPECT_NEAR(result2[0][0][0][0].real(), -1.425, prec);
    EXPECT_NEAR(result2[0][0][0][0].imag(), 0.338, prec);
    EXPECT_NEAR(result2[0][0][0][1].real(), 1.420, prec);
    EXPECT_NEAR(result2[0][0][0][1].imag(), -0.337, prec);
    EXPECT_NEAR(result2[1][0][0][0].real(), -4.261, prec);
    EXPECT_NEAR(result2[1][0][0][0].imag(), 1.011, prec);
    EXPECT_NEAR(result2[1][0][0][1].real(), 4.245, prec);
    EXPECT_NEAR(result2[1][0][0][1].imag(), -1.007, prec);
}

TEST_F(JsonModelTest, L1690)
{
    /*
    L1690XLineshape: (0.049,0.002)
    L1690Xlineshapeold: (0.350,0.014) 1.690 0.070
    L1690FF 0.9210.152
    Creating DecayChain with k=2, two_j=3
    L1690 k 2 Lineshape (0.049,0.002)3/2- helicity -10 parity 010
    (-1.590,-0.454)
    Single amp tensor:
    (-0.242,0.199)	(-0.243,0.200)
    (0.368,-0.303)	(0.369,-0.304)
    L1690XLineshape: (0.049,0.002)
    L1690Xlineshapeold: (0.350,0.014) 1.690 0.070
    L1690FF 0.9210.152
    Creating DecayChain with k=2, two_j=3
    L1690 k 2 Lineshape (0.049,0.002)3/2- helicity 10 parity 010
    (-11.254,-1.457)
    Single amp tensor:
    (-0.369,0.304)	(0.368,-0.303)
    (-0.243,0.200)	(0.242,-0.199)
    */

    MandelstamTuple σs = {1, 3, 2};
    // Create Breit-Wigner lineshape function for L1690
    auto originalBreitWigner = createBW(1.690, 0.070);
    complex formFactor1 = 0.921;
    complex formFactor2 = 0.152;
    auto Xlineshape = [originalBreitWigner, formFactor1, formFactor2](double s) -> complex
    {
        return originalBreitWigner(s) * formFactor1 * formFactor2;
    };
    auto L16901 = createDecayChainCoupling(
        2,                                              // k-value
        Xlineshape,                                     // Lineshape function
        "3/2-",                                         // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                     // ThreeBodySystem
        RecouplingType::NoRecoupling, {-1, 0}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, false // ParityRecoupling for Hij
    );
    Tensor4Dcomp result1 = decays.amplitude4dcomp(*L16901, σs, 1);
    double prec = 1e-3;
    EXPECT_NEAR(result1[0][0][0][0].real(), -0.242, prec);
    EXPECT_NEAR(result1[0][0][0][0].imag(), 0.199, prec);
    EXPECT_NEAR(result1[0][0][0][1].real(), -0.243, prec);
    EXPECT_NEAR(result1[0][0][0][1].imag(), 0.200, prec);
    EXPECT_NEAR(result1[1][0][0][0].real(), 0.368, prec);
    EXPECT_NEAR(result1[1][0][0][0].imag(), -0.303, prec);
    EXPECT_NEAR(result1[1][0][0][1].real(), 0.369, prec);
    EXPECT_NEAR(result1[1][0][0][1].imag(), -0.304, prec);
    auto L16902 = createDecayChainCoupling(
        2,                                              // k-value
        Xlineshape,                                     // Lineshape function
        "3/2-",                                         // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                     // ThreeBodySystem
        RecouplingType::NoRecoupling, {1, 0}, false,    // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, false // ParityRecoupling for Hij
    );
    Tensor4Dcomp result2 = decays.amplitude4dcomp(*L16902, σs, 1);
    EXPECT_NEAR(result2[0][0][0][0].real(), -0.369, prec);
    EXPECT_NEAR(result2[0][0][0][0].imag(), 0.304, prec);
    EXPECT_NEAR(result2[0][0][0][1].real(), 0.368, prec);
    EXPECT_NEAR(result2[0][0][0][1].imag(), -0.303, prec);
    EXPECT_NEAR(result2[1][0][0][0].real(), -0.243, prec);
    EXPECT_NEAR(result2[1][0][0][0].imag(), 0.200, prec);
    EXPECT_NEAR(result2[1][0][0][1].real(), 0.242, prec);
    EXPECT_NEAR(result2[1][0][0][1].imag(), -0.199, prec);
}

TEST_F(JsonModelTest, L2000)
{
    /*
    L2000XLineshape: (0.251,0.023)
    L2000Xlineshapeold: (0.251,0.023) 1.988 0.179
    L2000FF 1.0001.000
    Creating DecayChain with k=2, two_j=1
    L2000 k 2 Lineshape (0.251,0.023)1/2+ helicity 10 parity 011
    (-3.066,-2.684)
    Single amp tensor:
    (0.293,0.109)	(-0.292,-0.109)
    (0.875,0.327)	(-0.871,-0.326)
    L2000XLineshape: (0.251,0.023)
    L2000Xlineshapeold: (0.251,0.023) 1.988 0.179
    L2000FF 1.0001.000
    Creating DecayChain with k=2, two_j=1
    L2000 k 2 Lineshape (0.251,0.023)1/2+ helicity -10 parity 011
    (-5.667,-5.384)
    Single amp tensor:
    (-0.871,-0.326)	(-0.875,-0.327)
    (0.292,0.109)	(0.293,0.109)
    */

    MandelstamTuple σs = {1, 3, 2};
    // Create Breit-Wigner lineshape function for L2000
    auto originalBreitWigner = createBW(1.988, 0.179);
    complex formFactor1 = 1.000;
    complex formFactor2 = 1.000;
    auto Xlineshape = [originalBreitWigner, formFactor1, formFactor2](double s) -> complex
    {
        return originalBreitWigner(s) * formFactor1 * formFactor2;
    };
    auto L20001 = createDecayChainCoupling(
        2,                                             // k-value
        Xlineshape,                                    // Lineshape function
        "1/2+",                                        // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {1, 0}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, true // ParityRecoupling for Hij
    );
    Tensor4Dcomp result1 = decays.amplitude4dcomp(*L20001, σs, 1);
    double prec = 1e-3;
    EXPECT_NEAR(result1[0][0][0][0].real(), 0.293, prec);
    EXPECT_NEAR(result1[0][0][0][0].imag(), 0.109, prec);
    EXPECT_NEAR(result1[0][0][0][1].real(), -0.292, prec);
    EXPECT_NEAR(result1[0][0][0][1].imag(), -0.109, prec);
    EXPECT_NEAR(result1[1][0][0][0].real(), 0.875, prec);
    EXPECT_NEAR(result1[1][0][0][0].imag(), 0.327, prec);
    EXPECT_NEAR(result1[1][0][0][1].real(), -0.871, 1e-2);
    EXPECT_NEAR(result1[1][0][0][1].imag(), -0.326, prec);
    auto L20002 = createDecayChainCoupling(
        2,                                             // k-value
        Xlineshape,                                    // Lineshape function
        "1/2+",                                        // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {-1, 0}, false,  // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, true // ParityRecoupling for Hij
    );
    Tensor4Dcomp result2 = decays.amplitude4dcomp(*L20002, σs, 1);
    EXPECT_NEAR(result2[0][0][0][0].real(), -0.871, 1e-2);
    EXPECT_NEAR(result2[0][0][0][0].imag(), -0.326, prec);
    EXPECT_NEAR(result2[0][0][0][1].real(), -0.875, prec);
    EXPECT_NEAR(result2[0][0][0][1].imag(), -0.327, prec);
    EXPECT_NEAR(result2[1][0][0][0].real(), 0.292, prec);
    EXPECT_NEAR(result2[1][0][0][0].imag(), 0.109, prec);
    EXPECT_NEAR(result2[1][0][0][1].real(), 0.293, prec);
    EXPECT_NEAR(result2[1][0][0][1].imag(), 0.109, prec);
}

TEST_F(JsonModelTest, D1232)
{
    /*
    D1232XLineshape: (0.304,0.029)
    D1232Xlineshapeold: (0.653,0.062) 1.232 0.117
    D1232FF 0.9450.493
    Creating DecayChain with k=3, two_j=3
    D1232 k 3 Lineshape (0.304,0.029)3/2+ helicity -10 parity 101
    (-10.917,4.915)
    Single amp tensor:
    (0.585,-0.175)	(-0.762,0.228)
    (-0.259,0.078)	(0.338,-0.101)
    D1232XLineshape: (0.304,0.029)
    D1232Xlineshapeold: (0.653,0.062) 1.232 0.117
    D1232FF 0.9450.493
    Creating DecayChain with k=3, two_j=3
    D1232 k 3 Lineshape (0.304,0.029)3/2+ helicity 10 parity 101
    (-20.916,7.293)
    Single amp tensor:
    (0.338,-0.101)	(0.259,-0.078)
    (0.762,-0.228)	(0.585,-0.175)
    */

    MandelstamTuple σs = {1, 3, 2};
    // Create Breit-Wigner lineshape function for D1232
    auto originalBreitWigner = createBW(1.232, 0.117);
    complex formFactor1 = 0.945;
    complex formFactor2 = 0.493;
    auto Xlineshape = [originalBreitWigner, formFactor1, formFactor2](double s) -> complex
    {
        return originalBreitWigner(s) * formFactor1 * formFactor2;
    };
    auto D12321 = createDecayChainCoupling(
        3,                                             // k-value
        Xlineshape,                                    // Lineshape function
        "3/2+",                                        // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {-1, 0}, false,  // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {1, 0}, true // ParityRecoupling for Hij
    );
    Tensor4Dcomp result1 = decays.amplitude4dcomp(*D12321, σs, 1);
    double prec = 1e-3;
    EXPECT_NEAR(result1[0][0][0][0].real(), 0.585, prec);
    EXPECT_NEAR(result1[0][0][0][0].imag(), -0.175, prec);
    EXPECT_NEAR(result1[0][0][0][1].real(), -0.762, prec);
    EXPECT_NEAR(result1[0][0][0][1].imag(), 0.228, prec);
    EXPECT_NEAR(result1[1][0][0][0].real(), -0.259, prec);
    EXPECT_NEAR(result1[1][0][0][0].imag(), 0.078, prec);
    EXPECT_NEAR(result1[1][0][0][1].real(), 0.338, prec);
    EXPECT_NEAR(result1[1][0][0][1].imag(), -0.101, prec);

    auto D12322 = createDecayChainCoupling(
        3,                                             // k-value
        Xlineshape,                                    // Lineshape function
        "3/2+",                                        // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {1, 0}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {1, 0}, true // ParityRecoupling for Hij
    );
    Tensor4Dcomp result2 = decays.amplitude4dcomp(*D12322, σs, 1);
    EXPECT_NEAR(result2[0][0][0][0].real(), 0.338, prec);
    EXPECT_NEAR(result2[0][0][0][0].imag(), -0.101, prec);
    EXPECT_NEAR(result2[0][0][0][1].real(), 0.259, prec);
    EXPECT_NEAR(result2[0][0][0][1].imag(), -0.078, prec);
    EXPECT_NEAR(result2[1][0][0][0].real(), 0.762, prec);
    EXPECT_NEAR(result2[1][0][0][0].imag(), -0.228, prec);
    EXPECT_NEAR(result2[1][0][0][1].real(), 0.585, prec);
    EXPECT_NEAR(result2[1][0][0][1].imag(), -0.175, prec);
}

TEST_F(JsonModelTest, D1600)
{
    /*
    D1600XLineshape: (0.168,0.031)
    D1600Xlineshapeold: (0.360,0.066) 1.640 0.300
    D1600FF 0.9450.493
    Creating DecayChain with k=3, two_j=3
    D1600 k 3 Lineshape (0.168,0.031)3/2+ helicity -10 parity 101
    (10.394,-2.849)
    Single amp tensor:
    (-0.295,-0.211)	(0.384,0.274)
    (0.131,0.093)	(-0.170,-0.122)
    D1600XLineshape: (0.168,0.031)
    D1600Xlineshapeold: (0.360,0.066) 1.640 0.300
    D1600FF 0.9450.493
    Creating DecayChain with k=3, two_j=3
    D1600 k 3 Lineshape (0.168,0.031)3/2+ helicity 10 parity 101
    (6.135,-0.914)
    Single amp tensor:
    (-0.170,-0.122)	(-0.131,-0.093)
    (-0.384,-0.274)	(-0.295,-0.211)
    */

    MandelstamTuple σs = {1, 3, 2};
    // Create Breit-Wigner lineshape function for D1600
    auto originalBreitWigner = createBW(1.640, 0.300);
    complex formFactor1 = 0.945;
    complex formFactor2 = 0.493;
    auto Xlineshape = [originalBreitWigner, formFactor1, formFactor2](double s) -> complex
    {
        return originalBreitWigner(s) * formFactor1 * formFactor2;
    };
    auto D16001 = createDecayChainCoupling(
        3,                                             // k-value
        Xlineshape,                                    // Lineshape function
        "3/2+",                                        // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {-1, 0}, false,  // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {1, 0}, true // ParityRecoupling for Hij
    );
    Tensor4Dcomp result1 = decays.amplitude4dcomp(*D16001, σs, 1);
    double prec = 1e-3;
    EXPECT_NEAR(result1[0][0][0][0].real(), -0.295, prec);
    EXPECT_NEAR(result1[0][0][0][0].imag(), -0.211, prec);
    EXPECT_NEAR(result1[0][0][0][1].real(), 0.384, prec);
    EXPECT_NEAR(result1[0][0][0][1].imag(), 0.274, prec);
    EXPECT_NEAR(result1[1][0][0][0].real(), 0.131, prec);
    EXPECT_NEAR(result1[1][0][0][0].imag(), 0.093, prec);
    EXPECT_NEAR(result1[1][0][0][1].real(), -0.170, prec);
    EXPECT_NEAR(result1[1][0][0][1].imag(), -0.122, prec);
    auto D16002 = createDecayChainCoupling(
        3,                                             // k-value
        Xlineshape,                                    // Lineshape function
        "3/2+",                                        // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {1, 0}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {1, 0}, true // ParityRecoupling for Hij
    );
    Tensor4Dcomp result2 = decays.amplitude4dcomp(*D16002, σs, 1);
    EXPECT_NEAR(result2[0][0][0][0].real(), -0.170, prec);
    EXPECT_NEAR(result2[0][0][0][0].imag(), -0.122, prec);
    EXPECT_NEAR(result2[0][0][0][1].real(), -0.131, prec);
    EXPECT_NEAR(result2[0][0][0][1].imag(), -0.093, prec);
    EXPECT_NEAR(result2[1][0][0][0].real(), -0.384, prec);
    EXPECT_NEAR(result2[1][0][0][0].imag(), -0.274, prec);
    EXPECT_NEAR(result2[1][0][0][1].real(), -0.295, prec);
    EXPECT_NEAR(result2[1][0][0][1].imag(), -0.211, prec);
}

TEST_F(JsonModelTest, D1700)
{
    /*
    D1700XLineshape: (0.032,0.007)
    D1700Xlineshapeold: (0.333,0.075) 1.690 0.380
    D1700FF 0.9450.101
    Creating DecayChain with k=3, two_j=3
    D1700 k 3 Lineshape (0.032,0.007)3/2- helicity -10 parity 100
    (-29.253,-4.044)
    Single amp tensor:
    (0.044,0.033)	(-0.057,-0.043)
    (-0.027,-0.020)	(0.036,0.027)
    D1700XLineshape: (0.032,0.007)
    D1700Xlineshapeold: (0.333,0.075) 1.690 0.380
    D1700FF 0.9450.101
    Creating DecayChain with k=3, two_j=3
    D1700 k 3 Lineshape (0.032,0.007)3/2- helicity 10 parity 100
    (-36.287,-5.935)
    Single amp tensor:
    (-0.036,-0.027)	(-0.027,-0.020)
    (-0.057,-0.043)	(-0.044,-0.033)
    */

    MandelstamTuple σs = {1, 3, 2};
    // Create Breit-Wigner lineshape function for D1700
    auto originalBreitWigner = createBW(1.690, 0.380);
    complex formFactor1 = 0.945;
    complex formFactor2 = 0.101;
    auto Xlineshape = [originalBreitWigner, formFactor1, formFactor2](double s) -> complex
    {
        return originalBreitWigner(s) * formFactor1 * formFactor2;
    };
    auto D17001 = createDecayChainCoupling(
        3,                                              // k-value
        Xlineshape,                                     // Lineshape function
        "3/2-",                                         // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                     // ThreeBodySystem
        RecouplingType::NoRecoupling, {-1, 0}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {1, 0}, false // ParityRecoupling for Hij
    );
    Tensor4Dcomp result1 = decays.amplitude4dcomp(*D17001, σs, 1);
    double prec = 1e-3;
    EXPECT_NEAR(result1[0][0][0][0].real(), 0.044, prec);
    EXPECT_NEAR(result1[0][0][0][0].imag(), 0.033, prec);
    EXPECT_NEAR(result1[0][0][0][1].real(), -0.057, prec);
    EXPECT_NEAR(result1[0][0][0][1].imag(), -0.043, prec);
    EXPECT_NEAR(result1[1][0][0][0].real(), -0.027, prec);
    EXPECT_NEAR(result1[1][0][0][0].imag(), -0.020, prec);
    EXPECT_NEAR(result1[1][0][0][1].real(), 0.036, prec);
    EXPECT_NEAR(result1[1][0][0][1].imag(), 0.027, prec);
    auto D17002 = createDecayChainCoupling(
        3,                                              // k-value
        Xlineshape,                                     // Lineshape function
        "3/2-",                                         // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                     // ThreeBodySystem
        RecouplingType::NoRecoupling, {1, 0}, false,    // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {1, 0}, false // ParityRecoupling for Hij
    );
    Tensor4Dcomp result2 = decays.amplitude4dcomp(*D17002, σs, 1);
    EXPECT_NEAR(result2[0][0][0][0].real(), -0.036, prec);
    EXPECT_NEAR(result2[0][0][0][0].imag(), -0.027, prec);
    EXPECT_NEAR(result2[0][0][0][1].real(), -0.027, prec);
    EXPECT_NEAR(result2[0][0][0][1].imag(), -0.020, prec);
    EXPECT_NEAR(result2[1][0][0][0].real(), -0.057, prec);
    EXPECT_NEAR(result2[1][0][0][0].imag(), -0.043, prec);
    EXPECT_NEAR(result2[1][0][0][1].real(), -0.044, prec);
    EXPECT_NEAR(result2[1][0][0][1].imag(), -0.033, prec);
}

// Test for decay chain with ParityRecoupling
TEST_F(JsonModelTest, K892_Test)
{
    /*
    K892XLineshape: (0.593,0.031)
    K892Xlineshapeold: (1.244,0.066) 0.895 0.047
    K892FF 1.0000.477
    Creating DecayChain with k=1, two_j=2
    K892 k 1 Lineshape (0.593,0.031)1+ helicity -2-1 parity 001
    (1.722,-1.481)
    Single amp tensor:
    (2.765,-0.591)	(0.000,0.000)
    (0.000,0.000)	(0.000,0.000)
    K892XLineshape: (0.593,0.031)
    K892Xlineshapeold: (1.244,0.066) 0.895 0.047
    K892FF 1.0000.477
    Creating DecayChain with k=1, two_j=2
    K892 k 1 Lineshape (0.593,0.031)1+ helicity 01 parity 001
    (-1.050,-6.000)
    Single amp tensor:
    (0.000,0.000)	(0.000,0.000)
    (0.788,-0.169)	(0.000,0.000)
    K892XLineshape: (0.593,0.031)
    K892Xlineshapeold: (1.244,0.066) 0.895 0.047
    K892FF 1.0000.477
    Creating DecayChain with k=1, two_j=2
    K892 k 1 Lineshape (0.593,0.031)1+ helicity 0-1 parity 001
    (1.444,0.000)
    Single amp tensor:
    (0.000,0.000)	(-0.788,0.169)
    (0.000,0.000)	(0.000,0.000)
    K892XLineshape: (0.593,0.031)
    K892Xlineshapeold: (1.244,0.066) 0.895 0.047
    K892FF 1.0000.477
    Creating DecayChain with k=1, two_j=2
    K892 k 1 Lineshape (0.593,0.031)1+ helicity 21 parity 001
    (-4.537,-4.756)
    Single amp tensor:
    (0.000,0.000)	(0.000,0.000)
    (0.000,0.000)	(2.765,-0.591)
    */

    MandelstamTuple σs = {1, 3, 2};
    // Create Breit-Wigner lineshape function for K892
    auto originalBreitWigner = createBW(0.895, 0.047);
    complex formFactor1 = 1.000;
    complex formFactor2 = 0.477;
    auto Xlineshape = [originalBreitWigner, formFactor1, formFactor2](double s) -> complex
    {
        return originalBreitWigner(s) * formFactor1 * formFactor2;
    };

    auto K892_1 = createDecayChainCoupling(
        1,                                             // k-value
        Xlineshape,                                    // Lineshape function
        "1+",                                          // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {-2, -1}, false, // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 0}, true // ParityRecoupling for Hij
    );
    Tensor4Dcomp result1 = decays.amplitude4dcomp(*K892_1, σs, 1);
    double prec = 1e-3;
    EXPECT_NEAR(result1[0][0][0][0].real(), 2.755, prec);
    EXPECT_NEAR(result1[0][0][0][0].imag(), -0.582, prec);
    EXPECT_NEAR(result1[0][0][0][1].real(), 0.000, prec);
    EXPECT_NEAR(result1[0][0][0][1].imag(), 0.000, prec);
    EXPECT_NEAR(result1[1][0][0][0].real(), 0.000, prec);
    EXPECT_NEAR(result1[1][0][0][0].imag(), 0.000, prec);
    EXPECT_NEAR(result1[1][0][0][1].real(), 0.000, prec);
    EXPECT_NEAR(result1[1][0][0][1].imag(), 0.000, prec);

    auto K892_2 = createDecayChainCoupling(
        1,                                             // k-value
        Xlineshape,                                    // Lineshape function
        "1+",                                          // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {0, 1}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 0}, true // ParityRecoupling for Hij
    );
    Tensor4Dcomp result2 = decays.amplitude4dcomp(*K892_2, σs, 1);
    EXPECT_NEAR(result2[0][0][0][0].real(), 0.000, prec);
    EXPECT_NEAR(result2[0][0][0][0].imag(), 0.000, prec);
    EXPECT_NEAR(result2[0][0][0][1].real(), 0.000, prec);
    EXPECT_NEAR(result2[0][0][0][1].imag(), 0.000, prec);
    EXPECT_NEAR(result2[1][0][0][0].real(), 0.785, prec);
    EXPECT_NEAR(result2[1][0][0][0].imag(), -0.166, prec);
    EXPECT_NEAR(result2[1][0][0][1].real(), 0.000, prec);
    EXPECT_NEAR(result2[1][0][0][1].imag(), 0.000, prec);

    auto K892_3 = createDecayChainCoupling(
        1,                                             // k-value
        Xlineshape,                                    // Lineshape function
        "1+",                                          // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {0, -1}, false,  // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 0}, true // ParityRecoupling for Hij
    );
    Tensor4Dcomp result3 = decays.amplitude4dcomp(*K892_3, σs, 1);
    EXPECT_NEAR(result3[0][0][0][0].real(), 0.000, prec);
    EXPECT_NEAR(result3[0][0][0][0].imag(), 0.000, prec);
    EXPECT_NEAR(result3[0][0][0][1].real(), -0.785, prec);
    EXPECT_NEAR(result3[0][0][0][1].imag(), 0.166, prec);
    EXPECT_NEAR(result3[1][0][0][0].real(), 0.000, prec);
    EXPECT_NEAR(result3[1][0][0][0].imag(), 0.000, prec);
    EXPECT_NEAR(result3[1][0][0][1].real(), 0.000, prec);
    EXPECT_NEAR(result3[1][0][0][1].imag(), 0.000, prec);

    auto K892_4 = createDecayChainCoupling(
        1,                                             // k-value
        Xlineshape,                                    // Lineshape function
        "1+",                                          // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {2, 1}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 0}, true // ParityRecoupling for Hij
    );
    Tensor4Dcomp result4 = decays.amplitude4dcomp(*K892_4, σs, 1);
    EXPECT_NEAR(result4[0][0][0][0].real(), 0.000, prec);
    EXPECT_NEAR(result4[0][0][0][0].imag(), 0.000, prec);
    EXPECT_NEAR(result4[0][0][0][1].real(), 0.000, prec);
    EXPECT_NEAR(result4[0][0][0][1].imag(), 0.000, prec);
    EXPECT_NEAR(result4[1][0][0][0].real(), 0.0, prec);
    EXPECT_NEAR(result4[1][0][0][0].imag(), 0.0, prec);
    EXPECT_NEAR(result4[1][0][0][1].real(), 2.755, prec);
    EXPECT_NEAR(result4[1][0][0][1].imag(), -0.582, prec);
}

TEST_F(JsonModelTest, SumTest)
{
    /*
    const Lcpkpi = ThreeBodyDecay(
    ["Λ15201", "Λ15202", "Λ16001", "Λ16002",
     "Λ16701", "Λ16702", "Λ16901", "Λ16902",
     "Λ20001", "Λ20002", "D123201", "D123202",
     "D160001", "D160002", "D170001", "D170002",
     "K89201", "K89202", "K89203", "K89204"] .=> zip(
        [3.605+0.544im, -1.970+18.380im,
         10.063-1.216im, -6.987-4.450im,
         -0.240-0.102im, -0.404+0.715im,
         -1.590-0.454im, -11.254-1.457im,
         -3.066-2.684im, -5.667-5.384im,
         -10.917+4.915im, -20.916+7.293im,
         10.394-2.849im, 6.135-0.914im,
         -29.253-4.044im, -36.287-5.935im,
         1.722-1.481im, -1.050-6.000im,
         1.444+0im, -4.537-4.756im],
        [Λ15201, Λ15202, Λ16001, Λ16002,
         Λ16701, Λ16702, Λ16901, Λ16902,
         Λ20001, Λ20002, D123201, D123202,
         D160001, D160002, D170001, D170002,
         K89201, K89202, K89203, K89204],
    ),
)
    */

    MandelstamTuple σs = {1, 3, 2};

    // Create Breit-Wigner lineshape function for L1520
    auto originalBreitWigner = createBW(1.518, 0.015);
    complex formFactor1 = 0.921;
    complex formFactor2 = 0.152;
    auto Xlineshape = [originalBreitWigner, formFactor1, formFactor2](double s) -> complex
    {
        return originalBreitWigner(s) * formFactor1 * formFactor2;
    };

    // std::cout << "L1520XLineshape: " << Xlineshape(0.0) << std::endl;
    //  Create decay chain with ParityRecoupling for HRk
    auto L15201 = createDecayChainCoupling(
        2,                                              // k-value
        Xlineshape,                                     // Lineshape function
        "3/2-",                                         // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                     // ThreeBodySystem
        RecouplingType::NoRecoupling, {-1, 0}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, false // ParityRecoupling for Hij
    );
    // std::cout << "L1520Xlineshapeold: " << Xlineshape(0.0) << std::endl;

    // Calculate amplitude tensor
    Tensor4Dcomp resultL15201 = decays.amplitude4dcomp(*L15201, σs, 1);

    // Create decay chain with ParityRecoupling for HRk
    auto L15202 = createDecayChainCoupling(
        2,                                              // k-value
        Xlineshape,                                     // Lineshape function
        "3/2-",                                         // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                     // ThreeBodySystem
        RecouplingType::NoRecoupling, {1, 0}, false,    // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, false // ParityRecoupling for Hij
    );

    // Calculate amplitude tensor
    Tensor4Dcomp resultL15202 = decays.amplitude4dcomp(*L15202, σs, 1);

    // L1600

    auto originalBreitWignerL1600 = createBW(1.630, 0.250);
    complex formFactor1L1600 = 1.000;
    complex formFactor2L1600 = 0.577;
    auto XlineshapeL1600 = [originalBreitWignerL1600, formFactor1L1600, formFactor2L1600](double s) -> complex
    {
        return originalBreitWignerL1600(s) * formFactor1L1600 * formFactor2L1600;
    };

    auto L16001 = createDecayChainCoupling(
        2,                                              // k-value
        XlineshapeL1600,                                // Lineshape function
        "1/2-",                                         // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                     // ThreeBodySystem
        RecouplingType::NoRecoupling, {1, 0}, false,    // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, false // ParityRecoupling for Hij
    );

    Tensor4Dcomp resultL16001 = decays.amplitude4dcomp(*L16001, σs, 1);

    auto L16002 = createDecayChainCoupling(
        2,                                              // k-value
        XlineshapeL1600,                                // Lineshape function
        "1/2-",                                         // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                     // ThreeBodySystem
        RecouplingType::NoRecoupling, {-1, 0}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, false // ParityRecoupling for Hij
    );
    Tensor4Dcomp resultL16002 = decays.amplitude4dcomp(*L16002, σs, 1);

    // L1670
    // Create Breit-Wigner lineshape function for L1670
    auto originalBreitWignerL1670 = createBW(1.670, 0.030);
    complex formFactor1L1670 = 1.000;
    complex formFactor2L1670 = 1.000;
    auto XlineshapeL1670 = [originalBreitWignerL1670, formFactor1L1670, formFactor2L1670](double s) -> complex
    {
        return originalBreitWignerL1670(s) * formFactor1L1670 * formFactor2L1670;
    };

    auto L16701 = createDecayChainCoupling(
        2,                                             // k-value
        XlineshapeL1670,                               // Lineshape function
        "1/2+",                                        // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {-1, 0}, false,  // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, true // ParityRecoupling for Hij
    );

    Tensor4Dcomp resultL16701 = decays.amplitude4dcomp(*L16701, σs, 1);
    auto L16702 = createDecayChainCoupling(
        2,                                             // k-value
        XlineshapeL1670,                               // Lineshape function
        "1/2+",                                        // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {1, 0}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, true // ParityRecoupling for Hij
    );
    Tensor4Dcomp resultL16702 = decays.amplitude4dcomp(*L16702, σs, 1);

    // L1690
    auto originalBreitWignerL1690 = createBW(1.690, 0.070);
    complex formFactor1L1690 = 0.921;
    complex formFactor2L1690 = 0.152;
    auto XlineshapeL1690 = [originalBreitWignerL1690, formFactor1L1690, formFactor2L1690](double s) -> complex
    {
        return originalBreitWignerL1690(s) * formFactor1L1690 * formFactor2L1690;
    };
    auto L16901 = createDecayChainCoupling(
        2,                                              // k-value
        XlineshapeL1690,                                // Lineshape function
        "3/2-",                                         // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                     // ThreeBodySystem
        RecouplingType::NoRecoupling, {-1, 0}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, false // ParityRecoupling for Hij
    );
    Tensor4Dcomp resultL16901 = decays.amplitude4dcomp(*L16901, σs, 1);
    auto L16902 = createDecayChainCoupling(
        2,                                              // k-value
        XlineshapeL1690,                                // Lineshape function
        "3/2-",                                         // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                     // ThreeBodySystem
        RecouplingType::NoRecoupling, {1, 0}, false,    // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, false // ParityRecoupling for Hij
    );
    Tensor4Dcomp resultL16902 = decays.amplitude4dcomp(*L16902, σs, 1);

    // L2000
    auto originalBreitWignerL2000 = createBW(1.988, 0.179);
    complex formFactor1L2000 = 1.000;
    complex formFactor2L2000 = 1.000;
    auto XlineshapeL2000 = [originalBreitWignerL2000, formFactor1L2000, formFactor2L2000](double s) -> complex
    {
        return originalBreitWignerL2000(s) * formFactor1L2000 * formFactor2L2000;
    };
    auto L20001 = createDecayChainCoupling(
        2,                                             // k-value
        XlineshapeL2000,                               // Lineshape function
        "1/2+",                                        // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {1, 0}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, true // ParityRecoupling for Hij
    );
    Tensor4Dcomp resultL20001 = decays.amplitude4dcomp(*L20001, σs, 1);
    auto L20002 = createDecayChainCoupling(
        2,                                             // k-value
        XlineshapeL2000,                               // Lineshape function
        "1/2+",                                        // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {-1, 0}, false,  // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, true // ParityRecoupling for Hij
    );
    Tensor4Dcomp resultL20002 = decays.amplitude4dcomp(*L20002, σs, 1);

    // D1232
    auto originalBreitWignerD1232 = createBW(1.232, 0.117);
    complex formFactor1D1232 = 0.945;
    complex formFactor2D1232 = 0.493;
    auto XlineshapeD1232 = [originalBreitWignerD1232, formFactor1D1232, formFactor2D1232](double s) -> complex
    {
        return originalBreitWignerD1232(s) * formFactor1D1232 * formFactor2D1232;
    };
    auto D12321 = createDecayChainCoupling(
        3,                                             // k-value
        XlineshapeD1232,                               // Lineshape function
        "3/2+",                                        // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {-1, 0}, false,  // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {1, 0}, true // ParityRecoupling for Hij
    );
    Tensor4Dcomp resultD12321 = decays.amplitude4dcomp(*D12321, σs, 1);
    auto D12322 = createDecayChainCoupling(
        3,                                             // k-value
        XlineshapeD1232,                               // Lineshape function
        "3/2+",                                        // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {1, 0}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {1, 0}, true // ParityRecoupling for Hij
    );
    Tensor4Dcomp resultD12322 = decays.amplitude4dcomp(*D12322, σs, 1);

    // D1600
    auto originalBreitWignerD1600 = createBW(1.640, 0.300);
    complex formFactor1D1600 = 0.945;
    complex formFactor2D1600 = 0.493;
    auto XlineshapeD1600 = [originalBreitWignerD1600, formFactor1D1600, formFactor2D1600](double s) -> complex
    {
        return originalBreitWignerD1600(s) * formFactor1D1600 * formFactor2D1600;
    };
    auto D16001 = createDecayChainCoupling(
        3,                                             // k-value
        XlineshapeD1600,                               // Lineshape function
        "3/2+",                                        // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {-1, 0}, false,  // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {1, 0}, true // ParityRecoupling for Hij
    );
    Tensor4Dcomp resultD16001 = decays.amplitude4dcomp(*D16001, σs, 1);

    auto D16002 = createDecayChainCoupling(
        3,                                             // k-value
        XlineshapeD1600,                               // Lineshape function
        "3/2+",                                        // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {1, 0}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {1, 0}, true // ParityRecoupling for Hij
    );
    Tensor4Dcomp resultD16002 = decays.amplitude4dcomp(*D16002, σs, 1);

    auto originalBreitWignerD1700 = createBW(1.690, 0.380);
    complex formFactor1D1700 = 0.945;
    complex formFactor2D1700 = 0.101;
    auto XlineshapeD1700 = [originalBreitWignerD1700, formFactor1D1700, formFactor2D1700](double s) -> complex
    {
        return originalBreitWignerD1700(s) * formFactor1D1700 * formFactor2D1700;
    };
    auto D17001 = createDecayChainCoupling(
        3,                                              // k-value
        XlineshapeD1700,                                // Lineshape function
        "3/2-",                                         // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                     // ThreeBodySystem
        RecouplingType::NoRecoupling, {-1, 0}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {1, 0}, false // ParityRecoupling for Hij
    );
    Tensor4Dcomp resultD17001 = decays.amplitude4dcomp(*D17001, σs, 1);
    auto D17002 = createDecayChainCoupling(
        3,                                              // k-value
        XlineshapeD1700,                                // Lineshape function
        "3/2-",                                         // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                     // ThreeBodySystem
        RecouplingType::NoRecoupling, {1, 0}, false,    // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {1, 0}, false // ParityRecoupling for Hij
    );
    Tensor4Dcomp resultD17002 = decays.amplitude4dcomp(*D17002, σs, 1);

    // K892
    auto originalBreitWignerK892 = createBW(0.895, 0.047);
    complex formFactor1K892 = 1.000;
    complex formFactor2K892 = 0.477;
    auto XlineshapeK892 = [originalBreitWignerK892, formFactor1K892, formFactor2K892](double s) -> complex
    {
        return originalBreitWignerK892(s) * formFactor1K892 * formFactor2K892;
    };

    auto K892_1 = createDecayChainCoupling(
        1,                                             // k-value
        XlineshapeK892,                                // Lineshape function
        "1+",                                          // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {-2, -1}, false, // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 0}, true // ParityRecoupling for Hij
    );

    auto K892_2 = createDecayChainCoupling(
        1,                                             // k-value
        XlineshapeK892,                                // Lineshape function
        "1+",                                          // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {0, 1}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 0}, true // ParityRecoupling for Hij
    );

    auto K892_3 = createDecayChainCoupling(
        1,                                             // k-value
        XlineshapeK892,                                // Lineshape function
        "1+",                                          // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {0, -1}, false,  // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 0}, true // ParityRecoupling for Hij
    );

    auto K892_4 = createDecayChainCoupling(
        1,                                             // k-value
        XlineshapeK892,                                // Lineshape function
        "1+",                                          // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {2, 1}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 0}, true // ParityRecoupling for Hij
    );

    /*
    const Lcpkpi = ThreeBodyDecay(
    ["Λ15201", "Λ15202", "Λ16001", "Λ16002",
     "Λ16701", "Λ16702", "Λ16901", "Λ16902",
     "Λ20001", "Λ20002", "D123201", "D123202",
     "D160001", "D160002", "D170001", "D170002",
     "K89201", "K89202", "K89203", "K89204"] .=> zip(
        [3.605+0.544im, -1.970+18.380im,
         10.063-1.216im, -6.987-4.450im,
         -0.240-0.102im, -0.404+0.715im,
         -1.590-0.454im, -11.254-1.457im,
         -3.066-2.684im, -5.667-5.384im,
         -10.917+4.915im, -20.916+7.293im,
         10.394-2.849im, 6.135-0.914im,
         -29.253-4.044im, -36.287-5.935im,
         1.722-1.481im, -1.050-6.000im,
         1.444+0im, -4.537-4.756im],
        [Λ15201, Λ15202, Λ16001, Λ16002,
         Λ16701, Λ16702, Λ16901, Λ16902,
         Λ20001, Λ20002, D123201, D123202,
         D160001, D160002, D170001, D170002,
         K89201, K89202, K89203, K89204],
    ),
)
    */
    // Sum all results
    ThreeBodyAmplitudeModel model;
    model.add(L15201, "L1520", complex(3.605, 0.544));
    model.add(L15202, "L1520", complex(-1.970, 18.380));
    model.add(L16001, "L1600", complex(10.063, -1.216));
    model.add(L16002, "L1600", complex(-6.987, -4.450));
    model.add(L16701, "L1670", complex(-0.240, -0.102));
    model.add(L16702, "L1670", complex(-0.404, 0.715));
    model.add(L16901, "L1690", complex(-1.590, -0.454));
    model.add(L16902, "L1690", complex(-11.254, -1.457));
    model.add(L20001, "L2000", complex(-3.066, -2.684));
    model.add(L20002, "L2000", complex(-5.667, -5.384));
    model.add(D12321, "D1232", complex(-10.917, 4.915));
    model.add(D12322, "D1232", complex(-20.916, 7.293));
    model.add(D16001, "D1600", complex(10.394, -2.849));
    model.add(D16002, "D1600", complex(6.135, -0.914));
    model.add(D17001, "D1700", complex(-29.253, -4.044));
    model.add(D17002, "D1700", complex(-36.287, -5.935));
    model.add(K892_1, "K892", complex(1.722, -1.481));
    model.add(K892_2, "K892", complex(-1.050, -6.000));
    model.add(K892_3, "K892", complex(1.444, 0.000));
    model.add(K892_4, "K892", complex(-4.537, -4.756));

    Tensor4Dcomp result = model.amplitude4d(σs, 1);

    // ComplexF64[3.5276468872563234 - 8.266206050243916im; -6.597152023966293 - 3.9520803173505357im;;;; 4.597204011832391 + 21.223154258348426im; -25.919027836579062 + 0.31176559801143355im]
    double prec = 1e-4;
    EXPECT_NEAR(result[0][0][0][0].real(), 3.5276, prec);
    EXPECT_NEAR(result[0][0][0][0].imag(), -8.2662, prec);
    EXPECT_NEAR(result[1][0][0][0].real(), -6.5972, prec);
    EXPECT_NEAR(result[1][0][0][0].imag(), -3.9521, prec);
    EXPECT_NEAR(result[0][0][0][1].real(), 4.5972, prec);
    EXPECT_NEAR(result[0][0][0][1].imag(), 21.2232, prec);
    EXPECT_NEAR(result[1][0][0][1].real(), -25.9190, prec);
    EXPECT_NEAR(result[1][0][0][1].imag(), 0.3118, prec);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
