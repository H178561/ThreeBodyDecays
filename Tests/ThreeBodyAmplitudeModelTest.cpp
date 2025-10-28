#include <gtest/gtest.h>
#include "include/ThreeBodyAmplitudeModel.hh"
#include "include/ThreeBodyDecays.hh"
#include "include/ClebschGordan.hh"
#include <iostream>
#include <iomanip>
#include <complex>
#include <memory>
#include <cmath>

// Test fixture for ThreeBodyAmplitudeModel tests
class ThreeBodyAmplitudeModelTest : public ::testing::Test
{
protected:
    ThreeBodyAmplitudeModelTest()
        : tbs(nullptr),
          conservingParities('-', '+', '-', '+') // Initialize this in constructor
    {
    }

    void SetUp() override
    {
        // Constants from the tutorial example (Λb → J/ψ p K)
        constants["mJpsi"] = 3.09;
        constants["mp"] = 0.938;
        constants["mK"] = 0.49367;
        constants["mLb"] = 5.62;

        // Create masses and spins matching the tutorial
        ms = {
            constants["mJpsi"],
            constants["mp"],
            constants["mK"],
            constants["mLb"]};

        // Two times the spins (doubled representation)
        spins = {2, 1, 0, 1}; // J/ψ (spin 1), p (spin 1/2), K (spin 0), Λb (spin 1/2)

        // Create the system
        tbs = std::make_shared<ThreeBodySystem>(ms, spins);

        // Create ThreeBodyDecays for calculations
        decays = std::make_shared<ThreeBodyDecays>();
    }

    // Create Breit-Wigner lineshape function
    std::function<complex(double)> createBW(double mass, double width)
    {
        return [mass, width](double sigma) -> complex
        {
            return complex(1.0, 0.0) / complex(mass * mass - sigma, -mass * width);
        };
    }

    // Helper for comparing complex values with tolerance
    void expectNearComplex(const complex &expected, const complex &actual, double tol = 1e-5)
    {
        EXPECT_NEAR(std::real(expected), std::real(actual), tol);
        EXPECT_NEAR(std::imag(expected), std::imag(actual), tol);
    }

    std::map<std::string, double> constants;
    ThreeBodyMasses ms;
    ThreeBodySpins spins;
    std::shared_ptr<ThreeBodySystem> tbs; // Changed to shared_ptr
    ThreeBodyParities conservingParities;
    std::shared_ptr<ThreeBodyDecays> decays;
};

// Test case for Lambda resonance chains (2+3)
TEST_F(ThreeBodyAmplitudeModelTest, LambdaResonanceChains)
{
    // Create Lambda resonance chains matching the tutorial
    auto Lambda1520 = createDecayChainLS(
        1,                        // k-value
        createBW(1.5195, 0.0156), // Lineshape function
        "3/2-",                   // jp (Spin-parity)
        conservingParities,       // Parities
        *tbs                      // Dereference ThreeBodySystem
    );

    auto Lambda1690 = createDecayChainLS(
        1,                      // k-value
        createBW(1.685, 0.050), // Lineshape function
        "1/2+",                 // jp (Spin-parity)
        conservingParities,     // Parities
        *tbs                    // Dereference ThreeBodySystem
    );

    auto Lambda1810 = createDecayChainLS(
        1,                     // k-value
        createBW(1.80, 0.090), // Lineshape function
        "5/2+",                // jp (Spin-parity)
        conservingParities,    // Parities
        *tbs                   // Dereference ThreeBodySystem
    );

    // Verify the chains were created successfully
    ASSERT_NE(Lambda1520, nullptr);
    ASSERT_NE(Lambda1690, nullptr);
    ASSERT_NE(Lambda1810, nullptr);

    // Check basic properties of the chains
    EXPECT_EQ(Lambda1520->k, 1);
    EXPECT_EQ(Lambda1690->k, 1);
    EXPECT_EQ(Lambda1810->k, 1);

    // Test lineshape values at a specific point
    double testSigma = 2.5;
    expectNearComplex(complex(1.0, 0.0) / complex(1.5195 * 1.5195 - testSigma, -1.5195 * 0.0156),
                      Lambda1520->Xlineshape(testSigma));
}

// Test case for Pentaquark resonance chains (1+2)
TEST_F(ThreeBodyAmplitudeModelTest, PentaquarkResonanceChains)
{
    // Create Pentaquark resonance chains matching the tutorial
    auto Pc4312 = createDecayChainLS(
        3,                      // k-value
        createBW(4.312, 0.015), // Lineshape function
        "1/2+",                 // jp (Spin-parity)
        conservingParities,     // Parities
        *tbs                    // Dereference ThreeBodySystem
    );

    auto Pc4440 = createDecayChainLS(
        3,                      // k-value
        createBW(4.440, 0.010), // Lineshape function
        "1/2+",                 // jp (Spin-parity)
        conservingParities,     // Parities
        *tbs                    // Dereference ThreeBodySystem
    );

    auto Pc4457 = createDecayChainLS(
        3,                      // k-value
        createBW(4.457, 0.020), // Lineshape function
        "3/2+",                 // jp (Spin-parity)
        conservingParities,     // Parities
        *tbs                    // Dereference ThreeBodySystem
    );

    // Verify the chains were created successfully
    ASSERT_NE(Pc4312, nullptr);
    ASSERT_NE(Pc4440, nullptr);
    ASSERT_NE(Pc4457, nullptr);

    // Check basic properties of the chains
    EXPECT_EQ(Pc4312->k, 3);
    EXPECT_EQ(Pc4440->k, 3);
    EXPECT_EQ(Pc4457->k, 3);
}

// Test case for amplitude model creation and intensity calculation
TEST_F(ThreeBodyAmplitudeModelTest, AmplitudeModelIntensity)
{
    // Create all resonance chains
    auto Lambda1520 = createDecayChainLS(
        1, createBW(1.5195, 0.0156), "3/2+", conservingParities, *tbs);
    auto Lambda1690 = createDecayChainLS(
        1, createBW(1.685, 0.050), "1/2+", conservingParities, *tbs);
    auto Lambda1810 = createDecayChainLS(
        1, createBW(1.80, 0.090), "5/2+", conservingParities, *tbs);
    auto Pc4312 = createDecayChainLS(
        3, createBW(4.312, 0.015), "1/2+", conservingParities, *tbs);
    auto Pc4440 = createDecayChainLS(
        3, createBW(4.440, 0.010), "1/2+", conservingParities, *tbs);
    auto Pc4457 = createDecayChainLS(
        3, createBW(4.457, 0.020), "3/2+", conservingParities, *tbs);

    // Create amplitude model and add chains with coefficients
    ThreeBodyAmplitudeModel model;
    model.add(Lambda1520, "Lambda1520", complex(2.0, 0.0));
    model.add(Lambda1690, "Lambda1690", complex(2.1, 0.0));
    model.add(Lambda1810, "Lambda1810", complex(0.0, 1.4));
    model.add(Pc4312, "Pc4312", complex(0.4, 0.0));
    model.add(Pc4440, "Pc4440", complex(0.0, 0.3));
    model.add(Pc4457, "Pc4457", complex(0.0, -0.8));

    // Verify model size
    EXPECT_EQ(model.size(), 6);

    // Generate a test point in Dalitz plot
    MandelstamTuple σs = decays->x2σs({0.3, 0.3}, ms, 1);

    // Calculate intensity at this point
    double intensity = model.intensity(σs, 1);

    // The exact value depends on implementation details, but should be positive
    EXPECT_GT(intensity, 0.0);
    // Calculate individual component intensities
    auto componentIntensities = model.component_intensities(σs, 1);
    EXPECT_EQ(componentIntensities.size(), 6);

    // Sum of individual intensities should be different from total due to interference
    double sumIndividual = 0.0;
    for (const auto &value : componentIntensities)
    {
        sumIndividual += value;
    }

    // With interference terms, total intensity is not equal to sum of individual intensities
    // unless the phases are chosen so that there's no interference
    EXPECT_NE(intensity, sumIndividual);
}



// Add a new test case to match the 10-tutorial.jl result
TEST_F(ThreeBodyAmplitudeModelTest, MatchJuliaTutorialResult)
{
    // Create decay chains as in the tutorial
    // Since we don't have the exact parameters from the tutorial, we'll use approximate values
    auto Lambda1520 = createDecayChainLS(
        1, createBW(1.5195, 0.0156), "3/2+", conservingParities, *tbs);
    auto Lambda1690 = createDecayChainLS(
        1, createBW(1.685, 0.050), "1/2+", conservingParities, *tbs);
    auto Pc4312 = createDecayChainLS(
        3, createBW(4.312, 0.015), "1/2+", conservingParities, *tbs);

    // Create the model with the same coefficients
    ThreeBodyAmplitudeModel model;
    model.add(Lambda1520, "Lambda1520", complex(2.0, 0.0));
    model.add(Lambda1690, "Lambda1690", complex(2.1, 0.0)); // Using a phase to get imaginary component
    model.add(Pc4312, "Pc4312", complex(0.4, 0.0));         // Pure imaginary coefficient

    // Create the same point in phase space
    MandelstamTuple σs = decays->x2σs({0.3, 0.3}, ms, 1);
    // σ1 = 4.501751135564419, σ2 = 21.495750373843414, σ3 = 16.258552559492166
    MandelstamTuple σs2 = {4.501751135564419, 21.495750373843414, 16.258552559492166};

    // Use the same helicity configuration
    std::vector<int> two_λs = {2, 1, 0, 1}; // Doubled representation

    // Calculate amplitude
    complex amp = model.amplitude(σs2, two_λs, 1);
    Tensor4Dcomp amp2 = model.amplitude4d(σs2, 1);


    double intensity_valued = model.intensity(σs2, 1);

    ThreeBodyDecays tbd;
    complex res = tbd.amplitude(*Lambda1520, σs2, two_λs, 1);


    Tensor4Dcomp resultdalignedccomp = tbd.aligned_amplitude4dcomp(*Lambda1520, σs2);
    // Verify tensor dimensions
    std::cout << "Tensor4D resultdalignedccomp : " << std::endl;
    for (int i = 0; i < resultdalignedccomp.size(); ++i)
    {
        for (int j = 0; j < resultdalignedccomp[0].size(); ++j)
        {
            for (int k = 0; k < resultdalignedccomp[0][0].size(); ++k)
            {
                for (int z = 0; z < resultdalignedccomp[0][0][0].size(); ++z)
                {
                    std::cout << resultdalignedccomp[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }

    Tensor4Dcomp resultdccomp = tbd.amplitude4dcomp(*Lambda1520, σs2, 1);

    Tensor4Dcomp resultl1690 = tbd.amplitude4dcomp(*Lambda1690, σs2, 1);
    Tensor4Dcomp resultpc4312 = tbd.amplitude4dcomp(*Pc4312, σs2, 1);
    // Verify tensor dimensions
    std::cout << "\nTensor4D Lambda1520 : " << std::endl;
    for (int i = 0; i < resultdccomp.size(); ++i)
    {
        for (int j = 0; j < resultdccomp[0].size(); ++j)
        {
            for (int k = 0; k < resultdccomp[0][0].size(); ++k)
            {
                for (int z = 0; z < resultdccomp[0][0][0].size(); ++z)
                {
                    std::cout << resultdccomp[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }

    std::cout << "\nTensor4D Lambda1690 : " << std::endl;
    for (int i = 0; i < resultl1690.size(); ++i)
    {
        for (int j = 0; j < resultl1690[0].size(); ++j)
        {
            for (int k = 0; k < resultl1690[0][0].size(); ++k)
            {
                for (int z = 0; z < resultl1690[0][0][0].size(); ++z)
                {
                    std::cout << resultl1690[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }
    std::cout << "\nTensor4D Pc4312 : " << std::endl;
    for (int i = 0; i < resultpc4312.size(); ++i)
    {
        for (int j = 0; j < resultpc4312[0].size(); ++j)
        {
            for (int k = 0; k < resultpc4312[0][0].size(); ++k)
            {
                for (int z = 0; z < resultpc4312[0][0][0].size(); ++z)
                {
                    std::cout << resultpc4312[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }

    std::cout << "\nTensor4D ampall : " << std::endl;
    for (int i = 0; i < amp2.size(); ++i)
    {
        for (int j = 0; j < amp2[0].size(); ++j)
        {
            for (int k = 0; k < amp2[0][0].size(); ++k)
            {
                for (int z = 0; z < amp2[0][0][0].size(); ++z)
                {
                    std::cout << amp2[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }

    // Print result for comparison
    std::cout << "C++ amplitude: " << amp << std::endl;
    std::cout << "Julia amplitude (expected): 0.00276588 + 0.000205176im" << std::endl;
    std::cout << "C++ intensity: " << intensity_valued << std::endl;
    std::cout << "Julia intensity (expected): 2.5191201498154108" << std::endl;

    // Compare with expected Julia result
    complex expected(0.00276588, 0.000205176);
    expectNearComplex(expected, amp, 0.0001); // Allow for small numerical differences

    auto cg = clebschgordan(2, 1, 0, 1, 2, 1);
    std::cout << "C++ Clebsch-Gordan coefficient: " << cg << std::endl;

    SpinParity j = SpinParity("3/2+");
    ThreeBodySpins spins = {4, 2, 0, 2}; // J/ψ (spin 1), p (spin 1/2), K (spin 0), Λb (spin 1/2)
    std::vector<LSCoupling> two_lsLS = possible_lsLS(j, spins, conservingParities, 1);
    std::cout << "two_lsLS: [" << two_lsLS.size() << "]" << std::endl;
    for (size_t idx = 0; idx < two_lsLS.size(); ++idx)
    {
        const auto &two_ls = two_lsLS[idx].two_ls;
        const auto &two_LS = two_lsLS[idx].two_LS;
        std::cout << "LS coupling " << idx << ": " << two_ls[0]
                  << " " << two_ls[1] << ", L coupling: "
                  << two_LS[0] << " " << two_LS[1] << std::endl;
    }

}

// Add a test case to verify we match the unpolarized intensity from Julia
TEST_F(ThreeBodyAmplitudeModelTest, UnpolarizedIntensityJuliaMatch)
{
    // Create a model with multiple decay chains as in the tutorial
    auto Lambda1520 = createDecayChainLS(
        1, createBW(1.5195, 0.0156), "3/2+", conservingParities, *tbs);
    auto Lambda1690 = createDecayChainLS(
        1, createBW(1.685, 0.050), "1/2+", conservingParities, *tbs);
    auto Pc4312 = createDecayChainLS(
        3, createBW(4.312, 0.015), "1/2+", conservingParities, *tbs);

    // Create the model using values from the tutorial
    ThreeBodyAmplitudeModel model;
    model.add(Lambda1520, "Lambda1520", complex(2.0, 0.0));
    model.add(Lambda1690, "Lambda1690", complex(2.1, 0.0));
    model.add(Pc4312, "Pc4312", complex(0.4, 0.0));

    // Use the standard test point
    MandelstamTuple σs = {4.501751135564419, 21.495750373843414, 16.258552559492166};


    // Calculate amplitude tensor

    // Calculate intensity
    double intensity_value = model.intensity(σs, 1);

    // Compare with the expected value from Julia (10.029598796534886)
    std::cout << "C++ unpolarized intensity: " << intensity_value << std::endl;
    std::cout << "Julia unpolarized intensity: 2.5191201498154108" << std::endl;

    // Allow a small tolerance for numerical differences
    EXPECT_NEAR(intensity_value, 2.5191201498154108, 1e-6);
    ThreeBodySystem tbs2(ms, spins);
    ThreeBodySpins sp = tbs2.get_two_js();
    std::cout << sp[0] << " " << sp[1] << " " << sp[2] << " " << sp[3] << std::endl;
}

TEST_F(ThreeBodyAmplitudeModelTest, Lambda1520aligned_amplitude)
{
    // Test the Clebsch-Gordan coefficient calculation
    ThreeBodyDecays tbd;
    int k = 1;
    MandelstamTuple σs2 = {4.501751135564419, 21.495750373843414, 16.258552559492166};
    ThreeBodySpins spins = {2, 1, 0, 1}; // J/ψ (spin 1), p (spin 1/2), K (spin 0), Λb (spin 1/2)
    auto Lambda1520 = createDecayChainLS(
        1, createBW(1.5195, 0.0156), "3/2+", conservingParities, *tbs);

    Tensor4Dcomp resultdalignedccomp = tbd.aligned_amplitude4dcomp(*Lambda1520, σs2);
    // Verify tensor dimensions

    std::cout << "resultdalignedccomp dimensions: "
              << resultdalignedccomp.size() << " x "
              << (resultdalignedccomp.empty() ? 0 : resultdalignedccomp[0].size()) << " x "
              << (resultdalignedccomp.empty() || resultdalignedccomp[0].empty() ? 0 : resultdalignedccomp[0][0].size()) << " x "
              << (resultdalignedccomp.empty() || resultdalignedccomp[0].empty() || resultdalignedccomp[0][0].empty() ? 0 : resultdalignedccomp[0][0][0].size())
              << std::endl;

    std::cout << "Tensor4D resultdalignedccomp : " << std::endl;
    for (int i = 0; i < resultdalignedccomp.size(); ++i)
    {
        for (int j = 0; j < resultdalignedccomp[0].size(); ++j)
        {
            for (int k = 0; k < resultdalignedccomp[0][0].size(); ++k)
            {
                for (int z = 0; z < resultdalignedccomp[0][0][0].size(); ++z)
                {
                    std::cout << i << j << k << z << resultdalignedccomp[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }

    double prec = 1e-5;

    /*
    juliares dimensions: 3 x 2 x 1 x 2
    Tensor4D juliares :
    0000 Λ1520 aligned amplitude: -0.0016689984410686244 + 1.8041314291331854e-5im
    0001 Λ1520 aligned amplitude: -0.020276773169024702 + 0.00021918512837087678im
    0100 Λ1520 aligned amplitude: -0.03046088324379168 + 0.0003292719481756626im
    0101 Λ1520 aligned amplitude: -0.18420114996408743 + 0.0019911527521853517im
    1000 Λ1520 aligned amplitude: 0.028675687616797618 - 0.00030997458121258183im
    1001 Λ1520 aligned amplitude: 0.2604997644839328 - 0.0028159152268970385im
    1100 Λ1520 aligned amplitude: 0.2604997644839328 - 0.0028159152268970385im
    1101 Λ1520 aligned amplitude: -0.028675687616797618 + 0.00030997458121258183im
    2000 Λ1520 aligned amplitude: -0.18420114996408743 + 0.0019911527521853517im
    2001 Λ1520 aligned amplitude: 0.03046088324379168 - 0.0003292719481756626im
    2100 Λ1520 aligned amplitude: 0.020276773169024702 - 0.00021918512837087678im
    2101 Λ1520 aligned amplitude: -0.0016689984410686244 + 1.8041314291331854e-5im
    */
    // Expected values from the Julia tutorial
    Tensor4Dcomp juliares(3, std::vector<std::vector<std::vector<complex>>>(
                                 2, std::vector<std::vector<complex>>(
                                        1, std::vector<complex>(2, complex(0.0, 0.0)))));
    juliares[0][0][0][0] = complex(-0.0016689984410686244, 1.8041314291331854e-5);
    juliares[0][0][0][1] = complex(-0.020276773169024702, 0.00021918512837087678);
    juliares[0][1][0][0] = complex(-0.03046088324379168, 0.0003292719481756626);
    juliares[0][1][0][1] = complex(-0.18420114996408743, 0.0019911527521853517);

    juliares[1][0][0][0] = complex(0.028675687616797618, -0.00030997458121258183);
    juliares[1][0][0][1] = complex(0.2604997644839328, -0.0028159152268970385);
    juliares[1][1][0][0] = complex(0.2604997644839328, -0.0028159152268970385);
    juliares[1][1][0][1] = complex(-0.028675687616797618, 0.00030997458121258183);

    juliares[2][0][0][0] = complex(-0.18420114996408743, 0.0019911527521853517);
    juliares[2][0][0][1] = complex(0.03046088324379168, -0.0003292719481756626);
    juliares[2][1][0][0] = complex(0.020276773169024702, -0.00021918512837087678);
    juliares[2][1][0][1] = complex(-0.0016689984410686244, 1.8041314291331854e-5);

    // Compare the results with the expected values
    for (int i = 0; i < resultdalignedccomp.size(); ++i)
    {
        for (int j = 0; j < resultdalignedccomp[0].size(); ++j)
        {
            for (int k = 0; k < resultdalignedccomp[0][0].size(); ++k)
            {
                for (int z = 0; z < resultdalignedccomp[0][0][0].size(); ++z)
                {
                    EXPECT_NEAR(resultdalignedccomp[i][j][k][z].real(), juliares[i][j][k][z].real(), prec);
                    EXPECT_NEAR(resultdalignedccomp[i][j][k][z].imag(), juliares[i][j][k][z].imag(), prec);
                }
            }
        }
    }
    // Print result for comparison
}
TEST_F(ThreeBodyAmplitudeModelTest, Lambda1520amplitude)
{
    ThreeBodyDecays tbd;
    int k = 1;
    MandelstamTuple σs2 = {4.501751135564419, 21.495750373843414, 16.258552559492166};
    ThreeBodySpins spins = {2, 1, 0, 1}; // J/ψ (spin 1), p (spin 1/2), K (spin 0), Λb (spin 1/2)
    auto Lambda1520 = createDecayChainLS(
        1, createBW(1.5195, 0.0156), "3/2+", conservingParities, *tbs);

    Tensor4Dcomp resultdccomp = tbd.amplitude4dcomp(*Lambda1520, σs2, 1);
    // Verify tensor dimensions
    std::cout << "Tensor4D resultdccomp : " << std::endl;
    for (int i = 0; i < resultdccomp.size(); ++i)
    {
        for (int j = 0; j < resultdccomp[0].size(); ++j)
        {
            for (int k = 0; k < resultdccomp[0][0].size(); ++k)
            {
                for (int z = 0; z < resultdccomp[0][0][0].size(); ++z)
                {
                    std::cout << i << j << k << z << resultdccomp[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }

    /*
    Tensor4D resultdccomp :
    0000(-0.001669,1.80413e-05)     0001(-0.0202768,0.000219185)    0100(-0.0304609,0.000329272)    0101(-0.184201,0.00199115)
    1000(0.0286757,-0.000309975)    1001(0.2605,-0.00281592)        1100(0.2605,-0.00281592)        1101(-0.0286757,0.000309975)
    2000(-0.184201,0.00199115)      2001(0.0304609,-0.000329272)    2100(0.0202768,-0.000219185)    2101(-0.001669,1.80413e-05) */

    EXPECT_NEAR(resultdccomp[0][0][0][0].real(), -0.001669, 1e-6);
    EXPECT_NEAR(resultdccomp[0][0][0][0].imag(), 1.80413e-05, 1e-8);
    EXPECT_NEAR(resultdccomp[0][0][0][1].real(), -0.0202768, 1e-6);
    EXPECT_NEAR(resultdccomp[0][0][0][1].imag(), 0.000219185, 1e-8);
    EXPECT_NEAR(resultdccomp[0][1][0][0].real(), -0.0304609, 1e-6);
    EXPECT_NEAR(resultdccomp[0][1][0][0].imag(), 0.000329272, 1e-8);
    EXPECT_NEAR(resultdccomp[0][1][0][1].real(), -0.184201, 1e-6);
    EXPECT_NEAR(resultdccomp[0][1][0][1].imag(), 0.00199115, 1e-8);
    EXPECT_NEAR(resultdccomp[1][0][0][0].real(), 0.0286757, 1e-6);
    EXPECT_NEAR(resultdccomp[1][0][0][0].imag(), -0.000309975, 1e-8);
    EXPECT_NEAR(resultdccomp[1][0][0][1].real(), 0.2605, 1e-6);
    EXPECT_NEAR(resultdccomp[1][0][0][1].imag(), -0.00281592, 1e-8);
    EXPECT_NEAR(resultdccomp[1][1][0][0].real(), 0.2605, 1e-6);
    EXPECT_NEAR(resultdccomp[1][1][0][0].imag(), -0.00281592, 1e-8);
    EXPECT_NEAR(resultdccomp[1][1][0][1].real(), -0.0286757, 1e-6);
    EXPECT_NEAR(resultdccomp[1][1][0][1].imag(), 0.000309975, 1e-8);
    EXPECT_NEAR(resultdccomp[2][0][0][0].real(), -0.184201, 1e-6);
    EXPECT_NEAR(resultdccomp[2][0][0][0].imag(), 0.00199115, 1e-8);
    EXPECT_NEAR(resultdccomp[2][0][0][1].real(), 0.0304609, 1e-6);
    EXPECT_NEAR(resultdccomp[2][0][0][1].imag(), -0.000329272, 1e-8);
    EXPECT_NEAR(resultdccomp[2][1][0][0].real(), 0.0202768, 1e-6);
    EXPECT_NEAR(resultdccomp[2][1][0][0].imag(), -0.000219185, 1e-8);
    EXPECT_NEAR(resultdccomp[2][1][0][1].real(), -0.001669, 1e-6);
    EXPECT_NEAR(resultdccomp[2][1][0][1].imag(), 1.8041314291331854e-05, 1e-8);
    // Print result for comparison
}

TEST_F(ThreeBodyAmplitudeModelTest, Lambda1690aligned_amplitude)
{
    // Test the Clebsch-Gordan coefficient calculation
    ThreeBodyDecays tbd;
    int k = 1;
    MandelstamTuple σs2 = {4.501751135564419, 21.495750373843414, 16.258552559492166};
    ThreeBodySpins spins = {2, 1, 0, 1}; // J/ψ (spin 1), p (spin 1/2), K (spin 0), Λb (spin 1/2)
    auto Lambda1690 = createDecayChainLS(
        1, createBW(1.685, 0.050), "1/2+", conservingParities, *tbs);

    Tensor4Dcomp resultdalignedccomp = tbd.aligned_amplitude4dcomp(*Lambda1690, σs2);
    // Verify tensor dimensions

    std::cout << "resultdalignedccomp dimensions: "
              << resultdalignedccomp.size() << " x "
              << (resultdalignedccomp.empty() ? 0 : resultdalignedccomp[0].size()) << " x "
              << (resultdalignedccomp.empty() || resultdalignedccomp[0].empty() ? 0 : resultdalignedccomp[0][0].size()) << " x "
              << (resultdalignedccomp.empty() || resultdalignedccomp[0].empty() || resultdalignedccomp[0][0].empty() ? 0 : resultdalignedccomp[0][0][0].size())
              << std::endl;

    std::cout << "Tensor4D resultdalignedccomp : " << std::endl;
    for (int i = 0; i < resultdalignedccomp.size(); ++i)
    {
        for (int j = 0; j < resultdalignedccomp[0].size(); ++j)
        {
            for (int k = 0; k < resultdalignedccomp[0][0].size(); ++k)
            {
                for (int z = 0; z < resultdalignedccomp[0][0][0].size(); ++z)
                {
                    std::cout << i << j << k << z << resultdalignedccomp[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }

    /*
    0000 Λ1690 aligned amplitude: -0.0 + 0.0im
    0001 Λ1690 aligned amplitude: -0.018950450900958035 + 0.0009603310614203887im
    0100 Λ1690 aligned amplitude: -0.0 + 0.0im
    0101 Λ1690 aligned amplitude: 0.3458645964592334 - 0.017526997999221142im
    1000 Λ1690 aligned amplitude: -0.013399992338610148 + 0.0006790566057144318im
    1001 Λ1690 aligned amplitude: 0.2445632015286727 - 0.01239345913909232im
    1100 Λ1690 aligned amplitude: 0.2445632015286727 - 0.01239345913909232im
    1101 Λ1690 aligned amplitude: 0.013399992338610148 - 0.0006790566057144318im
    2000 Λ1690 aligned amplitude: 0.3458645964592334 - 0.017526997999221142im
    2001 Λ1690 aligned amplitude: -0.0 + 0.0im
    2100 Λ1690 aligned amplitude: 0.018950450900958035 - 0.0009603310614203887im
    2101 Λ1690 aligned amplitude: 0.0 - 0.0im
    */

    double prec = 1e-6;
    EXPECT_NEAR(resultdalignedccomp[0][0][0][0].real(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[0][0][0][0].imag(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[0][0][0][1].real(), -0.018950, prec);
    EXPECT_NEAR(resultdalignedccomp[0][0][0][1].imag(), 0.000960331, prec);
    EXPECT_NEAR(resultdalignedccomp[0][1][0][0].real(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[0][1][0][0].imag(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[0][1][0][1].real(), 0.345865, prec);
    EXPECT_NEAR(resultdalignedccomp[0][1][0][1].imag(), -0.017527, prec);
    EXPECT_NEAR(resultdalignedccomp[1][0][0][0].real(), -0.013400, prec);
    EXPECT_NEAR(resultdalignedccomp[1][0][0][0].imag(), 0.000679056, prec);
    EXPECT_NEAR(resultdalignedccomp[1][0][0][1].real(), 0.244563, prec);
    EXPECT_NEAR(resultdalignedccomp[1][0][0][1].imag(), -0.012393, prec);
    EXPECT_NEAR(resultdalignedccomp[1][1][0][0].real(), 0.244563, prec);
    EXPECT_NEAR(resultdalignedccomp[1][1][0][0].imag(), -0.012393, prec);
    EXPECT_NEAR(resultdalignedccomp[1][1][0][1].real(), 0.013400, prec);
    EXPECT_NEAR(resultdalignedccomp[1][1][0][1].imag(), -0.000679056, prec);
    EXPECT_NEAR(resultdalignedccomp[2][0][0][0].real(), 0.345865, prec);
    EXPECT_NEAR(resultdalignedccomp[2][0][0][0].imag(), -0.017527, prec);
    EXPECT_NEAR(resultdalignedccomp[2][0][0][1].real(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[2][0][0][1].imag(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[2][1][0][0].real(), 0.018950, prec);
    EXPECT_NEAR(resultdalignedccomp[2][1][0][0].imag(), -0.000960331, prec);
    EXPECT_NEAR(resultdalignedccomp[2][1][0][1].real(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[2][1][0][1].imag(), 0.0, prec);
}

TEST_F(ThreeBodyAmplitudeModelTest, Lambda1690amplitude)
{
    // Test the Clebsch-Gordan coefficient calculation
    ThreeBodyDecays tbd;
    int k = 1;
    MandelstamTuple σs2 = {4.501751135564419, 21.495750373843414, 16.258552559492166};
    ThreeBodySpins spins = {2, 1, 0, 1}; // J/ψ (spin 1), p (spin 1/2), K (spin 0), Λb (spin 1/2)
    auto Lambda1690 = createDecayChainLS(
        1, createBW(1.685, 0.050), "1/2+", conservingParities, *tbs);

    Tensor4Dcomp resultdalignedccomp = tbd.amplitude4dcomp(*Lambda1690, σs2, 1);
    // Verify tensor dimensions

    std::cout << "amplitude dimensions: "
              << resultdalignedccomp.size() << " x "
              << (resultdalignedccomp.empty() ? 0 : resultdalignedccomp[0].size()) << " x "
              << (resultdalignedccomp.empty() || resultdalignedccomp[0].empty() ? 0 : resultdalignedccomp[0][0].size()) << " x "
              << (resultdalignedccomp.empty() || resultdalignedccomp[0].empty() || resultdalignedccomp[0][0].empty() ? 0 : resultdalignedccomp[0][0][0].size())
              << std::endl;

    std::cout << "Tensor4D amplitude4dcomp : " << std::endl;
    for (int i = 0; i < resultdalignedccomp.size(); ++i)
    {
        for (int j = 0; j < resultdalignedccomp[0].size(); ++j)
        {
            for (int k = 0; k < resultdalignedccomp[0][0].size(); ++k)
            {
                for (int z = 0; z < resultdalignedccomp[0][0][0].size(); ++z)
                {
                    std::cout << i << j << k << z << resultdalignedccomp[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }

    /*
    0000(0,0)       0001(-0.0189505,0.000960331)    0100(0,0)       0101(0.345865,-0.017527)
1000(-0.0134,0.000679057)       1001(0.244563,-0.0123935)       1100(0.244563,-0.0123935)       1101(0.0134,-0.000679057)
2000(0.345865,-0.017527)        2001(0,0)       2100(0.0189505,-0.000960331)    2101(0,0)
    */

    double prec = 1e-6;
    EXPECT_NEAR(resultdalignedccomp[0][0][0][0].real(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[0][0][0][0].imag(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[0][0][0][1].real(), -0.0189505, prec);
    EXPECT_NEAR(resultdalignedccomp[0][0][0][1].imag(), 0.000960331, prec);
    EXPECT_NEAR(resultdalignedccomp[0][1][0][0].real(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[0][1][0][0].imag(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[0][1][0][1].real(), 0.345865, prec);
    EXPECT_NEAR(resultdalignedccomp[0][1][0][1].imag(), -0.017527, prec);
    EXPECT_NEAR(resultdalignedccomp[1][0][0][0].real(), -0.0134, prec);
    EXPECT_NEAR(resultdalignedccomp[1][0][0][0].imag(), 0.000679057, prec);
    EXPECT_NEAR(resultdalignedccomp[1][0][0][1].real(), 0.244563, prec);
    EXPECT_NEAR(resultdalignedccomp[1][0][0][1].imag(), -0.0123935, prec);
    EXPECT_NEAR(resultdalignedccomp[1][1][0][0].real(), 0.244563, prec);
    EXPECT_NEAR(resultdalignedccomp[1][1][0][0].imag(), -0.0123935, prec);
    EXPECT_NEAR(resultdalignedccomp[1][1][0][1].real(), 0.0134, prec);
    EXPECT_NEAR(resultdalignedccomp[1][1][0][1].imag(), -0.000679057, prec);
    EXPECT_NEAR(resultdalignedccomp[2][0][0][0].real(), 0.345865, prec);
    EXPECT_NEAR(resultdalignedccomp[2][0][0][0].imag(), -0.017527, prec);
    EXPECT_NEAR(resultdalignedccomp[2][0][0][1].real(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[2][0][0][1].imag(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[2][1][0][0].real(), 0.0189505, prec);
    EXPECT_NEAR(resultdalignedccomp[2][1][0][0].imag(), -0.000960331, prec);
    EXPECT_NEAR(resultdalignedccomp[2][1][0][1].real(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[2][1][0][1].imag(), 0.0, prec);
}

TEST_F(ThreeBodyAmplitudeModelTest, Pc4312aligned_amplitude)
{
    // Test the Clebsch-Gordan coefficient calculation
    ThreeBodyDecays tbd;
    int k = 1;
    MandelstamTuple σs2 = {4.501751135564419, 21.495750373843414, 16.258552559492166};
    ThreeBodySpins spins = {2, 1, 0, 1}; // J/ψ (spin 1), p (spin 1/2), K (spin 0), Λb (spin 1/2)
    auto Pc4312 = createDecayChainLS(
        3, createBW(4.312, 0.015), "1/2+", conservingParities, *tbs);

    Tensor4Dcomp resultdalignedccomp = tbd.aligned_amplitude4dcomp(*Pc4312, σs2);
    // Verify tensor dimensions

    std::cout << "resultdalignedccomp dimensions: "
              << resultdalignedccomp.size() << " x "
              << (resultdalignedccomp.empty() ? 0 : resultdalignedccomp[0].size()) << " x "
              << (resultdalignedccomp.empty() || resultdalignedccomp[0].empty() ? 0 : resultdalignedccomp[0][0].size()) << " x "
              << (resultdalignedccomp.empty() || resultdalignedccomp[0].empty() || resultdalignedccomp[0][0].empty() ? 0 : resultdalignedccomp[0][0][0].size())
              << std::endl;

    std::cout << "Tensor4D resultdalignedccomp : " << std::endl;
    for (int i = 0; i < resultdalignedccomp.size(); ++i)
    {
        for (int j = 0; j < resultdalignedccomp[0].size(); ++j)
        {
            for (int k = 0; k < resultdalignedccomp[0][0].size(); ++k)
            {
                for (int z = 0; z < resultdalignedccomp[0][0][0].size(); ++z)
                {
                    std::cout << i << j << k << z << resultdalignedccomp[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }

    /*
    0000 Pc4312 aligned amplitude: 0.21454237236042922 + 0.005943400512576101im
    0001 Pc4312 aligned amplitude: 0.1225800405096077 + 0.0033957966791401807im
    0100 Pc4312 aligned amplitude: 0.0 + 0.0im
    0101 Pc4312 aligned amplitude: -0.0 + 0.0im
    1000 Pc4312 aligned amplitude: -0.0866771778824653 - 0.0024011908593507803im
    1001 Pc4312 aligned amplitude: 0.15170436634790885 + 0.004202618805750164im
    1100 Pc4312 aligned amplitude: 0.15170436634790885 + 0.004202618805750164im
    1101 Pc4312 aligned amplitude: 0.0866771778824653 + 0.0024011908593507803im
    2000 Pc4312 aligned amplitude: -0.0 + 0.0im
    2001 Pc4312 aligned amplitude: 0.0 + 0.0im
    2100 Pc4312 aligned amplitude: -0.1225800405096077 - 0.0033957966791401807im
    2101 Pc4312 aligned amplitude: 0.21454237236042922 + 0.005943400512576101im
        */

    double prec = 1e-6;
    EXPECT_NEAR(resultdalignedccomp[0][0][0][0].real(), 0.21454237236042922, prec);
    EXPECT_NEAR(resultdalignedccomp[0][0][0][0].imag(), 0.005943400512576101, prec);
    EXPECT_NEAR(resultdalignedccomp[0][0][0][1].real(), 0.1225800405096077, prec);
    EXPECT_NEAR(resultdalignedccomp[0][0][0][1].imag(), 0.0033957966791401807, prec);
    EXPECT_NEAR(resultdalignedccomp[0][1][0][0].real(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[0][1][0][0].imag(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[0][1][0][1].real(), -0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[0][1][0][1].imag(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[1][0][0][0].real(), -0.0866771778824653, prec);
    EXPECT_NEAR(resultdalignedccomp[1][0][0][0].imag(), -0.0024011908593507803, prec);
    EXPECT_NEAR(resultdalignedccomp[1][0][0][1].real(), 0.15170436634790885, prec);
    EXPECT_NEAR(resultdalignedccomp[1][0][0][1].imag(), 0.004202618805750164, prec);
    EXPECT_NEAR(resultdalignedccomp[1][1][0][0].real(), 0.15170436634790885, prec);
    EXPECT_NEAR(resultdalignedccomp[1][1][0][0].imag(), 0.004202618805750164, prec);
    EXPECT_NEAR(resultdalignedccomp[1][1][0][1].real(), 0.0866771778824653, prec);
    EXPECT_NEAR(resultdalignedccomp[1][1][0][1].imag(), 0.0024011908593507803, prec);
    EXPECT_NEAR(resultdalignedccomp[2][0][0][0].real(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[2][0][0][0].imag(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[2][0][0][1].real(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[2][0][0][1].imag(), 0.0, prec);
    EXPECT_NEAR(resultdalignedccomp[2][1][0][0].real(), -0.1225800405096077, prec);
    EXPECT_NEAR(resultdalignedccomp[2][1][0][0].imag(), -0.0033957966791401807, prec);
    EXPECT_NEAR(resultdalignedccomp[2][1][0][1].real(), 0.21454237236042922, prec);
    EXPECT_NEAR(resultdalignedccomp[2][1][0][1].imag(), 0.005943400512576101, prec);
}

TEST_F(ThreeBodyAmplitudeModelTest, Pc4312amplitude)
{
    // Test the Clebsch-Gordan coefficient calculation
    ThreeBodyDecays tbd;
    int k = 1;
    MandelstamTuple σs2 = {4.501751135564419, 21.495750373843414, 16.258552559492166};
    ThreeBodySpins spins = {2, 1, 0, 1}; // J/ψ (spin 1), p (spin 1/2), K (spin 0), Λb (spin 1/2)
    auto Pc4312 = createDecayChainLS(
        3, createBW(4.312, 0.015), "1/2+", conservingParities, *tbs);

    Tensor4Dcomp resultdalignedccomp = tbd.amplitude4dcomp(*Pc4312, σs2, 1);
    // Verify tensor dimensions

    std::cout << "amplitude dimensions: "
              << resultdalignedccomp.size() << " x "
              << (resultdalignedccomp.empty() ? 0 : resultdalignedccomp[0].size()) << " x "
              << (resultdalignedccomp.empty() || resultdalignedccomp[0].empty() ? 0 : resultdalignedccomp[0][0].size()) << " x "
              << (resultdalignedccomp.empty() || resultdalignedccomp[0].empty() || resultdalignedccomp[0][0].empty() ? 0 : resultdalignedccomp[0][0][0].size())
              << std::endl;

    std::cout << "Tensor4D amplitude4dcomp : " << std::endl;
    for (int i = 0; i < resultdalignedccomp.size(); ++i)
    {
        for (int j = 0; j < resultdalignedccomp[0].size(); ++j)
        {
            for (int k = 0; k < resultdalignedccomp[0][0].size(); ++k)
            {
                for (int z = 0; z < resultdalignedccomp[0][0][0].size(); ++z)
                {
                    std::cout << i << j << k << z << resultdalignedccomp[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }

    /*
    0000(0.0152597,0.000422735)     0001(-0.00902943,-0.00025014)   0100(0.212104,0.00587586)       0101(-0.125506,-0.00347685)
1000(-0.143596,-0.00397798)     1001(0.0995363,0.00275742)      1100(0.0995363,0.00275742)      1101(0.143596,0.00397798)
2000(-0.125506,-0.00347685)     2001(-0.212104,-0.00587586)     2100(0.00902943,0.00025014)     2101(0.0152597,0.000422735)
    */

    double prec = 1e-6;
    EXPECT_NEAR(resultdalignedccomp[0][0][0][0].real(), 0.0152597, prec);
    EXPECT_NEAR(resultdalignedccomp[0][0][0][0].imag(), 0.000422735, prec);
    EXPECT_NEAR(resultdalignedccomp[0][0][0][1].real(), -0.00902943, prec);
    EXPECT_NEAR(resultdalignedccomp[0][0][0][1].imag(), -0.00025014, prec);
    EXPECT_NEAR(resultdalignedccomp[0][1][0][0].real(), 0.212104, prec);
    EXPECT_NEAR(resultdalignedccomp[0][1][0][0].imag(), 0.00587586, prec);
    EXPECT_NEAR(resultdalignedccomp[0][1][0][1].real(), -0.125506, prec);
    EXPECT_NEAR(resultdalignedccomp[0][1][0][1].imag(), -0.00347685, prec);
    EXPECT_NEAR(resultdalignedccomp[1][0][0][0].real(), -0.143596, prec);
    EXPECT_NEAR(resultdalignedccomp[1][0][0][0].imag(), -0.00397798, prec);
    EXPECT_NEAR(resultdalignedccomp[1][0][0][1].real(), 0.0995363, prec);
    EXPECT_NEAR(resultdalignedccomp[1][0][0][1].imag(), 0.00275742, prec);
    EXPECT_NEAR(resultdalignedccomp[1][1][0][0].real(), 0.0995363, prec);
    EXPECT_NEAR(resultdalignedccomp[1][1][0][0].imag(), 0.00275742, prec);
    EXPECT_NEAR(resultdalignedccomp[1][1][0][1].real(), 0.143596, prec);
    EXPECT_NEAR(resultdalignedccomp[1][1][0][1].imag(), 0.00397798, prec);
    EXPECT_NEAR(resultdalignedccomp[2][0][0][0].real(), -0.125506, prec);
    EXPECT_NEAR(resultdalignedccomp[2][0][0][0].imag(), -0.00347685, prec);
    EXPECT_NEAR(resultdalignedccomp[2][0][0][1].real(), -0.212104, prec);
    EXPECT_NEAR(resultdalignedccomp[2][0][0][1].imag(), -0.00587586, prec);
    EXPECT_NEAR(resultdalignedccomp[2][1][0][0].real(), 0.00902943, prec);
    EXPECT_NEAR(resultdalignedccomp[2][1][0][0].imag(), 0.00025014, prec);
    EXPECT_NEAR(resultdalignedccomp[2][1][0][1].real(), 0.0152597, prec);
    EXPECT_NEAR(resultdalignedccomp[2][1][0][1].imag(), 0.000422735, prec);
}

TEST_F(ThreeBodyAmplitudeModelTest, Lambda1520lambdacpkpi)
{
    ThreeBodyDecays tbd;
    int k = 1;
    MandelstamTuple σs2 = {4.501751135564419, 21.495750373843414, 16.258552559492166};
    ThreeBodySpins spins = {1, 0, 0, 1}; // p (spin 1/2), K (spin 0), pi (spin 0), Λc (spin 1/2)
    ThreeBodyMasses ms = {0.938, 0.494, 0.140, 2.286};
    ThreeBodySystem tbslc = {ms, spins};
    ThreeBodyParities parities('-', '-', '-', '-');

    auto Lambda1520 = createDecayChainLS(
        1, createBW(1.5195, 0.0156), "3/2+", conservingParities, tbslc);

    Tensor4Dcomp resultdccomp = tbd.amplitude4dcomp(*Lambda1520, σs2, 1);
    // Verify tensor dimensions
    std::cout << "Tensor4D resultdccomp : " << std::endl;
    for (int i = 0; i < resultdccomp.size(); ++i)
    {
        for (int j = 0; j < resultdccomp[0].size(); ++j)
        {
            for (int k = 0; k < resultdccomp[0][0].size(); ++k)
            {
                for (int z = 0; z < resultdccomp[0][0][0].size(); ++z)
                {
                    std::cout << i << j << k << z << resultdccomp[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }
}

#ifndef RUNNING_COMBINED_TESTS
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
#endif
