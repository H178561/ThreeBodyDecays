#include <gtest/gtest.h>
#include "../ThreeBodyAmplitudeModel.hh"
#include "../ThreeBodyDecays.hh"
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
        "3/2+",                   // jp (Spin-parity)
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
    double intensity = model.intensity(σs);

    // The exact value depends on implementation details, but should be positive
    EXPECT_GT(intensity, 0.0);

    // Calculate individual component intensities
    auto componentIntensities = model.component_intensities(σs);
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

// Test case for amplitude calculation with specific helicity values
TEST_F(ThreeBodyAmplitudeModelTest, SpecificHelicityAmplitude)
{
    // Create simplified model with just two chains
    auto Lambda1520 = createDecayChainLS(
        1, createBW(1.5195, 0.0156), "3/2+", conservingParities, *tbs);
    auto Pc4312 = createDecayChainLS(
        3, createBW(4.312, 0.015), "1/2+", conservingParities, *tbs);

    ThreeBodyAmplitudeModel model;
    model.add(Lambda1520, "Lambda1520", complex(1.0, 0.0));
    model.add(Pc4312, "Pc4312", complex(0.5, 0.0));

    // Calculate amplitude for specific helicity configuration
    MandelstamTuple σs = decays->x2σs({0.3, 0.3}, ms, 1);
    std::vector<int> two_λs = {2, 1, 0, 1}; // For 1, 1/2, 0, 1/2 helicities
    std::vector<int> refζs = {1, 2, 3, 1};  // Reference frames

    ThreeBodyDecays tbd;
    complex amp = model.amplitude(σs, two_λs, refζs);
    auto amp4d = tbd.amplitude4d(*Lambda1520, σs, refζs);

    std::cout << "Calculated amplitude: " << amp << std::endl;

    std::cout << "Tensor4D resultdc3 : " << std::endl;
    for (int i = 0; i < amp4d.size(); ++i)
    {
        for (int j = 0; j < amp4d[0].size(); ++j)
        {
            for (int k = 0; k < amp4d[0][0].size(); ++k)
            {
                for (int z = 0; z < amp4d[0][0][0].size(); ++z)
                {
                    std::cout << amp4d[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }

    // The amplitude should be a complex number
    // We can't predict exact values without detailed physics, but we can check it's not zero
    EXPECT_NE(std::abs(amp), 0.0);

    // Add more specific tests if possible
    std::cout << "Specific helicity amplitude: " << amp << std::endl;
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
    std::vector<int> refζs = {1, 2, 3, 1};  // Reference frames



    // Calculate amplitude
    complex amp = model.amplitude(σs, two_λs, refζs);
    Tensor4D amp2 = model.amplitude4d(σs, refζs);



    double intensity_valued = model.intensity(σs);

    ThreeBodyDecays tbd;
    Tensor4D resultdc3 = tbd.amplitude4d(*Lambda1520, σs, refζs);
    complex res = tbd.amplitude(*Lambda1520, σs, two_λs, refζs);
    std::cout << "Tensor4D resultdc3 : " << res << std::endl;
    // Verify tensor dimensions
    std::cout << "Tensor4D resultdc3 : " << std::endl;
        for (int i = 0; i < resultdc3.size(); ++i) {
            for (int j = 0; j < resultdc3[0].size(); ++j) {
                for (int k = 0; k < resultdc3[0][0].size(); ++k) {
                    for (int z = 0; z < resultdc3[0][0][0].size(); ++z) {
                        std::cout << resultdc3[i][j][k][z] << "\t";  // Tab für schöne Ausrichtung
                    }
                }
            }
            std::cout << "\n";
        }

    std::cout << "Tensor4D amp2 : " << std::endl;
    for (int i = 0; i < amp2.size(); ++i) {
        for (int j = 0; j < amp2[0].size(); ++j) {
            for (int k = 0; k < amp2[0][0].size(); ++k) {
                for (int z = 0; z < amp2[0][0][0].size(); ++z) {
                    std::cout << amp2[i][j][k][z] << "\t";  // Tab für schöne Ausrichtung
                }
            }
        }
        std::cout << "\n";
    }

    // Print result for comparison
    std::cout << "C++ amplitude: " << amp << std::endl;
    std::cout << "Julia amplitude (expected): 0.00498673748367451 + 0.0004851614900810134im" << std::endl;
    std::cout << "C++ intensity: " << intensity_valued << std::endl;
    std::cout << "Julia intensity (expected): 2.5191201498154108" << std::endl;



    // Compare with expected Julia result
    complex expected(0.00498673748367451, 0.0004851614900810134);
    expectNearComplex(expected, amp, 0.0001); // Allow for small numerical differences
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
    model.add(Lambda1520, "Lambda1520", complex(1.0, 0.0));
    model.add(Lambda1690, "Lambda1690", complex(0.5, 0.0));
    model.add(Pc4312, "Pc4312", complex(0.0, 1.0));

    // Use the standard test point
    MandelstamTuple σs = decays->x2σs({0.3, 0.3}, ms, 1);

    std::vector<int> refζs = {1, 2, 3, 1};

    // Calculate amplitude tensor


    // Calculate intensity
    double intensity_value = model.intensity(σs);

    // Compare with the expected value from Julia (10.029598796534886)
    std::cout << "C++ unpolarized intensity: " << intensity_value << std::endl;
    std::cout << "Julia unpolarized intensity: 10.029598796534886" << std::endl;

    // Allow a small tolerance for numerical differences
    EXPECT_NEAR(intensity_value, 10.029598796534886, 0.1);
}

/*
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}*/
