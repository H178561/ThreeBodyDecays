#include "include/FormFactors.hh"
#include <gtest/gtest.h>

#include <iostream>
#include <iomanip>
#include <vector>

using namespace FormFactors;

// Test fixture for Clebsch-Gordan coefficient tests
class BWTest : public ::testing::Test
{
protected:
    // Tolerance for floating point comparisons
    const double epsilon = 1e-10;

    // Helper function to compare complex values
    bool complex_approx_equal(complex a, complex b)
    {
        return std::abs(a - b) < epsilon;
    }

    // Test parameters for K*(892)
    const double K892_MASS = 0.8955;
    const double K892_WIDTH = 0.047299999999999995;
    const double K892_MA = 0.13957018; // Pion mass
    const double K892_MB = 0.493677;   // Kaon mass
    const int K892_L = 1;              // Angular momentum
    const double K892_D = 1.5;         // Blatt-Weisskopf parameter
};

// Test the multichannel Breit-Wigner implementation
TEST_F(BWTest, testMultichannelBW)
{
    std::cout << "\n=== Multichannel Breit-Wigner Test ===" << std::endl;

    // Create a single-channel MultichannelBreitWignerMulti using factory function
    auto mbw = make_multichannel_bw_single(K892_MASS, K892_WIDTH, K892_MA, K892_MB, K892_L, K892_D);

    // Test at different invariant masses
    std::vector<double> test_masses = {
        K892_MASS - 0.1, K892_MASS - 0.05, K892_MASS, K892_MASS + 0.05, K892_MASS + 0.1};

    double s_test = 1.5;
    std::cout << "\nInvariant Mass (GeV)  |  Magnitude  |  Real Part  |  Imaginary Part" << std::endl;
    std::cout << "--------------------|------------|-------------|----------------" << std::endl;
    double msub = sqrt(s_test);
    complex val = mbw(s_test);

    std::cout << std::fixed << std::setprecision(5)
              << std::setw(20) << s_test
              << " | " << std::setw(10) << std::abs(val)
              << " | " << std::setw(11) << val.real()
              << " | " << std::setw(14) << val.imag() << std::endl;

    double q = breakup(msub, K892_MA, K892_MB);

    // Impuls bei Resonanzmasse
    // double q0 = breakup(mass, m1, m2);

    double formFactor = BlattWeisskopf(q, K892_L, K892_D);
    val *= formFactor; // Simulate some variation in the value

    // Jula values: -0.8363631758621111 + 0.14486178164234825im

    std::cout << "After times " << formFactor << ": " << std::endl;
    std::cout << std::fixed << std::setprecision(5)
              << "Magnitude: " << std::abs(val) << ", "
              << "Real Part: " << val.real() << ", "
              << "Imaginary Part: " << val.imag() << std::endl;

    EXPECT_NEAR(val.real(), -0.83636, 1e-5);
    EXPECT_NEAR(val.imag(), 0.14486, 1e-5);
}

TEST_F(BWTest, testFlatteLineshape)
{
    /*

    "functions": [
    {
      "name": "L1405_Flatte",
      "type": "MultichannelBreitWigner",
      "x": "m_31_sq",
      "mass": 1.4051,
      "channels": [
        {
          "gsq": 0.23395150538434703,
          "ma": 0.938272046,
          "mb": 0.493677,
          "l": 0,
          "d": 0
        },
        {
          "gsq": 0.23395150538434703,
          "ma": 1.18937,
          "mb": 0.13957018,
          "l": 0,
          "d": 0
        }
      ]
    },
    */
    std::cout << "\n=== Flatte Lineshape Test ===" << std::endl;

    // Create a Flatte lineshape function
    auto flatte = make_flatte(
        1.4051, // mass
        0.328725260215546, 0.938272046, 0.493677, 0, 0, // first channel
        0.328725260215546, 1.18937, 0.13957018, 0, 0 // second channel
    );

    // Test at different invariant masses
    double test_mass = 3.2;

    std::cout << "\nInvariant Mass (GeV)  |  Value" << std::endl;
    std::cout << "--------------------|------------" << std::endl;

    // L1405:2 Lineshape(s=3.2): -0.7480602620708737 + 0.22521439829964177im
    complex val = flatte(test_mass);

    EXPECT_NEAR(val.real(), -0.74806, 1e-5);
    EXPECT_NEAR(val.imag(), 0.22521, 1e-5);
}

#ifndef RUNNING_COMBINED_TESTS
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
#endif
