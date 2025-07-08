#include <gtest/gtest.h>
#include "include/WignerFunctions.hh"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// Test fixture for Wigner matrix tests
class WignerFunctionsTest : public ::testing::Test
{
protected:
    const double epsilon = 1e-10;

    // Helper function to compare matrices with tolerance
    void expectMatrixNear(const Matrix2D &expected, const Matrix2D &actual, double tol = 1e-10)
    {
        ASSERT_EQ(expected.size(), actual.size()) << "Matrices have different row counts";

        for (size_t i = 0; i < expected.size(); ++i)
        {
            ASSERT_EQ(expected[i].size(), actual[i].size())
                << "Row " << i << " has different column counts";

            for (size_t j = 0; j < expected[i].size(); ++j)
            {
                EXPECT_NEAR(expected[i][j], actual[i][j], tol)
                    << "Difference at position (" << i << ", " << j << ")";
            }
        }
    }

    // Helper to print matrices for debugging
    void printMatrix(const Matrix2D &matrix)
    {
        for (const auto &row : matrix)
        {
            for (double val : row)
            {
                std::cout << std::setw(12) << std::setprecision(6) << val;
            }
            std::cout << std::endl;
        }
    }
};

// Test wignerd_doublearg_sign function
TEST_F(WignerFunctionsTest, WignerD_DoubleArg_Sign)
{
    int two_j = 4;         // j = 2
    double cosTheta = 0.7; // cos(θ)
    auto res = wignerd_doublearg_sign(two_j, cosTheta, true);

    // Print for debugging
    std::cout << "Wigner d-matrix for j=2, cosθ=0.7, positive sign:" << std::endl;
    printMatrix(res);

    /*
        0.722500     0.607021     0.312310     0.107121     0.022500
   -0.607021     0.340000     0.612250     0.360000     0.107121
    0.312310    -0.612250     0.235000     0.612250     0.312310
   -0.107121     0.360000    -0.612250     0.340000     0.607021
    0.022500    -0.107121     0.312310    -0.607021     0.722500
    */

    // Expected values based on analytical calculations
    Matrix2D expected = {
        {0.722500, 0.607021, 0.312310, 0.107121, 0.022500},
        {-0.607021, 0.340000, 0.612250, 0.360000, 0.107121},
        {0.312310, -0.612250, 0.235000, 0.612250, 0.312310},
        {-0.107121, 0.360000, -0.612250, 0.340000, 0.607021},
        {0.022500, -0.107121, 0.312310, -0.607021, 0.722500}};

    // Compare the result with expected values
    expectMatrixNear(expected, res, 1e-6);
}

// Test specific known values for Wigner d-functions
TEST_F(WignerFunctionsTest, WignerD_KnownValues)
{
    // Test d^2_00(θ) = (3cos²θ - 1)/2
    double cosTheta = 0.5;
    double expected = (3 * cosTheta * cosTheta - 1) / 2;
    double actual = wignerd_doublearg(4, 0, 0, cosTheta);
    EXPECT_NEAR(expected, actual, epsilon);

    // Test d^1_00(θ) = cosθ
    actual = wignerd_doublearg(2, 0, 0, cosTheta);
    EXPECT_NEAR(cosTheta, actual, epsilon);

    // Test d^1_11(θ) = (1 + cosθ)/2
    actual = wignerd_doublearg(2, 2, 2, cosTheta);
    EXPECT_NEAR((1 + cosTheta) / 2, actual, epsilon);
}

// Test symmetry properties of Wigner d-functions
TEST_F(WignerFunctionsTest, WignerD_Symmetry)
{
    double cosTheta = 0.6;

    // Test d^j_mn(θ) = (-1)^(m-n) d^j_-m-n(θ)
    int two_j = 4;  // j = 2
    int two_m = 2;  // m = 1
    int two_n = -2; // n = -1

    double val1 = wignerd_doublearg(two_j, two_m, two_n, cosTheta);
    double val2 = wignerd_doublearg(two_j, -two_m, -two_n, cosTheta);

    int sign = (std::abs(two_m - two_n) / 2) % 2 == 0 ? 1 : -1;
    EXPECT_NEAR(sign * val2, val1, epsilon);
}

// Test Wigner matrix at special angles
TEST_F(WignerFunctionsTest, WignerD_SpecialAngles)
{
    // At θ = 0 (cosθ = 1), d^j_mn(0) = δ_mn
    auto matrix_at_0 = wignerd_doublearg_sign(4, 1.0, true);

    // Should be identity matrix
    for (size_t i = 0; i < matrix_at_0.size(); ++i)
    {
        for (size_t j = 0; j < matrix_at_0[i].size(); ++j)
        {
            if (i == j)
            {
                EXPECT_NEAR(matrix_at_0[i][j], 1.0, epsilon);
            }
            else
            {
                EXPECT_NEAR(matrix_at_0[i][j], 0.0, epsilon);
            }
        }
    }
}

#ifndef RUNNING_COMBINED_TESTS
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
#endif
