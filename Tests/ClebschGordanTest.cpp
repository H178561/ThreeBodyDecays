#include <gtest/gtest.h>
#include "include/ClebschGordan.hh"
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>

// Test fixture for Clebsch-Gordan coefficient tests
class ClebschGordanTest : public ::testing::Test {
protected:
    // Tolerance for floating point comparisons
    const double epsilon = 1e-10;

    // Helper function to compare complex values
    bool complex_approx_equal(complex a, complex b) {
        return std::abs(a - b) < epsilon ;
    }
};

// Test integer angular momenta cases
TEST_F(ClebschGordanTest, IntegerAngularMomenta) {
    // j₁=1, j₂=1, j=2, m₁=1, m₂=1, m=2 should be -1.0
    EXPECT_NEAR(clebschgordan(1, 1, 1, 1, 2, 2), -1.0, epsilon);

    // j₁=1, j₂=1, j=0, m₁=1, m₂=-1, m=0 should be -0.57735026918962584
    EXPECT_NEAR(clebschgordan(1, 1, 1, -1, 0, 0), -0.57735026918962584, epsilon);

    // j₁=2, j₂=1, j=3, m₁=2, m₂=1, m=3 should be -1.0
    EXPECT_NEAR(clebschgordan(2, 2, 1, 1, 3, 3), -1.0, epsilon);

    // j₁=2, j₂=1, j=1, m₁=0, m₂=0, m=0 should be 0.63245553203367599
    EXPECT_NEAR(clebschgordan(2, 0, 1, 0, 1, 0), 0.63245553203367599, epsilon);
}

// Test half-integer angular momenta cases
TEST_F(ClebschGordanTest, HalfIntegerAngularMomenta) {
    // j₁=3/2, j₂=1/2, j=2, m₁=3/2, m₂=1/2, m=2 should be -1.0
    // Using doubled representation: two_j₁=3, two_j₂=1, two_j=4
    EXPECT_NEAR(clebschgordan_doublearg(3, 3, 1, 1, 4, 4), -1.0, epsilon);

    // j₁=1/2, j₂=1/2, j=0, m₁=1/2, m₂=-1/2, m=0 should be -0.70710678118654746
    // Using doubled representation: two_j₁=1, two_j₂=1, two_j=0
    EXPECT_NEAR(clebschgordan_doublearg(1, 1, 1, -1, 0, 0), -0.70710678118654746, epsilon);

    // j₁=3/2, j₂=1/2, j=1, m₁=1/2, m₂=1/2, m=1 should be 0.5
    // Using doubled representation: two_j₁=3, two_j₂=1, two_j=2
    EXPECT_NEAR(clebschgordan_doublearg(3, 1, 1, 1, 2, 2), 0.5, epsilon);
}

// Test symmetry properties of CG coefficients
TEST_F(ClebschGordanTest, SymmetryProperties) {
    // Test symmetry under exchange of particles:
    // <j₁m₁;j₂m₂|JM> = (-1)^(j₁+j₂-J) <j₂m₂;j₁m₁|JM>
    int j1 = 2, m1 = 1, j2 = 1, m2 = 0, J = 2, M = 1;
    complex direct = clebschgordan(j1, m1, j2, m2, J, M);
    complex exchanged = clebschgordan(j2, m2, j1, m1, J, M);
    complex expected = std::pow(-1.0, j1 + j2 - J) * exchanged;
    EXPECT_TRUE(complex_approx_equal(direct, expected));

    // Test that <j₁m₁;j₂m₂|JM> = 0 if m₁+m₂≠M
    EXPECT_NEAR(clebschgordan(2, 1, 1, 1, 3, 1), 0.0, epsilon);

    // Test that <j₁m₁;j₂m₂|JM> = 0 if J is outside triangle inequality
    EXPECT_NEAR(clebschgordan(2, 0, 1, 0, 4, 0), 0.0, epsilon);
    EXPECT_NEAR(clebschgordan(2, 0, 1, 0, 0, 0), 0.0, epsilon);
}

// Test orthogonality relations
TEST_F(ClebschGordanTest, OrthogonalityRelations) {
    // Test orthogonality relation:
    // ∑_{m₁,m₂} <j₁m₁;j₂m₂|JM><j₁m₁;j₂m₂|J'M'> = δ_{JJ'}δ_{MM'}

    int j1 = 1, j2 = 1;
    int J = 2, Jprime = 2;
    int M = 1, Mprime = 1;

    double sum = 0.0;
    for (int m1 = -j1; m1 <= j1; m1++) {
        for (int m2 = -j2; m2 <= j2; m2++) {
            if (m1 + m2 == M && m1 + m2 == Mprime) {
                double cg1 = clebschgordan(j1, m1, j2, m2, J, M);
                double cg2 = clebschgordan(j1, m1, j2, m2, Jprime, Mprime);
                sum += (cg1 * cg2);
            }
        }
    }

    // Should be 1.0 since J=J' and M=M'
    EXPECT_NEAR(sum, 1.0, epsilon);

    // Change J' to test orthogonality
    Jprime = 1;
    sum = 0.0;
    for (int m1 = -j1; m1 <= j1; m1++) {
        for (int m2 = -j2; m2 <= j2; m2++) {
            if (m1 + m2 == M && m1 + m2 == Mprime) {
                double cg1 = clebschgordan(j1, m1, j2, m2, J, M);
                double cg2 = clebschgordan(j1, m1, j2, m2, Jprime, Mprime);
                sum += (cg1 * cg2);
            }
        }
    }

    // Should be 0.0 since J≠J'
    EXPECT_NEAR(sum, 0.0, epsilon);
}

// Compare CG implementations
TEST_F(ClebschGordanTest, CompareImplementations) {
    // Test cases with different angular momenta values
    std::vector<std::tuple<int,int,int,int,int,int>> testCases = {
        {1, 0, 1, 0, 0, 0},
        {1, 1, 1, -1, 0, 0},
        {2, 1, 1, 0, 1, 1},
        {2, 2, 2, -2, 0, 0},
        {3, 1, 2, 1, 3, 2}
    };

    for (const auto& testCase : testCases) {
        int j1, m1, j2, m2, j, m;
        std::tie(j1, m1, j2, m2, j, m) = testCase;

        double resultNormal = clebschgordan(j1, m1, j2, m2, j, m);
        double resultDouble = clebschgordan_doublearg(2*j1, 2*m1, 2*j2, 2*m2, 2*j, 2*m);
        double resultDouble2 = CG_doublearg(2*j1, 2*m1, 2*j2, 2*m2, 2*j, 2*m);
        std::cout << "(" << j1 << "," << m1 << ","
                 << j2 << "," << m2 << ","
                 << j << "," << m << ") = "
                 << resultNormal << std::endl;

        EXPECT_TRUE(complex_approx_equal(resultNormal, resultDouble));

        std::cout << "CG_doublearg(" << j1 << "," << m1 << ","
                 << j2 << "," << m2 << ","
                 << j << "," << m << ") = "
                 << resultDouble << std::endl;
        std::cout << "CG_doublearg2(" << j1 << "," << m1 << ","
                 << j2 << "," << m2 << ","
                 << j << "," << m << ") = "
                 << resultDouble2 << std::endl;
    }


}

#ifndef RUNNING_COMBINED_TESTS
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
#endif
