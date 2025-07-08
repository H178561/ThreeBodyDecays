#include "include/WignerFunctions.hh"
#include "include/jacobi.hh"
#include <cmath>
#include <iostream>

extern bool debug; // Reference the global debug flag

// Function to get log factorial (moved from ThreeBodyDecays.cpp)
double getLogFactorial_wf(int n)
{
    static const int MAX_FACT = 100;
    static std::vector<double> logfact = []()
    {
        std::vector<double> vec(MAX_FACT + 1);
        vec[0] = 0.0;
        for (int i = 1; i <= MAX_FACT; ++i)
        {
            vec[i] = vec[i - 1] + std::log(i);
        }
        return vec;
    }();

    if (n < 0 || n > MAX_FACT)
    {
        throw std::runtime_error("Factorial out of range");
    }
    return logfact[n];
}

// Wigner d-hat function with doubled arguments
double wignerd_hat_doublearg(int two_j, int two_m1, int two_m2, double z)
{
    // Check valid range of input parameters
    if (std::abs(two_m1) > two_j || std::abs(two_m2) > two_j)
    {
        return 0.0;
    }

    // Determine sign factor
    int am1m2 = std::abs(two_m1 - two_m2);
    int ap1p2 = std::abs(two_m1 + two_m2);
    double factor = 1.0;
    if (((am1m2 + two_m1 - two_m2) % 8) == 4)
    {
        factor = -1.0;
    }

    // Compute j_minus_M, j_plus_M, j_minus_N, j_plus_N
    int two_am1 = std::abs(two_m1);
    int two_am2 = std::abs(two_m2);

    int two_M, two_N;
    if (two_am1 > two_am2)
    {
        two_M = two_am1;
        two_N = two_am2;
    }
    else
    {
        two_M = two_am2;
        two_N = two_am1;
    }

    int j_minus_M = (two_j - two_M) / 2;
    int j_plus_M = (two_j + two_M) / 2;
    int j_minus_N = (two_j - two_N) / 2;
    int j_plus_N = (two_j + two_N) / 2;

    // Compute logarithm of gamma ratio
    double log_gamma_ratio = 0.5 * (getLogFactorial_wf(j_minus_M) +
                                    getLogFactorial_wf(j_plus_M) -
                                    getLogFactorial_wf(j_minus_N) -
                                    getLogFactorial_wf(j_plus_N));

    // Compute the Jacobi polynomial
    double jacobi_val = boost::math::jacobi<double>(j_minus_M, am1m2 / 2, ap1p2 / 2, z);

    // Return the final result
    return factor * std::pow(2.0, -two_M / 2.0) * std::exp(log_gamma_ratio) * jacobi_val;
}

// Implementation of wignerd_doublearg
double wignerd_doublearg(int two_j, int two_m1, int two_m2, double z)
{
    // Special cases
    if (std::abs(z - 1.0) < 1e-10)
    {
        return (two_m1 == two_m2) ? 1.0 : 0.0;
    }

    if (std::abs(z + 1.0) < 1e-10)
    {
        if (two_m1 == -two_m2)
        {
            return ((two_j - two_m2) % 4 == 0) ? 1.0 : -1.0;
        }
        else
        {
            return 0.0;
        }
    }

    double hat = wignerd_hat_doublearg(two_j, two_m1, two_m2, z);
    double xi = std::pow(1.0 - z, std::abs(two_m1 - two_m2) / 4.0) *
                std::pow(1.0 + z, std::abs(two_m1 + two_m2) / 4.0);

    if (debug)
        std::cout << two_j << " " << two_m1 << " " << two_m2 << " " << z << " " << hat << " " << xi << " res " << hat * xi << std::endl;

    return hat * xi;
}

// Function to compute Wigner-d with sign adjustment for a single element
double wignerd_doublearg_sign_element(int two_j, int two_λ1, int two_λ2, double cosθ, bool ispositive)
{
    double phase_value = 0.0;
    if (ispositive)
    {
        phase_value = 1.0;
    }
    else
    {
        int abs_diff = std::abs(two_λ1 - two_λ2);
        phase_value = (abs_diff % 4 == 2) ? -1.0 : 1.0;
    }

    if (debug)
        std::cout << two_j << " " << two_λ1 << " " << two_λ2 << " " << cosθ << " " << phase_value << std::endl;

    return phase_value * wignerd_doublearg(two_j, two_λ1, two_λ2, cosθ);
}

// Function to compute Wigner-d matrix with sign adjustment
Matrix2D wignerd_doublearg_sign(int two_j, double cosθ, bool ispositive)
{
    // Create a range from -two_j to two_j with step size 2
    std::vector<int> range;
    for (int val = -two_j; val <= two_j; val += 2)
    {
        range.push_back(val);
    }

    // Initialize the result matrix
    size_t size = range.size();
    Matrix2D result(size, std::vector<double>(size, 0.0));

    // Compute the Wigner-D matrix with sign adjustment
    for (size_t i = 0; i < size; ++i)
    {
        for (size_t j = 0; j < size; ++j)
        {
            int two_λ1 = range[i];
            int two_λ2 = range[j];
            result[i][j] = wignerd_doublearg_sign_element(two_j, two_λ1, two_λ2, cosθ, ispositive);
            if (debug)
                std::cout << i << " " << j << " " << result[i][j] << std::endl;
        }
    }

    return result;
}
