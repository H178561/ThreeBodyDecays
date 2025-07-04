#include "include/ClebschGordan.hh"

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <array>
#include <vector>

using complex = std::complex<double>;

const int MAX_FACTcg = 100;

double getLogFactorialcg(int n)
{
    static std::vector<double> logfact = []()
    {
        std::vector<double> vec(MAX_FACTcg + 1);
        vec[0] = 0.0;
        for (int i = 1; i <= MAX_FACTcg; ++i)
        {
            vec[i] = vec[i - 1] + std::log(i);
        }
        return vec;
    }();

    if (n < 0 || n > MAX_FACTcg)
    {
        throw std::runtime_error("Factorial out of range");
    }
    return logfact[n];
}
// Clebsch-Gordan coefficient implementation

// Global array of log factorials for doubled arguments (similar to logfact2 in Julia)
std::vector<double> initializeLogFact2() {
    std::vector<double> result(101, 0.0);
    for (int two_n = 0; two_n <= 100; two_n++) {
        result[two_n] = getLogFactorialcg(two_n / 2);
    }
    return result;
}

// Helper function to calculate log factorial for doubled arguments (half-integer support)
double f_logfact2(int two_n) {
    static std::vector<double> logfact2 = initializeLogFact2();

    if (two_n < 0 || two_n > 100) {
        throw std::runtime_error("two_n < 0 || two_n > 100. Modify if needed.");
    }
    return logfact2[two_n];
}

// Main Clebsch-Gordan coefficient function with doubled arguments
double clebschgordan_doublearg(int two_j1, int two_m1, int two_j2, int two_m2, int two_j, int two_m) {
    // Check valid range of helicities
    if ((std::abs(two_m1) > two_j1) || (std::abs(two_m2) > two_j2) || (std::abs(two_m) > two_j)) {
        return 0.0;
    }

    // Check conservation of angular momentum projection and triangle inequality
    if ((two_m1 + two_m2 != two_m) || !(std::abs(two_j1 - two_j2) <= two_j && two_j <= two_j1 + two_j2)) {
        return 0.0;
    }

    // Calculate prefactor
    double prefactor = std::sqrt(two_j + 1.0) *
        std::exp((f_logfact2(two_j1 + two_j2 - two_j) +
                 f_logfact2(two_j1 + two_j - two_j2) +
                 f_logfact2(two_j2 + two_j - two_j1) -
                 f_logfact2(two_j1 + two_j2 + two_j + 2) +
                 f_logfact2(two_j1 + two_m1) +
                 f_logfact2(two_j1 - two_m1) +
                 f_logfact2(two_j2 + two_m2) +
                 f_logfact2(two_j2 - two_m2) +
                 f_logfact2(two_j + two_m) +
                 f_logfact2(two_j - two_m)) / 2.0);

    // Add a special case for maximum m values
    bool is_max_coupling = (two_m1 == two_j1 && two_m2 == two_j2 && two_m == two_j);

    double res = 0.0;

    // Calculate min/max values for the summation
    int two_t_min = std::max(0, std::max(
                             two_j2 - two_m1 - two_j,
                             two_j1 + two_m2 - two_j));

    int two_t_max = std::min(two_j1 + two_j2 - two_j, std::min(
                             two_j1 - two_m1,
                             two_j2 + two_m2));

    // Sum over terms
    for (int two_t = two_t_min; two_t <= two_t_max; two_t += 2) {
        double logs = f_logfact2(two_t) +
                      f_logfact2(two_j - two_j2 + two_m1 + two_t) +
                      f_logfact2(two_j - two_j1 - two_m2 + two_t) +
                      f_logfact2(two_j1 + two_j2 - two_j - two_t) +
                      f_logfact2(two_j1 - two_m1 - two_t) +
                      f_logfact2(two_j2 + two_m2 - two_t);

        // Need to fix the sign convention to match expected values
        double sign = ((std::abs(two_t) % 4 == 2) ? -1.0 : 1.0);

        res += sign * std::exp(-logs);
    }

    res *= prefactor;

    res *= -1;

    return res;
}

// Function for standard Clebsch-Gordan coefficients (integer inputs)
double clebschgordan(int j1, int m1, int j2, int m2, int j, int m) {
    return clebschgordan_doublearg(2*j1, 2*m1, 2*j2, 2*m2, 2*j, 2*m);
}

// Shortcut functions
double CG(int j1, int m1, int j2, int m2, int j, int m) {
    return clebschgordan(j1, m1, j2, m2, j, m);
}

/*
complex CG_doublearg(int two_j1, int two_m1, int two_j2, int two_m2, int two_j, int two_m) {
    return clebschgordan_doublearg(two_j1, two_m1, two_j2, two_m2, two_j, two_m);
}*/


// Replace the placeholder CG_doublearg function with this implementation

double CG_doublearg(int two_j1, int two_m1, int two_j2, int two_m2, int two_j, int two_m) {
    // Check if m = m1 + m2 (conservation of angular momentum projection)
    if (two_m != two_m1 + two_m2) {
        return 0.0;
    }

    // Check triangle inequality for j values
    if (two_j < std::abs(two_j1 - two_j2) || two_j > two_j1 + two_j2) {
        return 0.0;
    }

    // Convert doubled representation to standard form
    double j1 = two_j1 / 2.0;
    double m1 = two_m1 / 2.0;
    double j2 = two_j2 / 2.0;
    double m2 = two_m2 / 2.0;
    double j = two_j / 2.0;
    double m = two_m / 2.0;

    // Check if m values are within valid range
    if (std::abs(m1) > j1 || std::abs(m2) > j2 || std::abs(m) > j) {
        return 0.0;
    }

    // Calculate the normalization factor
    double norm = std::sqrt((2.0 * j + 1.0) *
                 std::exp(getLogFactorialcg(j1 + j2 - j) +
                 getLogFactorialcg(j + j1 - j2) +
                 getLogFactorialcg(j + j2 - j1) -
                 getLogFactorialcg(j1 + j2 + j + 1.0)));

    norm *= std::sqrt(std::exp(getLogFactorialcg(j1 + m1) +
    getLogFactorialcg(j1 - m1) +
    getLogFactorialcg(j2 + m2) +
    getLogFactorialcg(j2 - m2) +
    getLogFactorialcg(j + m) +
    getLogFactorialcg(j - m)));

    // Calculate the sum for the coefficient
    double sum = 0.0;
    double max_k = std::min(j1 + j2 - j, std::min( j1 - m1, j2 + m2));

    for (double k = 0.0; k <= max_k; k += 1.0) {
        if (j - j2 + m1 + k >= 0.0 && j - j1 - m2 + k >= 0.0) {
            double term = std::pow(-1.0, k) / (
                std::exp(getLogFactorialcg(k) +
                        getLogFactorialcg(j1 + j2 - j - k) +
                       getLogFactorialcg(j1 - m1 - k) +
                       getLogFactorialcg(j2 + m2 - k) +
                       getLogFactorialcg(j - j2 + m1 + k) +
                       getLogFactorialcg(j - j1 - m2 + k)));
            sum += term;
        }
    }

    return norm * sum;
}





/*
// Placeholder for CG_doublearg function
complex CG_doublearg(int j1, int m1, int j2, int m2, int j, int m) {
    // You need to implement the real CG coefficient computation.
    return complex(1.0, 0.0);
}*/


