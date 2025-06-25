#include "ThreeBodyDecays.hh"
// #include "ThreeBodyUtilities.hh"
#include "ClebschGordan.hh"
#include "FormFactors.hh"
#include "jacobi.hh"
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <array>
#include <unistd.h> // Für getcwd()

using complex = std::complex<double>;
bool debug = false;
bool debugls = false;
bool div2 = false;

// Kallen function implementation
double Kallen(double x, double y, double z)
{
    return x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * z * x;
}

// Check if two doubles are approximately equal
bool approx_equal(double a, double b, double epsilon)
{
    return std::abs(a - b) < epsilon;
}

double phase(int value)
{
    int phase = 1;
    // Return +1 if value is even, -1 if odd
    if (std::abs(value % 4) == 2)
    {
        phase = -1;
    }
    return (value > 0) ? -phase : phase;
}

double phase2(int value)
{
    return (value % 4 == 2) ? -1. : 1.;
}

double phaseabs(int value)
{
    // Return +1 if value is even, -1 if odd
    return (std::abs(value % 4) == 2) ? -1. : 1.;
}

// Helper function for padding indices
int pad(int index, int size)
{
    if (index < 0)
        return 0;
    if (index >= size)
        return size - 1;
    return index;
}

// Global array of log factorials (precomputed for performance)
// std::vector<double> logfact;
const int MAX_FACT = 100;
/*
// Initialize the log factorial table
void init_logfact() {
    logfact.resize(MAX_FACT + 1);
    logfact[0] = 0.0;
    for (int i = 1; i <= MAX_FACT; ++i) {
        logfact[i] = logfact[i-1] + std::log(i);
    }
}*/

double getLogFactorial(int n)
{
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

double jacobi_pols(int n, int a, int b, double z)
{
    return boost::math::jacobi<double>(n, a, b, z);
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
    double log_gamma_ratio = 0.5 * (getLogFactorial(j_minus_M) +
                                    getLogFactorial(j_plus_M) -
                                    getLogFactorial(j_minus_N) -
                                    getLogFactorial(j_plus_N));

    // Compute the Jacobi polynomial
    double jacobi_val = jacobi_pols(j_minus_M, am1m2 / 2, ap1p2 / 2, z);
    double norm_factor = std::pow(0.5, two_M / 2.0); // Oder: 1.0 / std::pow(2.0, two_M / 2.0)

    // Return the final result
    return factor * std::pow(2.0, -two_M / 2.0) * std::exp(log_gamma_ratio) * jacobi_val;
}

// Implementation of wignerd_doublearg (as defined earlier)
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

    // For a complete implementation, the general case would need to be computed
    // using either a recursive formula or the Jacobi polynomial approach
    // This is a placeholder for the actual computation
    double hat = wignerd_hat_doublearg(two_j, two_m1, two_m2, z);
    double xi = std::pow(1.0 - z, std::abs(two_m1 - two_m2) / 4.0) *
                std::pow(1.0 + z, std::abs(two_m1 + two_m2) / 4.0);

    if (debug)
        std::cout << two_j << " " << two_m1 << " " << two_m2 << " " << z << " " << hat << " " << xi << " res " << hat * xi << std::endl;
    return hat * xi;
}

/*
// Function to compute Wigner-d with sign adjustment for a single element
double wignerd_doublearg_sign_element(int two_j, int two_λ1, int two_λ2, double cosθ, bool ispositive)
{
    int sign = ispositive ? 1 : phase((two_λ1 - two_λ2) / 2);
    return sign * wignerd_doublearg(two_j, two_λ1, two_λ2, cosθ);
}*/

double wignerd_doublearg_sign_element(int two_j, int two_λ1, int two_λ2, double cosθ, bool ispositive)
{
    // double phase = ispositive ? 1.0 : std::pow(-1.0, (two_λ1 - two_λ2) / 2);

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

double scalar_product(std::vector<double> a, std::vector<double> b)
{
    double product = 0;
    for (int i = 0; i <= a.size() - 1; i++)
    {
        product = product + (a[i]) * (b[i]);
    }
    return product;
}

// Minkowski scalar product
double min_scalar_product(std::vector<double> v1, std::vector<double> v2)
{
    return v1[3] * v2[3] - (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
}

// Helper function to get indices i and j from k
std::pair<int, int> ij_from_k(int k)
{
    // std::cout << k << std::endl;
    switch (k)
    {
    case 1:
        return {2, 3};
    case 2:
        return {3, 1};
    case 3:
        return {1, 2};
    default:
        throw std::invalid_argument(
            "Invalid value for k. Must be 1, 2, or 3.");
    }
}

// Function to compute i, j, k indices from Wigner rotation
std::tuple<int, int, int> ijk(const AbstractWignerRotation &wr)
{
    int k = wr.get_k();

    // if (debug)
    //     std::cout << "ijk" << k << std::endl;
    auto [i, j] = ij_from_k(k);
    return {i, j, k};
}

// Helper function to determine if a system is sequential
bool issequential(int i, int j)
{
    int diff = j - i;
    return (diff == 1 || diff == -2);
}

// Helper function to create Wigner rotation
std::unique_ptr<AbstractWignerRotation> wr(int system_a, int reference_b,
                                           int particle_c)
{
    if (debug)
        std::cout << "Creating Wigner rotation for system_a: " << system_a << ", reference_b: " << reference_b << ", particle_c: " << particle_c << std::endl;
    if (system_a == reference_b)
    {
        return std::make_unique<TrivialWignerRotation>(particle_c);
    }

    bool S = issequential(system_a, reference_b);
    int A = S ? system_a : reference_b;
    int B = S ? reference_b : system_a;

    if (particle_c == 0)
    {
        return std::make_unique<Arg0WignerRotation>(A, S);
        // return std::make_unique<Arg0WignerRotation>( particle_c, S );
    }

    if (particle_c != system_a && particle_c != reference_b)
    {
        return std::make_unique<Arg3WignerRotation>(particle_c, S);
    }

    bool T = (particle_c == A);
    return std::make_unique<Arg2WignerRotation>(particle_c, !S, T);
}

// Explicit Wigner rotation functions
double cosζ21_for1(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(2, 1, 1))(0, σs, ms2);
}

double cosζ21_for2(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(2, 1, 2))(0, σs, ms2);
}

double cosζ13_for1(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(1, 3, 1))(0, σs, ms2);
}

double cosζ13_for3(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(1, 3, 3))(0, σs, ms2);
}

double cosζ32_for3(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(3, 2, 3))(0, σs, ms2);
}

double cosζ32_for2(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(3, 2, 2))(0, σs, ms2);
}

double cosζ12_for3(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(1, 2, 3))(0, σs, ms2);
}

double cosζ23_for1(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(2, 3, 1))(0, σs, ms2);
}

double cosζ31_for2(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(3, 1, 2))(0, σs, ms2);
}

double cosζ12_for0(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(1, 2, 0))(0, σs, ms2);
}

double cosζ23_for0(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(2, 3, 0))(0, σs, ms2);
}

double cosζ31_for0(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(3, 1, 0))(0, σs, ms2);
}

// Function to compute cosζ for Wigner rotation
double cosζ(const AbstractWignerRotation &wr, const std::array<double, 3> &sigma,
            const std::array<double, 4> &ms2)
{
    // Get the indices i, j, k from the Wigner rotation
    auto [i, j, k] = ijk(wr);

    // Calculate the energy and momentum terms
    double s = ms2[3]; // Parent mass squared
    double EE4σ = (sigma[0] + ms2[1] - ms2[2]) * (s - sigma[0] - ms2[0]);
    double pp4σ = std::sqrt(Kallen(sigma[0], ms2[1], ms2[2]) *
                            Kallen(s, sigma[0], ms2[0]));
    pp4σ = (pp4σ < 0) ? 0.0 : pp4σ; // Handle numerical errors

    // Calculate the rest term
    double rest = sigma[1] - ms2[0] - ms2[1];

    // Calculate cosζ using the same formula as Julia
    return (2 * sigma[0] * rest - EE4σ) / pp4σ;
}

// Function to compute cosθij
double ThreeBodyDecays::cosθij(const std::array<double, 3> &σs,
                               const std::array<double, 4> &msq, int k)
{
    auto [i, j] = ij_from_k(k);

    // Adjust indices for 0-based indexing in C++ arrays
    i--;
    j--;
    k--;

    double s = msq[3]; // Parent mass squared
    double EE4σ = (σs[k] + msq[i] - msq[j]) * (s - σs[k] - msq[k]);
    double pp4σ = std::sqrt(Kallen(σs[k], msq[i], msq[j]) *
                            Kallen(s, σs[k], msq[k]));
    double rest = σs[j] - msq[k] - msq[i];

    return (2 * σs[k] * rest - EE4σ) / pp4σ;
}

// Convenience functions for specific angles
double ThreeBodyDecays::cosθ23(MandelstamTuple σs, ThreeBodyMasses msq)
{
    return cosθij(σs, msq, 1);
}

double ThreeBodyDecays::cosθ31(MandelstamTuple σs, ThreeBodyMasses msq)
{
    return cosθij(σs, msq, 2);
}

double ThreeBodyDecays::cosθ12(MandelstamTuple σs, ThreeBodyMasses msq)
{
    return cosθij(σs, msq, 3);
}

// Utility function to perform modular arithmetic with 1-based indexing
int mod1(int k, int n)
{
    return ((k - 1) % n + n) % n + 1;
}

// circleorigin function: reorders the tuple elements based on k
std::array<double, 3> circleorigin(int k, const std::array<double, 3> &t)
{
    return {
        t[mod1(k - 1, 3) - 1], // First element
        t[mod1(k, 3) - 1],     // Second element
        t[mod1(k + 1, 3) - 1]  // Third element
    };
}

// Linear interpolation function
double fitin(double y, const std::pair<double, double> &limits)
{
    double a = limits.first;
    double b = limits.second;
    return a + y * (b - a);
}

// Calculate limits for a mass configuration
std::pair<double, double> lims(const std::array<double, 4> &ms, int k)
{
    auto [i, j] = ij_from_k(k);
    i--;
    j--;
    k--; // Convert to 0-based indexing

    double min_limit = std::pow(ms[i] + ms[j], 2);
    double max_limit = std::pow(ms[3] - ms[k], 2);

    return {min_limit, max_limit};
}

double ThreeBodyDecays::σjofk(double cos_theta, double σk,
                              const std::array<double, 4> &ms_squared, int k)
{
    auto [i, j] = ij_from_k(k);
    i--;
    j--;
    k--; // Convert to 0-based indexing

    double s = ms_squared[3]; // Parent mass squared

    // Calculate energy and momentum terms using the same formula as Julia
    double EE4σ = (σk + ms_squared[j] - ms_squared[i]) *
                  (σk + s - ms_squared[k]);
    double p2q24σ = Kallen(σk, ms_squared[i], ms_squared[j]) *
                    Kallen(s, σk, ms_squared[k]);
    p2q24σ = (p2q24σ < 0) ? 0.0
                          : p2q24σ; // Handle numerical errors like Julia

    // Calculate σi using the same formula as Julia
    double σi = s + ms_squared[j] -
                (EE4σ - std::sqrt(p2q24σ) * cos_theta) / (2 * σk);

    return σi;
}

// Convert from random numbers x to Mandelstam variables
MandelstamTuple ThreeBodyDecays::x2σs(const std::array<double, 2> &x,
                                      ThreeBodyMasses ms, int k)
{
    // Calculate σk through linear interpolation
    double σk = fitin(x[0], lims(ms, k));

    // Calculate squared masses
    std::array<double, 4> ms_squared;
    for (int i = 0; i < 4; i++)
    {
        ms_squared[i] = ms[i] * ms[i];
    }

    // Calculate σj based on the scattering angle
    double σj = σjofk(2 * x[1] - 1, σk, ms_squared, k);

    // Calculate σi through energy-momentum conservation
    double σi = 0.0;
    for (const auto &m_squared : ms_squared)
    {
        σi += m_squared;
    }
    σi -= σk + σj;

    // Create array of σ values and apply circleorigin
    std::array<double, 3> σt = {σi, σj, σk};
    return circleorigin(k, σt);
}

complex parseComplex(const std::string &str)
{
    std::istringstream iss(str);
    double real, imag;
    char plus_minus, i;
    iss >> real >> plus_minus >> imag >> i;
    if (plus_minus == '-')
    {
        imag = -imag;
    }
    return complex(real, imag);
}

double parseFraction(const std::string &str)
{
    size_t pos = str.find('/');
    if (pos != std::string::npos)
    {
        // Es ist ein Bruch "a/b"
        int numerator = stoi(str.substr(0, pos));    // Zähler
        int denominator = stoi(str.substr(pos + 1)); // Nenner
        return static_cast<double>(numerator) / denominator;
    }
    // Es ist eine ganze Zahl
    return std::stod(str);
}


// aligned_amplitude implementation that returns a Tensor4D
// Korrigierte aligned_amplitude4d-Funktion
Tensor4D ThreeBodyDecays::aligned_amplitude4d(const DecayChain &dc, const MandelstamTuple &σs)
{
    int k = dc.k;
    const auto &tbs = dc.tbs;
    int two_j = dc.two_j;
    const auto &two_js = tbs.two_js;

    // Get indices i, j from k (1-based indexing in result)
    auto [i, j] = ij_from_k(k);
    // Konvertieren zu 0-basierter Indexierung für Arrays
    int i_idx = i - 1;
    int j_idx = j - 1;
    int k_idx = k - 1;

    // Überprüfen der Grenzen für sicheren Array-Zugriff
    if (i_idx < 0 || i_idx >= 3 || j_idx < 0 || j_idx >= 3 || k_idx < 0 || k_idx >= 3)
    {
        std::cout << "Invalid indices in aligned_amplitude4d: i="
                  << i << ", j=" << j << ", k=" << k << std::endl;
        // Rückgabe eines leeren Tensors bei Fehler
        return Tensor4D();
    }

    // Quadrierte Massen berechnen
    std::array<double, 4> msq = {tbs.ms[0] * tbs.ms[0],
                                 tbs.ms[1] * tbs.ms[1],
                                 tbs.ms[2] * tbs.ms[2],
                                 tbs.ms[3] * tbs.ms[3]};

    // cosθ und Wigner d-Matrix berechnen
    double cosθ = cosθij(σs, msq, k);
    // std::cout << "cosθ: " << cosθ << std::endl;
    Matrix2D d_θ = wignerd_doublearg_sign(two_j, cosθ, true);
    // std::cout << "cosθ: " << cosθ << std::endl;
    /*
    // Wigner d Test ok
    std::cout << "d_θ dimensions: "
            << d_θ.size() << " x "
            << (d_θ.empty() ? 0 : d_θ[0].size()) << std::endl;

    for (int i = 0; i < d_θ.size(); ++i) {
        for (int j = 0; j < d_θ[0].size(); ++j) {
            std::cout << d_θ[i][j] << "\t";  // Tab für schöne Ausrichtung
        }
        std::cout << "\n";
    }
        */

    // Spins
    std::array<int, 3> two_js_Hij = {two_j, two_js[i_idx], two_js[j_idx]};
    std::array<int, 3> two_js_HRk = {two_js[3], two_j, two_js[k_idx]};
    if (debug)
        std::cout << "two_js_Hij: " << two_js_Hij[0] << " " << two_js_Hij[1] << " " << two_js_Hij[2] << std::endl;
    if (debug)
        std::cout << "two_js_HRk: " << two_js_HRk[0] << " " << two_js_HRk[1] << " " << two_js_HRk[2] << std::endl;
    // VRk-Matrix berechnen
    std::vector<std::vector<double>> VRk;
    int vrk_dim1 = two_j + 1;
    int vrk_dim2 = two_js[k_idx] + 1;
    // VRk-Matrix initialisieren
    VRk.resize(vrk_dim1, std::vector<double>(vrk_dim2, 0.0));

    for (int m1_idx = 0; m1_idx < vrk_dim1; m1_idx++)
    {
        int two_m1 = -two_j + 2 * m1_idx;
        for (int m2_idx = 0; m2_idx < vrk_dim2; m2_idx++)
        {
            int two_m2 = -two_js[k_idx] + 2 * m2_idx;

            std::array<int, 2> two_ms = {two_m1, two_m2};
            complex amp = amplitude_recoupling(dc.HRk, two_ms, two_js_HRk);
            double phase_value = (((two_js[k_idx] - two_m2) % 4 == 0) ? 1.0 : -1.0);

            VRk[m1_idx][m2_idx] = amp.real() * phase_value;
        }
    }

    if (debug)
        std::cout << "VRk dimensions: "
                  << VRk.size() << " x "
                  << (VRk.empty() ? 0 : VRk[0].size()) << std::endl;
    for (int i = 0; i < VRk.size(); ++i)
    {
        for (int j = 0; j < VRk[0].size(); ++j)
        {
            if (debug)
                std::cout << VRk[i][j] << "\t"; // Tab für schöne Ausrichtung
        }
        if (debug)
            std::cout << "\n";
    }

    // Vij-Matrix berechnen
    std::vector<std::vector<double>> Vij;
    int vij_dim1 = two_js[i_idx] + 1;
    int vij_dim2 = two_js[j_idx] + 1;

    // Vij-Matrix initialisieren
    Vij.resize(vij_dim1, std::vector<double>(vij_dim2, 0.0));

    for (int m1_idx = 0; m1_idx < vij_dim1; m1_idx++)
    {
        int two_m1 = -two_js[i_idx] + 2 * m1_idx;
        for (int m2_idx = 0; m2_idx < vij_dim2; m2_idx++)
        {
            int two_m2 = -two_js[j_idx] + 2 * m2_idx;

            std::array<int, 2> two_ms = {two_m1, two_m2};
            complex amp = amplitude_recoupling(dc.Hij, two_ms, two_js_Hij);
            double phase_value = phase(two_js[k_idx] - two_m2);
            Vij[m1_idx][m2_idx] = amp.real() * phase_value;
        }
    }


    // Δ-Verschiebungen berechnen
    int Δ_zk = (two_j - two_js[3] - two_js[k_idx]) / 2;
    int Δ_ij = (two_j - two_js[i_idx] + two_js[j_idx]) / 2;

    // Lineshape berechnen
    double lineshape = dc.Xlineshape(σs[k_idx]).real();

    // Dimensionen für den Ergebnistensor
    std::vector<int> dims = {two_js[0] + 1, two_js[1] + 1, two_js[2] + 1, two_js[3] + 1};

    // F0 mit Nullen initialisieren (Ergebnistensor)
    Tensor4D F0(dims[0], std::vector<std::vector<std::vector<double>>>(
                             dims[1], std::vector<std::vector<double>>(
                                          dims[2], std::vector<double>(dims[3], 0.0))));

    // Dimensionen für Berechnungstensor
    Tensor4D F(dims[i_idx], std::vector<std::vector<std::vector<double>>>(
                                dims[j_idx], std::vector<std::vector<double>>(
                                                 dims[k_idx], std::vector<double>(dims[3], 0.0))));

    // Hilfsfunktion zum Begrenzen von Indizes
    auto pad = [](int idx, int max_size)
    {
        return std::max(0, std::min(idx, max_size - 1));
    };

    // F berechnen (ähnlich @tullio in Julia)
    for (int _i = 0; _i < dims[i_idx]; _i++)
    {
        for (int _j = 0; _j < dims[j_idx]; _j++)
        {
            for (int _k = 0; _k < dims[k_idx]; _k++)
            {
                for (int _z = 0; _z < dims[3]; _z++)
                {
                    int vrk_idx1 = pad(_z + _k + Δ_zk, two_j + 1);
                    int vrk_idx2 = _k;
                    int d_idx1 = pad(_z + _k + Δ_zk, two_j + 1);
                    int d_idx2 = pad(_i - _j + Δ_ij, two_j + 1);

                    // Sicherstellung, dass alle Indizes innerhalb der Grenzen sind
                    if (vrk_idx1 < (int)VRk.size() && vrk_idx2 < (int)VRk[0].size() &&
                        d_idx1 < (int)d_θ.size() && d_idx2 < (int)d_θ[0].size() &&
                        _i < (int)Vij.size() && _j < (int)Vij[0].size())
                    {

                        F[_i][_j][_k][_z] +=
                            VRk[vrk_idx1][vrk_idx2] *
                            d_θ[d_idx1][d_idx2] *
                            Vij[_i][_j];
                    }
                }
            }
        }
    }

    // Normalisierung und Lineshape anwenden
    double d_norm = std::sqrt(two_j + 1.0);
    for (auto &dim1 : F)
    {
        for (auto &dim2 : dim1)
        {
            for (auto &dim3 : dim2)
            {
                for (auto &val : dim3)
                {
                    val *= d_norm * lineshape;
                }
            }
        }
    }
    if (debug)
        std::cout << d_norm << " " << lineshape << std::endl;

    // Zurück zur ursprünglichen Reihenfolge permutieren
    for (int _0 = 0; _0 < dims[0]; _0++)
    {
        for (int _1 = 0; _1 < dims[1]; _1++)
        {
            for (int _2 = 0; _2 < dims[2]; _2++)
            {
                for (int _3 = 0; _3 < dims[3]; _3++)
                {
                    // Indizes gemäß Permutation zuordnen
                    int perm_i = (i_idx == 0) ? _0 : (i_idx == 1) ? _1
                                                 : (i_idx == 2)   ? _2
                                                                  : _3;
                    int perm_j = (j_idx == 0) ? _0 : (j_idx == 1) ? _1
                                                 : (j_idx == 2)   ? _2
                                                                  : _3;
                    int perm_k = (k_idx == 0) ? _0 : (k_idx == 1) ? _1
                                                 : (k_idx == 2)   ? _2
                                                                  : _3;
                    int perm_z = 3; // Parent index is always last (3)

                    // Überprüfen der Grenzen für sicheren Array-Zugriff
                    if (perm_i < dims[i_idx] && perm_j < dims[j_idx] &&
                        perm_k < dims[k_idx] && _3 < dims[3])
                    {

                        F0[_0][_1][_2][_3] = F[perm_i][perm_j][perm_k][_3];
                    }
                }
            }
        }
    }

    return F0;
}

complex ThreeBodyDecays::amplitude_recoupling(
    const RecouplingLS &recoupling,
    const std::array<int, 2> &two_ms,
    const std::array<int, 3> &two_js)
{

    // std::function mit std::bind erlaubt kein dynamic_cast, daher müssen wir einen
    // anderen Ansatz verwenden. Wir können versuchen, verschiedene Recoupling-Typen
    // durch Funktionspointer zu identifizieren

    // Versuche, den gespeicherten void* Pointer zurückzugewinnen und zu prüfen
    using FunctionType = complex (*)(const std::array<int, 2> &, const std::array<int, 3> &);
    // std::cout << two_ms[0] << " " << two_ms[1] << " " << two_js[0] << " " << two_js[1] << " " << two_js[2] << std::endl;
    //  Direkte Verwendung der Funktion
    return recoupling(two_ms, two_js);
}

// Implementation of amplitude function that returns a Tensor4D
Tensor4D ThreeBodyDecays::amplitude4d(const DecayChain &dc,
                                      const MandelstamTuple &σs,
                                      const std::vector<int> &refζs)
{
    int k = dc.k;
    const auto &tbs = dc.tbs;
    int two_j = dc.two_j;
    const auto &two_js = tbs.two_js;
    const auto &ms = tbs.ms;
    // std::array<int, 4> two_js = {two_jst[0]*2, two_jst[1]*2, two_jst[2]*2, two_jst[3]*2};

    if (debug)
        std::cout << two_js[0] << " " << two_js[1] << " " << two_js[2] << " " << two_js[3] << std::endl;

    // Calculate squared masses
    std::array<double, 4> msq = {ms[0] * ms[0],
                                 ms[1] * ms[1],
                                 ms[2] * ms[2],
                                 ms[3] * ms[3]};

    // Get aligned amplitude
    auto F0 = aligned_amplitude4d(dc, σs);


    // Calculate alignment rotations
    std::vector<Matrix2D> d_ζs(4);
    for (size_t l = 0; l < 4; ++l)
    {
        int _two_j = two_js[l];
        int _refζ = refζs[l];
        std::unique_ptr<AbstractWignerRotation> _w = wr(k, _refζ, l % 4);
        double _cosζ = _w->cos_zeta(σs, msq);
        d_ζs[l] = wignerd_doublearg_sign(_two_j, _cosζ, _w->is_positive());
    }

    // Alias the d_ζs for clarity
    const auto &D1 = d_ζs[0];
    const auto &D2 = d_ζs[1];
    const auto &D3 = d_ζs[2];
    const auto &D0 = d_ζs[3];

    // Create a new tensor with the same dimensions as F0
    std::vector<int> dims;
    for (int i = 0; i < 4; i++)
    {
        dims.push_back(two_js[i] + 1);
    }

    Tensor4D F(dims[0], std::vector<std::vector<std::vector<double>>>(
                            dims[1], std::vector<std::vector<double>>(
                                         dims[2], std::vector<double>(dims[3], 0.0))));

    // Calculate the tensor contraction (equivalent to @tullio in Julia)
    // @tullio F[_i, _j, _k, _z] = D0[_z, _z′] * F0[_i′, _j′, _k′, _z′] * D1[_i′, _i] * D2[_j′, _j] * D3[_k′, _k]
    for (int _i = 0; _i < dims[0]; ++_i)
    {
        for (int _j = 0; _j < dims[1]; ++_j)
        {
            for (int _k = 0; _k < dims[2]; ++_k)
            {
                for (int _z = 0; _z < dims[3]; ++_z)
                {
                    // Sum over all primed indices
                    for (int _i_prime = 0; _i_prime < dims[0]; ++_i_prime)
                    {
                        if (_i_prime >= (int)D1.size() || _i >= (int)D1[0].size())
                            continue;

                        for (int _j_prime = 0; _j_prime < dims[1]; ++_j_prime)
                        {
                            if (_j_prime >= (int)D2.size() || _j >= (int)D2[0].size())
                                continue;

                            for (int _k_prime = 0; _k_prime < dims[2]; ++_k_prime)
                            {
                                if (_k_prime >= (int)D3.size() || _k >= (int)D3[0].size())
                                    continue;

                                for (int _z_prime = 0; _z_prime < dims[3]; ++_z_prime)
                                {
                                    if (_z_prime >= (int)D0.size() || _z >= (int)D0[0].size())
                                        continue;

                                    F[_i][_j][_k][_z] +=
                                        D0[_z][_z_prime] *
                                        F0[_i_prime][_j_prime][_k_prime][_z_prime] *
                                        D1[_i_prime][_i] *
                                        D2[_j_prime][_j] *
                                        D3[_k_prime][_k];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return F;
}

Tensor4Dcomp ThreeBodyDecays::aligned_amplitude4dcomp(const DecayChain &dc, const MandelstamTuple &σs)
{
    int k = dc.k;
    const auto &tbs = dc.tbs;
    int two_j = dc.two_j;
    // const auto &two_js = tbs.two_js;
    auto two_js = tbs.two_js;
    // std::array<int, 4> two_js = {two_jst[0]*2, two_jst[1]*2, two_jst[2]*2, two_jst[3]*2}
    if (debug)
        std::cout << "two_j" << two_j << " " << two_js[0] << " " << two_js[1] << " " << two_js[2] << " " << two_js[3] << std::endl;

    // div 2
    if (div2)
    {
        two_js[0] = two_js[0] / 2;
        two_js[1] = two_js[1] / 2;
        two_js[2] = two_js[2] / 2;
        two_js[3] = two_js[3] / 2;
    }
    // Get indices i, j from k (1-based indexing in result)
    auto [i, j] = ij_from_k(k);
    // Konvertieren zu 0-basierter Indexierung für Arrays
    int i_idx = i - 1;
    int j_idx = j - 1;
    int k_idx = k - 1;

    if (debug)
        std::cout << "i_idx" << i_idx << " " << j_idx << " " << k_idx << std::endl;

    // Überprüfen der Grenzen für sicheren Array-Zugriff
    if (i_idx < 0 || i_idx >= 3 || j_idx < 0 || j_idx >= 3 || k_idx < 0 || k_idx >= 3)
    {
        std::cout << "Invalid indices in aligned_amplitude4d: i="
                  << i << ", j=" << j << ", k=" << k << std::endl;
        // Rückgabe eines leeren Tensors bei Fehler
        return Tensor4Dcomp();
    }

    // Quadrierte Massen berechnen
    std::array<double, 4> msq = {tbs.ms[0] * tbs.ms[0],
                                 tbs.ms[1] * tbs.ms[1],
                                 tbs.ms[2] * tbs.ms[2],
                                 tbs.ms[3] * tbs.ms[3]};

    // cosθ und Wigner d-Matrix berechnen
    double cosθ = cosθij(σs, msq, k);
    // std::cout << "cosθ: " << cosθ << std::endl;
    if (debug)
        std::cout << "two_j " << two_j << " " << "cosθ " << cosθ << std::endl;
    Matrix2D d_θ = wignerd_doublearg_sign(two_j, cosθ, true);
    if (debug)
    {
        for (int i = 0; i < d_θ.size(); ++i)
        {
            for (int j = 0; j < d_θ[0].size(); ++j)
            {
                std::cout << d_θ[i][j] << "\t"; // Tab für schöne Ausrichtung
            }
            std::cout << "\n";
        }
    }

    // Spins
    std::array<int, 3> two_js_Hij = {two_j, two_js[i_idx], two_js[j_idx]};
    std::array<int, 3> two_js_HRk = {two_js[3], two_j, two_js[k_idx]};
    if (debug)
        std::cout << "two_js_Hij: " << two_js_Hij[0] << " " << two_js_Hij[1] << " " << two_js_Hij[2] << std::endl;
    if (debug)
        std::cout << "two_js_HRk: " << two_js_HRk[0] << " " << two_js_HRk[1] << " " << two_js_HRk[2] << std::endl;
    // VRk-Matrix berechnen
    std::vector<std::vector<double>> VRk;
    int vrk_dim1 = two_j + 1;
    int vrk_dim2 = two_js[k_idx] + 1;
    // VRk-Matrix initialisieren
    bool debugvrkvij = false;
    VRk.resize(vrk_dim1, std::vector<double>(vrk_dim2, 0.0));

    // std::cout << "ijk" << i << j << k << std::endl;
    for (int m1_idx = 0; m1_idx < vrk_dim1; m1_idx++)
    {
        int two_m1 = -two_j + 2 * m1_idx;
        for (int m2_idx = 0; m2_idx < vrk_dim2; m2_idx++)
        {
            int two_m2 = -two_js[k_idx] + 2 * m2_idx;

            std::array<int, 2> two_ms = {two_m1, two_m2};
            if (debug)
                std::cout << "VRk values" << two_m1 << " " << two_m2 << " " << two_js_HRk[0] << " " << two_js_HRk[1] << " " << two_js_HRk[2] << std::endl;

            complex amp = amplitude_recoupling(dc.HRk, two_ms, two_js_HRk);
            // complex phase_value = (((two_js[k_idx] - two_m2) % 4 == 0) ? complex(1.0, 0.0) : complex(-1.0, 0.0));
            //  In aligned_amplitude4dcomp für VRk
            double phase_value = phase2(two_js[k_idx] - two_m2);
            // phase_value = phaseabs(two_js[k_idx] - two_m2);

            double phase_value_old = (((two_js[k_idx] - two_m2) % 2 == 0) ? 1.0 : -1.0);

            // std::cout << (two_js[k_idx] - two_m2) << phase_value << phase_value_old << std::endl;
            if (debugvrkvij)
                std::cout << "VRk" << two_m1 << two_js[k_idx] << two_m2 << " " << amp << " " << phase_value << amp * phase_value << std::endl;
            VRk[m1_idx][m2_idx] = amp.real() * phase_value;
        }
    }

    // Nach der Berechnung von VRk
    if (debugvrkvij)
    {
        std::cout << "VRk Matrix berechnet: " << std::endl;
        for (int i = 0; i < VRk.size(); ++i)
        {
            for (int j = 0; j < VRk[0].size(); ++j)
            {

                std::cout << VRk[i][j] << "\t";
            }

            std::cout << "\n";
        }
    }

    // Vij-Matrix berechnen
    std::vector<std::vector<double>> Vij;
    int vij_dim1 = two_js[i_idx] + 1;
    int vij_dim2 = two_js[j_idx] + 1;

    // Vij-Matrix initialisieren
    Vij.resize(vij_dim1, std::vector<double>(vij_dim2, 0.0));
    if (debug)
        std::cout << "Vij dimensions: " << vij_dim1 << " x " << vij_dim2 << std::endl;
    for (int m1_idx = 0; m1_idx < vij_dim1; m1_idx++)
    {
        int two_m1 = -two_js[i_idx] + 2 * m1_idx;
        for (int m2_idx = 0; m2_idx < vij_dim2; m2_idx++)
        {
            int two_m2 = -two_js[j_idx] + 2 * m2_idx;

            std::array<int, 2> two_ms = {two_m1, two_m2};
            if (debug)
                std::cout << "Vij values" << two_m1 << " " << two_m2 << " " << two_js_Hij[0] << " " << two_js_Hij[1] << " " << two_js_Hij[2] << std::endl;
            complex amp = amplitude_recoupling(dc.Hij, two_ms, two_js_Hij);
            double phase_value = phase2(two_js[j_idx] - two_m2);
            // phase_value = phaseabs(two_js[k_idx] - two_m2);

            // phase_value = -phase2(two_js[k_idx] - two_m2);
            double phase_value_old = (((two_js[k_idx] - two_m2) % 2 == 0) ? 1.0 : -1.0);
            // std::cout << " Vij " << (two_js[k_idx] - two_m2) << phase_value << phase_value_old << std::endl;

            if (debugvrkvij)
                std::cout << "Vij" << two_m1 << two_js[k_idx] << two_m2 << amp << " " << phase_value << std::endl;
            //  In aligned_amplitude4dcomp für Vij
            // complex phase_value = (((two_js[j_idx] - two_m2) % 2 == 0) ? complex(1, 1) : complex(-1, -1));
            Vij[m1_idx][m2_idx] = amp.real() * phase_value;
            // Vij[m1_idx][m2_idx] = complex(1, 1);
        }
    }

    // Nach der Berechnung von Vij
    if (debugvrkvij)
    {
        std::cout << "Vij Matrix berechnet: " << std::endl;
        for (int i = 0; i < Vij.size(); ++i)
        {
            for (int j = 0; j < Vij[0].size(); ++j)
            {

                std::cout << Vij[i][j] << "\t";
            }

            std::cout << "\n";
        }
    }



    // Δ-Verschiebungen berechnen
    int Δ_zk = ((two_j - two_js[3] - two_js[k_idx]) / 2);
    int Δ_ij = ((two_j - two_js[i_idx] + two_js[j_idx]) / 2);
    if (debug)
        std::cout << "Δ_zk: " << Δ_zk << " Δ_ij: " << Δ_ij << std::endl;
    // try div 2
    // int Δ_zk = (two_j - two_js[3] / 2 - two_js[k_idx] / 2) / 2;
    // int Δ_ij = (two_j - two_js[i_idx] + two_js[j_idx] / 2) / 2;

    // Lineshape berechnen
    complex lineshape = dc.Xlineshape(σs[k_idx]);

    // Dimensionen für den Ergebnistensor
    std::vector<int> dims = {two_js[0] + 1, two_js[1] + 1, two_js[2] + 1, two_js[3] + 1};

    // F0 mit Nullen initialisieren (Ergebnistensor)
    Tensor4Dcomp F0(dims[0], std::vector<std::vector<std::vector<complex>>>(
                                 dims[1], std::vector<std::vector<complex>>(
                                              dims[2], std::vector<complex>(dims[3], 0.0))));

    // Dimensionen für Berechnungstensor
    Tensor4Dcomp F(dims[i_idx], std::vector<std::vector<std::vector<complex>>>(
                                    dims[j_idx], std::vector<std::vector<complex>>(
                                                     dims[k_idx], std::vector<complex>(dims[3], 0.0))));

    // Hilfsfunktion zum Begrenzen von Indizes
    auto pad = [](int idx, int max_size)
    {
        return std::max(0, std::min(idx, max_size - 1));
    };

    for (int _i = 0; _i < F.size(); _i++)
    {
        for (int _j = 0; _j < F[0].size(); _j++)
        {
            for (int _k = 0; _k < F[0][0].size(); _k++)
            {
                for (int _z = 0; _z < F[0][0][0].size(); _z++)
                {

                    // Verwenden Sie genau die gleiche Indexberechnung wie in Julia
                    int d2_idx1 = _z + _k + Δ_zk;
                    int d2_idx2 = _i - _j + Δ_ij;
                    int vrk_idx1 = pad(_z + _k + Δ_zk, two_j + 1);
                    int d_idx1 = pad(_z + _k + Δ_zk, two_j + 1);
                    int d_idx2 = pad(_i - _j + Δ_ij, two_j + 1);
                    int idx1 = _z + _k + Δ_zk;
                    int idx2 = _i - _j + Δ_ij;
                    vrk_idx1 = (d2_idx1 < 0 || d2_idx1 > two_j + 1) ? -1 : d2_idx1;
                    d_idx1 = vrk_idx1;
                    d_idx2 = (d2_idx2 < 0 || d2_idx2 > two_j + 1) ? -1 : d2_idx2;

                    if (debug)
                        std::cout << "result " << _i << _j << _k << _z << " " << vrk_idx1 << " " << " " << d_idx1 << " " << d_idx2 << " " << std::endl;

                    // Stellen Sie sicher, dass die Indizes innerhalb der Grenzen sind
                    if (d_idx1 >= 0 && d_idx2 >= 0 &&
                        vrk_idx1 < VRk.size() && _k < VRk[0].size() &&
                        d_idx1 < d_θ.size() && d_idx2 < d_θ[0].size() &&
                        _i < Vij.size() && _j < Vij[0].size())
                    {
                        auto res = VRk[vrk_idx1][_k] *
                                   d_θ[d_idx1][d_idx2] *
                                   Vij[_i][_j];
                        if (debug)
                            std::cout << "result " << _i << _j << _k << _z << " " << vrk_idx1 << " " << " " << d_idx1 << " " << d_idx2 << " " << res << std::endl;
                        if (debug)
                            std::cout << "VRk " << VRk[vrk_idx1][_k] << "d_θ " << d_θ[d_idx1][d_idx2] << "Vij " << Vij[_i][_j] << std::endl;
                        // KRITISCH: Exakt die gleiche Berechnung wie in Julia verwenden
                        F[_i][_j][_k][_z] +=
                            VRk[vrk_idx1][_k] *
                            d_θ[d_idx1][d_idx2] *
                            Vij[_i][_j];
                    }
                    // else
                    //{
                    //     F[_i][_j][_k][_z] = 0.0;
                    // }
                }
            }
        }
    }

    // print F
    if (debug)
    {
        std::cout << "Tensor4D F : " << std::endl;

        for (int i = 0; i < F.size(); ++i)
        {
            for (int j = 0; j < F[0].size(); ++j)
            {
                for (int k = 0; k < F[0][0].size(); ++k)
                {
                    for (int z = 0; z < F[0][0][0].size(); ++z)
                    {
                        std::cout << F[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                    }
                }
            }
            std::cout << "\n";
        }
    }

    // Normalisierung und Lineshape anwenden
    double d_norm = std::sqrt(two_j + 1.0);
    if (debug)
        std::cout << d_norm << " " << lineshape << std::endl;

    for (auto &dim1 : F)
    {
        for (auto &dim2 : dim1)
        {
            for (auto &dim3 : dim2)
            {
                for (auto &val : dim3)
                {
                    val *= d_norm * lineshape;
                }
            }
        }
    }
    if (debug)
        std::cout << d_norm << " " << lineshape << std::endl;

    // Zurück zur ursprünglichen Reihenfolge permutieren
    for (int _0 = 0; _0 < dims[0]; _0++)
    {
        for (int _1 = 0; _1 < dims[1]; _1++)
        {
            for (int _2 = 0; _2 < dims[2]; _2++)
            {
                for (int _3 = 0; _3 < dims[3]; _3++)
                {
                    // Indizes gemäß Permutation zuordnen
                    int perm_i = (i_idx == 0) ? _0 : (i_idx == 1) ? _1
                                                 : (i_idx == 2)   ? _2
                                                                  : _3;
                    int perm_j = (j_idx == 0) ? _0 : (j_idx == 1) ? _1
                                                 : (j_idx == 2)   ? _2
                                                                  : _3;
                    int perm_k = (k_idx == 0) ? _0 : (k_idx == 1) ? _1
                                                 : (k_idx == 2)   ? _2
                                                                  : _3;
                    int perm_z = 3; // Parent index is always last (3)

                    // Überprüfen der Grenzen für sicheren Array-Zugriff
                    if (perm_i < dims[i_idx] && perm_j < dims[j_idx] &&
                        perm_k < dims[k_idx] && _3 < dims[3])
                    {

                        F0[_0][_1][_2][_3] = F[perm_i][perm_j][perm_k][_3]; // * d_norm * lineshape;
                    }
                }
            }
        }
    }

    if (debug)
        std::cout << "Tensor4D F0 : " << std::endl;
    if (debug)
    {
        for (int i = 0; i < F0.size(); ++i)
        {
            for (int j = 0; j < F0[0].size(); ++j)
            {
                for (int k = 0; k < F0[0][0].size(); ++k)
                {
                    for (int z = 0; z < F0[0][0][0].size(); ++z)
                    {
                        std::cout << F0[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                    }
                }
            }
            std::cout << "\n";
        }
    }

    return F0;
}

// Implementation of amplitude function that returns a Tensor4D
Tensor4Dcomp ThreeBodyDecays::amplitude4dcomp(const DecayChain &dc,
                                              const MandelstamTuple &σs,
                                              const int &k_amp,
                                              std::vector<int> refζ)
{
    int k = dc.k;
    const auto &tbs = dc.tbs;
    int two_j = dc.two_j;
    // const auto &two_js = tbs.two_js;
    const auto &ms = tbs.ms;

    if (refζ == std::vector<int>{-1, -1, -1, -1})
    {
        // if (k_amp == 0)
        //{
        //  Setze die Referenz-Spin-Quantenzahlen auf 0
        //    refζ = {0, 0, 0, 0};
        //}
        if (k_amp == 1)
        {
            // Setze die Referenz-Spin-Quantenzahlen auf 1
            refζ = {1, 1, 1, 1};
        }
        else if (k_amp == 2)
        {
            // Setze die Referenz-Spin-Quantenzahlen auf 2
            refζ = {2, 2, 2, 2};
        }
        else if (k_amp == 3)
        {
            // Setze die Referenz-Spin-Quantenzahlen auf 3
            refζ = {3, 3, 3, 3};
        }
        else
        {
            std::cout << "Invalid k_amp value: " << k_amp << std::endl;
            return Tensor4Dcomp();
        }
    }
    const std::vector<int> refζs = refζ;

    // div 2
    auto two_js = tbs.two_js;

    if (div2)
    {
        two_js[0] = two_js[0] / 2;
        two_js[1] = two_js[1] / 2;
        two_js[2] = two_js[2] / 2;
        two_js[3] = two_js[3] / 2;
    }
    // std::array<int, 4> two_js = {two_jst[0]*2, two_jst[1]*2, two_jst[2]*2, two_jst[3]*2};

    if (debug)
        std::cout << two_js[0] << " " << two_js[1] << " " << two_js[2] << " " << two_js[3] << std::endl;

    // Calculate squared masses
    std::array<double, 4> msq = {ms[0] * ms[0],
                                 ms[1] * ms[1],
                                 ms[2] * ms[2],
                                 ms[3] * ms[3]};

    // Get aligned amplitude
    auto F0 = aligned_amplitude4dcomp(dc, σs);

    // Print the dimensions of F0

    if (debug)
    {
        std::cout << "F0 dimensions: "
                  << F0.size() << " x "
                  << (F0.empty() ? 0 : F0[0].size()) << " x "
                  << (F0.empty() || F0[0].empty() ? 0 : F0[0][0].size()) << " x "
                  << (F0.empty() || F0[0].empty() || F0[0][0].empty() ? 0 : F0[0][0][0].size())
                  << std::endl;

        for (int i = 0; i < F0.size(); ++i)
        {
            for (int j = 0; j < F0[0].size(); ++j)
            {
                for (int k = 0; k < F0[0][0].size(); ++k)
                {
                    for (int z = 0; z < F0[0][0][0].size(); ++z)
                    {
                        std::cout << F0[i][j][k][z] << "\t"; // Tab für schöne Ausrichtung
                    }
                }
            }
            std::cout << "\n";
        }
    }

    if (debug)
        std::cout
            << msq[0] << " " << msq[1] << " " << msq[2] << " " << msq[3] << std::endl;
    // Calculate alignment rotations
    std::vector<Matrix2D> d_ζs(4);
    for (size_t l = 0; l < 4; ++l)
    {
        int _two_j = two_js[l];
        int _refζ = refζs[l];
        int julia_mod_result = (l + 1) % 4; // Dies funktioniert für l=0,1,2,3
        if (julia_mod_result == 0)
            julia_mod_result = 4; // Für den Fall l=3
        if (julia_mod_result == 4)
            julia_mod_result = 0; // Für den Fall l=3

        std::unique_ptr<AbstractWignerRotation> _w = wr(k, _refζ, julia_mod_result);
        double _cosζ = _w->cos_zeta(σs, msq);
        if (debug)
            std::cout << _cosζ << std::endl;
        d_ζs[l] = wignerd_doublearg_sign(_two_j, _cosζ, _w->is_positive());
    }

    for (size_t l = 0; l < 4; ++l)
    {
        if (debug)
            std::cout << "d_ζs[" << l << "]" << std::endl;
        for (int i = 0; i < d_ζs[l].size(); ++i)
        {
            for (int j = 0; j < d_ζs[l][0].size(); ++j)
            {
                if (debug)
                    std::cout << d_ζs[l][i][j] << "\t"; // Tab für schöne Ausrichtung
            }
            if (debug)
                std::cout << "\n";
        }
    }

    // Alias the d_ζs for clarity
    const auto &D1 = d_ζs[0];
    const auto &D2 = d_ζs[1];
    const auto &D3 = d_ζs[2];
    const auto &D0 = d_ζs[3];

    // Create a new tensor with the same dimensions as F0
    std::vector<int> dims;
    for (int i = 0; i < 4; i++)
    {
        dims.push_back(two_js[i] + 1);
    }

    Tensor4Dcomp F(dims[0], std::vector<std::vector<std::vector<complex>>>(
                                dims[1], std::vector<std::vector<complex>>(
                                             dims[2], std::vector<complex>(dims[3], 0.0))));

    // Calculate the tensor contraction (equivalent to @tullio in Julia)
    // @tullio F[_i, _j, _k, _z] = D0[_z, _z′] * F0[_i′, _j′, _k′, _z′] * D1[_i′, _i] * D2[_j′, _j] * D3[_k′, _k]
    for (int _i = 0; _i < dims[0]; ++_i)
    {
        for (int _j = 0; _j < dims[1]; ++_j)
        {
            for (int _k = 0; _k < dims[2]; ++_k)
            {
                for (int _z = 0; _z < dims[3]; ++_z)
                {
                    // Sum over all primed indices
                    for (int _i_prime = 0; _i_prime < dims[0]; ++_i_prime)
                    {
                        if (_i_prime >= (int)D1.size() || _i >= (int)D1[0].size())
                            continue;

                        for (int _j_prime = 0; _j_prime < dims[1]; ++_j_prime)
                        {
                            if (_j_prime >= (int)D2.size() || _j >= (int)D2[0].size())
                                continue;

                            for (int _k_prime = 0; _k_prime < dims[2]; ++_k_prime)
                            {
                                if (_k_prime >= (int)D3.size() || _k >= (int)D3[0].size())
                                    continue;

                                for (int _z_prime = 0; _z_prime < dims[3]; ++_z_prime)
                                {
                                    if (_z_prime >= (int)D0.size() || _z >= (int)D0[0].size())
                                        continue;

                                    if (debug)
                                        std::cout << _i << " " << _j << " " << _k << " " << _z << " " << _i_prime << " " << _j_prime << " " << _k_prime << " " << _z_prime << std::endl;

                                    F[_i][_j][_k][_z] +=
                                        D0[_z][_z_prime] *
                                        F0[_i_prime][_j_prime][_k_prime][_z_prime] *
                                        D1[_i_prime][_i] *
                                        D2[_j_prime][_j] *
                                        D3[_k_prime][_k];
                                    if (debug)
                                        std::cout << D0[_z][_z_prime] << " " << F0[_i_prime][_j_prime][_k_prime][_z_prime] << " " << D1[_i_prime][_i] << " " << D2[_j_prime][_j] << " " << D3[_k_prime][_k] << std::endl;
                                    if (debug)
                                        std::cout << F[_i][_j][_k][_z] << std::endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return F;
}

// Calculate amplitude for specific helicity values
complex ThreeBodyDecays::amplitude(const DecayChain &dc,
                                   const MandelstamTuple &σs,
                                   const std::vector<int> &two_λs,
                                   const int &k_amp,
                                   const std::vector<int> &refζs)
{
    // Get full 4D tensor of amplitudes for all helicity combinations
    Tensor4Dcomp F0 = amplitude4dcomp(dc, σs, k_amp, refζs);

    // Calculate indices from helicity values
    std::vector<int> indices(4);
    for (int i = 0; i < 4; ++i)
    {
        // Match Julia's div(_two_j + _two_λ, 2) + 1, but adjust for 0-based indexing
        indices[i] = (dc.tbs.two_js[i] + two_λs[i]) / 2;
    }

    // Check if indices are within bounds
    if (indices[0] >= 0 && indices[0] < (int)F0.size() &&
        indices[1] >= 0 && indices[1] < (int)F0[0].size() &&
        indices[2] >= 0 && indices[2] < (int)F0[0][0].size() &&
        indices[3] >= 0 && indices[3] < (int)F0[0][0][0].size())
    {

        return F0[indices[0]][indices[1]][indices[2]][indices[3]];
    }

    // Return 0 if indices are out of bounds
    return 0.0;
}

double ThreeBodyDecays::intensity(const DecayChain &dc, const MandelstamTuple &σs, const int &k_amp, const complex weight, const std::vector<int> refζs)
{
    // Get the 4D amplitude tensor
    auto amp = amplitude4dcomp(dc, σs, k_amp, refζs);

    // Sum squared amplitudes (similar to Julia's sum(abs2, ...))
    double total_intensity = 0.0;

    for (const auto &dim1 : amp)
    {
        for (const auto &dim2 : dim1)
        {
            for (const auto &dim3 : dim2)
            {
                for (const auto &val : dim3)
                {
                    total_intensity += val.real() * val.real() * weight.real() +
                                       val.imag() * val.imag() * weight.imag();
                }
            }
        }
    }

    return total_intensity;
}

// Implementation of cos_zeta for TrivialWignerRotation
double TrivialWignerRotation::cos_zeta(const std::array<double, 3> &sigma,
                                       const std::array<double, 4> &ms2) const
{
    return 1.0; // Trivial rotation always returns 1
}

// Implementation of cos_zeta for Arg0WignerRotation
double Arg0WignerRotation::cos_zeta(const std::array<double, 3> &sigma,
                                    const std::array<double, 4> &ms2) const
{
    auto [i, j, k] = ijk(*this);
    // Keep 1-based indexing like Julia
    i--;
    j--;
    k--; // Convert to 0-based indexing for array access

    double s = ms2[3]; // Parent mass squared
    double EE4s = (s + ms2[i] - sigma[i]) * (s + ms2[k] - sigma[k]);
    double pp4s = std::sqrt(Kallen(s, ms2[i], sigma[i]) *
                            Kallen(s, ms2[k], sigma[k]));
    pp4s = (pp4s < 0) ? 0.0 : pp4s; // Handle numerical errors
    double rest = sigma[j] - ms2[i] - ms2[k];
    if (debug)
        std::cout << "wr0" << std::endl;

    return (EE4s - 2 * s * rest) / pp4s;
}

// Implementation of cos_zeta for Arg2WignerRotation
double Arg2WignerRotation::cos_zeta(const std::array<double, 3> &sigma,
                                    const std::array<double, 4> &ms2) const
{
    auto [i, j, k] = ijk(*this);
    if (debug)
        std::cout << "wr2" << std::endl;
    // Keep 1-based indexing like Julia
    i--;
    j--;
    k--; // Convert to 0-based indexing for array access
    if (debug)
        std::cout << i << " " << j << " " << k << std::endl;
    // Swap i and j if not even (like Julia)
    if (!is_even())
    {
        std::swap(i, j);
    }
    if (debug)
        std::cout << i << " " << j << std::endl;

    // Handle massless case
    if (std::abs(ms2[k]) < 1e-10)
    {
        return 1.0;
    }

    double s = ms2[3]; // Parent mass squared
    double EE4mksq = (s + ms2[k] - sigma[k]) * (sigma[i] - ms2[k] - ms2[j]);
    double pp4mksq = std::sqrt(Kallen(s, ms2[k], sigma[k]) *
                               Kallen(ms2[k], ms2[j], sigma[i]));
    pp4mksq = (pp4mksq < 0) ? 0.0 : pp4mksq; // Handle numerical errors
    double rest = sigma[j] - s - ms2[j];
    if (debug)
        std::cout << s << " " << EE4mksq << " " << pp4mksq << " " << rest << " " << (2 * ms2[k] * rest + EE4mksq) / pp4mksq << std::endl;
    return (2 * ms2[k] * rest + EE4mksq) / pp4mksq;
}

// Implementation of cos_zeta for Arg3WignerRotation
double Arg3WignerRotation::cos_zeta(const std::array<double, 3> &sigma,
                                    const std::array<double, 4> &ms2) const
{
    auto [i, j, k] = ijk(*this);
    if (debug)
        std::cout << "wr3" << std::endl;
    // Keep 1-based indexing like Julia
    i--;
    j--;
    k--; // Convert to 0-based indexing for array access

    // Handle massless case
    if (std::abs(ms2[k]) < 1e-10)
    {
        return 1.0;
    }

    double s = ms2[3]; // Parent mass squared
    double EE4m1sq = (sigma[i] - ms2[j] - ms2[k]) *
                     (sigma[j] - ms2[k] - ms2[i]);
    double pp4m1sq = std::sqrt(Kallen(sigma[i], ms2[j], ms2[k]) *
                               Kallen(sigma[j], ms2[k], ms2[i]));
    pp4m1sq = (pp4m1sq < 0) ? 0.0 : pp4m1sq; // Handle numerical errors
    double rest = ms2[i] + ms2[j] - sigma[k];
    return (2 * ms2[k] * rest + EE4m1sq) / pp4m1sq;
}

// Implementierung des Konstruktors für ThreeBodyParities
ThreeBodyParities::ThreeBodyParities(char P1, char P2, char P3, char P0) : P1_(P1), P2_(P2), P3_(P3), P0_(P0)
{
}

// Implementierung des Paritätsoperators
char ThreeBodyParities::operator*(const ThreeBodyParities &other) const
{
    char result = (P0_ == other.P0_) ? '+' : '-';
    return result;
}

// Implementation of SpinParity constructor
SpinParity::SpinParity(const std::string &jp)
{
    if (jp.empty())
    {
        throw std::invalid_argument("jp string cannot be empty");
    }

    // Get the parity character (last character)
    p_ = jp.back();

    // Parse the spin value (rest of the string)
    std::string spin_str = jp.substr(0, jp.size() - 1);
    str = spin_str; // Store the original string for parsing

    // Check if spin string is empty
    if (spin_str.empty())
    {
        throw std::invalid_argument("Spin value cannot be empty");
    }

    // Check if spin is a fraction
    size_t slash_pos = spin_str.find('/');
    if (slash_pos != std::string::npos)
    {
        // Handle fraction
        std::string num_str = spin_str.substr(0, slash_pos);
        std::string denom_str = spin_str.substr(slash_pos + 1);

        if (num_str.empty() || denom_str.empty())
        {
            throw std::invalid_argument("Invalid fraction format");
        }

        try
        {
            int numerator = std::stoi(num_str);
            int denominator = std::stoi(denom_str);
            if (denominator == 0)
            {
                throw std::invalid_argument("Division by zero");
            }

            // For fractions, double the numerator
            two_j_ = (2 * numerator) / denominator;
        }
        catch (const std::exception &e)
        {
            throw std::invalid_argument("Invalid fraction: " + spin_str);
        }
    }
    else
    {
        // Integer spin - multiply by 2 for doubled representation
        try
        {
            int spin = std::stoi(spin_str);
            two_j_ = spin * 2; // This doubles the integer spin
        }
        catch (const std::exception &e)
        {
            throw std::invalid_argument("Invalid spin: " + spin_str);
        }
    }
}

std::vector<std::array<int, 2>> possible_ls(
    const SpinParity &jp1,
    const SpinParity &jp2,
    const SpinParity &jp)
{
    std::vector<std::array<int, 2>> two_ls;
    if (debugls)
        std::cout << "possible_ls: " << jp1.get_two_j() << " " << jp2.get_two_j() << " " << jp.get_two_j() << std::endl;
    // Loop through possible s values with step 2
    for (int two_s = std::abs(jp1.get_two_j() - jp2.get_two_j());
         two_s <= jp1.get_two_j() + jp2.get_two_j();
         two_s += 2)
    {

        // Loop through possible l values with step 2
        for (int two_l = std::abs(jp.get_two_j() - two_s);
             two_l <= jp.get_two_j() + two_s;
             two_l += 2)
        {

            // Check parity condition - match the Julia logic exactly
            int negative_count = 0;
            if (jp1.get_p() == '-')
                negative_count++;
            if (jp2.get_p() == '-')
                negative_count++;
            if (jp.get_p() == '-')
                negative_count++;

            char combined_parity = (negative_count % 2 == 0) ? '+' : '-';
            bool odd_l = ((two_l / 2) % 2) != 0;
            char expected_parity = odd_l ? '-' : '+';

            if (combined_parity == expected_parity)
            {
                two_ls.push_back({two_l, two_s});
            }
        }
    }

    // Sort by the first element (l value)
    std::sort(two_ls.begin(), two_ls.end(),
              [](const std::array<int, 2> &a, const std::array<int, 2> &b)
              {
                  return a[0] < b[0];
              });

    return two_ls;
}

// Calculate possible ls couplings for i,j particles
std::vector<std::array<int, 2>> possible_ls_ij(
    const SpinParity &jp,
    const std::array<int, 4> &two_js,
    const ThreeBodyParities &Ps,
    int k)
{
    auto [i, j] = ij_from_k(k);
    i--;
    j--; // Convert to 0-based indexing for array access

    // Select parity based on i and j indices
    char p_i, p_j;
    switch (i)
    {
    case 0:
        p_i = Ps.get_P1();
        break;
    case 1:
        p_i = Ps.get_P2();
        break;
    case 2:
        p_i = Ps.get_P3();
        break;
    default:
        p_i = '+'; // Default value
    }

    switch (j)
    {
    case 0:
        p_j = Ps.get_P1();
        break;
    case 1:
        p_j = Ps.get_P2();
        break;
    case 2:
        p_j = Ps.get_P3();
        break;
    default:
        p_j = '+'; // Default value
    }

    // Create SpinParity objects for particles i and j div 2
    SpinParity jp_i(std::to_string(two_js[i]) + "/2" + p_i);
    SpinParity jp_j(std::to_string(two_js[j]) + "/2" + p_j);

    return possible_ls(jp_i, jp_j, jp);
}

// Calculate possible LS couplings for resonance and k particle
std::vector<std::array<int, 2>> possible_ls_Rk(
    const SpinParity &jp,
    const std::array<int, 4> &two_js,
    const ThreeBodyParities &Ps,
    int k)
{
    k--; // Convert to 0-based indexing for array access

    // Select parity for k
    char p_k;
    switch (k)
    {
    case 0:
        p_k = Ps.get_P1();
        break;
    case 1:
        p_k = Ps.get_P2();
        break;
    case 2:
        p_k = Ps.get_P3();
        break;
    default:
        p_k = '+'; // Default value
    }

    if (debugls)
        std::cout << "possible_ls_Rk: " << two_js[k] << " " << p_k << two_js[3] << " " << Ps.get_P0() << std::endl;

    // Create SpinParity objects for particle k and the parent div 2
    SpinParity jp_k(std::to_string(two_js[k]) + "/2" + p_k);
    SpinParity jp_0(std::to_string(two_js[3]) + "/2" + Ps.get_P0());

    return possible_ls(jp, jp_k, jp_0);
}

// Implementation of possible_lsLS function matching Julia's implementation
std::vector<LSCoupling> possible_lsLS(
    const SpinParity &jp,
    const std::array<int, 4> &two_js,
    const ThreeBodyParities &Ps,
    int k)
{
    // Get the list of possible (l,s) values
    std::vector<std::array<int, 2>> lsv = possible_ls_ij(jp, two_js, Ps, k);
    if (debugls)
        std::cout << "lsv size: " << lsv.size() << " with " << two_js[0] << " " << two_js[1] << " " << two_js[2] << " " << two_js[3] << " " << k << std::endl;
    for (const auto &ls : lsv)
    {
        if (debugls)
            std::cout << "  ls: " << ls[0] << ", " << ls[1] << std::endl;
    }

    // Get the list of possible (L,S) values
    std::vector<std::array<int, 2>> LSv = possible_ls_Rk(jp, two_js, Ps, k);
    if (debugls)
        std::cout << "LSv size: " << LSv.size() << std::endl;
    for (const auto &LS : LSv)
    {
        if (debugls)
            std::cout << "  LS: " << LS[0] << ", " << LS[1] << std::endl;
    }
    // Create the Cartesian product (like Julia's Iterators.product)
    std::vector<LSCoupling> result;
    for (const auto &ls : lsv)
    {
        for (const auto &LS : LSv)
        {
            result.push_back(LSCoupling(ls, LS));
        }
    }

    return result;
}

// Updated createDecayChainLS implementation that uses the LSCoupling structure
std::shared_ptr<DecayChain> createDecayChainLS(
    int k,
    std::function<std::complex<double>(double)> Xlineshape,
    const std::string &jp,
    const ThreeBodyParities &Ps,
    const ThreeBodySystem &tbs)
{
    // Parse spin-parity
    SpinParity SP(jp);
    if (debug)
        std::cout << "Creating DecayChain with spin-parity: " << jp << std::endl;
    int two_j = SP.get_two_j();
    if (debug)
        std::cout << "Creating DecayChain with k=" << k << ", two_j=" << two_j << std::endl;
    RecouplingLS Hij;
    RecouplingLS HRk;

    // Calculate possible LS couplings
    std::vector<LSCoupling> two_lsLS = possible_lsLS(SP, tbs.two_js, Ps, k);

    if (two_lsLS.empty())
    {
        throw std::runtime_error("No possible LS couplings found for the given configuration");
    }

    // Sort by two_LS[0] (first element of LS)
    std::sort(two_lsLS.begin(), two_lsLS.end(),
              [](const LSCoupling &a, const LSCoupling &b)
              {
                  return a.two_LS[0] < b.two_LS[0];
              });

    // Use the first (smallest L) coupling
    const auto &two_ls = two_lsLS[0].two_ls;
    const auto &two_LS = two_lsLS[0].two_LS;

    if (debugls)
        std::cout << "LS coupling: " << two_ls[0] << " " << two_ls[1] << std::endl;
    if (debugls)
        std::cout << "L coupling: " << two_LS[0] << " " << two_LS[1] << std::endl;

    HRk = createRecouplingFunction(RecouplingType::LSRecoupling, two_LS, false);
    Hij = createRecouplingFunction(RecouplingType::LSRecoupling, two_ls, false);
    // Create and return the DecayChain
    return std::make_shared<DecayChain>(
        k,          // k-value
        two_j,      // two_j
        Xlineshape, // Lineshape function
        HRk,        // HRk recoupling
        Hij,        // Hij recoupling
        tbs         // ThreeBodySystem
    );
}

complex cg_double(int j1, int j2, int j3, int m1, int m2, int m3)
{
    // Placeholder for the actual Clebsch-Gordan coefficient calculation
    // This should be replaced with the actual implementation
    return complex(0.0, 0.0);
}

// Helper functions
int one(int /*x*/) { return 1; }
int zero(int /*x*/) { return 0; }

double jls_coupling(int two_j1, int two_λ1, int two_j2, int two_λ2, int two_j, int two_l, int two_s)
{
    int T1 = one(two_λ1);
    double sqrt_factor = std::sqrt(static_cast<double>(two_l * T1 + 1) / (two_j * T1 + 1));
    // two_j1 = two_j1 / 2;
    // two_j2 = two_j2 / 2;
    //   two_j = two_j / 2;
    //   two_λ1 = two_λ1 / 2;
    //   two_λ2 = two_λ2 / 2;
    //   two_l = two_l / 2;
    //   two_s = two_s / 2;

    double cg1 = -clebschgordan_doublearg(two_j1, two_λ1, two_j2, -two_λ2, two_s, two_λ1 - two_λ2);
    double cg2 = -clebschgordan_doublearg(two_l, 0, two_s, two_λ1 - two_λ2, two_j, two_λ1 - two_λ2);
    // std::cout << " " << two_j1 << ", " << two_λ1 << ", " << two_j2 << ", " << -two_λ2 << ", " << two_s << ", " << two_λ1 - two_λ2 << std::endl;
    // std::cout << "cg1: " << cg1 << std::endl;
    // std::cout << " " << two_l << ", " << 0 << ", " << two_s << ", " << two_λ1 - two_λ2 << ", " << two_j << ", " << two_λ1 - two_λ2 << std::endl;
    // std::cout << "cg2: " << cg2 << std::endl;
    // std::cout << sqrt_factor << std::endl;
    // std::cout << sqrt_factor * cg1 * cg2 << std::endl;
    double res = sqrt_factor * cg1 * cg2;
    if (debug)
        std::cout << "jls_coupling: " << two_j1 << ", " << two_λ1 << ", " << two_j2 << ", " << two_λ2 << ", " << two_j << ", " << two_l << ", " << two_s << " res: " << res << std::endl;

    return res;
}

RecouplingLS createRecouplingFunction(
    RecouplingType type,
    const std::array<int, 2> &params,
    bool parityPhase) // Parameter für ParityRecoupling
{
    switch (type)
    {
    case RecouplingType::NoRecoupling:
    {
        // NoRecoupling verwenden
        NoRecoupling noRecoupling(params[0], params[1]);

        return [noRecoupling](const std::array<int, 2> &two_ms, const std::array<int, 3> &two_js) -> double
        {
            //std::cout << noRecoupling.get_two_λa() << noRecoupling.get_two_λb() << " " << two_ms[0] << two_ms[1] << std::endl;

            // Exakte Implementierung der Julia-Logik: (cs.two_λa == two_λa) * (cs.two_λb == two_λb)
            if (noRecoupling.get_two_λa() == two_ms[0] && noRecoupling.get_two_λb() == two_ms[1])
            {
                return 1.0;
            }
            return 0.0;
        };
    }
    case RecouplingType::ParityRecoupling:
    {
        // Parameter extrahieren
        int two_λa = params[0];
        int two_λb = params[1];
        bool ηηηphaseisplus = parityPhase;

        return [two_λa, two_λb, ηηηphaseisplus](
                   const std::array<int, 2> &two_ms,
                   const std::array<int, 3> &two_js) -> double
        {
            // Erste Bedingung: cs.two_λa == two_λa and cs.two_λb == two_λb
            if (two_λa == two_ms[0] && two_λb == two_ms[1])
            {
                return 1.0;//*ηηηphaseisplus;
            }

            // Zweite Bedingung: cs.two_λa == -two_λa and cs.two_λb == -two_λb
            // Dies prüft, ob der Parameter gleich dem negativen des übergebenen Werts ist
            if (two_λa == -two_ms[0] && two_λb == -two_ms[1])
            {
                // Exakt die Julia-Logik: 2 * ηηηphaseisplus - 1 (= 1 wenn true, -1 wenn false)
                return 2 * ηηηphaseisplus - 1;
            }

            return 0.0;
        };
    }

    case RecouplingType::LSRecoupling:
    {
        int two_l = 0; // You may need to pass two_l into the params or as part of two_js
        int two_s = 0; // Same for two_s
        two_l = params[0];
        two_s = params[1];

        return [two_l, two_s](const std::array<int, 2> &two_ms, const std::array<int, 3> &two_js) -> complex
        {
            int two_j1 = two_js[1];
            int two_j2 = two_js[2];
            int two_j = two_js[0];
            int two_λ1 = two_ms[0];
            int two_λ2 = two_ms[1];

            // Call the jls_coupling function that matches the Julia implementation
            double result = jls_coupling(two_j1, two_λ1, two_j2, two_λ2, two_j, two_l, two_s);
            return result;
        };
    }
    // Weitere Recoupling-Typen hier hinzufügen
    default:
    {
        // Standardfall ist NoRecoupling mit (0,0)
        return [](const std::array<int, 2> &two_ms, const std::array<int, 3> &two_js) -> double
        {
            if (two_ms[0] == 0 && two_ms[1] == 0)
            {
                return 1.0;
            }
            return 0.0;
        };
    }
    }
}

std::shared_ptr<DecayChain> createDecayChainCoupling(
    int k,
    std::function<std::complex<double>(double)> Xlineshape,
    const std::string &jp,
    const ThreeBodySystem &tbs,
    RecouplingType HRkType,
    const std::array<int, 2> &HRkParams,
    bool HRkParityPhase,
    RecouplingType HijType,
    const std::array<int, 2> &HijParams,
    bool HijParityPhase)
{
    // Parse spin-parity
    SpinParity SP(jp);
    if (debug)
        std::cout << "Creating DecayChain with spin-parity: " << jp << std::endl;
    int two_j = SP.get_two_j();
    // std::cout << "Creating DecayChain with k=" << k << ", two_j=" << two_j << std::endl;
    RecouplingLS Hij;
    RecouplingLS HRk;

    // Auto-calculate LS couplings if requested

    // Use provided recoupling parameters
    Hij = createRecouplingFunction(HijType, HijParams, HijParityPhase);
    HRk = createRecouplingFunction(HRkType, HRkParams, HRkParityPhase);

    // Create and return the DecayChain
    return std::make_shared<DecayChain>(
        k,          // k-value
        two_j,      // two_j
        Xlineshape, // Lineshape function
        HRk,        // HRk recoupling
        Hij,        // Hij recoupling
        tbs         // ThreeBodySystem
    );
}

std::vector<std::shared_ptr<DecayChain>> createDecayChainsLS(
    int k,
    std::function<complex(double)> Xlineshape,
    const std::string &jp,
    const ThreeBodyParities &Ps,
    const ThreeBodySystem &tbs)
{
    // Parse spin-parity
    SpinParity SP(jp);
    int two_j = SP.get_two_j();

    // Calculate possible LS couplings
    std::vector<LSCoupling> two_lsLS = possible_lsLS(SP, tbs.two_js, Ps, k);

    if (two_lsLS.empty())
    {
        throw std::runtime_error("No possible LS couplings found for the given configuration");
    }

    // Create a decay chain for each possible LS coupling
    std::vector<std::shared_ptr<DecayChain>> result;
    for (const auto &coupling : two_lsLS)
    {
        const auto &two_ls = coupling.two_ls;
        const auto &two_LS = coupling.two_LS;

        RecouplingLS HRk = createRecouplingFunction(RecouplingType::LSRecoupling, two_LS, false);
        RecouplingLS Hij = createRecouplingFunction(RecouplingType::LSRecoupling, two_ls, false);

        // Create the DecayChain
        result.push_back(std::make_shared<DecayChain>(
            k,          // k-value
            two_j,      // two_j
            Xlineshape, // Lineshape function
            HRk,        // HRk recoupling
            Hij,        // Hij recoupling
            tbs         // ThreeBodySystem
            ));
    }

    return result;
}
