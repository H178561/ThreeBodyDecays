#include "include/ThreeBodyDecays.hh"
// #include "ThreeBodyUtilities.hh"
#include "include/ClebschGordan.hh"
#include "include/FormFactors.hh"
#include "include/jacobi.hh"
#include "include/WignerFunctions.hh" // Include for Wigner D-functions
#include "include/SpinParity.hh"      // Include for spin-parity definitions
#include "include/WignerRotation.hh"  // Include for Wigner rotation operations

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <array>
#include <unistd.h> // For getcwd()

using complex = std::complex<double>;
bool debug = false;
bool debugls = false;
bool div2 = false;

/**
 * @brief Calculates the Kallen (triangle) function
 *
 * Computes λ(x,y,z) = x² + y² + z² - 2xy - 2yz - 2zx
 * This function appears in relativistic kinematics for three-body decays.
 *
 * @param x First squared mass
 * @param y Second squared mass
 * @param z Third squared mass
 * @return The value of the Kallen function
 */
double Kallen(double x, double y, double z)
{
    return x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * z * x;
}

/**
 * @brief Checks if two double values are approximately equal
 *
 * @param a First value
 * @param b Second value
 * @param epsilon Maximum allowed difference
 * @return true if values are within epsilon of each other, false otherwise
 */
bool approx_equal(double a, double b, double epsilon)
{
    return std::abs(a - b) < epsilon;
}

/**
 * @brief Calculates a phase factor based on integer value
 *
 * Returns a phase that depends on the parity of the value:
 * If |value % 4| == 2, phase is -1, otherwise +1
 * The sign is inverted if value > 0
 *
 * @param value Integer input
 * @return Phase factor (+1 or -1)
 */
double phase(int value)
{
    int phase = 1;
    if (std::abs(value % 4) == 2)
    {
        phase = -1;
    }
    return (value > 0) ? -phase : phase;
}

/**
 * @brief Simplified phase factor calculation
 *
 * Returns -1 if value % 4 == 2, otherwise returns 1
 *
 * @param value Integer input
 * @return Phase factor (+1 or -1)
 */
double phase2(int value)
{
    return (value % 4 == 2) ? -1. : 1.;
}

/**
 * @brief Phase factor based on absolute value
 *
 * Returns -1 if |value| % 4 == 2, otherwise returns 1
 *
 * @param value Integer input
 * @return Phase factor (+1 or -1)
 */
double phaseabs(int value)
{
    return (std::abs(value % 4) == 2) ? -1. : 1.;
}

/**
 * @brief Constrains an index within bounds
 *
 * @param index The index to constrain
 * @param size Upper bound
 * @return Index clamped to [0, size-1]
 */
int pad(int index, int size)
{
    if (index < 0)
        return 0;
    if (index >= size)
        return size - 1;
    return index;
}

// Maximum factorial value for precomputation
const int MAX_FACT = 100;

/**
 * @brief Returns the logarithm of n factorial
 *
 * Uses precomputed values for efficiency. First call initializes
 * the lookup table using an immediately-invoked lambda.
 *
 * @param n Input value
 * @return log(n!)
 * @throws std::runtime_error if n is out of range
 */
double getLogFactorial(int n)
{
    // Static initialization of logarithm factorial table
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

/**
 * @brief Evaluates Jacobi polynomial at point z
 *
 * Computes the Jacobi polynomial P_n^(a,b)(z) which is orthogonal with
 * respect to the weight function (1-x)^a(1+x)^b on the interval [-1,1]
 *
 * @param n Order of the polynomial
 * @param a First parameter of the Jacobi polynomial
 * @param b Second parameter of the Jacobi polynomial
 * @param z Point at which to evaluate the polynomial
 * @return Value of the Jacobi polynomial
 */
double jacobi_pols(int n, int a, int b, double z)
{
    return boost::math::jacobi<double>(n, a, b, z);
}

/**
 * @brief Calculates the standard Euclidean scalar product of two vectors
 *
 * @param a First vector
 * @param b Second vector
 * @return Scalar product of a and b
 */
double scalar_product(std::vector<double> a, std::vector<double> b)
{
    double product = 0;
    for (int i = 0; i <= a.size() - 1; i++)
    {
        product = product + (a[i]) * (b[i]);
    }
    return product;
}

/**
 * @brief Calculates Minkowski scalar product with metric (+,-,-,-)
 *
 * The time component is assumed to be at index 3 of the vectors.
 *
 * @param v1 First 4-vector
 * @param v2 Second 4-vector
 * @return Minkowski scalar product v1·v2
 */
double min_scalar_product(std::vector<double> v1, std::vector<double> v2)
{
    return v1[3] * v2[3] - (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
}

/**
 * @brief Computes the cosine of the helicity angle θij
 *
 * Calculates the cosine of the angle between particles i and j
 * in the rest frame of the decay, given Mandelstam variables
 * and particle masses.
 *
 * @param σs Array of Mandelstam variables (s12, s23, s31)
 * @param msq Array of squared masses (m1², m2², m3², M²)
 * @param k Index specifying which particle is spectator (1, 2, or 3)
 * @return Cosine of the helicity angle
 */
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

/**
 * @brief Convenience function to calculate cos(θ23)
 *
 * @param σs Mandelstam variables
 * @param msq Squared masses
 * @return Cosine of θ23
 */
double ThreeBodyDecays::cosθ23(MandelstamTuple σs, ThreeBodyMasses msq)
{
    return cosθij(σs, msq, 1);
}

/**
 * @brief Convenience function to calculate cos(θ31)
 *
 * @param σs Mandelstam variables
 * @param msq Squared masses
 * @return Cosine of θ31
 */
double ThreeBodyDecays::cosθ31(MandelstamTuple σs, ThreeBodyMasses msq)
{
    return cosθij(σs, msq, 2);
}

/**
 * @brief Convenience function to calculate cos(θ12)
 *
 * @param σs Mandelstam variables
 * @param msq Squared masses
 * @return Cosine of θ12
 */
double ThreeBodyDecays::cosθ12(MandelstamTuple σs, ThreeBodyMasses msq)
{
    return cosθij(σs, msq, 3);
}

/**
 * @brief Performs modular arithmetic with 1-based indexing
 *
 * This helper function handles the circular nature of three-body indices.
 *
 * @param k Index to transform
 * @param n Modulus
 * @return Result of modular arithmetic with 1-based indexing
 */
int mod1(int k, int n)
{
    return ((k - 1) % n + n) % n + 1;
}

/**
 * @brief Reorders tuple elements based on spectator index k
 *
 * @param k Spectator particle index (1, 2, or 3)
 * @param t Original tuple
 * @return Reordered tuple with elements arranged according to k
 */
std::array<double, 3> circleorigin(int k, const std::array<double, 3> &t)
{
    return {
        t[mod1(k - 1, 3) - 1], // First element
        t[mod1(k, 3) - 1],     // Second element
        t[mod1(k + 1, 3) - 1]  // Third element
    };
}

/**
 * @brief Performs linear interpolation between limits
 *
 * @param y Value in range [0,1]
 * @param limits Pair containing minimum and maximum limits
 * @return Interpolated value
 */
double fitin(double y, const std::pair<double, double> &limits)
{
    double a = limits.first;
    double b = limits.second;
    return a + y * (b - a);
}

/**
 * @brief Calculates kinematic limits for Mandelstam variable σk
 *
 * @param ms Array of particle masses
 * @param k Spectator particle index
 * @return Pair containing (min, max) allowed values for σk
 */
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

/**
 * @brief Calculates σj given σk and cosθ
 *
 * @param cos_theta Cosine of the helicity angle
 * @param σk Mandelstam variable σk
 * @param ms_squared Squared masses of all particles
 * @param k Spectator particle index
 * @return The value of σj
 */
double ThreeBodyDecays::σjofk(double cos_theta, double σk,
                              const std::array<double, 4> &ms_squared, int k)
{
    auto [i, j] = ij_from_k(k);
    i--;
    j--;
    k--; // Convert to 0-based indexing

    double s = ms_squared[3]; // Parent mass squared

    // Calculate energy and momentum terms
    double EE4σ = (σk + ms_squared[j] - ms_squared[i]) *
                  (σk + s - ms_squared[k]);
    double p2q24σ = Kallen(σk, ms_squared[i], ms_squared[j]) *
                    Kallen(s, σk, ms_squared[k]);
    p2q24σ = (p2q24σ < 0) ? 0.0
                          : p2q24σ; // Handle numerical errors

    // Calculate σi
    double σi = s + ms_squared[j] -
                (EE4σ - std::sqrt(p2q24σ) * cos_theta) / (2 * σk);

    return σi;
}

/**
 * @brief Converts random numbers to physical Mandelstam variables
 *
 * This function is useful for Monte Carlo integration over phase space.
 *
 * @param x Array of two random numbers in range [0,1]
 * @param ms Masses of all particles
 * @param k Spectator particle index
 * @return Tuple of Mandelstam variables (s12, s23, s31)
 */
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

/**
 * @brief Parses a complex number from string representation
 *
 * @param str String in the format "a+bi" or "a-bi"
 * @return Complex number
 */
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

/**
 * @brief Parses a fraction from string representation
 *
 * @param str String in the format "a/b" or a decimal number
 * @return Double value
 */
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

/**
 * @brief Applies a recoupling scheme to calculate amplitude components
 *
 * @param recoupling Function defining the recoupling scheme
 * @param two_ms Array of helicity values (doubled)
 * @param two_js Array of spin values (doubled)
 * @return Complex amplitude component
 */
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


/**
 * @brief Calculates the aligned 4D amplitude tensor for the decay chain
 *
 * This function computes the 4D amplitude tensor for the given decay chain,
 * applying the necessary Wigner rotations and recouplings. It also handles
 * the alignment of spins and the calculation of the lineshape.
 *
 * @param dc The decay chain information
 * @param σs The Mandelstam variables
 * @return The aligned 4D amplitude tensor
 */
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

    // dimensions for F0 in order (j1, j2, j3, j4)
    std::vector<int> dims = {two_js[0] + 1, two_js[1] + 1, two_js[2] + 1, two_js[3] + 1};

    // initialise F0 with dimensions (2j1+1, 2j2+1, 2j3+1, 2j4+1)
    Tensor4Dcomp F0(dims[0], std::vector<std::vector<std::vector<complex>>>(
                                 dims[1], std::vector<std::vector<complex>>(
                                              dims[2], std::vector<complex>(dims[3], 0.0))));

    // dimensions for F in order (i, j, k, 3)
    Tensor4Dcomp F(dims[i_idx], std::vector<std::vector<std::vector<complex>>>(
                                    dims[j_idx], std::vector<std::vector<complex>>(
                                                     dims[k_idx], std::vector<complex>(dims[3], 0.0))));

    // Helper function to pad indices
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

                    // Calculate shifted indices
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

                    // Ensure indices are within bounds
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
        int mod_result = (l + 1) % 4; // Dies funktioniert für l=0,1,2,3
        if (mod_result == 0)
            mod_result = 4; // Für den Fall l=3
        if (mod_result == 4)
            mod_result = 0; // Für den Fall l=3

        std::unique_ptr<AbstractWignerRotation> _w = wr(k, _refζ, mod_result);
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

    // Calculate the tensor contraction
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
        // div(_two_j + _two_λ, 2) + 1, but adjust for 0-based indexing
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

/**
 * @brief Calculates the total intensity from the decay chain
 *
 * Sums the squared amplitudes over all helicity states, weighted by the
 * given complex weight. This function corresponds to the intensity
 * calculation in the decay process.
 *
 * @param dc The decay chain information
 * @param σs The Mandelstam variables
 * @param k_amp The index of the amplitude to calculate
 * @param weight The weight to apply to the amplitude
 * @param refζs The reference spin-quantum numbers
 * @return The total intensity
 */
double ThreeBodyDecays::intensity(const DecayChain &dc, const MandelstamTuple &σs, const int &k_amp, const complex weight, const std::vector<int> refζs)
{
    // Get the 4D amplitude tensor
    auto amp = amplitude4dcomp(dc, σs, k_amp, refζs);

    // Sum squared amplitudes
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

            // Check parity condition
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

// Implementation of possible_lsLS function
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
    // Create the Cartesian product
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

/**
 * @brief Creates a recoupling function for the given type and parameters
 *
 * This function returns a callable recoupling function (e.g., LS, Parity)
 * based on the specified type and parameters. The recoupling function can
 * then be used to compute recoupled amplitudes.
 *
 * @param type The type of recoupling to perform
 * @param params The parameters for the recoupling (e.g., L and S values)
 * @param parityPhase Optional parameter for parity recoupling
 * @return A recoupling function
 */
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
            // std::cout << noRecoupling.get_two_λa() << noRecoupling.get_two_λb() << " " << two_ms[0] << two_ms[1] << std::endl;

            //  (cs.two_λa == two_λa) * (cs.two_λb == two_λb)
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
                return 1.0; //*ηηηphaseisplus;
            }

            // Zweite Bedingung: cs.two_λa == -two_λa and cs.two_λb == -two_λb
            // Dies prüft, ob der Parameter gleich dem negativen des übergebenen Werts ist
            if (two_λa == -two_ms[0] && two_λb == -two_ms[1])
            {
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

            // Call the jls_coupling function
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
