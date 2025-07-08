#ifndef THREEBODYDECAYS_HH
#define THREEBODYDECAYS_HH

#include <complex>
#include <map>
#include <memory> // for std::unique_ptr
#include <string>
#include <tuple>
#include <vector>
#include <array>
#include <functional> // for std::function

#include "SpinParity.hh"     // For spin and parity representation
#include "WignerRotation.hh" // For Wigner rotation operations

using complex = std::complex<double>;
// Type definitions for commonly used structures
using MandelstamTuple = std::array<double, 3>;    // (s12, s23, s31)
using MassTuple = std::array<double, 4>;          // (m1, m2, m3, M)
using SpinTuple = std::array<int, 4>;             // (j1, j2, j3, J)
using RecouplingLS =
    std::function<complex(const std::array<int, 2> &, const std::array<int, 3> &)>;
using LineshapeFunction = std::function<complex(double)>;
using HelicityFunction =
    std::function<complex(const std::array<int, 2> &, const std::array<int, 3> &)>;
using ThreeBodyMasses = std::array<double, 4>;
using ThreeBodySpins = std::array<int, 4>;
using Matrix2D = std::vector<std::vector<double>>;
using Tensor4Dcomp = std::vector<std::vector<std::vector<std::vector<complex>>>>;
using Tensor4D = std::vector<std::vector<std::vector<std::vector<double>>>>;

/**
 * @brief Calculates Minkowski scalar product with metric (+,-,-,-)
 *
 * @param v1 First 4-vector
 * @param v2 Second 4-vector
 * @return Minkowski scalar product v1·v2
 */
double min_scalar_product(std::vector<double> v1, std::vector<double> v2);

/**
 * @brief Checks if two double values are approximately equal
 *
 * @param a First value
 * @param b Second value
 * @param epsilon Maximum allowed difference
 * @return true if values are within epsilon of each other, false otherwise
 */
bool approx_equal(double a, double b, double epsilon = 1e-6);

/**
 * @brief Represents a three-body system with masses and spins
 *
 * Contains the masses and spin quantum numbers for three decay particles
 * and their parent particle.
 */
class ThreeBodySystem
{
public:
    std::array<double, 4> ms;     // Masses (m1, m2, m3, M)
    std::array<int, 4> two_js;    // Doubled spin values (2*j1, 2*j2, 2*j3, 2*J)

    /**
     * @brief Constructs a three-body system
     *
     * @param ms_ Array of masses (m1, m2, m3, M)
     * @param spins Array of spins (j1, j2, j3, J)
     */
    ThreeBodySystem(const std::array<double, 4> &ms_,
                    std::array<int, 4> &spins) : ms(ms_)
    {
        // Initialize two_js with doubled spin values
        for (int i = 0; i < 4; ++i)
        {
            two_js[i] = spins[i];
        }
    }

    /**
     * @brief Gets the array of masses
     * @return Reference to the masses array
     */
    const std::array<double, 4> &get_ms() const { return ms; }

    /**
     * @brief Gets the array of doubled spin values
     * @return Reference to the doubled spins array
     */
    const std::array<int, 4> &get_two_js() const { return two_js; }
};

/**
 * @brief Represents the intrinsic parities of particles in a three-body decay
 */
class ThreeBodyParities
{
public:
    /**
     * @brief Constructs parity representation for a three-body system
     *
     * @param P1 Parity of first particle ('+' or '-')
     * @param P2 Parity of second particle ('+' or '-')
     * @param P3 Parity of third particle ('+' or '-')
     * @param P0 Parity of parent particle ('+' or '-')
     */
    ThreeBodyParities(char P1, char P2, char P3, char P0);

    // Getters
    char get_P1() const { return P1_; }
    char get_P2() const { return P2_; }
    char get_P3() const { return P3_; }
    char get_P0() const { return P0_; }

    /**
     * @brief Parity product operator
     *
     * @param other Another ThreeBodyParities object
     * @return '+' if parities of parent particles match, '-' otherwise
     */
    char operator*(const ThreeBodyParities &other) const;

private:
    char P1_, P2_, P3_, P0_;  // Parity values for each particle
};

/**
 * @brief Represents a decay chain in three-body decay
 *
 * Contains all the information needed to calculate the amplitude
 * for a specific decay chain through an intermediate resonance.
 */
class DecayChain
{
public:
    int k;                      // Spectator particle index (1, 2, or 3)
    int two_j;                  // Doubled spin of the resonance
    LineshapeFunction Xlineshape; // Lineshape function for the resonance
    HelicityFunction HRk;       // Helicity coupling function for resonance-spectator
    HelicityFunction Hij;       // Helicity coupling function for i-j pair
    ThreeBodySystem tbs;        // Three-body system information

    /**
     * @brief Constructs a decay chain
     *
     * @param k_ Spectator particle index
     * @param two_j_ Doubled spin of the resonance
     * @param Xlineshape_ Lineshape function for the resonance
     * @param HRk_ Helicity coupling function for resonance-spectator
     * @param Hij_ Helicity coupling function for i-j pair
     * @param tbs_ Three-body system
     */
    DecayChain(int k_, int two_j_, LineshapeFunction Xlineshape_,
               HelicityFunction HRk_, HelicityFunction Hij_,
               ThreeBodySystem tbs_) :
                                       k(k_),
                                       two_j(two_j_),
                                       Xlineshape(Xlineshape_),
                                       HRk(HRk_),
                                       Hij(Hij_),
                                       tbs(tbs_)
    {
    }
};

/**
 * @brief Implements the "no recoupling" scheme for helicity amplitudes
 *
 * Returns 1.0 only when the helicities match the specified values,
 * and 0.0 otherwise.
 */
class NoRecoupling
{
private:
    int two_λa_;  // Doubled helicity value for first particle
    int two_λb_;  // Doubled helicity value for second particle

public:
    /**
     * @brief Constructs a NoRecoupling object
     *
     * @param two_λa Doubled helicity for first particle
     * @param two_λb Doubled helicity for second particle
     */
    NoRecoupling(int two_λa, int two_λb) : two_λa_(two_λa), two_λb_(two_λb) {}

    int get_two_λa() const { return two_λa_; }
    int get_two_λb() const { return two_λb_; }

    /**
     * @brief Function call operator implementing the no-recoupling logic
     *
     * @param two_ms Array of helicity values to test against
     * @param two_js Array of spin values (not used in this implementation)
     * @return 1.0 if helicities match, 0.0 otherwise
     */
    complex operator()(const std::array<int, 2> &two_ms, const std::array<int, 3> &two_js) const
    {
        return (two_λa_ == two_ms[0] && two_λb_ == two_ms[1]) ? complex(1.0, 0.0) : complex(0.0, 0.0);
    }
};

/**
 * @brief Implements parity-based recoupling for helicity amplitudes
 *
 * Returns values based on whether helicities match directly or
 * with opposite signs, accounting for parity properties.
 */
class ParityRecoupling
{
private:
    int two_λa_;          // Doubled helicity value for first particle
    int two_λb_;          // Doubled helicity value for second particle
    bool ηηηphaseisplus_; // Phase parameter for parity recoupling

public:
    /**
     * @brief Constructs a ParityRecoupling object
     *
     * @param two_λa Doubled helicity for first particle
     * @param two_λb Doubled helicity for second particle
     * @param ηηηphaseisplus Phase parameter (true for +1, false for -1)
     */
    ParityRecoupling(int two_λa, int two_λb, bool ηηηphaseisplus) :
        two_λa_(two_λa), two_λb_(two_λb), ηηηphaseisplus_(ηηηphaseisplus) {}

    int get_two_λa() const { return two_λa_; }
    int get_two_λb() const { return two_λb_; }
    bool is_ηηηphaseisplus() const { return ηηηphaseisplus_; }

    /**
     * @brief Function call operator implementing parity recoupling logic
     *
     * @param two_ms Array of helicity values to test against
     * @param two_js Array of spin values (not used in this implementation)
     * @return Complex amplitude according to parity recoupling rules
     */
    complex operator()(const std::array<int, 2> &two_ms, const std::array<int, 3> &two_js) const
    {
        // Direct match
        if (two_λa_ == two_ms[0] && two_λb_ == two_ms[1])
            return complex(1.0, 0.0);

        // Match with opposite signs (parity-related)
        if (two_λa_ == -two_ms[0] && two_λb_ == -two_ms[1])
            return complex(2 * ηηηphaseisplus_ - 1, 0.0);

        return complex(0.0, 0.0);
    }
};

/**
 * @brief Enumeration of supported recoupling schemes
 */
enum class RecouplingType
{
    NoRecoupling,     // No recoupling (delta function)
    ParityRecoupling, // Parity-based recoupling
    LSRecoupling      // LS coupling scheme
};

/**
 * @brief Creates a decay chain with automatic LS coupling
 *
 * @param k Spectator particle index
 * @param Xlineshape Lineshape function for the resonance
 * @param jp Spin-parity string for the resonance (e.g., "1-", "2+")
 * @param Ps Three-body parities
 * @param tbs Three-body system
 * @return Shared pointer to created DecayChain
 */
std::shared_ptr<DecayChain> createDecayChainLS(
    int k,
    std::function<std::complex<double>(double)> Xlineshape,
    const std::string &jp,
    const ThreeBodyParities &Ps,
    const ThreeBodySystem &tbs);

/**
 * @brief Creates multiple decay chains with all possible LS couplings
 *
 * @param k Spectator particle index
 * @param Xlineshape Lineshape function for the resonance
 * @param jp Spin-parity string for the resonance
 * @param Ps Three-body parities
 * @param tbs Three-body system
 * @return Vector of shared pointers to DecayChain objects
 */
std::vector<std::shared_ptr<DecayChain>> createDecayChainsLS(
    int k,
    std::function<std::complex<double>(double)> Xlineshape,
    const std::string &jp,
    const ThreeBodyParities &Ps,
    const ThreeBodySystem &tbs);

/**
 * @brief Creates a decay chain with explicit coupling parameters
 *
 * @param k Spectator particle index
 * @param Xlineshape Lineshape function for the resonance
 * @param jp Spin-parity string for the resonance
 * @param tbs Three-body system
 * @param HRkType Recoupling type for resonance-spectator coupling
 * @param HRkParams Parameters for resonance-spectator coupling
 * @param HRkParityPhase Phase parameter for ParityRecoupling (if used)
 * @param HijType Recoupling type for i-j pair coupling
 * @param HijParams Parameters for i-j pair coupling
 * @param HijParityPhase Phase parameter for ParityRecoupling (if used)
 * @return Shared pointer to created DecayChain
 */
std::shared_ptr<DecayChain> createDecayChainCoupling(
    int k,
    std::function<std::complex<double>(double)> Xlineshape,
    const std::string &jp,
    const ThreeBodySystem &tbs,
    RecouplingType HRkType = RecouplingType::NoRecoupling,
    const std::array<int, 2> &HRkParams = {0, 0},
    bool HRkParityPhase = false,
    RecouplingType HijType = RecouplingType::NoRecoupling,
    const std::array<int, 2> &HijParams = {0, 0},
    bool HijParityPhase = false);

/**
 * @brief Creates a recoupling function from parameters
 *
 * @param type Type of recoupling scheme to create
 * @param params Parameters for the recoupling (interpretation depends on type)
 * @param parityPhase Phase parameter for ParityRecoupling
 * @return Function object implementing the requested recoupling
 */
RecouplingLS createRecouplingFunction(
    RecouplingType type,
    const std::array<int, 2> &params,
    bool parityPhase = true);

/**
 * @brief Structure to hold LS coupling information
 *
 * Contains both (l,s) values for i-j coupling and (L,S) values for
 * resonance-spectator coupling.
 */
struct LSCoupling
{
    std::array<int, 2> two_ls;  // Doubled (l,s) values for i-j coupling
    std::array<int, 2> two_LS;  // Doubled (L,S) values for resonance-spectator coupling

    /**
     * @brief Constructs an LSCoupling object
     *
     * @param ls Doubled (l,s) values
     * @param LS Doubled (L,S) values
     */
    LSCoupling(const std::array<int, 2> &ls, const std::array<int, 2> &LS) : two_ls(ls), two_LS(LS) {}
};

/**
 * @brief Calculates all possible LS couplings for a resonance
 *
 * @param jp Spin-parity of the resonance
 * @param two_js Doubled spins of all particles
 * @param Ps Parities of all particles
 * @param k Spectator particle index
 * @return Vector of all possible LS couplings
 */
std::vector<LSCoupling> possible_lsLS(
    const SpinParity &jp,
    const std::array<int, 4> &two_js,
    const ThreeBodyParities &Ps,
    int k);

/**
 * @brief Calculates possible (l,s) values for two particles
 *
 * @param jp1 Spin-parity of first particle
 * @param jp2 Spin-parity of second particle
 * @param jp Spin-parity of the combined system
 * @return Vector of possible (l,s) values
 */
std::vector<std::array<int, 2>> possible_ls(
    const SpinParity &jp1,
    const SpinParity &jp2,
    const SpinParity &jp);

/**
 * @brief Calculates possible (l,s) values for i-j pair
 *
 * @param jp Spin-parity of the resonance
 * @param two_js Doubled spins of all particles
 * @param Ps Parities of all particles
 * @param k Spectator particle index
 * @return Vector of possible (l,s) values
 */
std::vector<std::array<int, 2>> possible_ls_ij(
    const SpinParity &jp,
    const std::array<int, 4> &two_js,
    const ThreeBodyParities &Ps,
    int k);

/**
 * @brief Calculates possible (L,S) values for resonance-spectator coupling
 *
 * @param jp Spin-parity of the resonance
 * @param two_js Doubled spins of all particles
 * @param Ps Parities of all particles
 * @param k Spectator particle index
 * @return Vector of possible (L,S) values
 */
std::vector<std::array<int, 2>> possible_ls_Rk(
    const SpinParity &jp,
    const std::array<int, 4> &two_js,
    const ThreeBodyParities &Ps,
    int k);

/**
 * @brief Main class implementing three-body decay calculations
 *
 * Provides methods to calculate decay amplitudes and related quantities
 * for three-body decays with intermediate resonances.
 */
class ThreeBodyDecays
{
public:
    /**
     * @brief Calculates aligned amplitude in matrix form
     *
     * @param dc Decay chain
     * @param sigma Mandelstam variables
     * @return 2D matrix of complex amplitudes
     */
    std::vector<std::vector<complex>> aligned_amplitude(const DecayChain &dc, const MandelstamTuple &sigma);

    /**
     * @brief Calculates amplitude for specific helicity values
     *
     * @param dc Decay chain
     * @param σs Mandelstam variables
     * @param two_λs Doubled helicity values for all particles
     * @param k_amp Amplitude type index
     * @param refζs Reference angles for Wigner rotations
     * @return Complex amplitude
     */
    complex amplitude(const DecayChain &dc, const MandelstamTuple &σs, const std::vector<int> &two_λs, const int &k_amp, const std::vector<int> &refζs = {-1, -1, -1, -1});

    /**
     * @brief Calculates 4D tensor of aligned amplitudes
     *
     * @param dc Decay chain
     * @param σs Mandelstam variables
     * @return 4D tensor of complex amplitudes
     */
    Tensor4Dcomp aligned_amplitude4dcomp(const DecayChain &dc, const MandelstamTuple &σs);

    /**
     * @brief Calculates 4D tensor of amplitudes with rotations
     *
     * @param dc Decay chain
     * @param σs Mandelstam variables
     * @param k_amp Amplitude type index
     * @param refζ Reference angles for Wigner rotations
     * @return 4D tensor of complex amplitudes
     */
    Tensor4Dcomp amplitude4dcomp(const DecayChain &dc, const MandelstamTuple &σs, const int &k_amp, std::vector<int> refζ = {-1, -1, -1, -1});

    /**
     * @brief Calculates total decay intensity
     *
     * @param dc Decay chain
     * @param σs Mandelstam variables
     * @param k_amp Amplitude type index
     * @param weight Complex weight for the amplitude
     * @param refζs Reference angles for Wigner rotations
     * @return Total intensity (proportional to decay rate)
     */
    double intensity(const DecayChain &dc, const MandelstamTuple &σs, const int &k_amp, const complex weight = complex(1,1), const std::vector<int> refζs = {-1, -1, -1, -1});

    /**
     * @brief Applies recoupling scheme to calculate amplitude components
     *
     * @param recoupling Recoupling function
     * @param two_ms Doubled helicity values
     * @param two_js Doubled spin values
     * @return Complex amplitude component
     */
    complex amplitude_recoupling(const RecouplingLS &recoupling, const std::array<int, 2> &two_ms, const std::array<int, 3> &two_js);

    /**
     * @brief Parses a fraction from string representation
     *
     * @param str String in format "a/b" or decimal number
     * @return Double value
     */
    double parseFraction(const std::string &str);

    /**
     * @brief Parses a complex number from string representation
     *
     * @param str String in format "a+bi" or "a-bi"
     * @return Complex number
     */
    complex parseComplex(const std::string &str);

    /**
     * @brief Converts random numbers to Mandelstam variables
     *
     * @param x Two random numbers in range [0,1]
     * @param ms Masses of all particles
     * @param k Spectator particle index
     * @return Tuple of Mandelstam variables
     */
    MandelstamTuple x2σs(const std::array<double, 2> &x, ThreeBodyMasses ms, int k);

    /**
     * @brief Calculates σj given σk and cosθ
     *
     * @param cos_theta Cosine of helicity angle
     * @param σk Mandelstam variable σk
     * @param ms_squared Squared masses
     * @param k Spectator particle index
     * @return Value of σj
     */
    double σjofk(double cos_theta, double σk, const std::array<double, 4> &ms_squared, int k);

    /**
     * @brief Calculates cosine of helicity angle
     *
     * @param σs Mandelstam variables
     * @param msq Squared masses
     * @param k Spectator particle index
     * @return Cosine of helicity angle
     */
    double cosθij(const std::array<double, 3> &σs, const std::array<double, 4> &msq, int k);

    /**
     * @brief Convenience function for cos(θ23)
     *
     * @param σs Mandelstam variables
     * @param msq Squared masses
     * @return Cosine of θ23
     */
    double cosθ23(MandelstamTuple σs, ThreeBodyMasses msq);

    /**
     * @brief Convenience function for cos(θ31)
     *
     * @param σs Mandelstam variables
     * @param msq Squared masses
     * @return Cosine of θ31
     */
    double cosθ31(MandelstamTuple σs, ThreeBodyMasses msq);

    /**
     * @brief Convenience function for cos(θ12)
     *
     * @param σs Mandelstam variables
     * @param msq Squared masses
     * @return Cosine of θ12
     */
    double cosθ12(MandelstamTuple σs, ThreeBodyMasses msq);

protected:
    void A_test();

private:
};

#endif
