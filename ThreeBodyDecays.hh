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

using complex = std::complex<double>;
using MandelstamTuple = std::array<double, 3>;
using MassTuple = std::array<double, 4>;
using SpinTuple = std::array<int, 4>;
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

// Forward declarations
class AbstractWignerRotation;

// Kallen function declaration
double Kallen(double x, double y, double z);

// Function to compute i, j, k indices from Wigner rotation
std::tuple<int, int, int> ijk(const AbstractWignerRotation &wr);

// Helper function to determine if a system is sequential
bool issequential(int system_a, int reference_b);

// In ThreeBodyDecays.hh
double min_scalar_product(std::vector<double> v1, std::vector<double> v2);

// In ThreeBodyDecays.hh
bool approx_equal(double a, double b, double epsilon = 1e-6);

// Helper function to create Wigner rotation
std::unique_ptr<AbstractWignerRotation> wr(int system_a, int reference_b,
                                           int particle_c = 0);

// ThreeBodySystem class
// Korrigierte ThreeBodySystem Klasse
class ThreeBodySystem
{
public:
    std::array<double, 4> ms;
    std::array<int, 4> two_js;

    ThreeBodySystem(const std::array<double, 4> &ms_,
                    std::array<int, 4> &spins) : ms(ms_)
    // two_js(two_js_)
    // spins_(spins)
    {
        // Initialisiere spins_ mit verdoppelten Werten
        for (int i = 0; i < 4; ++i)
        {
            two_js[i] = spins[i];//2;
        }
    }

    const std::array<double, 4> &get_ms() const { return ms; }
    const std::array<int, 4> &get_two_js() const { return two_js; }
};

class ThreeBodyParities
{
public:
    ThreeBodyParities(char P1, char P2, char P3, char P0);

    // Getters
    char get_P1() const { return P1_; }
    char get_P2() const { return P2_; }
    char get_P3() const { return P3_; }
    char get_P0() const { return P0_; }

    // Parity product operator
    char operator*(const ThreeBodyParities &other) const;

private:
    char P1_, P2_, P3_, P0_;
};

class SpinParity
{
public:
    SpinParity(const std::string &jp);

    // Getters
    int get_two_j() const {
        return two_j_; // Now properly returns the doubled value
    }
    int get_j() const {
        return two_j_ / 2; // Returns the original j value
    }
    char get_p() const { return p_; }

private:
    int two_j_; // 2 * j
    char p_;    // '+' or '-'
    std::string str; // Store the original string for parsing
};

class DecayChain
{
public:
    int k;
    // std::array<double, 4> ms;
    // std::array<int, 4> two_js;
    int two_j;
    LineshapeFunction Xlineshape; // Lineshape function for intermediate state
    HelicityFunction HRk;         // Helicity function for intermediate state
    HelicityFunction Hij;
    ThreeBodySystem tbs;

    /*DecayChain( int k_, const std::array<double, 4>& ms_,
                const std::array<int, 4>& two_js_, int two_j_,
                LineshapeFunction Xlineshape_, HelicityFunction Hij_,
                HelicityFunction HRk_ , ThreeBodySystem tbs_) :*/
    DecayChain(int k_, int two_j_, LineshapeFunction Xlineshape_,
               HelicityFunction HRk_, HelicityFunction Hij_,
               ThreeBodySystem tbs_) :

                                       k(k_),
                                       // ms( ms_ ),
                                       // two_js(  ),
                                       two_j(two_j_),
                                       Xlineshape(Xlineshape_),
                                       HRk(HRk_),
                                       Hij(Hij_),
                                       tbs(tbs_)
    {
        // for (int i = 0; i < 4; ++i) {
        //     two_js[i] = two_js_[i] * 2;
        // }
    }
};

// Abstract base class for Wigner rotations
class AbstractWignerRotation
{
public:
    virtual ~AbstractWignerRotation() = default;
    virtual double operator()(double w, const std::array<double, 3> &sigma,
                              const std::array<double, 4> &ms2) const = 0;
    virtual int get_k() const = 0;
    virtual bool is_positive() const = 0;
    virtual bool is_even() const = 0;
    virtual double cos_zeta(const std::array<double, 3> &sigma,
                            const std::array<double, 4> &ms2) const = 0;
};

// Trivial Wigner rotation
class TrivialWignerRotation : public AbstractWignerRotation
{
private:
    int k;

public:
    TrivialWignerRotation(int k_) : k(k_) {}
    double operator()(double w, const std::array<double, 3> &sigma,
                      const std::array<double, 4> &ms2) const override
    {
        return 1.0; // Trivial rotation always returns 1
    }
    int get_k() const override { return k; }
    bool is_positive() const override { return true; }
    bool is_even() const override { return true; }
    double cos_zeta(const std::array<double, 3> &sigma,
                    const std::array<double, 4> &ms2) const override;
};

// Base Wigner rotation class
class WignerRotation : public AbstractWignerRotation
{
protected:
    int k;
    bool is_positive_;
    bool is_even_;

public:
    WignerRotation(int k_, bool is_positive_, bool is_even_) : k(k_), is_positive_(is_positive_), is_even_(is_even_)
    {
    }
    virtual double operator()(double w, const std::array<double, 3> &sigma,
                              const std::array<double, 4> &ms2) const = 0;
    int get_k() const override { return k; }
    bool is_positive() const override { return is_positive_; }
    bool is_even() const override { return is_even_; }
    virtual double cos_zeta(const std::array<double, 3> &sigma,
                            const std::array<double, 4> &ms2) const = 0;
};

// Arg0 Wigner rotation
class Arg0WignerRotation : public WignerRotation
{
public:
    Arg0WignerRotation(int k_, bool is_positive_) : WignerRotation(k_, is_positive_, true)
    {
    }
    double operator()(double w, const std::array<double, 3> &sigma,
                      const std::array<double, 4> &ms2) const override
    {
        auto [i, j, k] = ijk(*this);
        double s = ms2[3]; // Parent mass squared
        double EE4s = (s + ms2[i] - sigma[i]) * (s + ms2[k] - sigma[k]);
        double pp4s = std::sqrt(Kallen(s, ms2[i], sigma[i]) *
                                Kallen(s, ms2[k], sigma[k]));
        double rest = sigma[j] - ms2[i] - ms2[k];
        return (EE4s - 2 * s * rest) / pp4s;
    }
    double cos_zeta(const std::array<double, 3> &sigma,
                    const std::array<double, 4> &ms2) const override;
};

// Arg2 Wigner rotation
class Arg2WignerRotation : public WignerRotation
{
public:
    Arg2WignerRotation(int k_, bool is_positive_, bool is_even_) : WignerRotation(k_, is_positive_, is_even_)
    {
    }
    double operator()(double w, const std::array<double, 3> &sigma,
                      const std::array<double, 4> &ms2) const override
    {
        auto [i, j, k] = ijk(*this);
        if (!is_even())
        {
            std::swap(i, j);
        }
        if (std::abs(ms2[k]) < 1e-10)
            return 1.0; // Handle massless case

        double s = ms2[3];
        double EE4mksq = (s + ms2[k] - sigma[k]) *
                         (sigma[i] - ms2[k] - ms2[j]);
        double pp4mksq = std::sqrt(Kallen(s, ms2[k], sigma[k]) *
                                   Kallen(ms2[k], ms2[j], sigma[i]));
        double rest = sigma[j] - s - ms2[j];
        return (2 * ms2[k] * rest + EE4mksq) / pp4mksq;
    }
    double cos_zeta(const std::array<double, 3> &sigma,
                    const std::array<double, 4> &ms2) const override;
};

// Arg3 Wigner rotation
class Arg3WignerRotation : public WignerRotation
{
public:
    Arg3WignerRotation(int k_, bool is_positive_) : WignerRotation(k_, is_positive_, true)
    {
    }
    double operator()(double w, const std::array<double, 3> &sigma,
                      const std::array<double, 4> &ms2) const override
    {
        auto [i, j, k] = ijk(*this);
        if (std::abs(ms2[k]) < 1e-10)
            return 1.0; // Handle massless case

        double EE4m1sq = (sigma[i] - ms2[j] - ms2[k]) *
                         (sigma[j] - ms2[k] - ms2[i]);
        double pp4m1sq = std::sqrt(Kallen(sigma[i], ms2[j], ms2[k]) *
                                   Kallen(sigma[j], ms2[k], ms2[i]));
        double rest = ms2[i] + ms2[j] - sigma[k];
        return (2 * ms2[k] * rest + EE4m1sq) / pp4m1sq;
    }
    double cos_zeta(const std::array<double, 3> &sigma,
                    const std::array<double, 4> &ms2) const override;
};


// Entferne das 'override' Schlüsselwort in der NoRecoupling Klasse
// In ThreeBodyDecays.hh
// In ThreeBodyDecays.hh
class NoRecoupling
{
private:
    int two_λa_;
    int two_λb_;

public:
    NoRecoupling(int two_λa, int two_λb) : two_λa_(two_λa), two_λb_(two_λb) {}

    int get_two_λa() const { return two_λa_; }
    int get_two_λb() const { return two_λb_; }

    // Diese Methode implementiert exakt die Julia-Logik
    complex operator()(const std::array<int, 2> &two_ms, const std::array<int, 3> &two_js) const
    {
        // (cs.two_λa == two_λa) * (cs.two_λb == two_λb)
        return (two_λa_ == two_ms[0] && two_λb_ == two_ms[1]) ? complex(1.0, 0.0) : complex(0.0, 0.0);
    }
};

class ParityRecoupling
{
private:
    int two_λa_;
    int two_λb_;
    bool ηηηphaseisplus_;

public:
    ParityRecoupling(int two_λa, int two_λb, bool ηηηphaseisplus) : two_λa_(two_λa), two_λb_(two_λb), ηηηphaseisplus_(ηηηphaseisplus) {}

    int get_two_λa() const { return two_λa_; }
    int get_two_λb() const { return two_λb_; }
    bool is_ηηηphaseisplus() const { return ηηηphaseisplus_; }

    // Diese Methode implementiert exakt die Julia-Logik
    complex operator()(const std::array<int, 2> &two_ms, const std::array<int, 3> &two_js) const
    {
        // (cs.two_λa == two_λa) * (cs.two_λb == two_λb) && return 1
        if (two_λa_ == two_ms[0] && two_λb_ == two_ms[1])
            return complex(1.0, 0.0);

        // (cs.two_λa == -two_λa) * (cs.two_λb == -two_λb) && return 2 * cs.ηηηphaseisplus - 1
        if (two_λa_ == -two_ms[0] && two_λb_ == -two_ms[1])
            return complex(2 * ηηηphaseisplus_ - 1, 0.0);

        return complex(0.0, 0.0);
    }
};

using LineshapeFunction = std::function<complex(double)>;

enum class RecouplingType
{
    NoRecoupling,
    ParityRecoupling,
    LSRecoupling
    // Weitere Typen können bei Bedarf hinzugefügt werden
};

// std::shared_ptr<DecayChain> createDecayChainLS(
//     int k, std::function<complex( double )> Xlineshape, const std::string& jp,
//     const ThreeBodyParities& Ps, std::shared_ptr<ThreeBodySystem> tbs , RecouplingType recouplingType);
//  In ThreeBodyDecays.hh: Deklaration der Funktion anpassen
//  In ThreeBodyDecays.hh
//  In ThreeBodyDecays.hh
/**/

// Simple version with automatic LS coupling calculation
std::shared_ptr<DecayChain> createDecayChainLS(
    int k,
    std::function<std::complex<double>(double)> Xlineshape,
    const std::string &jp,
    const ThreeBodyParities &Ps,
    const ThreeBodySystem &tbs);

std::vector<std::shared_ptr<DecayChain>> createDecayChainsLS(
    int k,
    std::function<std::complex<double>(double)> Xlineshape,
    const std::string &jp,
    const ThreeBodyParities &Ps,
    const ThreeBodySystem &tbs);

// Version with explicit recoupling parameters (always disables autoCalculateLS)
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

RecouplingLS createRecouplingFunction(
    RecouplingType type,
    const std::array<int, 2> &params,
    bool parityPhase = true); // Neuer Parameter für ParityRecoupling

// Structure to hold LS coupling information, similar to Julia's implementation
struct LSCoupling
{
    std::array<int, 2> two_ls;
    std::array<int, 2> two_LS;

    // Constructor for convenience
    LSCoupling(const std::array<int, 2> &ls, const std::array<int, 2> &LS) : two_ls(ls), two_LS(LS) {}
};

// Function to calculate possible LS couplings
std::vector<LSCoupling> possible_lsLS(
    const SpinParity &jp,
    const std::array<int, 4> &two_js,
    const ThreeBodyParities &Ps,
    int k);

// Helper function for possible_ls_ij and possible_ls_Rk
std::vector<std::array<int, 2>> possible_ls(
    const SpinParity &jp1,
    const SpinParity &jp2,
    const SpinParity &jp);

// Calculate possible ls couplings for i,j particles
std::vector<std::array<int, 2>> possible_ls_ij(
    const SpinParity &jp,
    const std::array<int, 4> &two_js,
    const ThreeBodyParities &Ps,
    int k);

// Calculate possible LS couplings for resonance and k particle
std::vector<std::array<int, 2>> possible_ls_Rk(
    const SpinParity &jp,
    const std::array<int, 4> &two_js,
    const ThreeBodyParities &Ps,
    int k);

class ThreeBodyDecays
{
public:
    std::vector<std::vector<complex>> aligned_amplitude(const DecayChain &dc, const MandelstamTuple &sigma);
    complex amplitude(const DecayChain &dc, const MandelstamTuple &σs, const std::vector<int> &two_λs, const int &k_amp, const std::vector<int> &refζs = {-1, -1, -1, -1});
    Tensor4D aligned_amplitude4d(const DecayChain &dc, const MandelstamTuple &σs);
    Tensor4D amplitude4d(const DecayChain &dc, const MandelstamTuple &σs, const std::vector<int> &refζs);
    Tensor4Dcomp aligned_amplitude4dcomp(const DecayChain &dc, const MandelstamTuple &σs);
    Tensor4Dcomp amplitude4dcomp(const DecayChain &dc, const MandelstamTuple &σs, const int &k_amp, std::vector<int> refζ = {-1, -1, -1, -1});
    double intensity(const DecayChain &dc, const MandelstamTuple &σs, const int &k_amp, const complex weight = complex(1,1), const std::vector<int> refζs = {-1, -1, -1, -1});

    complex amplitude_recoupling(const RecouplingLS &recoupling, const std::array<int, 2> &two_ms, const std::array<int, 3> &two_js);
    double parseFraction(const std::string &str);
    complex parseComplex(const std::string &str);
    MandelstamTuple x2σs(const std::array<double, 2> &x, ThreeBodyMasses ms,
                         int k);
    double σjofk(double cos_theta, double σk,
                 const std::array<double, 4> &ms_squared, int k);

    double cosθij(const std::array<double, 3> &σs, const std::array<double, 4> &msq, int k);
    double cosθ23(MandelstamTuple σs, ThreeBodyMasses msq);
    double cosθ31(MandelstamTuple σs, ThreeBodyMasses msq);
    double cosθ12(MandelstamTuple σs, ThreeBodyMasses msq);

protected:
    void A_test();

private:
};

// Explicit Wigner rotation functions
double cosζ21_for1(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2);
double cosζ21_for2(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2);
double cosζ13_for1(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2);
double cosζ13_for3(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2);
double cosζ32_for3(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2);
double cosζ32_for2(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2);
double cosζ12_for3(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2);
double cosζ23_for1(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2);
double cosζ31_for2(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2);
double cosζ12_for0(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2);
double cosζ23_for0(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2);
double cosζ31_for0(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2);

#endif
