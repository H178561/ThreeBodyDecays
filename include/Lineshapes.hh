#ifndef LINESHAPES_HH
#define LINESHAPES_HH

#include <complex>
#include <functional>
#include <vector>
#include "FormFactors.hh"

using complex = std::complex<double>;


// Base class for lineshapes
class Lineshape
{
public:
    virtual ~Lineshape() = default;
    virtual complex operator()(double sigma) const = 0;
};

// BreitWigner class
class BreitWigner : public Lineshape
{
public:
    BreitWigner(double mass, double width) : m_(mass), gamma_(width) {}

    complex operator()(double sigma) const override
    {
        return complex(1.0, 0.0) / (complex(m_ * m_ - sigma, 0.0) -
                                    complex(0.0, m_ * gamma_));
    }

private:
    double m_;     // mass
    double gamma_; // width
};



// Multichannel Breit-Wigner class
class MultichannelBreitWigner : public Lineshape
{
public:
    MultichannelBreitWigner(double mass, double gsq, double ma, double mb, int l, double d)
        : m_(mass), gsq_(gsq), ma_(ma), mb_(mb), l_(l), d_(d)
    {
    }

    complex operator()(double sigma) const override
    {
        double p0 = FormFactors::breakup(std::sqrt(sigma), ma_, mb_);
        double FF = FormFactors::BlattWeisskopf(p0, l_, d_);
        double width = gsq_ * 2.0 * p0 / (std::sqrt(sigma) * (FF * FF));

        return complex(1.0, 0.0) / (complex(m_ * m_ - sigma, 0.0) -
                                    complex(0.0, width));
    }

private:
    double m_;   // mass
    double gsq_; // squared coupling constant
    double ma_;  // mass of first channel particle
    double mb_;  // mass of second channel particle
    int l_;      // orbital angular momentum
    double d_;   // Blatt-Weisskopf parameter
};

class BreitWignerExtended : public Lineshape
{
public:
    BreitWignerExtended(double mass, double width, double ma, double mb, int l, double d)
        : m_(mass), gamma_(width), ma_(ma), mb_(mb), l_(l), d_(d)
    {
    }

    complex operator()(double sigma) const override
    {
        double p0 = FormFactors::breakup(m_, ma_, mb_);
        double FF = FormFactors::BlattWeisskopf(p0, l_, d_);
        double gsq = m_ * gamma_ / (2.0 * p0) * m_ / (FF * FF);

        // use of multichannel Breit-Wigner
        MultichannelBreitWigner mbw(m_, gsq, ma_, mb_, l_, d_);
        return mbw(sigma);
    }

private:
    double m_;     // mass
    double gamma_; // width
    double ma_;    // mass of first channel particle
    double mb_;    // mass of second channel particle
    int l_;        // orbital angular momentum
    double d_;     // Blatt-Weisskopf parameter
};

// BuggBW class
class BuggBW : public Lineshape
{
public:
    BuggBW(double mass, double width, double s0, double s1, double s2) : mass_(mass), width_(width), s0_(s0), s1_(s1), s2_(s2)
    {
    }

    complex operator()(double sigma) const override
    {
        // Bugg BW formula: 1/(m2 - s - i * m * (s - sA)/(m2 - sA) * Γ * exp(-s2 * s))
        complex denominator = mass_ * mass_ - sigma;
        complex numerator = mass_ * (sigma - s0_) / (mass_ * mass_ - s0_);
        complex exponential = std::exp(-s2_ * sigma);
        complex width_term = width_ * exponential;

        return 1.0 /
               (denominator - complex(0.0, 1.0) * numerator * width_term);
    }

private:
    double mass_;
    double width_;
    double s0_; // Adler zero
    double s1_; // threshold
    double s2_; // exponential factor
};


// Helper functions to create lineshape functions
inline std::function<complex(double)> make_breit_wigner(double mass,
                                                        double width)
{
    BreitWigner bw(mass, width);
    return [bw](double sigma)
    { return bw(sigma); };
}

inline std::function<complex(double)> make_bugg_bw(double mass,
                                                   double width, double s0,
                                                   double s1, double s2)
{
    return BuggBW(mass, width, s0, s1, s2);
}

// Struktur zur Repräsentation eines Kanals
struct Channel
{
    double gsq; // quadrierte Kopplungskonstante
    double ma;  // Masse des ersten Teilchens
    double mb;  // Masse des zweiten Teilchens
    int l;      // Drehimpuls
    double d;   // Blatt-Weisskopf-Parameter
};

// Erweiterte MultichannelBreitWigner Klasse für mehrere Kanäle
class MultichannelBreitWignerMulti : public Lineshape
{
public:
    // Konstruktor mit einem Vektor von Kanälen
    MultichannelBreitWignerMulti(double mass, const std::vector<Channel> &channels)
        : m_(mass), channels_(channels)
    {
    }

    // Konstruktor für einen einzelnen Kanal (Kompatibilität mit vorhandenem Code)
    MultichannelBreitWignerMulti(double mass, double gsq, double ma, double mb, int l, double d)
        : m_(mass)
    {
        channels_.push_back({gsq, ma, mb, l, d});
    }

    complex operator()(double sigma) const override
    {
        // Berechnung der laufenden Breite über alle Kanäle
        double total_width = 0.0;

        for (const auto &channel : channels_)
        {
            double p = FormFactors::breakup(std::sqrt(sigma), channel.ma, channel.mb);
            double FF = FormFactors::BlattWeisskopf(p, channel.l, channel.d);

            // Beitrag dieses Kanals zur Gesamtbreite
            total_width += channel.gsq * 2.0 * p / std::sqrt(sigma) * (FF * FF);
        }

        // Breit-Wigner-Formel mit der berechneten Breite
        return complex(1.0, 0.0) /
               (complex(m_ * m_ - sigma, 0.0) - complex(0.0, m_ * total_width / m_));
    }

private:
    double m_;                      // Masse
    std::vector<Channel> channels_; // Vektor der Kanäle
};

// Factory-Funktion für mehrere Kanäle
inline std::function<complex(double)> make_multichannel_bw(
    double mass, const std::vector<Channel> &channels)
{
    MultichannelBreitWignerMulti mbw(mass, channels);
    return [mbw](double sigma)
    { return mbw(sigma); };
}

// Factory-Funktion für einen einzelnen Kanal
inline std::function<complex(double)> make_multichannel_bw_single(
    double mass, double width, double ma, double mb, int l, double d)
{
    // Berechnung von gsq wie in der Julia-Implementierung
    double p0 = FormFactors::breakup(mass, ma, mb);
    double FF = FormFactors::BlattWeisskopf(p0, l, d);
    double gsq = mass * width / (2.0 * p0) * mass / (FF * FF);

    MultichannelBreitWignerMulti mbw(mass, gsq, ma, mb, l, d);
    return [mbw](double sigma)
    { return mbw(sigma); };
}

// Flatte class
class Flatte : public Lineshape
{
public:
    Flatte(double mass, double gsq1, double ma1, double mb1, int l1, double d1,
           double gsq2, double ma2, double mb2, int l2, double d2)
        : m_(mass)
    {
        // Create two channels with l=0, d=1.0
        channels_.push_back(Channel{gsq1, ma1, mb1, l1, d1});
        channels_.push_back(Channel{gsq2, ma2, mb2, l2, d2});
    }

    complex operator()(double sigma) const override
    {
        MultichannelBreitWignerMulti mbw(m_, channels_);
        return mbw(sigma);
    }

private:
    double m_;                      // mass
    std::vector<Channel> channels_; // vector of channels
};

inline std::function<complex(double)> make_flatte(double mass,
                                                  double gsq1, double ma1, double mb1, int l1, double d1,
                                                  double gsq2, double ma2, double mb2, int l2, double d2)
{
    Flatte f(mass, gsq1, ma1, mb1, l1, d1, gsq2, ma2, mb2, l2, d2);
    return [f](double sigma)
    { return f(sigma); };
}



//////////////////////// New Lineshapes ////////////////////////




// BreitWignerMinL class - corresponds to Julia BreitWignerMinL
class BreitWignerMinL : public Lineshape
{
private:
    // Helper function to calculate F² like in Julia
    double F2_function(int l, double p, double p0, double d) const {
        if (l == 0) return 1.0;

        double pR = p * d;
        double p0R = p0 * d;
        double pR2 = pR * pR;
        double p0R2 = p0R * p0R;

        if (l == 1) {
            return (1.0 + p0R2) / (1.0 + pR2);
        } else if (l == 2) {
            return (9.0 + 3.0*p0R2 + p0R2*p0R2) / (9.0 + 3.0*pR2 + pR2*pR2);
        } else if (l == 3) {
            double p0R4 = p0R2 * p0R2;
            double p0R6 = p0R4 * p0R2;
            double pR4 = pR2 * pR2;
            double pR6 = pR4 * pR2;
            return (225.0 + 45.0*p0R2 + 6.0*p0R4 + p0R6) / (225.0 + 45.0*pR2 + 6.0*pR4 + pR6);
        }
        return 1.0; // fallback
    }

    double my_breakup(double s, double m1_sq, double m2_sq) const
    {
        if (s < 0.0) return 0.0;

        // Kallen function: λ(s, m1², m2²) = s² + m1⁴ + m2⁴ - 2s*m1² - 2s*m2² - 2m1²*m2²
        double lambda = s*s + m1_sq*m1_sq + m2_sq*m2_sq - 2.0*s*m1_sq - 2.0*s*m2_sq - 2.0*m1_sq*m2_sq;

        if (lambda < 0.0) return 0.0;

        return std::sqrt(lambda) / (2.0 * std::sqrt(s));
    }

public:
    BreitWignerMinL(double mass, double width, int l, int minL,
                    double m1, double m2, double mk, double m0,
                    double dR = 1.5, double dLambdac = 5.0)
        : m_(mass), gamma0_(width), l_(l), minL_(minL),
          m1_(m1), m2_(m2), mk_(mk), m0_(m0),
          dR_(dR), dLambdac_(dLambdac)
    {
    }

    complex operator()(double sigma) const override
    {
        // Calculate breakup momenta
        /*double p = FormFactors::breakup(std::sqrt(sigma), m1_, m2_);
        double p0 = FormFactors::breakup(m_, m1_, m2_);
        double q = FormFactors::breakup(m0_, std::sqrt(sigma), mk_);
        double q0 = FormFactors::breakup(m0_, m_, mk_);
        */
        double p = my_breakup(sigma, m1_*m1_, m2_*m2_);      // my_breakup(σ, m1², m2²)
        double p0 = my_breakup(m_*m_, m1_*m1_, m2_*m2_);     // my_breakup(m², m1², m2²)
        double q = my_breakup(m0_*m0_, sigma, mk_*mk_);      // my_breakup(m0², σ, mk²)
        double q0 = my_breakup(m0_*m0_, m_*m_, mk_*mk_);     // my_breakup(m0², m², mk²)


        // Handle edge cases
        if (p0 <= 0.0 || q0 <= 0.0) {
            return complex(0.0, 0.0);
        }

         // Calculate F² functions (Julia style)
        double F2_l = F2_function(l_, p, p0, dR_);             // F²(l, p, p0, dR)
        double F2_minL = F2_function(minL_, q, q0, dLambdac_); // F²(minL, q, q0, dΛc)

        // Calculate running width: Γ = Γ₀ * (p / p0)^(2l + 1) * m / sqrt(σ) * F²(l, p, p0, dR)
        double p_ratio_power = std::pow(p / p0, 2 * l_ + 1);
        double gamma = gamma0_ * p_ratio_power * m_ / std::sqrt(sigma) * F2_l;

        // Momentum factors
        double p_factor = std::pow(p / p0, l_);        // (p / p0)^l
        double q_factor = std::pow(q / q0, minL_);     // (q / q0)^minL

        // Julia formula: 1 / (m² - σ - 1im * m * Γ) * (p / p0)^l * (q / q0)^minL * sqrt(F²(l, p, p0, dR) * F²(minL, q, q0, dΛc))
        complex denominator = complex(m_ * m_ - sigma, -m_ * gamma);  // Note: -1im * m * Γ
        complex bw_part = complex(1.0, 0.0) / denominator;

        double momentum_ff_part = p_factor * q_factor * std::sqrt(F2_l * F2_minL);

        return bw_part * momentum_ff_part;
    }

private:
    double m_;          // mass
    double gamma0_;     // width at pole
    int l_;             // orbital angular momentum (decay)
    int minL_;          // minimal orbital angular momentum (production)
    double m1_, m2_;    // masses of decay products
    double mk_;         // mass of spectator particle
    double m0_;         // mass of parent particle
    double dR_;         // interaction radius for decay (default 1.5 GeV^-1)
    double dLambdac_;   // interaction radius for production (default 5.0 GeV^-1)
};

// BuggBreitWignerMinL class - corresponds to Julia BuggBreitWignerMinL
class BuggBreitWignerMinL : public Lineshape
{
private:
    // Constants from Julia code
    static constexpr double mK = 0.493677;   // Kaon mass
    static constexpr double mπ = 0.13957018; // Pion mass

public:
    BuggBreitWignerMinL(double mass, double width, double gamma, int l, int minL,
                       double m1, double m2, double mk, double m0)
        : m_(mass), gamma0_(width), gamma_(gamma), l_(l), minL_(minL),
          m1_(m1), m2_(m2), mk_(mk), m0_(m0)
    {
    }

    // Constructor with default gamma = 1.1 (like Julia default)
    BuggBreitWignerMinL(double mass, double width, int l, int minL,
                       double m1, double m2, double mk, double m0)
        : BuggBreitWignerMinL(mass, width, 1.1, l, minL, m1, m2, mk, m0)
    {
    }

    complex operator()(double sigma) const override
    {
        // Calculate σA = mK² - mπ²/2 (Adler zero)
        double sigma_A = mK * mK - (mπ * mπ) / 2.0;

        // Calculate running width: Γ = (σ - σA)/(m² - σA) * Γ₀ * exp(-γ * σ)
        double width_factor = (sigma - sigma_A) / (m_ * m_ - sigma_A);
        double exponential_factor = std::exp(-gamma_ * sigma);
        double gamma = width_factor * gamma0_ * exponential_factor;

        // Bugg BW formula: 1 / (m² - σ - 1im * m * Γ)
        complex denominator = complex(m_ * m_ - sigma, -m_ * gamma);

        return complex(1.0, 0.0) / denominator;
    }

private:
    double m_;          // mass
    double gamma0_;     // width at pole (Γ₀)
    double gamma_;      // exponential parameter (γ)
    int l_;             // orbital angular momentum (decay)
    int minL_;          // minimal orbital angular momentum (production)
    double m1_, m2_;    // masses of decay products
    double mk_;         // mass of spectator particle
    double m0_;         // mass of parent particle
};

// Factory function for BreitWignerMinL
inline std::function<complex(double)> make_breit_wigner_minl(
    double mass, double width, int l, int minL,
    double m1, double m2, double mk, double m0,
    double dR = 1.5, double dLambdac = 5.0)
{
    BreitWignerMinL bw(mass, width, l, minL, m1, m2, mk, m0, dR, dLambdac);
    return [bw](double sigma) { return bw(sigma); };
}

// Factory function for BuggBreitWignerMinL
inline std::function<complex(double)> make_bugg_breit_wigner_minl(
    double mass, double width, double gamma, int l, int minL,
    double m1, double m2, double mk, double m0)
{
    BuggBreitWignerMinL bw(mass, width, gamma, l, minL, m1, m2, mk, m0);
    return [bw](double sigma) { return bw(sigma); };
}

// Factory function with default gamma = 1.1
inline std::function<complex(double)> make_bugg_breit_wigner_minl(
    double mass, double width, int l, int minL,
    double m1, double m2, double mk, double m0)
{
    BuggBreitWignerMinL bw(mass, width, l, minL, m1, m2, mk, m0);
    return [bw](double sigma) { return bw(sigma); };
}


// Flatte1405 class - corresponds to Julia Flatte1405
class Flatte1405 : public Lineshape
{
private:
    // Constants for Λ(1405) → πΣ channels
    static constexpr double mπ = 0.13957018;  // Pion mass
    static constexpr double mΣ = 1.18937;       // Sigma mass (approximate)


    double my_breakup(double s, double m1_sq, double m2_sq) const
    {
        if (s < 0.0) return 0.0;

        // Kallen function: λ(s, m1², m2²) = s² + m1⁴ + m2⁴ - 2s*m1² - 2s*m2² - 2m1²*m2²
        double lambda = s*s + m1_sq*m1_sq + m2_sq*m2_sq - 2.0*s*m1_sq - 2.0*s*m2_sq - 2.0*m1_sq*m2_sq;

        if (lambda < 0.0) return 0.0;

        return std::sqrt(lambda) / (2.0 * std::sqrt(s));
    }

public:
    Flatte1405(double mass, double width, int l, int minL,
               double m1, double m2, double mk, double m0)
        : m_(mass), gamma0_(width), l_(l), minL_(minL),
          m1_(m1), m2_(m2), mk_(mk), m0_(m0)
    {
    }

    complex operator()(double sigma) const override
    {
        // Channel 1: m1, m2 (from constructor)
        double p = my_breakup(sigma, m1_*m1_, m2_*m2_);
        double p0 = my_breakup(m_*m_, mπ*mπ, mΣ*mΣ);

        // Channel 2: π, Σ (fixed)
        double p_prime = my_breakup(sigma, mπ*mπ, mΣ*mΣ);
        double p0_prime = my_breakup(m_*m_, mπ*mπ, mΣ*mΣ);

        // Handle edge cases
        if (p0 <= 0.0 || p0_prime <= 0.0) {
            return complex(0.0, 0.0);
        }

        // Calculate running widths
        double gamma1 = gamma0_ * (p / p0) * m_ / std::sqrt(sigma);
        double gamma2 = gamma0_ * (p_prime / p0_prime) * m_ / std::sqrt(sigma);
        double gamma = gamma1 + gamma2;

        // Flatte formula: 1 / (m² - σ - 1im * m * Γ)
        complex denominator = complex(m_*m_ - sigma, -m_ * gamma);

        return complex(1.0, 0.0) / denominator;
    }

private:
    double m_;          // mass
    double gamma0_;     // width at pole (Γ₀)
    int l_;             // orbital angular momentum (decay)
    int minL_;          // minimal orbital angular momentum (production)
    double m1_, m2_;    // masses of decay products
    double mk_;         // mass of spectator particle
    double m0_;         // mass of parent particle
};

// Factory function for Flatte1405
inline std::function<complex(double)> make_flatte1405(
    double mass, double width, int l, int minL,
    double m1, double m2, double mk, double m0)
{
    Flatte1405 f(mass, width, l, minL, m1, m2, mk, m0);
    return [f](double sigma) { return f(sigma); };
}

// L1670Flatte class - corresponds to Julia L1670Flatte
class L1670Flatte : public Lineshape
{
private:
    // Helper function: k(m, ma, mb) = breakup(m², ma², mb²)
    double k_function(double m, double ma, double mb) const {
        return my_breakup(m*m, ma*ma, mb*mb);
    }

    double my_breakup(double s, double m1_sq, double m2_sq) const
    {
        if (s < 0.0) return 0.0;

        // Kallen function: λ(s, m1², m2²) = s² + m1⁴ + m2⁴ - 2s*m1² - 2s*m2² - 2m1²*m2²
        double lambda = s*s + m1_sq*m1_sq + m2_sq*m2_sq - 2.0*s*m1_sq - 2.0*s*m2_sq - 2.0*m1_sq*m2_sq;

        if (lambda < 0.0) return 0.0;

        return std::sqrt(lambda) / (2.0 * std::sqrt(s));
    }

public:
    L1670Flatte(double mass, double width, int l, int minL,
                double m1, double m2, double mk, double m0)
        : mf_(mass), width_(width), l_(l), minL_(minL),
          m1_(m1), m2_(m2), mk_(mk), m0_(m0)
    {
    }

    complex operator()(double sigma) const override
    {
        // Constants from Julia implementation
        const double Gamma0 = 0.0272;
        const double g = 0.258;
        const double ma2 = 1.115683;
        const double mb2 = 0.547862;
        const complex iepsilon(0.0, 1e-8);  // iϵ = 1e-8im

        double m = std::sqrt(sigma);

        // Calculate k(m + iϵ, ma2, mb2)
        // For small imaginary part, we can approximate: k(m + iϵ) ≈ k(m) + iϵ * dk/dm
        double k_real = k_function(m, ma2, mb2);

        // Julia formula: D = mf - m - 0.5im * (Γ₀ + g * k(m + iϵ, ma2, mb2))
        complex D = complex(mf_ - m, -0.5 * (Gamma0 + g * k_real));

        // Return -1/D (note the minus sign as mentioned in Julia comment)
        return complex(-1.0, 0.0) / D;
    }

private:
    double mf_;         // fitted mass (mf parameter)
    double width_;      // width parameter (unused in this implementation)
    int l_;             // orbital angular momentum (decay)
    int minL_;          // minimal orbital angular momentum (production)
    double m1_, m2_;    // masses of decay products
    double mk_;         // mass of spectator particle
    double m0_;         // mass of parent particle
};

// Factory function for L1670Flatte
inline std::function<complex(double)> make_l1670_flatte(
    double mass, double width, int l, int minL,
    double m1, double m2, double mk, double m0)
{
    L1670Flatte f(mass, width, l, minL, m1, m2, mk, m0);
    return [f](double sigma) { return f(sigma); };
}

#endif
