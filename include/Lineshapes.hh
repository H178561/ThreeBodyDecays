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


// BreitWignerMinL class - corresponds to Julia BreitWignerMinL
class BreitWignerMinL : public Lineshape
{
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
        double p = FormFactors::breakup(std::sqrt(sigma), m1_, m2_);
        double p0 = FormFactors::breakup(m_, m1_, m2_);
        double q = FormFactors::breakup(m0_, std::sqrt(sigma), mk_);
        double q0 = FormFactors::breakup(m0_, m_, mk_);

        // Handle edge cases
        if (p0 <= 0.0 || q0 <= 0.0) {
            return complex(0.0, 0.0);
        }

        // Calculate form factors
        double FF_l_p = FormFactors::BlattWeisskopf(p, l_, dR_);
        double FF_l_p0 = FormFactors::BlattWeisskopf(p0, l_, dR_);
        double FF_minL_q = FormFactors::BlattWeisskopf(q, minL_, dLambdac_);
        double FF_minL_q0 = FormFactors::BlattWeisskopf(q0, minL_, dLambdac_);

        double F2_l = (FF_l_p * FF_l_p) / (FF_l_p0 * FF_l_p0);
        double F2_minL = (FF_minL_q * FF_minL_q) / (FF_minL_q0 * FF_minL_q0);

        // Calculate running width
        double p_ratio_power = std::pow(p / p0, 2 * l_ + 1);
        double gamma = gamma0_ * p_ratio_power * m_ / std::sqrt(sigma) * F2_l;

        // Calculate momentum factors
        double p_factor = std::pow(p / p0, l_);
        double q_factor = std::pow(q / q0, minL_);

        // Breit-Wigner denominator
        complex denominator = complex(m_ * m_ - sigma, -m_ * gamma);

        // Complete amplitude
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

// Factory function for BreitWignerMinL
inline std::function<complex(double)> make_breit_wigner_minl(
    double mass, double width, int l, int minL,
    double m1, double m2, double mk, double m0,
    double dR = 1.5, double dLambdac = 5.0)
{
    BreitWignerMinL bw(mass, width, l, minL, m1, m2, mk, m0, dR, dLambdac);
    return [bw](double sigma) { return bw(sigma); };
}

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

#endif
