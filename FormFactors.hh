#ifndef FORM_FACTORS_HH
#define FORM_FACTORS_HH

#include <complex>
#include <functional>
#include <vector>

using complex = std::complex<double>;

namespace FormFactors
{
    // Basic Blatt-Weisskopf form factor
    double BlattWeisskopf(double q, int L, double d);

    // Breakup momentum calculation
    double breakup(double M, double m1, double m2);

    /*

    // Enhanced Blatt-Weisskopf with threshold behavior
    double BlattWeisskopfEnhanced(double q, double q0, int L, double d);

    // Relativistic Breit-Wigner with Blatt-Weisskopf form factors
    std::complex<double> RelativisticBreitWigner(
        double s,
        double mass,
        double width,
        double m1,
        double m2,
        int L,
        double d
    );

    // Create enhanced lineshape function with form factors
    std::function<std::complex<double>(double)> createEnhancedBreitWigner(
        double mass,
        double width,
        double m1,
        double m2,
        int L,
        double d = 1.0
    );


    // Flatte lineshape (for coupled channels)
    complex Flatte(
        double s,
        double mass,
        double g1, double m1a, double m1b,
        double g2, double m2a, double m2b
    );



    // Create Flatte lineshape function
    std::function<complex(double)> createFlatte(
        double mass,
        double g1, double m1a, double m1b,
        double g2, double m2a, double m2b
    );
    */
}

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

/*
function (bw::BreitWigner)(σ::Number)
    @unpack m, ma, mb, l, d = bw
    _p0 = breakup(m, ma, mb)
    FF = BlattWeisskopf{l}(d)
    gsq = m * bw.Γ / (2_p0) * m / FF(_p0)^2
    mbw = MultichannelBreitWigner(m, SVector((; gsq, ma, mb, l, d)))
    mbw(σ)
end
*/

/*

BW(σ, m, Γ) = 1 / (m^2 - σ - 1im * m * Γ)

# ## MultichannelBreitWigner

@with_kw struct MultichannelBreitWigner{N} <: AbstractFlexFunc
    m::Float64
    channels::SVector{N,<:NamedTuple{(:gsq, :ma, :mb, :l, :d)}}
end

# MultichannelBreitWigner constructor from Vector
function MultichannelBreitWigner(
    m::Real,
    channels::Vector{<:NamedTuple{(:gsq, :ma, :mb, :l, :d)}},
)
    N = length(channels)
    return MultichannelBreitWigner(m, SVector{N}(channels...))
end

function (bw::MultichannelBreitWigner)(σ::Number)
    m0 = bw.m
    mΓ = sum(bw.channels) do channel
        @unpack gsq, ma, mb, l, d = channel
        FF = BlattWeisskopf{l}(d)
        _p = breakup(sqrt(σ), ma, mb)
        gsq * 2_p / sqrt(σ) * FF(_p)^2
    end
    BW(σ, m0, mΓ / m0)
end
(bw::MultichannelBreitWigner)(σ::Real) = (bw)(σ + 1im * eps())


function MultichannelBreitWigner(m::Real, Γ::Real, ma::Number, mb::Number, l::Int, d::Real)
    _p0 = breakup(m, ma, mb)
    FF = BlattWeisskopf{l}(d)
    gsq = m * Γ / (2_p0) * m / FF(_p0)^2
    return MultichannelBreitWigner(m, SVector((; gsq, ma, mb, l, d)))
end
MultichannelBreitWigner(m::Float64, Γ::Float64) =
    MultichannelBreitWigner{1}(m, Γ, 0.0, 0.0, 0, 1.0)

*/

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
                                    complex(0.0, m_ * width));
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

// Flatte class
class Flatte : public Lineshape
{
public:
    Flatte(double mass, double g1, double g2, double m1, double m2) : m_(mass), g1_(g1), g2_(g2), m1_(m1), m2_(m2)
    {
    }

    complex operator()(double sigma) const override
    {
        // Phase space factors
        double rho1 = std::sqrt(1.0 - 4.0 * m1_ * m1_ / sigma);
        double rho2 = std::sqrt(1.0 - 4.0 * m2_ * m2_ / sigma);

        // Running width
        complex running_width =
            complex(0.0, m_ * (g1_ * rho1 + g2_ * rho2));

        return complex(1.0, 0.0) /
               (complex(m_ * m_ - sigma, 0.0) - running_width);
    }

private:
    double m_;  // mass
    double g1_; // coupling to first channel
    double g2_; // coupling to second channel
    double m1_; // mass of first channel particle
    double m2_; // mass of second channel particle
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

inline std::function<complex(double)> make_flatte(double mass, double g1,
                                                  double g2, double m1,
                                                  double m2)
{
    Flatte f(mass, g1, g2, m1, m2);
    return [f](double sigma)
    { return f(sigma); };
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
#endif
