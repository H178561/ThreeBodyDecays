#ifndef FORM_FACTORS_HH
#define FORM_FACTORS_HH

#include <complex>
#include <functional>

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
        double p0 = FormFactors::breakup(m_, ma_, mb_);
        double FF = FormFactors::BlattWeisskopf(p0, l_, d_);
        double width = gsq_ * p0 * p0 / (m_ * m_ * FF * FF);

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
    double m1_; // mass of first channel particles
    double m2_; // mass of second channel particles
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

#endif
