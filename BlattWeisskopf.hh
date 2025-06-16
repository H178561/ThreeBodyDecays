#ifndef BLATT_WEISSKOPF_HH
#define BLATT_WEISSKOPF_HH

#include <complex>
#include <functional>

namespace ThreeBodyDecays
{
    // Basic Blatt-Weisskopf form factor
    double BlattWeisskopf(double q, int L, double d);

    // Breakup momentum calculation
    double breakup(double M, double m1, double m2);

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

    // Flatte lineshape (for coupled channels)
    std::complex<double> Flatte(
        double s,
        double mass,
        double g1, double m1a, double m1b,
        double g2, double m2a, double m2b
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

    // Create Flatte lineshape function
    std::function<std::complex<double>(double)> createFlatte(
        double mass,
        double g1, double m1a, double m1b,
        double g2, double m2a, double m2b
    );
}

#endif
