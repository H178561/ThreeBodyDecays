#ifndef FORM_FACTORS_HH
#define FORM_FACTORS_HH

#include <complex>
#include <functional>
#include <vector>

using complex = std::complex<double>;

namespace FormFactors
{
    // Basic Blatt-Weisskopf form factor
    /**
     * @brief Calculates the Blatt-Weisskopf barrier factor
     * @param q Momentum value
     * @param L Angular momentum
     * @param d Interaction radius (GeV^-1)
     * @return The barrier factor value
     */
    double BlattWeisskopf(double q, int L, double d);

    // Breakup momentum calculation
    /**
     * @brief Calculates the breakup momentum for a two-body decay
     * @param M Total invariant mass of the system
     * @param m1 Mass of the first particle
     * @param m2 Mass of the second particle
     * @return The breakup momentum, or 0 if the decay is not kinematically allowed
     */
    double breakup(double M, double m1, double m2);

    double MomentumPower(double p, int L);


}

#endif
