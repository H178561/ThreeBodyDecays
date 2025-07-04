#include <cmath>
#include <complex>
#include "include/FormFactors.hh"

using complex = std::complex<double>;

namespace FormFactors
{

double BlattWeisskopf(double q, int L, double d)
{
    // Calculate z²
    double z2 = pow((q * d),2);

    // Return appropriate factor based on angular momentum
    switch(L)
    {
        case 0:
            return 1.0;  // Always 1 for S-wave

        case 1:
            // P-wave: sqrt(z² / (1 + z²))
            return std::sqrt(z2 / (1.0 + z2));

        case 2:
            // D-wave: sqrt(z⁴ / (9 + 3z² + z⁴))
            return std::sqrt(z2 * z2 / (9.0 + 3.0 * z2 + z2 * z2));

        /*case 3:
            // F-wave: sqrt(z⁶ / (225 + 45z² + 6z⁴ + z⁶))
            {
                double z4 = z2 * z2;
                double z6 = z4 * z2;
                return std::sqrt(z6 / (225.0 + 45.0 * z2 + 6.0 * z4 + z6));
            }

        case 4:
            // G-wave: sqrt(z⁸ / (11025 + 1575z² + 135z⁴ + 10z⁶ + z⁸))
            {
                double z4 = z2 * z2;
                double z6 = z4 * z2;
                double z8 = z6 * z2;
                return std::sqrt(z8 / (11025.0 + 1575.0 * z2 + 135.0 * z4 + 10.0 * z6 + z8));
            }*/

        default:
            return 1.0;  // Return 1.0 for unsupported L values
    }
}

double breakup(double M, double m1, double m2)
{
    if (M < m1 + m2) return 0.0;
    double lambda = (M * M - (m1 + m2) * (m1 + m2)) * (M * M - (m1 - m2) * (m1 - m2));
    return 0.5 * std::sqrt(lambda) / M;
}

/*
double BlattWeisskopfEnhanced(double q, double q0, int L, double d)
{
    if (q0 <= 0.0) return BlattWeisskopf(q, L, d);

    double BL_q = BlattWeisskopf(q, L, d);
    double BL_q0 = BlattWeisskopf(q0, L, d);

    if (BL_q0 == 0.0) return BL_q;

    return BL_q / BL_q0;
}

std::complex<double> RelativisticBreitWigner(
    double s,
    double mass,
    double width,
    double m1,
    double m2,
    int L,
    double d)
{
    using namespace std::complex_literals;

    double M = std::sqrt(s);
    double q = breakup(M, m1, m2);
    double q0 = breakup(mass, m1, m2);

    if (q0 == 0.0) {
        // Fallback to simple Breit-Wigner
        return 1.0 / (mass * mass - s - 1i * mass * width);
    }

    // Form factor ratio
    double FF_ratio = BlattWeisskopfEnhanced(q, q0, L, d);

    // Mass-dependent width
    double rho = q / M;
    double rho0 = q0 / mass;
    double width_s = width * (rho / rho0) * FF_ratio * FF_ratio;

    return 1.0 / (mass * mass - s - 1i * std::sqrt(s) * width_s);
}



std::function<std::complex<double>(double)> createEnhancedBreitWigner(
    double mass,
    double width,
    double m1,
    double m2,
    int L,
    double d)
{
    return [mass, width, m1, m2, L, d](double s) -> std::complex<double> {
        return RelativisticBreitWigner(s, mass, width, m1, m2, L, d);
    };
}

*/


}
