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

        case 3:
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
            }
        case 5:
            // H-wave: sqrt(z¹⁰ / (1102500 + 157500z² + 13500z⁴ + 1000z⁶ + 45z⁸ + z¹⁰))
            {
                double z4 = z2 * z2;
                double z6 = z4 * z2;
                double z8 = z6 * z2;
                double z10 = z8 * z2;
                return std::sqrt(z10 / (1102500.0 + 157500.0 * z2 + 13500.0 * z4 + 1000.0 * z6 + 45.0 * z8 + z10));
            }

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

// In FormFactors.cpp - Implementation:
double MomentumPower(double p, int L) {
    if (L == 0) return 1.0;
    return std::pow(p, L);
}


double BlattWeisskopfRatio2(double p, double p0, int l, double radius) {
        double z = (p * radius) * (p * radius);
        double z0 = (p0 * radius) * (p0 * radius);

        switch(l) {
            case 0: return 1.0;
            case 1: return (1 + z0) / (1 + z);
            case 2: return (9 + 3*z0 + z0*z0) / (9 + 3*z + z*z);
            // weitere Fälle für höhere l-Werte
            default:
                throw std::invalid_argument("Unsupported angular momentum");
        }
    }

}
