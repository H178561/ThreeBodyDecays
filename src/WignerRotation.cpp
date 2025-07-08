#include "include/WignerRotation.hh"
#include <cmath>
#include <stdexcept>
#include <iostream>

extern bool debug; // External debug flag

// Helper function to get indices i and j from k
std::pair<int, int> ij_from_k(int k)
{
    switch (k)
    {
    case 1:
        return {2, 3};
    case 2:
        return {3, 1};
    case 3:
        return {1, 2};
    default:
        throw std::invalid_argument(
            "Invalid value for k. Must be 1, 2, or 3.");
    }
}

// Function to compute i, j, k indices from Wigner rotation
std::tuple<int, int, int> ijk(const AbstractWignerRotation &wr)
{
    int k = wr.get_k();

    auto [i, j] = ij_from_k(k);
    return {i, j, k};
}

// Helper function to determine if a system is sequential
bool issequential(int i, int j)
{
    int diff = j - i;
    return (diff == 1 || diff == -2);
}

// Helper function to create Wigner rotation
std::unique_ptr<AbstractWignerRotation> wr(int system_a, int reference_b,
                                           int particle_c)
{
    if (debug)
        std::cout << "Creating Wigner rotation for system_a: " << system_a << ", reference_b: " << reference_b << ", particle_c: " << particle_c << std::endl;
    if (system_a == reference_b)
    {
        return std::make_unique<TrivialWignerRotation>(particle_c);
    }

    bool S = issequential(system_a, reference_b);
    int A = S ? system_a : reference_b;
    int B = S ? reference_b : system_a;

    if (particle_c == 0)
    {
        return std::make_unique<Arg0WignerRotation>(A, S);
    }

    if (particle_c != system_a && particle_c != reference_b)
    {
        return std::make_unique<Arg3WignerRotation>(particle_c, S);
    }

    bool T = (particle_c == A);
    return std::make_unique<Arg2WignerRotation>(particle_c, !S, T);
}

// Implementation of cos_zeta for TrivialWignerRotation
double TrivialWignerRotation::cos_zeta(const std::array<double, 3> &sigma,
                                       const std::array<double, 4> &ms2) const
{
    return 1.0; // Trivial rotation always returns 1
}

// Implementation of cos_zeta for Arg0WignerRotation
double Arg0WignerRotation::cos_zeta(const std::array<double, 3> &sigma,
                                    const std::array<double, 4> &ms2) const
{
    auto [i, j, k] = ijk(*this);
    // Keep 1-based indexing like Julia
    i--;
    j--;
    k--; // Convert to 0-based indexing for array access

    double s = ms2[3]; // Parent mass squared
    double EE4s = (s + ms2[i] - sigma[i]) * (s + ms2[k] - sigma[k]);
    double pp4s = std::sqrt(Kallen(s, ms2[i], sigma[i]) *
                            Kallen(s, ms2[k], sigma[k]));
    pp4s = (pp4s < 0) ? 0.0 : pp4s; // Handle numerical errors
    double rest = sigma[j] - ms2[i] - ms2[k];
    if (debug)
        std::cout << "wr0" << std::endl;

    return (EE4s - 2 * s * rest) / pp4s;
}

// Implementation of operator() for Arg0WignerRotation
double Arg0WignerRotation::operator()(double w, const std::array<double, 3> &sigma,
                                     const std::array<double, 4> &ms2) const
{
    auto [i, j, k] = ijk(*this);
    i--;
    j--;
    k--; // Convert to 0-based indexing for array access

    double s = ms2[3]; // Parent mass squared
    double EE4s = (s + ms2[i] - sigma[i]) * (s + ms2[k] - sigma[k]);
    double pp4s = std::sqrt(Kallen(s, ms2[i], sigma[i]) *
                            Kallen(s, ms2[k], sigma[k]));
    pp4s = (pp4s < 0) ? 0.0 : pp4s; // Handle numerical errors
    double rest = sigma[j] - ms2[i] - ms2[k];

    return (EE4s - 2 * s * rest) / pp4s;
}

// Implementation of cos_zeta for Arg2WignerRotation
double Arg2WignerRotation::cos_zeta(const std::array<double, 3> &sigma,
                                    const std::array<double, 4> &ms2) const
{
    auto [i, j, k] = ijk(*this);
    if (debug)
        std::cout << "wr2" << std::endl;
    // Keep 1-based indexing like Julia
    i--;
    j--;
    k--; // Convert to 0-based indexing for array access
    if (debug)
        std::cout << i << " " << j << " " << k << std::endl;
    // Swap i and j if not even (like Julia)
    if (!is_even())
    {
        std::swap(i, j);
    }
    if (debug)
        std::cout << i << " " << j << std::endl;

    // Handle massless case
    if (std::abs(ms2[k]) < 1e-10)
    {
        return 1.0;
    }

    double s = ms2[3]; // Parent mass squared
    double EE4mksq = (s + ms2[k] - sigma[k]) * (sigma[i] - ms2[k] - ms2[j]);
    double pp4mksq = std::sqrt(Kallen(s, ms2[k], sigma[k]) *
                               Kallen(ms2[k], ms2[j], sigma[i]));
    pp4mksq = (pp4mksq < 0) ? 0.0 : pp4mksq; // Handle numerical errors
    double rest = sigma[j] - s - ms2[j];
    if (debug)
        std::cout << s << " " << EE4mksq << " " << pp4mksq << " " << rest << " " << (2 * ms2[k] * rest + EE4mksq) / pp4mksq << std::endl;
    return (2 * ms2[k] * rest + EE4mksq) / pp4mksq;
}

// Implementation of operator() for Arg2WignerRotation
double Arg2WignerRotation::operator()(double w, const std::array<double, 3> &sigma,
                                     const std::array<double, 4> &ms2) const
{
    auto [i, j, k] = ijk(*this);
    if (!is_even())
    {
        std::swap(i, j);
    }
    i--;
    j--;
    k--; // Convert to 0-based indexing for array access

    if (std::abs(ms2[k]) < 1e-10)
        return 1.0; // Handle massless case

    double s = ms2[3];
    double EE4mksq = (s + ms2[k] - sigma[k]) *
                     (sigma[i] - ms2[k] - ms2[j]);
    double pp4mksq = std::sqrt(Kallen(s, ms2[k], sigma[k]) *
                               Kallen(ms2[k], ms2[j], sigma[i]));
    pp4mksq = (pp4mksq < 0) ? 0.0 : pp4mksq; // Handle numerical errors
    double rest = sigma[j] - s - ms2[j];
    return (2 * ms2[k] * rest + EE4mksq) / pp4mksq;
}

// Implementation of cos_zeta for Arg3WignerRotation
double Arg3WignerRotation::cos_zeta(const std::array<double, 3> &sigma,
                                    const std::array<double, 4> &ms2) const
{
    auto [i, j, k] = ijk(*this);
    if (debug)
        std::cout << "wr3" << std::endl;
    // Keep 1-based indexing like Julia
    i--;
    j--;
    k--; // Convert to 0-based indexing for array access

    // Handle massless case
    if (std::abs(ms2[k]) < 1e-10)
    {
        return 1.0;
    }

    double s = ms2[3]; // Parent mass squared
    double EE4m1sq = (sigma[i] - ms2[j] - ms2[k]) *
                     (sigma[j] - ms2[k] - ms2[i]);
    double pp4m1sq = std::sqrt(Kallen(sigma[i], ms2[j], ms2[k]) *
                               Kallen(sigma[j], ms2[k], ms2[i]));
    pp4m1sq = (pp4m1sq < 0) ? 0.0 : pp4m1sq; // Handle numerical errors
    double rest = ms2[i] + ms2[j] - sigma[k];
    return (2 * ms2[k] * rest + EE4m1sq) / pp4m1sq;
}

// Implementation of operator() for Arg3WignerRotation
double Arg3WignerRotation::operator()(double w, const std::array<double, 3> &sigma,
                                     const std::array<double, 4> &ms2) const
{
    auto [i, j, k] = ijk(*this);
    i--;
    j--;
    k--; // Convert to 0-based indexing for array access

    if (std::abs(ms2[k]) < 1e-10)
        return 1.0; // Handle massless case

    double EE4m1sq = (sigma[i] - ms2[j] - ms2[k]) *
                     (sigma[j] - ms2[k] - ms2[i]);
    double pp4m1sq = std::sqrt(Kallen(sigma[i], ms2[j], ms2[k]) *
                               Kallen(sigma[j], ms2[k], ms2[i]));
    pp4m1sq = (pp4m1sq < 0) ? 0.0 : pp4m1sq; // Handle numerical errors
    double rest = ms2[i] + ms2[j] - sigma[k];
    return (2 * ms2[k] * rest + EE4m1sq) / pp4m1sq;
}

// Implementation of explicit Wigner rotation functions
double cosζ21_for1(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(2, 1, 1))(0, σs, ms2);
}

double cosζ21_for2(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(2, 1, 2))(0, σs, ms2);
}

double cosζ13_for1(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(1, 3, 1))(0, σs, ms2);
}

double cosζ13_for3(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(1, 3, 3))(0, σs, ms2);
}

double cosζ32_for3(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(3, 2, 3))(0, σs, ms2);
}

double cosζ32_for2(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(3, 2, 2))(0, σs, ms2);
}

double cosζ12_for3(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(1, 2, 3))(0, σs, ms2);
}

double cosζ23_for1(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(2, 3, 1))(0, σs, ms2);
}

double cosζ31_for2(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(3, 1, 2))(0, σs, ms2);
}

double cosζ12_for0(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(1, 2, 0))(0, σs, ms2);
}

double cosζ23_for0(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(2, 3, 0))(0, σs, ms2);
}

double cosζ31_for0(const std::array<double, 3> &σs,
                   const std::array<double, 4> &ms2)
{
    return (*wr(3, 1, 0))(0, σs, ms2);
}
