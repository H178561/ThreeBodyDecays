#ifndef WIGNER_ROTATION_HH
#define WIGNER_ROTATION_HH

#include <memory>
#include <array>
#include <tuple>

// Forward declarations
double Kallen(double x, double y, double z);

// Helper function to get indices i and j from k
std::pair<int, int> ij_from_k(int k);

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

// Function to compute i, j, k indices from Wigner rotation
std::tuple<int, int, int> ijk(const AbstractWignerRotation &wr);

// Helper function to determine if a system is sequential
bool issequential(int system_a, int reference_b);

// Helper function to create Wigner rotation
std::unique_ptr<AbstractWignerRotation> wr(int system_a, int reference_b,
                                           int particle_c = 0);

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
                      const std::array<double, 4> &ms2) const override;
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
                      const std::array<double, 4> &ms2) const override;
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
                      const std::array<double, 4> &ms2) const override;
    double cos_zeta(const std::array<double, 3> &sigma,
                    const std::array<double, 4> &ms2) const override;
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

#endif // WIGNER_ROTATION_HH
