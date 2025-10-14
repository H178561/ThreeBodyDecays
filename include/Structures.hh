#ifndef STRUCTURES_HH
#define STRUCTURES_HH


#include <complex>
#include <map>
#include <memory> // for std::unique_ptr
#include <string>
#include <tuple>
#include <vector>
#include <array>
#include <functional> // for std::function

using complex = std::complex<double>;
using MandelstamTuple = std::array<double, 3>;    // (s12, s23, s31)
using MassTuple = std::array<double, 4>;          // (m1, m2, m3, m0)
using SpinTuple = std::array<int, 4>;             // (j1, j2, j3, J)
using RecouplingLS =
    std::function<complex(const std::array<int, 2> &, const std::array<int, 3> &)>;
using LineshapeFunction = std::function<complex(double)>;
using HelicityFunction =
    std::function<complex(const std::array<int, 2> &, const std::array<int, 3> &)>;
using ThreeBodyMasses = std::array<double, 4>;
using ThreeBodySpins = std::array<int, 4>;
using Matrix2D = std::vector<std::vector<double>>;
using Tensor4Dcomp = std::vector<std::vector<std::vector<std::vector<complex>>>>;
using Tensor4D = std::vector<std::vector<std::vector<std::vector<double>>>>;


// Ersetzen Sie Zeile 20:
// using MassTuple = std::array<double, 4>;          // (m1, m2, m3, M)

struct MassTupleexpl {
    double m1, m2, m3;  // Erste drei Werte (positional)
    double m0;          // Muss IMMER als .m0 = value gesetzt werden

    // KEIN Default-Konstruktor - verhindert MassTuple{}
    MassTupleexpl() = delete;

    // KEIN Konstruktor mit 4 Parametern - verhindert MassTuple(a,b,c,d)
    MassTupleexpl(double, double, double, double) = delete;

    // NUR dieser Konstruktor ist erlaubt - für designated initializers
    MassTupleexpl(double m1_, double m2_, double m3_)
        : m1(m1_), m2(m2_), m3(m3_), m0(0.0) {}  // m0 wird später per designated initializer gesetzt

    // Array-Zugriff für Rückwärtskompatibilität
    double operator[](int i) const {
        switch(i) {
            case 0: return m1;
            case 1: return m2;
            case 2: return m3;
            case 3: return m0;
            default: throw std::out_of_range("Invalid index for MassTuple");
        }
    }

    double& operator[](int i) {
        switch(i) {
            case 0: return m1;
            case 1: return m2;
            case 2: return m3;
            case 3: return m0;
            default: throw std::out_of_range("Invalid index for MassTuple");
        }
    }

    // Konvertierung zu std::array
    std::array<double, 4> to_array() const {
        return {m1, m2, m3, m0};
    }

    // STL-Interface
    auto begin() const { return to_array().begin(); }
    auto end() const { return to_array().end(); }
    size_t size() const { return 4; }
};



// Ersetzen Sie diese Zeile:
// using MandelstamTuple = std::array<double, 3>;    // (s12, s23, s31)

struct MandelstamTupleexpl {
    double s12, s23, s31;  // ALLE müssen explizit benannt werden

    // ALLE Konstruktoren sind gelöscht - erzwingt designated initializers
    MandelstamTupleexpl() = delete;
    MandelstamTupleexpl(double, double, double) = delete;
    MandelstamTupleexpl(double, double) = delete;
    MandelstamTupleexpl(double) = delete;

    // Array-Zugriff für Rückwärtskompatibilität
    double operator[](int i) const {
        switch(i) {
            case 0: return s12;
            case 1: return s23;
            case 2: return s31;
            default: throw std::out_of_range("Invalid index for MandelstamTuple");
        }
    }

    double& operator[](int i) {
        switch(i) {
            case 0: return s12;
            case 1: return s23;
            case 2: return s31;
            default: throw std::out_of_range("Invalid index for MandelstamTuple");
        }
    }

    // Konvertierung zu std::array
    std::array<double, 3> to_array() const {
        return {s12, s23, s31};
    }

    // STL-Interface
    auto begin() const { return to_array().begin(); }
    auto end() const { return to_array().end(); }
    size_t size() const { return 3; }
};




#endif // STRUCTURES_HH
