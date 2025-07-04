#ifndef CLEBSCH_GORDON_HH
#define CLEBSCH_GORDON_HH

#include <complex>
#include <vector>

using complex = std::complex<double>;

// Define constant
extern const int MAX_FACTcg;

// Declarations only - no implementationsZZZZ
double getLogFactorialcg(int n);
std::vector<double> initializeLogFact2();
double f_logfact2(int two_n);
double clebschgordan_doublearg(int two_j1, int two_m1, int two_j2, int two_m2, int two_j, int two_m);
double clebschgordan(int j1, int m1, int j2, int m2, int j, int m);
double CG(int j1, int m1, int j2, int m2, int j, int m);
double CG_doublearg(int two_j1, int two_m1, int two_j2, int two_m2, int two_j, int two_m);
complex CG_0(int two_j1, int two_m1, int two_j2, int two_m2, int two_j, int two_m);

#endif // CLEBSCH_GORDON_HH
