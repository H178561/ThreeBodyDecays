#ifndef WIGNER_FUNCTIONS_HH
#define WIGNER_FUNCTIONS_HH

#include <vector>
#include <array>

// Define the Matrix2D type if not already defined elsewhere
using Matrix2D = std::vector<std::vector<double>>;

// Wigner d-hat function with doubled arguments
double wignerd_hat_doublearg(int two_j, int two_m1, int two_m2, double z);

// Wigner d-function with doubled arguments
double wignerd_doublearg(int two_j, int two_m1, int two_m2, double z);

// Function to compute Wigner-d with sign adjustment for a single element
double wignerd_doublearg_sign_element(int two_j, int two_λ1, int two_λ2, double cosθ, bool ispositive);

// Function to compute Wigner-d matrix with sign adjustment
Matrix2D wignerd_doublearg_sign(int two_j, double cosθ, bool ispositive);

#endif // WIGNER_FUNCTIONS_HH
