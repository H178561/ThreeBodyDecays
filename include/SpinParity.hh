#ifndef SPIN_PARITY_HH
#define SPIN_PARITY_HH

#include <string>

class SpinParity
{
public:
    SpinParity(const std::string &jp);

    // Getters
    int get_two_j() const {
        return two_j_; // Now properly returns the doubled value
    }
    int get_j() const {
        return two_j_ / 2; // Returns the original j value
    }
    char get_p() const { return p_; }

private:
    int two_j_; // 2 * j
    char p_;    // '+' or '-'
    std::string str; // Store the original string for parsing
};

#endif // SPIN_PARITY_HH
