#include "include/SpinParity.hh"
#include <stdexcept>

// Implementation of SpinParity constructor
SpinParity::SpinParity(const std::string &jp)
{
    if (jp.empty())
    {
        throw std::invalid_argument("jp string cannot be empty");
    }

    // Get the parity character (last character)
    p_ = jp.back();

    // Parse the spin value (rest of the string)
    std::string spin_str = jp.substr(0, jp.size() - 1);
    str = spin_str; // Store the original string for parsing

    // Check if spin string is empty
    if (spin_str.empty())
    {
        throw std::invalid_argument("Spin value cannot be empty");
    }

    // Check if spin is a fraction
    size_t slash_pos = spin_str.find('/');
    if (slash_pos != std::string::npos)
    {
        // Handle fraction
        std::string num_str = spin_str.substr(0, slash_pos);
        std::string denom_str = spin_str.substr(slash_pos + 1);

        if (num_str.empty() || denom_str.empty())
        {
            throw std::invalid_argument("Invalid fraction format");
        }

        try
        {
            int numerator = std::stoi(num_str);
            int denominator = std::stoi(denom_str);
            if (denominator == 0)
            {
                throw std::invalid_argument("Division by zero");
            }

            // For fractions, double the numerator
            two_j_ = (2 * numerator) / denominator;
        }
        catch (const std::exception &e)
        {
            throw std::invalid_argument("Invalid fraction: " + spin_str);
        }
    }
    else
    {
        // Integer spin - multiply by 2 for doubled representation
        try
        {
            int spin = std::stoi(spin_str);
            two_j_ = spin * 2; // This doubles the integer spin
        }
        catch (const std::exception &e)
        {
            throw std::invalid_argument("Invalid spin: " + spin_str);
        }
    }
}
