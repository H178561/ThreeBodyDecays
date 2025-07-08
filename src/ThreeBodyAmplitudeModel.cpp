#include "include/ThreeBodyAmplitudeModel.hh"
#include "include/ThreeBodyDecays.hh"
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <array>

/**
 * @brief Adds a decay chain to the amplitude model
 *
 * If no label is provided, automatically generates one using the index.
 * Each chain is associated with a complex coefficient that weights its
 * contribution to the overall amplitude.
 *
 * @param chain Shared pointer to the DecayChain object
 * @param label Optional name to identify the chain
 * @param coefficient Complex weight applied to this chain's amplitude
 */
void ThreeBodyAmplitudeModel::add(std::shared_ptr<DecayChain> chain,
                                  const std::string &label,
                                  complex coefficient)
{
    // Generate a default label if none provided
    std::string actual_label = label.empty() ? "chain_" + std::to_string(chains_.size()) : label;

    chains_.emplace_back(chain, actual_label, coefficient);
}

/**
 * @brief Removes a decay chain by its index
 *
 * @param index Zero-based index of the chain to remove
 * @throws std::out_of_range if index is invalid
 */
void ThreeBodyAmplitudeModel::remove(size_t index)
{
    if (index >= chains_.size())
    {
        throw std::out_of_range("Chain index out of range");
    }
    chains_.erase(chains_.begin() + index);
}

/**
 * @brief Removes a decay chain by its label
 *
 * @param label The label of the chain to remove
 * @note Silently does nothing if label is not found
 */
void ThreeBodyAmplitudeModel::remove(const std::string &label)
{
    auto it = std::find_if(chains_.begin(), chains_.end(),
                           [&label](const auto &item)
                           {
                               return std::get<1>(item) == label;
                           });

    if (it != chains_.end())
    {
        chains_.erase(it);
    }
}

/**
 * @brief Removes all decay chains from the model
 */
void ThreeBodyAmplitudeModel::clear()
{
    chains_.clear();
}

/**
 * @brief Calculates the 4D amplitude tensor for all helicity states
 *
 * Computes a combined amplitude tensor that includes contributions from
 * all decay chains, weighted by their coefficients. This is the core method
 * for amplitude calculations in the model.
 *
 * @param σs Mandelstam variables
 * @param k_amp Amplitude type index
 * @param refζs Reference angles for Wigner rotations
 * @return 4D tensor of complex amplitudes
 * @throws std::runtime_error if the model contains no decay chains
 */
Tensor4Dcomp ThreeBodyAmplitudeModel::amplitude4d(const MandelstamTuple &σs,
                                                  const int &k_amp,
                                                  const std::vector<int> &refζs) const
{
    if (chains_.empty())
    {
        throw std::runtime_error("No decay chains in the amplitude model");
    }

    ThreeBodyDecays tbd;

    // Get dimensions from the first chain
    auto &[first_chain, first_label, first_coef] = chains_[0];
    auto first_amp = tbd.amplitude4dcomp(*first_chain, σs, k_amp, refζs);

    // Optimization: Skip computation if there's only one chain
    if (chains_.size() == 1)
    {
        const complex &coef = std::get<2>(chains_[0]);

        // If coefficient is 1.0, return as is
        if (coef == complex(1.0, 0.0))
        {
            return first_amp;
        }

        // Apply coefficient to all elements
        Tensor4Dcomp result = first_amp;
        double coef_mag = std::abs(coef);

        // Apply coefficient to all elements (similar to Julia broadcasting)
        for (auto &dim1 : result)
        {
            for (auto &dim2 : dim1)
            {
                for (auto &dim3 : dim2)
                {
                    for (auto &val : dim3)
                    {
                        val *= coef;
                    }
                }
            }
        }
        return result;
    }

    // Initialize result tensor with zeros
    Tensor4Dcomp result(first_amp.size(),
                        std::vector<std::vector<std::vector<complex>>>(
                            first_amp[0].size(),
                            std::vector<std::vector<complex>>(
                                first_amp[0][0].size(),
                                std::vector<complex>(first_amp[0][0][0].size(), 0.0))));

    // Sum all amplitudes with coefficients (similar to Julia sum)
    for (const auto &[chain, label, coef] : chains_)
    {
        auto chain_amp = tbd.amplitude4dcomp(*chain, σs, k_amp, refζs);
        double coef_mag = std::abs(coef);

        // Add to result
        for (size_t i = 0; i < result.size(); ++i)
        {
            for (size_t j = 0; j < result[i].size(); ++j)
            {
                for (size_t k = 0; k < result[i][j].size(); ++k)
                {
                    for (size_t l = 0; l < result[i][j][k].size(); ++l)
                    {
                        result[i][j][k][l] += chain_amp[i][j][k][l] * coef;
                    }
                }
            }
        }
    }

    return result;
}

/**
 * @brief Calculates the combined amplitude for specific helicity values
 *
 * @param σs Mandelstam variables
 * @param two_λs Doubled helicity values for all particles
 * @param k_amp Amplitude type index
 * @param refζs Reference angles for Wigner rotations
 * @return Complex amplitude
 */
complex ThreeBodyAmplitudeModel::amplitude(const MandelstamTuple &σs,
                                           const std::vector<int> &two_λs,
                                           const int &k_amp,
                                           const std::vector<int> &refζs) const
{
    Tensor4Dcomp F0 = amplitude4d(σs, k_amp, refζs);

    // Calculate indices from helicity values
    std::vector<int> indices(4);
    for (int i = 0; i < 4; ++i)
    {
        // Match Julia's div(_two_j + _two_λ, 2) + 1, but adjust for 0-based indexing
        indices[i] = (two_λs[i]);
    }

    // Check if indices are within bounds
    if (indices[0] >= 0 && indices[0] < (int)F0.size() &&
        indices[1] >= 0 && indices[1] < (int)F0[0].size() &&
        indices[2] >= 0 && indices[2] < (int)F0[0][0].size() &&
        indices[3] >= 0 && indices[3] < (int)F0[0][0][0].size())
    {

        return F0[indices[0]][indices[1]][indices[2]][indices[3]];
    }

    // Return 0 if indices are out of bounds
    return 0.0;
}

/**
 * @brief Calculates the combined amplitude for specific helicity values (legacy method)
 *
 * This method computes amplitudes directly from each decay chain and sums them,
 * rather than using the pre-computed 4D tensor. It's maintained for backward
 * compatibility.
 *
 * @param σs Mandelstam variables
 * @param two_λs Doubled helicity values for all particles
 * @param k_amp Amplitude type index
 * @param refζs Reference angles for Wigner rotations
 * @return Complex amplitude
 * @throws std::runtime_error if the model contains no decay chains
 */
complex ThreeBodyAmplitudeModel::amplitudes(const MandelstamTuple &σs,
                                            const std::vector<int> &two_λs,
                                            const int &k_amp,
                                            const std::vector<int> &refζs) const
{
    if (chains_.empty())
    {
        throw std::runtime_error("No decay chains in the amplitude model");
    }

    ThreeBodyDecays tbd;
    complex result(0.0, 0.0);

    // Sum amplitudes with coefficients (like Julia's sum)
    for (const auto &[chain, label, coef] : chains_)
    {
        result += coef * tbd.amplitude(*chain, σs, two_λs, k_amp, refζs);
    }

    return result;
}

/**
 * @brief Calculates the total decay intensity
 *
 * Computes the sum over all helicity combinations of the squared magnitude
 * of the combined amplitude, representing the probability density.
 *
 * @param σs Mandelstam variables
 * @param k_amp Amplitude type index
 * @param refζs Reference angles for Wigner rotations
 * @return Total intensity (probability density)
 */
double ThreeBodyAmplitudeModel::intensity(const MandelstamTuple &σs, const int &k_amp, const std::vector<int> refζs) const
{
    if (chains_.empty())
    {
        return 0.0;
    }

    // Get the 4D amplitude tensor
    auto amp = amplitude4d(σs, k_amp, refζs);

    bool printout = false;
    if(printout){
    // debug print out amp
    for (const auto &dim1 : amp)
    {
        for (const auto &dim2 : dim1)
        {
            for (const auto &dim3 : dim2)
            {
                for (const auto &val : dim3)
                {
                    std::cout << val << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }}

    // Sum squared amplitudes (similar to Julia's sum(abs2, ...))
    double total_intensity = 0.0;

    for (const auto &dim1 : amp)
    {
        for (const auto &dim2 : dim1)
        {
            for (const auto &dim3 : dim2)
            {
                for (const auto &val : dim3)
                {
                    double curintens = std::norm(val);
                    total_intensity += curintens;
                }
            }
        }
    }

    return total_intensity;
}

/**
 * @brief Calculates individual intensity contributions from each chain
 *
 * Returns the intensity that would result from each chain individually,
 * without interference effects. This is useful for understanding the
 * relative contributions of different resonances.
 *
 * @param σs Mandelstam variables
 * @param k_amp Amplitude type index
 * @param refζs Reference angles for Wigner rotations
 * @return Vector of individual intensities, one for each chain
 */
std::vector<double> ThreeBodyAmplitudeModel::component_intensities(const MandelstamTuple &σs, const int &k_amp, const std::vector<int> refζs) const
{
    std::vector<double> result;
    result.reserve(chains_.size());

    // Default reference frames
    ThreeBodyDecays tbd;

    // Calculate intensity for each component separately
    for (const auto &[chain, label, coef] : chains_)
    {
        auto chain_amp = tbd.amplitude4dcomp(*chain, σs, k_amp, refζs);

        // Sum squared amplitudes for this chain
        double chain_intensity = 0.0;
        for (const auto &dim1 : chain_amp)
        {
            for (const auto &dim2 : dim1)
            {
                for (const auto &dim3 : dim2)
                {
                    for (const auto &val : dim3)
                    {
                        complex ampval = val * coef;
                        double curintens = std::norm(ampval);
                        chain_intensity += curintens;
                    }
                }
            }
        }

        result.push_back(chain_intensity);
    }

    return result;
}

/**
 * @brief Calculates interference terms between all chain pairs
 *
 * Returns a matrix where element (i,j) represents the interference
 * contribution between chains i and j to the total intensity.
 * Diagonal elements (i,i) represent the intensity of chain i alone.
 *
 * @param σs Mandelstam variables
 * @param k_amp Amplitude type index
 * @param refζs Reference angles for Wigner rotations
 * @return Matrix of interference terms
 */
std::vector<std::vector<double>> ThreeBodyAmplitudeModel::interference_terms(const MandelstamTuple &σs, const int &k_amp, const std::vector<int> refζs) const
{
    size_t n = chains_.size();
    std::vector<std::vector<double>> result(n, std::vector<double>(n, 0.0));

    if (n <= 1)
    {
        return result;
    }

    // Default reference frames
    ThreeBodyDecays tbd;

    // Calculate all amplitudes once
    std::vector<Tensor4Dcomp> amplitudes;
    for (const auto &[chain, label, coef] : chains_)
    {
        amplitudes.push_back(tbd.amplitude4dcomp(*chain, σs, k_amp, refζs));
    }

    // Calculate interference terms
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = i; j < n; ++j)
        {
            double interference_term = 0.0;

            // Get the coefficients
            complex coef_i = std::get<2>(chains_[i]);
            complex coef_j = std::get<2>(chains_[j]);

            // Calculate the interference contribution
            for (size_t d1 = 0; d1 < amplitudes[i].size(); ++d1)
            {
                for (size_t d2 = 0; d2 < amplitudes[i][d1].size(); ++d2)
                {
                    for (size_t d3 = 0; d3 < amplitudes[i][d1][d2].size(); ++d3)
                    {
                        for (size_t d4 = 0; d4 < amplitudes[i][d1][d2][d3].size(); ++d4)
                        {
                            if (i == j)
                            {
                                // Diagonal terms
                                interference_term += amplitudes[i][d1][d2][d3][d4].real() *
                                                     amplitudes[j][d1][d2][d3][d4].real() *
                                                     std::norm(coef_i);
                            }
                            else
                            {
                                // Off-diagonal terms
                                interference_term += 2 * amplitudes[i][d1][d2][d3][d4].real() *
                                                     amplitudes[j][d1][d2][d3][d4].real() *
                                                     std::abs(coef_i) * std::abs(coef_j);
                            }
                        }
                    }
                }
            }

            result[i][j] = interference_term;
            if (i != j)
            {
                result[j][i] = interference_term; // Matrix is symmetric
            }
        }
    }

    return result;
}
