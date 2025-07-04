#include "include/ThreeBodyAmplitudeModel.hh"
#include "include/ThreeBodyDecays.hh"
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <array>

void ThreeBodyAmplitudeModel::add(std::shared_ptr<DecayChain> chain,
                                  const std::string &label,
                                  complex coefficient)
{
    // Generate a default label if none provided
    std::string actual_label = label.empty() ? "chain_" + std::to_string(chains_.size()) : label;

    chains_.emplace_back(chain, actual_label, coefficient);
}

void ThreeBodyAmplitudeModel::remove(size_t index)
{
    if (index >= chains_.size())
    {
        throw std::out_of_range("Chain index out of range");
    }
    chains_.erase(chains_.begin() + index);
}

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

void ThreeBodyAmplitudeModel::clear()
{
    chains_.clear();
}

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

    // Skip computation if there's only one chain
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
        // std::cout << "Applying coefficient: " << coef << " on first amplitude tensor" << std::endl;

        // Add to result
        for (size_t i = 0; i < result.size(); ++i)
        {
            for (size_t j = 0; j < result[i].size(); ++j)
            {
                for (size_t k = 0; k < result[i][j].size(); ++k)
                {
                    for (size_t l = 0; l < result[i][j][k].size(); ++l)
                    {
                        // std::cout << i << " " << j << " " << k << " " << l << " " << chain_amp[i][j][k][l] << " " << coef_mag << " " << " endres " << chain_amp[i][j][k][l] * coef_mag << std::endl;
                        result[i][j][k][l] += chain_amp[i][j][k][l] * coef;
                        // result[i][j][k][l] += complex(chain_amp[i][j][k][l].real() * coef.real(),chain_amp[i][j][k][l].imag() * coef.imag());
                        // std::cout << "Adding " << label << " amplitude: " << chain_amp[i][j][k][l] << " with coefficient: " << coef << " to result at index [" << i << "][" << j << "][" << k << "][" << l << "]" << " = " << chain_amp[i][j][k][l] * coef << std::endl;
                    }
                }
            }
        }
    }

    return result;
}

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
                    //std::cout << curintens << std::endl;
                }
            }
        }
    }

    return total_intensity;
}

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
        // std::cout << "Calculating intensity for chain: " << label << " with coefficient: " << coef << std::endl;

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
                        // chain_intensity += (val.real()*val.real() * coef.real()  * coef.real()  +
                        //                    (val.imag() * val.imag() * coef.imag() * coef.imag();

                        complex ampval = complex(val.real() * coef.real(), val.imag() * coef.imag());
                        // ampval = val;
                        ampval = val * coef;
                        double curintens = std::norm(ampval);

                        chain_intensity += curintens;

                        // std::cout << "Adding intensity for chain: " << label << " with coefficient: " << coef
                        //           << " and value: " << val << " with weighted amp" << val*coef << " gives current intensity: " << curintens
                        //           << " total: " << chain_intensity << std::endl;
                    }
                }
            }
        }

        /*
        double realintensity;
                for ( int i = 0; i < chain_amp.size(); ++i ) {
                    for ( int j = 0; j < chain_amp[0].size(); ++j ) {
                        for ( int k = 0; k < chain_amp[0][0].size(); ++k ) {
                            for ( int z = 0; z < chain_amp[0][0][0].size(); ++z ) {
                                //realintensity += A_chain_values[i][j][k][z].real() * A_chain_values[i][j][k][z].real() * weight.real() * weight.real() +
                                //                + A_chain_values[i][j][k][z].imag() * A_chain_values[i][j][k][z].imag() * weight.imag() * weight.imag();
                                double amp_intensity = chain_amp[i][j][k][z].real() * chain_amp[i][j][k][z].real() * coef.real() * coef.real() +
                                                chain_amp[i][j][k][z].imag() * chain_amp[i][j][k][z].imag() * coef.imag() * coef.imag();
                                //if(deubevt) std::cout << "Amplitude[" << i << "][" << j << "][" << k << "][" << z << "] = " << realintensity << std::endl;

                                realintensity += amp_intensity;
                                std::cout << label << " Amplitudemodel[" << i << "][" << j << "][" << k << "][" << z << "] = "
                                          << chain_amp[i][j][k][z] << coef << amp_intensity << std::endl;
                            }
                        }
                    }
                }*/

        // result.push_back(chain_intensity);
        result.push_back(chain_intensity);
    }

    return result;
}

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
                            // For simplicity, we just use the magnitude
                            // In a full implementation, you'd use the phase too
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
