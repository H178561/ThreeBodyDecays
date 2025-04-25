#include "ThreeBodyAmplitudeModel.hh"
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <iostream>

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

Tensor4D ThreeBodyAmplitudeModel::amplitude4d(const MandelstamTuple &σs,
                                              const std::vector<int> &refζs) const
{
    if (chains_.empty())
    {
        throw std::runtime_error("No decay chains in the amplitude model");
    }

    ThreeBodyDecays tbd;

    // Get dimensions from the first chain
    auto &[first_chain, first_label, first_coef] = chains_[0];
    auto first_amp = tbd.amplitude4d(*first_chain, σs, refζs);

    // Skip computation if there's only one chain
    if (chains_.size() == 1)
    {
        // Always apply the coefficient (even if it's 1.0)
        Tensor4D result = first_amp;

        // Apply coefficient to all elements
        for (auto &dim1 : result)
        {
            for (auto &dim2 : dim1)
            {
                for (auto &dim3 : dim2)
                {
                    for (auto &val : dim3)
                    {
                        // Note: In this tensor we only store real values, so we
                        // only apply the magnitude of the coefficient here.
                        // Phase information is handled in the complex amplitude function.
                        val *= std::abs(first_coef);
                    }
                }
            }
        }
        return result;
    }

    // Initialize result tensor with zeros (same dimensions as first_amp)
    Tensor4D result(first_amp.size(),
                    std::vector<std::vector<std::vector<double>>>(
                        first_amp[0].size(),
                        std::vector<std::vector<double>>(
                            first_amp[0][0].size(),
                            std::vector<double>(first_amp[0][0][0].size(), 0.0))));

    // Sum all amplitudes with coefficient magnitudes
    for (const auto &[chain, label, coef] : chains_)
    {
        auto chain_amp = tbd.amplitude4d(*chain, σs, refζs);
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
                        result[i][j][k][l] += chain_amp[i][j][k][l] * coef_mag;
                    }
                }
            }
        }
    }

    return result;
}

complex ThreeBodyAmplitudeModel::amplitude(const MandelstamTuple &σs,
                                           const std::vector<int> &two_λs,
                                           const std::vector<int> &refζs) const
{
    if (chains_.empty())
    {
        return complex(0.0, 0.0);
    }

    ThreeBodyDecays tbd;
    complex result(0.0, 0.0);

    // Sum amplitudes with full complex coefficients
    for (const auto &[chain, label, coef] : chains_)
    {
        // Get amplitude for this chain with the specific helicities
        complex chain_amp = tbd.amplitude(*chain, σs, two_λs, refζs);

        // Apply the complex coefficient and add to total
        result += coef * chain_amp;

        // Debug output
        std::cout << "Chain: " << label
                  << ", Coef: " << coef
                  << ", Amp: " << chain_amp
                  << ", Contrib: " << (coef * chain_amp) << std::endl;
    }

    return result;
}

double ThreeBodyAmplitudeModel::intensity(const MandelstamTuple &σs) const
{
    if (chains_.empty())
    {
        return 0.0;
    }

    // Default reference frames
    std::vector<int> refζs = {1, 2, 3, 1};
    ThreeBodyDecays tbd;

    // Sum over all possible helicity configurations (summed_over_polarization in Julia)
    double total_intensity = 0.0;

    // Get two_js from the first chain
    const auto &[first_chain, _, __] = chains_[0];
    const auto &two_js = first_chain->tbs.two_js;

    // Generate all possible helicity combinations (like itr in Julia)
    std::vector<std::vector<int>> helicity_combos;
    generateHelicityCombinations(two_js, {}, helicity_combos);

    // Sum over all helicity combinations
    for (const auto &two_λs : helicity_combos)
    {
        complex amp = amplitude(σs, two_λs, refζs);
        total_intensity += std::norm(amp); // |amp|²
    }

    return total_intensity;
}

void ThreeBodyAmplitudeModel::generateHelicityCombinations(
    const std::array<int, 4> &two_js,
    std::vector<int> current,
    std::vector<std::vector<int>> &result) const
{

    if (current.size() == two_js.size())
    {
        result.push_back(current);
        return;
    }

    int idx = current.size();
    int two_j = two_js[idx];

    for (int two_λ = -two_j; two_λ <= two_j; two_λ += 2)
    {
        std::vector<int> next = current;
        next.push_back(two_λ);
        generateHelicityCombinations(two_js, next, result);
    }
}

std::vector<double> ThreeBodyAmplitudeModel::component_intensities(const MandelstamTuple &σs) const
{
    std::vector<double> result;
    result.reserve(chains_.size());

    // Default reference frames
    std::vector<int> refζs = {1, 2, 3, 1};
    ThreeBodyDecays tbd;

    // Calculate intensity for each component separately
    for (const auto &[chain, label, coef] : chains_)
    {
        const auto &two_js = chain->tbs.two_js;

        // Generate all possible helicity combinations
        std::vector<std::vector<int>> helicity_combos;
        generateHelicityCombinations(two_js, {}, helicity_combos);

        // Sum over all helicity combinations for this chain
        double chain_intensity = 0.0;
        for (const auto &two_λs : helicity_combos)
        {
            complex amp = tbd.amplitude(*chain, σs, two_λs, refζs);
            chain_intensity += std::norm(coef * amp); // |coef * amp|²
        }

        result.push_back(chain_intensity);
    }

    return result;
}

std::vector<std::vector<double>> ThreeBodyAmplitudeModel::interference_terms(const MandelstamTuple &σs) const
{
    size_t n = chains_.size();
    std::vector<std::vector<double>> result(n, std::vector<double>(n, 0.0));

    if (n <= 1)
    {
        return result;
    }

    // Default reference frames
    std::vector<int> refζs = {1, 2, 3, 1};
    ThreeBodyDecays tbd;

    // Calculate interference terms
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = i; j < n; ++j)
        {
            const auto &[chain_i, label_i, coef_i] = chains_[i];
            const auto &[chain_j, label_j, coef_j] = chains_[j];

            // Generate all possible helicity combinations (use the first chain's two_js)
            const auto &two_js = chain_i->tbs.two_js;
            std::vector<std::vector<int>> helicity_combos;
            generateHelicityCombinations(two_js, {}, helicity_combos);

            // Sum over all helicity combinations
            double interference_term = 0.0;
            for (const auto &two_λs : helicity_combos)
            {
                complex amp_i = tbd.amplitude(*chain_i, σs, two_λs, refζs);
                complex amp_j = tbd.amplitude(*chain_j, σs, two_λs, refζs);

                // When i == j, this is just the individual intensity
                // When i != j, this is the interference term
                if (i == j)
                {
                    interference_term += std::norm(coef_i * amp_i);
                }
                else
                {
                    // For interference, we need both Re(a*b) and Im(a*b)
                    complex prod = (coef_i * amp_i) * std::conj(coef_j * amp_j);
                    interference_term += 2.0 * std::real(prod);
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
