#ifndef THREEBODYAMPLITUDEMODEL_HH
#define THREEBODYAMPLITUDEMODEL_HH

#include "ThreeBodyDecays.hh"
#include <string>
#include <vector>
#include <tuple>

/**
 * @brief Model class for combining multiple decay chains with coefficients
 *
 * This class manages a collection of decay chains, each with an associated
 * label and complex coefficient. It enables the calculation of combined
 * amplitudes and intensities from all chains, including interference effects.
 */
class ThreeBodyAmplitudeModel
{
public:
    /**
     * @brief Adds a decay chain to the amplitude model
     *
     * @param chain Shared pointer to the DecayChain object
     * @param label Optional name to identify the chain (auto-generated if empty)
     * @param coefficient Complex weight applied to this chain's amplitude
     */
    void add(std::shared_ptr<DecayChain> chain,
             const std::string &label = "",
             complex coefficient = complex(1.0, 0.0));

    /**
     * @brief Removes a decay chain by its index
     *
     * @param index Zero-based index of the chain to remove
     * @throws std::out_of_range if index is invalid
     */
    void remove(size_t index);

    /**
     * @brief Removes a decay chain by its label
     *
     * @param label The label of the chain to remove
     * @note Silently does nothing if label is not found
     */
    void remove(const std::string &label);

    /**
     * @brief Removes all decay chains from the model
     */
    void clear();

    /**
     * @brief Returns the number of decay chains in the model
     *
     * @return Number of chains
     */
    size_t size() const { return chains_.size(); }

    /**
     * @brief Calculates the 4D amplitude tensor for all helicity states
     *
     * Computes a combined amplitude tensor that includes contributions from
     * all decay chains, weighted by their coefficients.
     *
     * @param σs Mandelstam variables
     * @param k_amp Amplitude type index
     * @param refζs Reference angles for Wigner rotations
     * @return 4D tensor of complex amplitudes
     */
    Tensor4Dcomp amplitude4d(const MandelstamTuple &σs,
                             const int &k_amp,
                             const std::vector<int> &refζs = {-1, -1, -1, -1}) const;

    /**
     * @brief Calculates the combined amplitude for specific helicity values
     *
     * Legacy method, use amplitude() instead.
     *
     * @param σs Mandelstam variables
     * @param two_λs Doubled helicity values for all particles
     * @param k_amp Amplitude type index
     * @param refζs Reference angles for Wigner rotations
     * @return Complex amplitude
     */
    complex amplitudes(const MandelstamTuple &σs,
                       const std::vector<int> &two_λs,
                       const int &k_amp,
                       const std::vector<int> &refζs = {-1, -1, -1, -1}) const;

    /**
     * @brief Calculates the combined amplitude for specific helicity values
     *
     * @param σs Mandelstam variables
     * @param two_λs Doubled helicity values for all particles
     * @param k_amp Amplitude type index
     * @param refζs Reference angles for Wigner rotations
     * @return Complex amplitude
     */
    complex amplitude(const MandelstamTuple &σs,
                      const std::vector<int> &two_λs,
                      const int &k_amp,
                      const std::vector<int> &refζs = {-1, -1, -1, -1}) const;

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
    double intensity(const MandelstamTuple &σs, const int &k_amp, const std::vector<int> refζs = {-1, -1, -1, -1}) const;

    /**
     * @brief Calculates individual intensity contributions from each chain
     *
     * Returns the intensity that would result from each chain individually,
     * without interference effects.
     *
     * @param σs Mandelstam variables
     * @param k_amp Amplitude type index
     * @param refζs Reference angles for Wigner rotations
     * @return Vector of individual intensities, one for each chain
     */
    std::vector<double> component_intensities(const MandelstamTuple &σs, const int &k_amp, const std::vector<int> refζs = {-1, -1, -1, -1}) const;

    /**
     * @brief Calculates interference terms between all chain pairs
     *
     * Returns a matrix where element (i,j) represents the interference
     * contribution between chains i and j to the total intensity.
     *
     * @param σs Mandelstam variables
     * @param k_amp Amplitude type index
     * @param refζs Reference angles for Wigner rotations
     * @return Matrix of interference terms
     */
    std::vector<std::vector<double>> interference_terms(const MandelstamTuple &σs, const int &k_amp, const std::vector<int> refζs = {-1, -1, -1, -1}) const;

private:
    /**
     * Stores decay chains with their labels and coefficients.
     * Each tuple contains:
     * - shared_ptr to DecayChain
     * - string label
     * - complex coefficient
     */
    std::vector<std::tuple<std::shared_ptr<DecayChain>, std::string, complex>> chains_;
};

#endif // THREEBODYAMPLITUDEMODEL_HH
