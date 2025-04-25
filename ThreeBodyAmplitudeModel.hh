#ifndef THREEBODYAMPLITUDEMODEL_HH
#define THREEBODYAMPLITUDEMODEL_HH

#include "ThreeBodyDecays.hh"
#include <string>
#include <vector>
#include <tuple>

// Class to manage multiple decay chains with coefficients
class ThreeBodyAmplitudeModel
{
public:
    // Add a decay chain with an optional name and coefficient
    void add(std::shared_ptr<DecayChain> chain,
             const std::string &label = "",
             complex coefficient = complex(1.0, 0.0));

    // Remove a chain by index or name
    void remove(size_t index);
    void remove(const std::string &label);

    // Clear all chains
    void clear();

    // Get number of chains
    size_t size() const { return chains_.size(); }

    // Calculate amplitude tensor (for all helicity states)
    Tensor4D amplitude4d(const MandelstamTuple &σs,
                         const std::vector<int> &refζs) const;

    // Calculate specific amplitude (for specific helicity values)
    complex amplitude(const MandelstamTuple &σs,
                      const std::vector<int> &two_λs,
                      const std::vector<int> &refζs) const;

    // Calculate total intensity (probability density)
    double intensity(const MandelstamTuple &σs) const;

    // Calculate individual component intensities
    std::vector<double> component_intensities(const MandelstamTuple &σs) const;

    // Calculate interference terms
    std::vector<std::vector<double>> interference_terms(const MandelstamTuple &σs) const;

private:
    // Store chains with their labels and coefficients
    std::vector<std::tuple<std::shared_ptr<DecayChain>, std::string, complex>> chains_;

    // Helper method to generate all helicity combinations
    void generateHelicityCombinations(
        const std::array<int, 4> &two_js,
        std::vector<int> current,
        std::vector<std::vector<int>> &result) const;
};

#endif // THREEBODYAMPLITUDEMODEL_HH
