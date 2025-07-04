#ifndef THREEBODYAMPLITUDEMODEL_HH
#define THREEBODYAMPLITUDEMODEL_HH

#include "include/ThreeBodyDecays.hh"
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
    Tensor4Dcomp amplitude4d(const MandelstamTuple &σs,
                             const int &k_amp,
                             const std::vector<int> &refζs = {-1, -1, -1, -1}) const;

    // Calculate specific amplitude (for specific helicity values)
    complex amplitudes(const MandelstamTuple &σs,
                       const std::vector<int> &two_λs,
                       const int &k_amp,
                       const std::vector<int> &refζs = {-1, -1, -1, -1}) const;
    complex amplitude(const MandelstamTuple &σs,
                      const std::vector<int> &two_λs,
                      const int &k_amp,
                      const std::vector<int> &refζs = {-1, -1, -1, -1}) const;

    // Calculate total intensity (probability density)
    double intensity(const MandelstamTuple &σs, const int &k_amp, const std::vector<int> refζs = {-1, -1, -1, -1}) const;

    // Calculate individual component intensities
    std::vector<double> component_intensities(const MandelstamTuple &σs, const int &k_amp, const std::vector<int> refζs = {-1, -1, -1, -1}) const;

    // Calculate interference terms
    std::vector<std::vector<double>> interference_terms(const MandelstamTuple &σs, const int &k_amp, const std::vector<int> refζs = {-1, -1, -1, -1}) const;

private:
    // Store chains with their labels and coefficients
    std::vector<std::tuple<std::shared_ptr<DecayChain>, std::string, complex>> chains_;
};

#endif // THREEBODYAMPLITUDEMODEL_HH
