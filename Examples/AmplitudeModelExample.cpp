#include "../ThreeBodyAmplitudeModel.hh"
#include <iostream>
#include <iomanip>

int main()
{
    // Create masses and spins for the system
    ThreeBodyMasses ms = {1.0, 2.0, 3.0, 15.0};
    ThreeBodySpins spins = {1, 0, 0, 1}; // h0=1 refers to parent particle spin

    auto tbs = ThreeBodySystem(ms, spins);
    ThreeBodyParities parities('+', '+', '+', '+');

    // Create lineshape functions with different resonance parameters
    auto breitWigner1 = [](double sigma) -> complex
    {
        double mass = 5.0;
        double width = 0.2;
        return complex(1.0, 0.0) / complex(mass * mass - sigma, -mass * width);
    };

    auto breitWigner2 = [](double sigma) -> complex
    {
        double mass = 6.5;
        double width = 0.4;
        return complex(1.0, 0.0) / complex(mass * mass - sigma, -mass * width);
    };

    // Create different decay chains in different channels
    auto dc1 = createDecayChainCoupling(
        1,            // k-value (channel 1)
        breitWigner1, // First resonance
        "2+",         // Spin-parity 2+
        //parities,     // Parities
        tbs,          // ThreeBodySystem
        RecouplingType::NoRecoupling, {0, 0});

    auto dc2 = createDecayChainCoupling(
        2,            // k-value (channel 2)
        breitWigner2, // Second resonance
        "1-",         // Spin-parity 1-
        //parities,     // Parities
        tbs,          // ThreeBodySystem
        RecouplingType::ParityRecoupling, {2, 0}, true);

    // Create a model and add components (similar to Julia's collection)
    ThreeBodyAmplitudeModel model;
    model.add(dc1, "1", complex(1.0, 0.0));  // First chain with weight 1.0
    model.add(dc2, "2", complex(0.5, 0.25)); // Second chain with weight 0.5+0.25i

    std::cout << "Created amplitude model with " << model.size() << " decay chains" << std::endl;

    // Demonstrate using the model
    // Create a kinematic point (Mandelstam variables)
    ThreeBodyDecays tbd;
    MandelstamTuple σs = tbd.x2σs({0.3, 0.3}, ms, 1);

    std::cout << "Mandelstam variables: "
              << σs[0] << ", " << σs[1] << ", " << σs[2] << std::endl;

    // Calculate intensity at this kinematic point
    double intensity = model.intensity(σs, 0);
    std::cout << "Intensity at this point: " << intensity << std::endl;

    // Calculate amplitude for specific helicity configuration
    std::vector<int> two_λs = {1, 1, 1, 1}; // For 1/2 helicity we use 1 in doubled representation
    std::vector<int> refζs = {1, 2, 3, 1};  // Reference frames

    complex amp = model.amplitude(σs, two_λs, 0, refζs);
    std::cout << "Amplitude for specific helicities: "
              << std::real(amp) << " + " << std::imag(amp) << "i" << std::endl;

    return 0;
}
