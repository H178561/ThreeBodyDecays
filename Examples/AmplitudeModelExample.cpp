#include "../ThreeBodyAmplitudeModel.hh"
#include <iostream>
#include <iomanip>

int main()
{
    // Create masses and spins for the system
    ThreeBodyMasses ms = {0.938272046, 0.13957018, 0.493677, 2.28646};

    ThreeBodySpins spins = {1, 0, 0, 1}; // h0=1 refers to parent particle spin

    auto tbs = ThreeBodySystem(ms, spins);
    ThreeBodyParities parities('+', '-', '-', '+');

    // Create lineshape functions with different resonance parameters
    auto breitWignerL1520 = [](double sigma) -> complex
    {
        double mass = 1.518;
        double width = 0.015;
        return complex(1.0, 0.0) / complex(mass * mass - sigma, -mass * width);
    };

    auto breitWignerD1232 = [](double sigma) -> complex
    {
        double mass = 1.232;
        double width = 0.117;
        return complex(1.0, 0.0) / complex(mass * mass - sigma, -mass * width);
    };

    // Create different decay chains in different channels
    auto L1520 = createDecayChainCoupling(
        2,                                              // k-value
        breitWignerL1520,                               // Lineshape function
        "3/2-",                                         // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                     // ThreeBodySystem
        RecouplingType::NoRecoupling, {-1, 0}, false,   // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {0, 1}, false // ParityRecoupling for Hij
    );

    auto D1232 = createDecayChainCoupling(
        3,                                             // k-value
        breitWignerD1232,                              // Lineshape function
        "3/2+",                                        // jp (Spin-parity)
        ThreeBodySystem(ms, spins),                    // ThreeBodySystem
        RecouplingType::NoRecoupling, {-1, 0}, false,  // NoRecoupling for HRk
        RecouplingType::ParityRecoupling, {1, 0}, true // ParityRecoupling for Hij
    );

    // Create a model and add components (similar to Julia's collection)
    ThreeBodyAmplitudeModel model;
    model.add(L1520, "L1520", complex(1.0, 0.0));  // First chain with weight 1.0
    model.add(D1232, "D1232", complex(0.5, 0.25)); // Second chain with weight 0.5+0.25i

    std::cout << "Created amplitude model with " << model.size() << " decay chains" << std::endl;

    // Demonstrate using the model
    // Create a kinematic point (Mandelstam variables)
    ThreeBodyDecays tbd;
    MandelstamTuple σs = tbd.x2σs({0.3, 0.3}, ms, 1);

    std::cout << "Mandelstam variables: "
              << σs[0] << ", " << σs[1] << ", " << σs[2] << std::endl;

    // Calculate intensity at this kinematic point
    double intensity = model.intensity(σs, 1);
    std::cout << "Intensity at this point: " << intensity << std::endl;

    // Calculate amplitude for specific helicity configuration
    std::vector<int> two_λs = {0, 0, 0, 1}; // For 1/2 helicity we use 1 in doubled representation
    std::vector<int> refζs = {1, 2, 3, 1};  // Reference frames

    complex amp = model.amplitude(σs, two_λs, 1);
    std::cout << "Amplitude for specific helicities: "
              << std::real(amp) << " + " << std::imag(amp) << "i" << std::endl;

    return 0;
}
