#include <gtest/gtest.h>
#include "../OrientationAngles.hh"
#include <iostream>
#include <vector>
#include <tuple>

using namespace OrientationAngles;

// Test fixture for ThreeBodyDecays tests
class OrientationAnglesTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Initialize the four vectors from the notebook
        four_vectors["Dst"] = FourVector(0.0570074, -0.026685, 0.0479813, 2.0085299);
        four_vectors["D"] = FourVector(-0.108349, -0.056907, -0.295222, 1.8967276);
        four_vectors["Pi"] = FourVector(0.1034199, 0.0873037, 0.2482869, 0.3153471);
    }

    std::map<std::string, FourVector> four_vectors;
};

// Test masses of the four-vectors
TEST_F(OrientationAnglesTest, TestMasses)
{
    // Test individual masses
    EXPECT_NEAR(four_vectors["Dst"].mass(), 2.0069699, 1e-6);
    EXPECT_NEAR(four_vectors["D"].mass(), 1.87, 0.01);
    EXPECT_NEAR(four_vectors["Pi"].mass(), 0.14, 0.01);

    // Test sqrt(s) - sum of all vectors
    FourVector total = four_vectors["Dst"] + four_vectors["D"] + four_vectors["Pi"];
    EXPECT_NEAR(total.mass(), 4.22028, 0.01); // Updated from 4.20 to 4.22
}

// Test pure boost to rest frame
TEST_F(OrientationAnglesTest, TestPureBoost)
{
    // Boost to rest frame
    auto four_vectors_rf = pure_B(four_vectors);

    // Check if total momentum is zero in rest frame
    FourVector total_rf{0, 0, 0, 0};
    for (const auto &[name, p] : four_vectors_rf)
    {
        total_rf += p;
    }

    EXPECT_NEAR(total_rf.px(), 0.0, 1e-10);
    EXPECT_NEAR(total_rf.py(), 0.0, 1e-10);
    EXPECT_NEAR(total_rf.pz(), 0.0, 1e-10);

    // Check if masses are preserved
    EXPECT_NEAR(four_vectors_rf["Dst"].mass(), four_vectors["Dst"].mass(), 1e-10);
    EXPECT_NEAR(four_vectors_rf["D"].mass(), four_vectors["D"].mass(), 1e-10);
    EXPECT_NEAR(four_vectors_rf["Pi"].mass(), four_vectors["Pi"].mass(), 1e-10);

    // Verify actual values of boosted vectors match expected values
    EXPECT_NEAR(four_vectors_rf["Dst"].px(), 0.0322264, 1e-6);
    EXPECT_NEAR(four_vectors_rf["Dst"].py(), -0.0284512, 1e-6);
    EXPECT_NEAR(four_vectors_rf["Dst"].pz(), 0.0474835, 1e-6);
    EXPECT_NEAR(four_vectors_rf["Dst"].E(), 2.00799, 1e-5);

    EXPECT_NEAR(four_vectors_rf["D"].px(), -0.131764, 1e-6);
    EXPECT_NEAR(four_vectors_rf["D"].py(), -0.0585758, 1e-6);
    EXPECT_NEAR(four_vectors_rf["D"].pz(), -0.295692, 1e-6);
    EXPECT_NEAR(four_vectors_rf["D"].E(), 1.89833, 1e-5);

    EXPECT_NEAR(four_vectors_rf["Pi"].px(), 0.0995372, 1e-6);
    EXPECT_NEAR(four_vectors_rf["Pi"].py(), 0.087027, 1e-6);
    EXPECT_NEAR(four_vectors_rf["Pi"].pz(), 0.248209, 1e-6);
    EXPECT_NEAR(four_vectors_rf["Pi"].E(), 0.313957, 1e-5);
}

// Test helicity angles
TEST_F(OrientationAnglesTest, TestHelicityAngles)
{
    // Boost to rest frame
    auto four_vectors_rf = pure_B(four_vectors);

    // Define three different topologies - matching the Julia notebook format
    auto topology1 = std::make_tuple(
        std::make_tuple("Pi", "D"), "Dst");
    std::tuple<std::tuple<std::string, std::string>, std::string> topology2{
        std::make_tuple("D", "Dst"), "Pi"};
    std::tuple<std::tuple<std::string, std::string>, std::string> topology3{
        std::make_tuple("Dst", "Pi"), "D"};

    // Calculate angles for each topology
    auto angles1 = helicity_angles(four_vectors_rf, topology1);
    auto angles2 = helicity_angles(four_vectors_rf, topology2);
    auto angles3 = helicity_angles(four_vectors_rf, topology3);

    std::cout << angles1.size() << " angles for topology 1" << angles1[0].decay << " " << angles1[0].theta << " " << angles1[0].phi << std::endl;
    // Verify specific angle values (based on exact values from the notebook)
    // Topology (Pi,D),Dst: from a84fbe38-6d26-4738-a9bc-aa479d32f22a in notebook
    ASSERT_EQ(angles1.size(), 2); // Should have 2 decay entries
    // Check first decay (child1_θ, child1_φ)
    EXPECT_NEAR(angles1[0].theta, 2.40584, 0.00001);
    EXPECT_NEAR(angles1[0].phi, 2.21136, 0.00001);
    // Check second decay
    EXPECT_NEAR(angles1[1].theta, 0.453726, 0.00001);
    EXPECT_NEAR(angles1[1].phi, 0.418314, 0.00001);

    // Topology (D,Dst),Pi: from a84fbe38-6d26-4738-a9bc-aa479d32f22a in notebook
    ASSERT_EQ(angles2.size(), 2); // Should have 2 decay entries
    // Check first decay
    EXPECT_NEAR(angles2[0].theta, 2.36202, 0.00001);
    EXPECT_NEAR(angles2[0].phi, -2.43617, 0.00001);
    // Check second decay
    EXPECT_NEAR(angles2[1].theta, 2.65214, 0.00001);
    EXPECT_NEAR(angles2[1].phi, -2.42316, 0.00001);

    // Topology (Dst,Pi),D: from a84fbe38-6d26-4738-a9bc-aa479d32f22a in notebook
    ASSERT_EQ(angles3.size(), 2); // Should have 2 decay entries
    // Check first decay
    EXPECT_NEAR(angles3[0].theta, 0.238360, 0.00001);
    EXPECT_NEAR(angles3[0].phi, -1.18248, 0.00001);
    // Check second decay
    EXPECT_NEAR(angles3[1].theta, 2.97652, 0.00001);
    EXPECT_NEAR(angles3[1].phi, -1.69129, 0.00001);

    // Print results in a format similar to the Julia DataFrame output
    std::cout << "\nHelicity angles for all topologies (from C++ implementation):" << std::endl;

    std::cout << "\nTopology (Pi,D),Dst:" << std::endl;
    for (const auto &angle : angles1)
    {
        std::cout << "  " << angle.decay << ": theta=" << angle.theta
                  << ", phi=" << angle.phi << std::endl;
    }

    std::cout << "\nTopology (D,Dst),Pi:" << std::endl;
    for (const auto &angle : angles2)
    {
        std::cout << "  " << angle.decay << ": theta=" << angle.theta
                  << ", phi=" << angle.phi << std::endl;
    }

    std::cout << "\nTopology (Dst,Pi),D:" << std::endl;
    for (const auto &angle : angles3)
    {
        std::cout << "  " << angle.decay << ": theta=" << angle.theta
                  << ", phi=" << angle.phi << std::endl;
    }

    // Also verify the total masses of each subsystem - from the notebook
    double sqrt_s = 4.22028; // Total mass from notebook

    // Print masses of each subsystem as in the notebook
    std::cout << "\nTotal mass (sqrt_s): " << sqrt_s << std::endl;

    // Masses of individual particles should match the notebook
    std::cout << "Masses of individual particles:" << std::endl;
    std::cout << "  Dst: " << four_vectors_rf["Dst"].mass() << std::endl;
    std::cout << "  D: " << four_vectors_rf["D"].mass() << std::endl;
    std::cout << "  Pi: " << four_vectors_rf["Pi"].mass() << std::endl;
}

#ifndef RUNNING_COMBINED_TESTS
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
#endif
