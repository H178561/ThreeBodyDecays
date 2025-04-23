#ifndef ORIENTATION_ANGLES_HH
#define ORIENTATION_ANGLES_HH

#include <array>
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <cmath>
#include <tuple>

namespace OrientationAngles
{

    // Simple 4-vector class
    class FourVector
    {
    public:
        FourVector(); // Add default constructor
        FourVector(double x, double y, double z, double e);

        double px() const { return m_p[0]; }
        double py() const { return m_p[1]; }
        double pz() const { return m_p[2]; }
        double E() const { return m_p[3]; }

        double mass() const;

        // Operations
        FourVector operator+(const FourVector &other) const;
        FourVector &operator+=(const FourVector &other);

        // Transformations
        FourVector boost(double bx, double by, double bz) const;
        FourVector boost_z(double beta) const;
        FourVector rotate_y(double angle) const;
        FourVector rotate_z(double angle) const;

    private:
        std::array<double, 4> m_p; // px, py, pz, E
    };

    // Helper struct for spherical coordinates
    struct SphericalCoordinates
    {
        double cosTheta;
        double phi;
    };

    SphericalCoordinates spherical_coordinates(const FourVector &p);
    double boost_gamma(const FourVector &p);

    // Transformation functions
    FourVector Rz(double phi, const FourVector &p);
    FourVector Ry(double theta, const FourVector &p);
    FourVector Bz(double gamma, const FourVector &p);

    // Pure boost to rest frame
    FourVector pure_B(const FourVector &p, const FourVector &p_ref);
    std::map<std::string, FourVector> pure_B(const std::map<std::string, FourVector> &system);

    // Decay tree representation
    class DecayNode
    {
    public:
        enum class NodeType
        {
            Particle,
            Decay
        };

        // Constructors for different types
        static DecayNode createParticle(const std::string &name);
        static DecayNode createDecay(const DecayNode &left, const DecayNode &right);
        static DecayNode createTopology(const std::tuple<std::tuple<std::string, std::string>, std::string> &topology);

        NodeType type() const { return m_type; }
        const std::string &name() const { return m_name; }
        const DecayNode *left() const { return m_left.get(); }
        const DecayNode *right() const { return m_right.get(); }

        // For debugging
        void printTopology(int indent = 0) const;

        DecayNode(NodeType type, const std::string &name);
        DecayNode(const DecayNode &other);            // Add copy constructor
        DecayNode &operator=(const DecayNode &other); // Add copy assignment operator

    private:
        NodeType m_type;
        std::string m_name;
        std::unique_ptr<DecayNode> m_left;
        std::unique_ptr<DecayNode> m_right;
    };

    class HelicityTransformation
    {
    public:
        // Apply transformation
        static DecayNode transform(const DecayNode &node, const std::map<std::string, FourVector> &momenta);
    };

    // Utility functions
    DecayNode add_indices_order(const DecayNode &node);
    DecayNode add_transform_through(const DecayNode &node, const std::map<std::string, FourVector> &momenta);

    // Compute decay angles
    struct DecayAngles
    {
        double theta;
        double phi;
        std::string decay;
    };

    std::vector<DecayAngles> decay_angles(const DecayNode &tree);
    std::vector<DecayAngles> helicity_angles(const std::map<std::string, FourVector> &four_vectors_rf,
                                             const std::tuple<std::tuple<std::string, std::string>, std::string> &topology);

    class OrientationAngles
    {
    public:
        // Calculate Euler angles from two sets of Cartesian coordinates
        static std::array<double, 3> calculateEulerAngles(
            const std::array<double, 3> &v1,
            const std::array<double, 3> &v2);

        // Calculate orientation angles from Mandelstam variables
        static std::array<double, 3> mandelstamToEulerAngles(
            const std::array<double, 3> &sigma,
            const std::array<double, 4> &masses_squared);

        // Normalize a vector
        static std::array<double, 3> normalize(const std::array<double, 3> &v);

        // Calculate the dot product of two vectors
        static double dotProduct(const std::array<double, 3> &v1, const std::array<double, 3> &v2);

        // Calculate the cross product of two vectors
        static std::array<double, 3> crossProduct(
            const std::array<double, 3> &v1,
            const std::array<double, 3> &v2);
    };

} // namespace ThreeBodyDecays

#endif // THREE_BODY_DECAYS_HH
