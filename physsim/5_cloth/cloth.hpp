#include <vislab/core/array_fwd.hpp>

#include <Eigen/Eigen>
#include <memory>
#include <set>

namespace vislab
{
    class TrimeshGeometry;
}

namespace physsim
{
    /**
     * @brief Helper class for the implementation of cloth simulations.
     */
    struct Cloth
    {
    public:
        /**
         * @brief Enumeration of spring topologies.
         */
        enum class ETopology
        {
            /**
             * @brief Structural springs, i.e., horizontal and vertical connections.
             */
            Structural,

            /**
             * @brief Structural and diagonal springs, i.e., diagonal connections.
             */
            Diagonal
        };

        /**
         * @brief Constructor for a new cloth.
         * @param resolution Resolution of the grid.
         * @param topology Connectivity of the springs.
         * @param origin Origin of the mesh grid.
         * @param axis1 Horizontal axis for the grid (x coordinate spacing).
         * @param axis2 Vertical axis for the grid (y coordinate spacing).
         */
        Cloth(const Eigen::Vector2i& resolution, const ETopology& topology, const Eigen::Vector3f& origin, const Eigen::Vector3f& axis1, const Eigen::Vector3f& axis2);

        /**
         * @brief Advance to next time step.
         * @param stepSize Numerical integration step size.
         */
        void advance(double stepSize);

        /**
         * @brief Resets the cloth back into its rest state.
         */
        void reset();

        /**
         * @brief Resolution of the cloth grid.
         * @return Grid resolution.
         */
        const Eigen::Vector2i& resolution() const;

        /**
         * @brief Pins a grid vertex.
         * @param gridIndex Index of vertex to hold fixed.
         */
        void pin(const Eigen::Vector2i& gridIndex);

        /**
         * @brief Unpins a grid vertex.
         * @param gridIndex Index of vertex to release.
         */
        void unpin(const Eigen::Vector2i& gridIndex);

        /**
         * @brief Mass of all bodies.
         */
        double mass;

        /**
         * @brief Stiffness coefficient of the springs.
         */
        double stiffness;

        /**
         * @brief Damping coefficient of the springs.
         */
        double damping;

        /**
         * @brief Gravitational acceleration.
         */
        Eigen::Vector3f gravity;

        /**
         * @brief Mesh that stores the geometry.
         */
        std::shared_ptr<vislab::TrimeshGeometry> mesh;

        /**
         * @brief Buffer that contains all spring velocities.
         */
        std::shared_ptr<vislab::Array3f> velocities;

        /**
         * @brief Buffer in which the acceleration is accumulated each frame.
         */
        std::shared_ptr<vislab::Array3f> accelerations;

    private:
        /**
         * @brief Helper function that applies a spring force from startPoint to endPoint.
         * @param startPoint Start of the spring.
         * @param endPoint End of the spring. This is where the force is added.
         * @param L Rest length of the spring.
         */
        void addSpringForce(const Eigen::Vector2i& startPoint, const Eigen::Vector2i& endPoint, const double& L);

        /**
         * @brief Linear indices of vertices that are fixed.
         */
        std::set<int> mFixed;

        /**
         * @brief Resolution of the grid.
         */
        Eigen::Vector2i mResolution;

        /**
         * @brief Connectivity of the springs.
         */
        ETopology mTopology;

        /**
         * @brief Origin of the mesh grid.
         */
        Eigen::Vector3f mOrigin;

        /**
         * @brief Horizontal axis for the grid (x coordinate spacing).
         */
        Eigen::Vector3f mAxis1;

        /**
         * @brief Vertical axis for the grid (y coordinate spacing).
         */
        Eigen::Vector3f mAxis2;
    };
}
