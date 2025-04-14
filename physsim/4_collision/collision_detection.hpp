#pragma once

#include "contact.hpp"
#include "rigid_body.hpp"

#include <Eigen/Eigen>
#include <set>
#include <vector>

namespace physsim
{
    struct TransformedMesh;

    /**
     * @brief Enumeration of broad phase collision detection methods.
     */
    enum class EBroadPhaseMethod
    {
        None,
        AABB,
        SweepAndPrune
    };

    /**
     * @brief Class to detect and resolve collisions between rigid bodies.
     */
    class CollisionDetection
    {
    public:
        CollisionDetection(const std::vector<std::shared_ptr<RigidBody>>& objects);

        /**
         * @brief Performs a collision detection.
         * @param broadPhaseMethod Broad phase detection method.
         * @param eps Impulse epsilon, typically in [0,1]. A value of 1 is energy conserving. Everything below adds dampening.
         * @param stepSize Numerical integration step size.
         */
        void computeCollisionDetection(EBroadPhaseMethod broadPhaseMethod, double eps, double stepSize);

    private:
        /**
         * @brief Perform the broad phase collision detection.
         * @param broadPhaseMethod Broad phase method to apply.
         */
        void computeBroadPhase(EBroadPhaseMethod broadPhaseMethod);

        /**
         * @brief Perform narrow phase collision detection.
         */
        void computeNarrowPhase();

        /**
         * @brief Apply impulses for the found contacts.
         * @param eps Impulse epsilon.
         * @param stepSize Numerical integration step size.
         */
        void applyImpulse(double eps, double stepSize);

        /**
         * @brief Clear all data structures.
         */
        void clearDataStructures();

        /**
         * @brief All rigid body objects in the scene.
         */
        std::vector<std::shared_ptr<RigidBody>> mObjects;

        /**
         * @brief Result of broadphase, pairs of objects with possible collisions.
         */
        std::vector<std::pair<std::size_t, std::size_t>> mOverlappingBodys;

        /**
         * @brief Set of vertex indices that penetrate a face, used to avoid duplicates.
         */
        std::set<int> mPenetratingVertices;

        /**
         * @brief Set of pairs of vertex indices that represent a penetrating edge, used to avoid duplicates.
         */
        std::set<std::pair<int, int>> mPenetratingEdges;

        /**
         * @brief Computed contact points.
         */
        std::vector<Contact> mContacts;
    };
}
