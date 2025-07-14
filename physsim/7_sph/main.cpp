#include <init_vislab.hpp>

#include "physsim_window.hpp"
#include "simulation.hpp"

#include <imgui.h>
#include <vislab/geometry/points.hpp>
#include <vislab/graphics/actor.hpp>
#include <vislab/graphics/colormap_texture.hpp>
#include <vislab/graphics/const_texture.hpp>
#include <vislab/graphics/diffuse_bsdf.hpp>
#include <vislab/graphics/perspective_camera.hpp>
#include <vislab/graphics/point_light.hpp>
#include <vislab/graphics/rectangle_geometry.hpp>
#include <vislab/graphics/scene.hpp>
#include <vislab/graphics/sphere_geometry.hpp>
#include <vislab/graphics/transform.hpp>

#include <vislab/core/array.hpp>

#include "nearest_neighbors.hpp"

#include <random>

using namespace vislab;

namespace physsim
{
    /*
     * \brief Cubic spline kernel from https://github.com/InteractiveComputerGraphics/SPlisHSPlasH/blob/master/SPlisHSPlasH/SPHKernels.h
     */
    class CubicKernel
    {
    protected:
        float m_radius;
        float m_k;
        float m_l;

    public:
        float getRadius() { return m_radius; }
        void setRadius(float val)
        {
            m_radius       = val;
            const float pi = static_cast<float>(EIGEN_PI);

            const float h3 = m_radius * m_radius * m_radius;
            m_k            = 8.f / (pi * h3);
            m_l            = 48.f / (pi * h3);
        }

    public:
        float W(const float r)
        {
            float res     = 0.0;
            const float q = r / m_radius;
            if (q <= 1.0)
            {
                if (q <= 0.5)
                {
                    const float q2 = q * q;
                    const float q3 = q2 * q;
                    res            = m_k * (6.f * q3 - 6.f * q2 + 1.f);
                }
                else
                {
                    res = m_k * (2.f * pow(1.f - q, 3.f));
                }
            }
            return res;
        }

        float W(const Eigen::Vector3f& r)
        {
            return W(r.norm());
        }

        Eigen::Vector3f gradW(const Eigen::Vector3f& r)
        {
            Eigen::Vector3f res;
            const float rl = r.norm();
            const float q  = rl / m_radius;
            if ((rl > 1.0e-9) && (q <= 1.0))
            {
                Eigen::Vector3f gradq = r / rl;
                gradq /= m_radius;
                if (q <= 0.5f)
                {
                    res = m_l * q * (3.f * q - 2.f) * gradq;
                }
                else
                {
                    const float factor = 1.f - q;
                    res                = m_l * (-factor * factor) * gradq;
                }
            }
            else
                res.setZero();

            return res;
        }
    };

    /**
     * @brief Implementation of weakly-coupled smooth particle hydrodynamics with Akinci boundary conditions.
     */
    class SmoothedParticleHydrodynamicsSimulation : public Simulation
    {
    public:
        /**
         * @brief Initializes the scene.
         */
        void init() override
        {
            // set extent of simulation domain
            mDomain = Eigen::AlignedBox3f(Eigen::Vector3f(-10, -5, 0), Eigen::Vector3f(10, 5, 10));

            // allocate buffers needed by the simulation
            mPositions         = std::make_shared<Array3f>("positions");
            mVelocities        = std::make_shared<Array3f>("velocities");
            mAccelerations     = std::make_shared<Array3f>("accelerations");
            mDensities         = std::make_shared<Array1f>("densities");
            mMasses            = std::make_shared<Array1f>("masses");
            mPressures         = std::make_shared<Array1f>("pressure");
            mBoundaryParticles = std::make_shared<Array3f>("boundaryParticles");
            mBoundaryMasses    = std::make_shared<Array1f>("boundaryMasses");

            // set simulation parameters
            mStepSize      = 1E-2;
            mGravity       = Eigen::Vector3f(0, 0, -9.8065f);
            mStiffness     = 50000;
            mViscosity     = 5;
            mExponent      = 7;
            mRho0          = 1000;
            mSupportRadius = 1.f;

            // create ground plane
            auto rectActor = std::make_shared<Actor>("ground");
            rectActor->components.add(std::make_shared<RectangleGeometry>());
            rectActor->components.add(std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(1, 1, 1))));
            rectActor->components.add(std::make_shared<Transform>(Eigen::Vector3d::Zero(), Eigen::Quaterniond::Identity(), Eigen::Vector3d(10., 5., 1.)));

            // create a point light
            auto lightActor = std::make_shared<Actor>("light");
            lightActor->components.add(std::make_shared<PointLight>(Spectrum(500.0, 500.0, 500.0)));
            lightActor->components.add(std::make_shared<Transform>(Eigen::Vector3d(5, 20, 10)));

            // create a camera
            auto cameraActor = std::make_shared<Actor>("camera");
            auto camera      = std::make_shared<PerspectiveCamera>();
            camera->setLookAt(Eigen::Vector3d(0, 0, 3));
            camera->setPosition(Eigen::Vector3d(2, 25, 6));
            camera->setUp(Eigen::Vector3d(0, 0, 1));
            camera->setNear(0.01);
            camera->setFar(100);
            cameraActor->components.add(camera);

            mSpheres                  = std::make_shared<Actor>("spheres");
            auto sphereGeometry       = std::make_shared<SphereGeometry>(0.25 * mSupportRadius);
            sphereGeometry->positions = mPositions;
            sphereGeometry->data      = mPressures;
            mSpheres->components.add(sphereGeometry);
            auto colormapTexture                       = std::make_shared<ColormapTexture>();
            colormapTexture->transferFunction.minValue = 0;
            colormapTexture->transferFunction.maxValue = 100000;
            colormapTexture->transferFunction.values.clear();
            colormapTexture->transferFunction.values.insert(std::make_pair(0., Eigen::Vector4d(103 / 255., 169 / 255., 207 / 255., 1)));
            colormapTexture->transferFunction.values.insert(std::make_pair(0.25, Eigen::Vector4d(209 / 255., 229 / 255., 240 / 255., 1)));
            colormapTexture->transferFunction.values.insert(std::make_pair(0.5, Eigen::Vector4d(247 / 255., 247 / 255., 247 / 255., 1)));
            colormapTexture->transferFunction.values.insert(std::make_pair(0.75, Eigen::Vector4d(253 / 255., 219 / 255., 199 / 255., 1)));
            colormapTexture->transferFunction.values.insert(std::make_pair(1., Eigen::Vector4d(239 / 255., 138 / 255., 98 / 255., 1)));
            mSpheres->components.add(std::make_shared<DiffuseBSDF>(colormapTexture));
            mSpheres->components.add(std::make_shared<Transform>(Eigen::Vector3d(0, 0, 0)));

            // add elements to scene
            scene->actors.push_back(rectActor);
            scene->actors.push_back(lightActor);
            scene->actors.push_back(cameraActor);
            scene->actors.push_back(mSpheres);
        }

        /**
         * @brief Restarts the simulation.
         */
        void restart() override
        {
            // fill a portion of the domain with particles (uniform distribution with stratified sampling)
            Eigen::AlignedBox3f domainToFill(Eigen::Vector3f(-8, -3, 2), Eigen::Vector3f(-2, 3, 9));
            Eigen::Vector3i initialRes = ((domainToFill.max() - domainToFill.min()) / mSupportRadius * 2).cast<int>();
            mPositions->setSize(initialRes.prod());
            mVelocities->setSize(initialRes.prod());
            mAccelerations->setSize(initialRes.prod());
            mDensities->setSize(initialRes.prod());
            mPressures->setSize(initialRes.prod());
            std::default_random_engine rng;
            std::uniform_real_distribution<float> rnd(0, 1);
            for (int iz = 0; iz < initialRes.z(); ++iz)
                for (int iy = 0; iy < initialRes.y(); ++iy)
                    for (int ix = 0; ix < initialRes.x(); ++ix)
                    {
                        Eigen::Index linearIndex = (iz * initialRes.y() + iy) * initialRes.x() + ix;
                        Eigen::Vector3f relpos   = Eigen::Vector3f(
                            (ix + rnd(rng)) / (float)initialRes.x(),
                            (iy + rnd(rng)) / (float)initialRes.y(),
                            (iz + rnd(rng)) / (float)initialRes.z());
                        Eigen::Vector3f pos = domainToFill.min() + relpos.cwiseProduct((domainToFill.max() - domainToFill.min()));
                        mPositions->setValue(linearIndex, pos);
                        mVelocities->setValue(linearIndex, Eigen::Vector3f::Zero());
                    }

            // place boundary particles on the walls
            Eigen::Vector3i boundaryRes = ((mDomain.max() - mDomain.min()) / mSupportRadius * 4).cast<int>();
            mBoundaryParticles->setSize(2 * boundaryRes.x() * boundaryRes.y() + 2 * boundaryRes.y() * boundaryRes.z() + 2 * boundaryRes.x() * boundaryRes.z());
            Eigen::Index linearIndex = 0;
            for (int iy = 0; iy < boundaryRes.y(); ++iy)
                for (int ix = 0; ix < boundaryRes.x(); ++ix)
                {
                    Eigen::Vector3f relpos = Eigen::Vector3f(
                        ix / (boundaryRes.x() - 1.f),
                        iy / (boundaryRes.y() - 1.f),
                        0);
                    Eigen::Vector3f pos = mDomain.min() + relpos.cwiseProduct((mDomain.max() - mDomain.min()));
                    mBoundaryParticles->setValue(linearIndex++, pos);

                    relpos = Eigen::Vector3f(
                        ix / (boundaryRes.x() - 1.f),
                        iy / (boundaryRes.y() - 1.f),
                        1);
                    pos = mDomain.min() + relpos.cwiseProduct((mDomain.max() - mDomain.min()));
                    mBoundaryParticles->setValue(linearIndex++, pos);
                }
            for (int iz = 0; iz < boundaryRes.z(); ++iz)
                for (int iy = 0; iy < boundaryRes.y(); ++iy)
                {
                    Eigen::Vector3f relpos = Eigen::Vector3f(
                        0,
                        iy / (boundaryRes.y() - 1.f),
                        iz / (boundaryRes.z() - 1.f));
                    Eigen::Vector3f pos = mDomain.min() + relpos.cwiseProduct((mDomain.max() - mDomain.min()));
                    mBoundaryParticles->setValue(linearIndex++, pos);

                    relpos = Eigen::Vector3f(
                        1,
                        iy / (boundaryRes.y() - 1.f),
                        iz / (boundaryRes.z() - 1.f));
                    pos = mDomain.min() + relpos.cwiseProduct((mDomain.max() - mDomain.min()));
                    mBoundaryParticles->setValue(linearIndex++, pos);
                }
            for (int iz = 0; iz < boundaryRes.z(); ++iz)
                for (int ix = 0; ix < boundaryRes.x(); ++ix)
                {
                    Eigen::Vector3f relpos = Eigen::Vector3f(
                        ix / (boundaryRes.x() - 1.f),
                        0,
                        iz / (boundaryRes.z() - 1.f));
                    Eigen::Vector3f pos = mDomain.min() + relpos.cwiseProduct((mDomain.max() - mDomain.min()));
                    mBoundaryParticles->setValue(linearIndex++, pos);

                    relpos = Eigen::Vector3f(
                        ix / (boundaryRes.x() - 1.f),
                        1,
                        iz / (boundaryRes.z() - 1.f));
                    pos = mDomain.min() + relpos.cwiseProduct((mDomain.max() - mDomain.min()));
                    mBoundaryParticles->setValue(linearIndex++, pos);
                }

            // create a cubic spline kernel
            CubicKernel W;
            W.setRadius(mSupportRadius);

            // build nearest neighbor data structure for the boundary particles
            mBoundaryKNN = std::make_shared<NearestNeighbors3f>();
            mBoundaryKNN->setPoints(mBoundaryParticles);

            // compute the mass of boundary particles
            mBoundaryMasses->setSize(mBoundaryParticles->getSize());
            for (Eigen::Index i = 0; i < mBoundaryParticles->getSize(); ++i)
            {
                Eigen::Vector3f xk = mBoundaryParticles->getValue(i);
                float volk         = 0;
                NearestNeighbors3f::RadiusResult nnl;
                if (mBoundaryKNN->closestRadius(xk, mSupportRadius, nnl) != 0)
                    for (auto& nl : nnl)
                    {
                        Eigen::Vector3f xl  = mBoundaryParticles->getValue(nl.first);
                        Eigen::Vector3f xkl = xk - xl;
                        volk += W.W(xkl);
                    }
                float massk = mRho0 / volk;
                mBoundaryMasses->setValue(i, massk);
            }

            // compute the mass of domain particles
#if 1
            mMasses->setSize(mPositions->getSize());
            auto KNN = std::make_shared<NearestNeighbors3f>();
            KNN->setPoints(mPositions);
            for (Eigen::Index i = 0; i < mPositions->getSize(); ++i)
            {
                Eigen::Vector3f xi = mPositions->getValue(i);
                float voli         = 0;
                NearestNeighbors3f::RadiusResult nn;
                if (KNN->closestRadius(xi, mSupportRadius, nn) != 0)
                    for (auto& nl : nn)
                    {
                        Eigen::Vector3f xk  = mPositions->getValue(nl.first);
                        Eigen::Vector3f xik = xi - xk;
                        voli += W.W(xik);
                    }
                NearestNeighbors3f::RadiusResult nnl;
                if (mBoundaryKNN->closestRadius(xi, mSupportRadius, nnl) != 0)
                    for (auto& nl : nnl)
                    {
                        Eigen::Vector3f xl  = mBoundaryParticles->getValue(nl.first);
                        Eigen::Vector3f xil = xi - xl;
                        voli += W.W(xil);
                    }
                float massi = mRho0 / voli;
                mMasses->setValue(i, massi);
            }
#else
            mMasses->setSize(mPositions->getSize());
            float volume = std::powf(mSupportRadius / 2, 3) * 0.8;
            for (Eigen::Index i = 0; i < mPositions->getSize(); ++i)
                mMasses->setValue(i, volume * mRho0);
#endif
        }

        /**
         * @brief Advances the simulation one time step forward.
         * @param elapsedTime Elapsed time in milliseconds during the last frame.
         * @param totalTime Total time in milliseconds since the beginning of the first frame.
         * @param timeStep Time step of the simulation. Restarts when resetting the simulation.
         */
        void advance(double elapsedTime, double totalTime, int64_t timeStep) override
        {
            double dt = mStepSize;
            float h2  = mSupportRadius * mSupportRadius;

            // build data structure for nearest neighbor search
            auto knn = std::make_shared<NearestNeighbors3f>();
            knn->setPoints(mPositions);
            CubicKernel W;
            W.setRadius(mSupportRadius);

            // density estimation
#ifndef _DEBUG
#pragma omp parallel for
#endif
            for (Eigen::Index i = 0; i < mPositions->getSize(); ++i)
            {
                // get current position and velocity
                Eigen::Vector3f xi = mPositions->getValue(i);
                float rhoi         = 0;

                // compute density from boundary
                NearestNeighbors3f::RadiusResult nn;
                if (mBoundaryKNN->closestRadius(xi, mSupportRadius, nn) != 0)
                {
                    for (auto& n : nn)
                    {
                        float massk         = mBoundaryMasses->getValue(n.first).x();
                        Eigen::Vector3f xk  = mBoundaryParticles->getValue(n.first);
                        Eigen::Vector3f xik = xi - xk;
                        rhoi += massk * W.W(xik);
                    }
                }

                // compute density from interior particles
                if (knn->closestRadius(xi, mSupportRadius, nn) != 0)
                {
                    for (auto& n : nn)
                    {
                        float massj         = mMasses->getValue(n.first).x();
                        Eigen::Vector3f xj  = mPositions->getValue(n.first);
                        Eigen::Vector3f xij = xi - xj;
                        rhoi += massj * W.W(xij);
                    }
                }
                rhoi = std::max(rhoi, mRho0); // clamp
                mDensities->setValue(i, rhoi);
            }

            // pressure estimation based on EOS equation (WCSPH)
#ifndef _DEBUG
#pragma omp parallel for
#endif
            for (Eigen::Index i = 0; i < mPositions->getSize(); ++i)
            {
                // TODO: compute pressure p_i from mDensities, using mStiffness (k), mRho0 (rest density) and mExponent (gamma).
                float rhoi = mDensities->getValue(i).x();
                float pi = mStiffness * (std::pow(rhoi / mRho0, mExponent) - 1.0f);
                mPressures->setValue(i, pi);
            }

// pressure and viscosity acceleration
#ifndef _DEBUG
#pragma omp parallel for
#endif
            for (Eigen::Index i = 0; i < mPositions->getSize(); ++i)
            {
                // get current position and velocity
                Eigen::Vector3f xi = mPositions->getValue(i);
                Eigen::Vector3f vi = mVelocities->getValue(i);
                float pi           = mPressures->getValue(i).x();
                float rhoi         = mDensities->getValue(i).x();

                // find closest points
                Eigen::Vector3f ai_p(0, 0, 0);
                NearestNeighbors3f::RadiusResult nn;
                if (knn->closestRadius(xi, mSupportRadius, nn) != 0)
                {
                    for (auto& n : nn)
                    {
                        Eigen::Vector3f xj = mPositions->getValue(n.first);
                        Eigen::Vector3f vj = mVelocities->getValue(n.first);
                        float pj           = mPressures->getValue(n.first).x();
                        float rhoj         = mDensities->getValue(n.first).x();
                        float massj        = mMasses->getValue(n.first).x();

                        // TODO: pressure acceleration
                        Eigen::Vector3f xij = xi - xj;
                        Eigen::Vector3f grad_W_ij = W.gradW(xij);
                        ai_p += -massj * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj)) * grad_W_ij;

                        // TODO: viscosity acceleration
                        Eigen::Vector3f vij = vi - vj;
                        float dot_product = xij.dot(vij);
                        float xij_norm_sq = xij.squaredNorm();
                        float eta = 0.01f * mSupportRadius * mSupportRadius; // small epsilon to avoid division by zero
                        float Pi_ij = mViscosity * dot_product / (rhoj) / (xij_norm_sq + eta);
                        ai_p += massj * Pi_ij * grad_W_ij;

                    }
                }

                // boundary handling
                if (mBoundaryKNN->closestRadius(xi, mSupportRadius, nn) != 0)
                {
                    for (auto& n : nn)
                    {
                        Eigen::Vector3f xk = mBoundaryParticles->getValue(n.first);
                        float pk           = pi;    // pressure-mirroring
                        float rhok         = mRho0; // boundary is liquid at rest
                        float massk        = mBoundaryMasses->getValue(n.first).x();

                        // TODO: pressure acceleration (Akinci)
                        Eigen::Vector3f xik = xi - xk;
                        Eigen::Vector3f grad_W_ik = W.gradW(xik);
                        ai_p += -massk * (pi / (rhoi * rhoi) + pi / (mRho0 * mRho0)) * grad_W_ik;

                        // TODO: viscosity acceleration (Akinci)
                        Eigen::Vector3f vik = vi; // boundary velocity is zero, so vi - 0 = vi
                        float dot_product = xik.dot(vik);
                        float xik_norm_sq = xik.squaredNorm();
                        float eta = 0.01f * mSupportRadius * mSupportRadius; // small epsilon to avoid division by zero
                        float Pi_ik = mViscosity * dot_product / (mRho0) / (xik_norm_sq + eta);
                        ai_p += massk * Pi_ik * grad_W_ik;
                    }
                }

                // store acceleration
                mAccelerations->setValue(i, ai_p);
            }

            // advance particles
#ifndef _DEBUG
#pragma omp parallel for
#endif
            for (Eigen::Index i = 0; i < mPositions->getSize(); ++i)
            {
                // get current position and velocity
                Eigen::Vector3f xi    = mPositions->getValue(i);
                Eigen::Vector3f vi    = mVelocities->getValue(i);
                Eigen::Vector3f ai_p  = mAccelerations->getValue(i);
                Eigen::Vector3f ai_np = mGravity;

                // symplectic euler update
                Eigen::Vector3f vi_new = vi + dt * (ai_p + ai_np);
                Eigen::Vector3f xi_new = xi + dt * vi_new;

                // write result
                mPositions->setValue(i, xi_new);
                mVelocities->setValue(i, vi_new);
            }

            // notify that the geometry has changed
            mSpheres->components.get<Geometry>()->markChanged();
        }

        /**
         * @brief Adds graphical user interface elements with imgui.
         */
        void gui() override
        {
            ImGui::PushItemWidth(100);

            double stepSizeMin = 1E-3, stepSizeMax = 1E-1;
            ImGui::SliderScalar("dt", ImGuiDataType_Double, &mStepSize, &stepSizeMin, &stepSizeMax);

            ImGui::SliderFloat("stiffness", &mStiffness, 0, 200000);
            ImGui::SliderFloat("exponent", &mExponent, 0, 10);
            ImGui::SliderFloat("viscosity", &mViscosity, 0, 50);

            ImGui::PopItemWidth();
        }

        ///**
        // * @brief Helper function to create a diffuse sphere.
        // * @param position Position of the sphere in world space.
        // * @param scale Scale of the sphere in world space.
        // * @param color Color of the sphere.
        // * @return Reference to the sphere shape that was created.
        // */
        // std::shared_ptr<Actor> addSphere(const Eigen::Vector3d& position, const double& scale, const Spectrum& color)
        //{
        //    // create a sphere and assign the bsdf
        //    auto actor = std::make_shared<Actor>();
        //    actor->components.add(std::make_shared<SphereGeometry>(scale));
        //    actor->components.add(std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(color)));
        //    actor->components.add(std::make_shared<Transform>(position));
        //    scene->actors.push_back(actor);
        //    return actor;
        //}

    private:
        /**
         * @brief Fluid simulation domain. We will assume that all particles will be constrained into a box.
         */
        Eigen::AlignedBox3f mDomain;

        /**
         * @brief Positions of all particles.
         */
        std::shared_ptr<vislab::Array3f> mPositions;

        /**
         * @brief Velocities of all particles.
         */
        std::shared_ptr<vislab::Array3f> mVelocities;

        /**
         * @brief Pressure accelerations of all particles.
         */
        std::shared_ptr<vislab::Array3f> mAccelerations;

        /**
         * @brief Densities of all particles.
         */
        std::shared_ptr<vislab::Array1f> mDensities;

        /**
         * @brief Masses of all particles.
         */
        std::shared_ptr<vislab::Array1f> mMasses;

        /**
         * @brief Pressures of all particles.
         */
        std::shared_ptr<vislab::Array1f> mPressures;

        /**
         * @brief Nearest neighbor search data structure for the static boundary particles.
         */
        std::shared_ptr<NearestNeighbors3f> mBoundaryKNN;

        /**
         * @brief Set of static boundary particles.
         */
        std::shared_ptr<vislab::Array3f> mBoundaryParticles;

        /**
         * @brief Masses of all boundary particles.
         */
        std::shared_ptr<vislab::Array1f> mBoundaryMasses;

        /**
         * @brief Sphere geometries that are used for rendering.
         */
        std::shared_ptr<Actor> mSpheres;

        /**
         * @brief Support radius of the SPH kernels.
         */
        float mSupportRadius;

        /**
         * @brief Rest density.
         */
        float mRho0;
        /**
         * @brief Exponent of the tait equation.
         */
        float mExponent;

        /**
         * @brief Stiffness of the density to pressure conversion.
         */
        float mStiffness;

        /**
         * @brief Viscosity of the material.
         */
        float mViscosity;

        /**
         * @brief Integration step size.
         */
        double mStepSize;

        /**
         * @brief Gravity vector.
         */
        Eigen::Vector3f mGravity;
    };
}

int main()
{
    vislab::Init();

    physsim::PhyssimWindow window(
        800,     // width
        600,     // height
        "7_sph", // title
        std::make_shared<physsim::SmoothedParticleHydrodynamicsSimulation>(),
        false // fullscreen
    );

    return window.run();
}
