#include <init_vislab.hpp>

#include "physsim_window.hpp"
#include "rigid_body.hpp"
#include "rigid_body_integrator.hpp"
#include "simulation.hpp"

#include <imgui.h>
#include <vislab/core/array.hpp>
#include <vislab/geometry/face_normals.hpp>
#include <vislab/geometry/surfaces.hpp>
#include <vislab/graphics/actor.hpp>
#include <vislab/graphics/const_texture.hpp>
#include <vislab/graphics/diffuse_bsdf.hpp>
#include <vislab/graphics/perspective_camera.hpp>
#include <vislab/graphics/point_light.hpp>
#include <vislab/graphics/rectangle_geometry.hpp>
#include <vislab/graphics/scene.hpp>
#include <vislab/graphics/transform.hpp>
#include <vislab/graphics/trimesh_geometry.hpp>

using namespace vislab;

namespace physsim
{
    /**
     * @brief Rigid body update with gyroscopic force.
     */
    class GyroscopeSimulation : public Simulation
    {
    public:
        /**
         * @brief Enumeration of integration methods.
         */
        enum EMethod
        {
            ExplicitEuler,
            SymplecticEuler,
            Implicit,
        };

        /**
         * @brief Initializes the scene.
         */
        void init() override
        {
            // initial simulation parameters
            mMethod   = EMethod::ExplicitEuler;
            mStepSize = 1E-2;
            mRigidBody.setMass(1.);
            mRigidBody.setInertiaBody(Eigen::Matrix3d::diag(Eigen::Vector3d(24.1449, 28.436, 118.812)).block(0, 0, 3, 3));

            // create ground plane
            auto rectActor = std::make_shared<Actor>("ground");
            rectActor->components.add(std::make_shared<RectangleGeometry>());
            rectActor->components.add(std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(1, 1, 1))));
            Eigen::Matrix4d transform;
            transform << 15, 0, 0, 0,
                0, 0, 15, -3.5,
                0, 15, 0, 0,
                0, 0, 0, 1;
            rectActor->components.add(std::make_shared<Transform>(transform));

            // create a T-shaped mesh
            auto surfaces = std::make_shared<Surfaces3f>();
            auto surface  = surfaces->createSurface();
            surface->positions->setValues({ { -1, -1, -1 },
                                            { 1, -1, -1 },
                                            { -1, 1, -1 },
                                            { 1, 1, -1 },
                                            { -1, -1, 1 },
                                            { 1, -1, 1 },
                                            { -1, 1, 1 },
                                            { 1, 1, 1 },
                                            { -3, -1, -1 },
                                            { 3, -1, -1 },
                                            { -3, 1, -1 },
                                            { 3, 1, -1 },
                                            { -3, -1, 1 },
                                            { 3, -1, 1 },
                                            { -3, 1, 1 },
                                            { 3, 1, 1 },
                                            { -1, -3, -1 },
                                            { 1, -3, -1 },
                                            { -1, -3, 1 },
                                            { 1, -3, 1 } });

            // bring center of mass to (0,0,0)
            auto& positions      = surface->positions->getData();
            Eigen::Vector3f mean = Eigen::Vector3f::Zero();
            for (const auto& p : positions)
                mean += p;
            mean /= positions.size();
            for (auto& p : positions)
                p -= mean;
            surfaces->recomputeBoundingBox();

            // create the index buffer
            surface->indices->setValues({
                0,
                2,
                1,
                2,
                3,
                1,
                9,
                11,
                13,
                11,
                15,
                13,
                2,
                6,
                3,
                6,
                7,
                3,
                5,
                7,
                4,
                7,
                6,
                4,
                12,
                14,
                8,
                14,
                10,
                8,
                18,
                16,
                17,
                18,
                17,
                19,
                8,
                10,
                0,
                10,
                2,
                0,
                3,
                9,
                1,
                3,
                11,
                9,
                12,
                6,
                14,
                12,
                4,
                6,
                7,
                5,
                15,
                15,
                5,
                13,
                14,
                2,
                10,
                2,
                14,
                6,
                8,
                0,
                12,
                4,
                12,
                0,
                3,
                7,
                11,
                7,
                15,
                11,
                5,
                1,
                9,
                5,
                9,
                13,
                0,
                1,
                17,
                0,
                17,
                16,
                19,
                5,
                18,
                5,
                4,
                18,
                0,
                16,
                4,
                4,
                16,
                18,
                1,
                5,
                19,
                17,
                1,
                19,
            });

            // compute face normals
            FaceNormals3f faceNormals;
            faceNormals.inputSurfaces.setData(surfaces);
            faceNormals.outputSurfaces.setData(std::make_shared<Surfaces3f>());
            auto info = faceNormals.update();
            assert(info.success());

            // create a triangle mesh
            auto meshActor = std::make_shared<Actor>();
            meshActor->components.add(std::make_shared<TrimeshGeometry>(faceNormals.outputSurfaces.getData()));
            meshActor->components.add(std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(1, 0, 0))));
            meshActor->components.add(std::make_shared<Transform>());
            mRigidBody.setActor(meshActor);

            // create a point light
            auto lightActor = std::make_shared<Actor>("light");
            lightActor->components.add(std::make_shared<PointLight>(Spectrum(100.0, 100.0, 100.0)));
            lightActor->components.add(std::make_shared<Transform>(Eigen::Vector3d(-4, 5, 4)));

            // create a camera
            auto cameraActor = std::make_shared<Actor>("camera");
            auto camera      = std::make_shared<PerspectiveCamera>();
            camera->setLookAt(Eigen::Vector3d(0, 0, 0));
            camera->setPosition(Eigen::Vector3d(-10, 0, 0));
            camera->setUp(Eigen::Vector3d(0, 1, 0));
            camera->setNear(0.01);
            camera->setFar(100);
            cameraActor->components.add(camera);

            // add elements to scene
            scene->actors.push_back(rectActor);
            scene->actors.push_back(meshActor);
            scene->actors.push_back(lightActor);
            scene->actors.push_back(cameraActor);
        }

        /**
         * @brief Restarts the simulation.
         */
        void restart() override
        {
            mRigidBody.setRotation(Eigen::Quaterniond::Identity());
            mRigidBody.setAngularVelocity(Eigen::Vector3d(0.1, 10, 0));
        }

        /**
         * @brief Advances the simulation one time step forward.
         * @param elapsedTime Elapsed time in milliseconds during the last frame.
         * @param totalTime Total time in milliseconds since the beginning of the first frame.
         * @param timeStep Time step of the simulation. Restarts when resetting the simulation.
         */
        void advance(double elapsedTime, double totalTime, int64_t timeStep) override
        {
            switch (mMethod)
            {
            case EMethod::ExplicitEuler:
            {
                explicitEuler(mRigidBody, mStepSize);
                break;
            }
            case EMethod::SymplecticEuler:
            {
                symplecticEuler(mRigidBody, mStepSize);
                break;
            }
            case EMethod::Implicit:
            {
                implicitEuler(mRigidBody, mStepSize);
                break;
            }
            }
        }

        /**
         * @brief Adds graphical user interface elements with imgui.
         */
        void gui() override
        {
            ImGui::PushItemWidth(100);

            ImGui::Combo("method", (int*)&mMethod, "explicit euler\0symplectic euler\0implicit euler\0\0");

            double stepSizeMin = 1E-3, stepSizeMax = 1E-1;
            ImGui::SliderScalar("dt", ImGuiDataType_Double, &mStepSize, &stepSizeMin, &stepSizeMax);

            ImGui::PopItemWidth();
        }

    private:
        /**
         * @brief Numerical integration method.
         */
        EMethod mMethod;

        /**
         * @brief Rigid body state.
         */
        RigidBody mRigidBody;

        /**
         * @brief Integration step size.
         */
        double mStepSize;
    };
}

int main()
{
    vislab::Init();

    physsim::PhyssimWindow window(
        800,           // width
        600,           // height
        "3_gyroscope", // title
        std::make_shared<physsim::GyroscopeSimulation>(),
        false // fullscreen
    );

    return window.run();
}
