#include <init_vislab.hpp>

#include "cloth.hpp"
#include "physsim_window.hpp"
#include "simulation.hpp"

#include <imgui.h>
#include <vislab/core/array.hpp>
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
     * @brief Cloth simulation.
     */
    class ClothSimulation : public Simulation
    {
    public:
        /**
         * @brief Initializes the scene.
         */
        void init() override
        {
            // initial simulation parameters
            mStepSize = 1E-3;
            mGravity << 0, 0, -9.81f;
            mMass      = 0.01;
            mStiffness = 3000.;
            mDamping   = 0.01;

            // compute rest lengths from the resolution
            Eigen::Vector2i resolution(40, 60);
            double Lx = 1. / (resolution.x() - 1);
            double Ly = 1.5 / (resolution.y() - 1);

            // allocate cloth simulations
            mCloth[0] = std::make_shared<Cloth>(resolution, Cloth::ETopology::Structural, Eigen::Vector3f(0.1, 0, 2.5), Eigen::Vector3f(Lx, 0, 0), Eigen::Vector3f(0, -Ly, 0));
            mCloth[1] = std::make_shared<Cloth>(resolution, Cloth::ETopology::Diagonal, Eigen::Vector3f(-1.1, 0, 2.5), Eigen::Vector3f(Lx, 0, 0), Eigen::Vector3f(0, -Ly, 0));
            mCloth[0]->pin(Eigen::Vector2i(0, 0));
            mCloth[1]->pin(Eigen::Vector2i(0, 0));
            mCloth[0]->pin(Eigen::Vector2i(resolution.x() - 1, 0));
            mCloth[1]->pin(Eigen::Vector2i(resolution.x() - 1, 0));

            // create shapes
            for (int i = 0; i < 2; ++i)
            {
                auto meshActor = std::make_shared<Actor>();
                meshActor->components.add(mCloth[i]->mesh);
                meshActor->components.add(std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(0.5 + i * 0.5, 1, 1 - i * 0.5))));
                meshActor->components.add(std::make_shared<Transform>());
                scene->actors.push_back(meshActor);
            }

            // create ground plane
            auto rectActor = std::make_shared<Actor>("ground");
            rectActor->components.add(std::make_shared<RectangleGeometry>());
            rectActor->components.add(std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(1, 1, 1))));
            rectActor->components.add(std::make_shared<Transform>(Eigen::Vector3d::Zero(), Eigen::Quaterniond::Identity(), Eigen::Vector3d(5., 5., 5.)));

            // create a point light
            auto lightActor = std::make_shared<Actor>("light");
            lightActor->components.add(std::make_shared<PointLight>(Spectrum(100.0, 100.0, 100.0)));
            lightActor->components.add(std::make_shared<Transform>(Eigen::Vector3d(1, -5, 5)));

            // create a camera
            auto cameraActor = std::make_shared<Actor>("camera");
            auto camera      = std::make_shared<PerspectiveCamera>();
            camera->setLookAt(Eigen::Vector3d(0, 0, 1.3));
            camera->setPosition(Eigen::Vector3d(0, -4, 2));
            camera->setUp(Eigen::Vector3d(0, 0, 1));
            camera->setNear(0.01);
            camera->setFar(100);
            cameraActor->components.add(camera);

            // add elements to scene
            scene->actors.push_back(rectActor);
            scene->actors.push_back(lightActor);
            scene->actors.push_back(cameraActor);
        }

        /**
         * @brief Restarts the simulation.
         */
        void restart() override
        {
            for (int i = 0; i < 2; ++i)
                mCloth[i]->reset();
        }

        /**
         * @brief Advances the simulation one time step forward.
         * @param elapsedTime Elapsed time in milliseconds during the last frame.
         * @param totalTime Total time in milliseconds since the beginning of the first frame.
         * @param timeStep Time step of the simulation. Restarts when resetting the simulation.
         */
        void advance(double elapsedTime, double totalTime, int64_t timeStep) override
        {
            // process both cloths
            for (int i = 0; i < 2; ++i)
            {
                // advance positions
                mCloth[i]->damping   = mDamping;
                mCloth[i]->stiffness = mStiffness;
                mCloth[i]->mass      = mMass;
                mCloth[i]->gravity   = mGravity;
                mCloth[i]->advance(mStepSize);
            }
        }

        /**
         * @brief Adds graphical user interface elements with imgui.
         */
        void gui() override
        {
            ImGui::PushItemWidth(100);

            double stepSizeMin = 1E-3, stepSizeMax = 3E-3;
            ImGui::SliderScalar("dt", ImGuiDataType_Double, &mStepSize, &stepSizeMin, &stepSizeMax);

            double dampingMin = 0, dampingMax = 5E-1;
            ImGui::SliderScalar("damping", ImGuiDataType_Double, &mDamping, &dampingMin, &dampingMax);

            double stiffnessMin = 0, stiffnessMax = 5000;
            ImGui::SliderScalar("stiffness", ImGuiDataType_Double, &mStiffness, &stiffnessMin, &stiffnessMax);

            double massMin = 1E-2, massMax = 1E-1;
            ImGui::SliderScalar("mass", ImGuiDataType_Double, &mMass, &massMin, &massMax);

            ImGui::PopItemWidth();
        }

    private:
        /**
         * @brief Cloth solvers.
         */
        std::shared_ptr<Cloth> mCloth[2];

        /**
         * @brief Integration step size.
         */
        double mStepSize;

        /**
         * @brief Stiffness coefficient of the springs.
         */
        double mStiffness;

        /**
         * @brief Damping coefficient of the springs.
         */
        double mDamping;

        /**
         * @brief Mass of all bodies.
         */
        double mMass;

        /**
         * @brief Gravitational acceleration.
         */
        Eigen::Vector3f mGravity;
    };
}

int main()
{
    vislab::Init();

    physsim::PhyssimWindow window(
        800,       // width
        600,       // height
        "5_cloth", // title
        std::make_shared<physsim::ClothSimulation>(),
        false // fullscreen
    );

    return window.run();
}
