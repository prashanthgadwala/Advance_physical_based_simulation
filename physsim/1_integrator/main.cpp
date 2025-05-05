#include <init_vislab.hpp>

#include "physsim_window.hpp"
#include "simulation.hpp"

#include <imgui.h>
#include <vislab/graphics/actor.hpp>
#include <vislab/graphics/const_texture.hpp>
#include <vislab/graphics/diffuse_bsdf.hpp>
#include <vislab/graphics/perspective_camera.hpp>
#include <vislab/graphics/point_light.hpp>
#include <vislab/graphics/rectangle_geometry.hpp>
#include <vislab/graphics/scene.hpp>
#include <vislab/graphics/sphere_geometry.hpp>
#include <vislab/graphics/transform.hpp>

using namespace vislab;

namespace physsim
{
    /**
     * @brief Comparison of explicit numerical integration methods.
     */
    class IntegratorSimulation : public Simulation
    {
    public:
        /**
         * @brief Enumeration of integration methods.
         */
        enum EMethod
        {
            Analytic,
            ExplicitEuler,
            SymplecticEuler
        };

        /**
         * @brief Initializes the scene.
         */
        void init() override
        {
            // initial simulation parameters
            mStepSize        = 2E-2;
            mInitialPosition = Eigen::Vector3d(0, 0, 0);
            mInitialVelocity = Eigen::Vector3d(0, 3, 4);
            mGravity         = Eigen::Vector3d(0, 0, -9.8065);

            // create ground plane
            auto rectActor = std::make_shared<Actor>("ground");
            rectActor->components.add(std::make_shared<RectangleGeometry>());
            rectActor->components.add(std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(1, 1, 1))));
            rectActor->components.add(std::make_shared<Transform>(Eigen::Vector3d(0, 1.5, -1), Eigen::Quaterniond::Identity(), Eigen::Vector3d(10, 10, 10)));

            // create a point light
            auto lightActor = std::make_shared<Actor>("light");
            lightActor->components.add(std::make_shared<PointLight>(Spectrum(100.0, 100.0, 100.0)));
            lightActor->components.add(std::make_shared<Transform>(Eigen::Vector3d(-5, -3, 3)));

            // create a camera
            auto cameraActor = std::make_shared<Actor>("camera");
            auto camera      = std::make_shared<PerspectiveCamera>();
            camera->setLookAt(Eigen::Vector3d(0, 1.5, 0));
            camera->setPosition(Eigen::Vector3d(-3.5, 1.5, 0));
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
            // delete all sphere shapes (except for the first three: ground, light, camera)
            scene->actors.resize(3);

            // reset the simulation
            mTime        = 0;
            mPosition[0] = mPosition[1] = mPosition[2] = mInitialPosition;
            mVelocity[0] = mVelocity[1] = mVelocity[2] = mInitialVelocity;
        }

        /**
         * @brief Advances the simulation one time step forward.
         * @param elapsedTime Elapsed time in milliseconds during the last frame.
         * @param totalTime Total time in milliseconds since the beginning of the first frame.
         * @param timeStep Time step of the simulation. Restarts when resetting the simulation.
         */
        void advance(double elapsedTime, double totalTime, int64_t timeStep) override
        {
            if (mTime > 1)
                return;

            // TODO: analytical solution
            // x(t) = x_0 + v_0*t + 0.5*a*t^2
            // v(t) = v_0 + a*t
            mPosition[EMethod::Analytic] = mInitialPosition + mInitialVelocity * mTime + 0.5 * mGravity * mTime * mTime; // <- adjust
            mVelocity[EMethod::Analytic] = mInitialVelocity + mTime * mGravity; // <- adjust
            addSphere(mPosition[EMethod::Analytic], 0.03, Spectrum(1, 0, 0));

            // TODO: explicit euler
            // x' = x + dt*v
            // v' = v + dt*a
            mPosition[EMethod::ExplicitEuler] = mPosition[EMethod::ExplicitEuler] + mStepSize * mVelocity[EMethod::ExplicitEuler];
            mVelocity[EMethod::ExplicitEuler] = mVelocity[EMethod::ExplicitEuler] + mStepSize * mGravity;
            addSphere(mPosition[EMethod::ExplicitEuler], 0.03, Spectrum(0, 1, 0));

            // TODO: symplectic euler
            // v' = v + dt*a
            // x' = x + dt*v'
            mVelocity[EMethod::SymplecticEuler] = mVelocity[EMethod::SymplecticEuler] + mStepSize * mGravity;
            mPosition[EMethod::SymplecticEuler] = mPosition[EMethod::SymplecticEuler] + mStepSize * mVelocity[EMethod::SymplecticEuler];
            addSphere(mPosition[EMethod::SymplecticEuler], 0.03, Spectrum(0, 0, 1));

            mTime += mStepSize;
        }

        /**
         * @brief Adds graphical user interface elements with imgui.
         */
        void gui() override
        {
            ImGui::PushItemWidth(100);

            double stepSizeMin = 1E-3, stepSizeMax = 1E-1;
            ImGui::SliderScalar("dt", ImGuiDataType_Double, &mStepSize, &stepSizeMin, &stepSizeMax);

            static ImVec4 color_analytic(1, 0, 0, 1);
            ImGui::ColorEdit4("color_analytic", (float*)&color_analytic, ImGuiColorEditFlags_NoPicker | ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel | ImGuiColorEditFlags_NoTooltip);
            ImGui::SameLine();
            ImGui::Text("analytic");

            static ImVec4 color_explicit(0, 1, 0, 1);
            ImGui::ColorEdit4("color_explicit", (float*)&color_explicit, ImGuiColorEditFlags_NoPicker | ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel | ImGuiColorEditFlags_NoTooltip);
            ImGui::SameLine();
            ImGui::Text("explicit euler");

            static ImVec4 color_symplectic(0, 0, 1, 1);
            ImGui::ColorEdit4("color_symplectic", (float*)&color_symplectic, ImGuiColorEditFlags_NoPicker | ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel | ImGuiColorEditFlags_NoTooltip);
            ImGui::SameLine();
            ImGui::Text("symplectic euler");

            ImGui::PopItemWidth();
        }

        /**
         * @brief Helper function to create a diffuse sphere.
         * @param position Position of the sphere in world space.
         * @param scale Scale of the sphere in world space.
         * @param color Color of the sphere.
         * @return Reference to the sphere shape that was created.
         */
        std::shared_ptr<Actor> addSphere(const Eigen::Vector3d& position, const double& scale, const Spectrum& color)
        {
            // create a sphere and assign the bsdf
            auto actor          = std::make_shared<Actor>();
            actor->components.add(std::make_shared<SphereGeometry>(scale));
            actor->components.add(std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(color)));
            actor->components.add(std::make_shared<Transform>(position));
            scene->actors.push_back(actor);
            return actor;
        }

    private:
        /**
         * @brief Position of the body.
         */
        Eigen::Vector3d mPosition[3];

        /**
         * @brief Linear velocity of the body.
         */
        Eigen::Vector3d mVelocity[3];

        /**
         * @brief Simulation time.
         */
        double mTime;

        /**
         * @brief Initial position of the body.
         */
        Eigen::Vector3d mInitialPosition;

        /**
         * @brief Initial linear velocity of the body.
         */
        Eigen::Vector3d mInitialVelocity;

        /**
         * @brief Integration step size.
         */
        double mStepSize;

        /**
         * @brief Gravity vector.
         */
        Eigen::Vector3d mGravity;
    };
}

int main()
{
    vislab::Init();

    physsim::PhyssimWindow window(
        800,            // width
        600,            // height
        "1_integrator", // title
        std::make_shared<physsim::IntegratorSimulation>(),
        false // fullscreen
    );

    return window.run();
}
