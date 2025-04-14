#include <init_vislab.hpp>

#include "ocean.hpp"
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

#include <random>

using namespace vislab;

namespace physsim
{
    /**
     * @brief Ocean simulation using a simple Gerstner model.
     */
    class OceanSimulation : public Simulation
    {
    public:
        /**
         * @brief Initializes the scene.
         */
        void init() override
        {
            // grid resolution for the ocean surface
            Eigen::Vector2i resolution(100, 100);

            // allocate ocean simulation
            mOcean = std::make_shared<Ocean>(resolution, Eigen::Vector3f(-1, 1, 0.2), Eigen::Vector3f(2, 0, 0), Eigen::Vector3f(0, -2, 0));

            // add random waves
            std::default_random_engine rng;
            std::uniform_real_distribution<double> rnd;
            for (int i = 0; i < 50; ++i)
            {
                double lambda    = rnd(rng) * 1 + 0.1;
                double amplitude = rnd(rng) * 0.007;
                double angle     = rnd(rng) * 0.8 - 0.4;
                double steepness = rnd(rng) * 0.8;
                mOcean->waves.push_back(Wave{
                    lambda,    // Wavelength lambda in [m]
                    amplitude, // Amplitude in [m]
                    angle,     // direction angle in [rad]
                    steepness, // Steepness Q in [0,1]
                    0.0        // Phase offset
                });
            }

            // create triangle shape
            auto meshActor = std::make_shared<Actor>();
            meshActor->components.add(mOcean->mesh);
            meshActor->components.add(std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(0.5, 1, 1))));
            meshActor->components.add(std::make_shared<Transform>());
            scene->actors.push_back(meshActor);

            // create ground plane
            auto rectActor = std::make_shared<Actor>("ground");
            rectActor->components.add(std::make_shared<RectangleGeometry>());
            rectActor->components.add(std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(1, 1, 1))));
            rectActor->components.add(std::make_shared<Transform>());

            // create a point light
            auto lightActor = std::make_shared<Actor>("light");
            lightActor->components.add(std::make_shared<PointLight>(Spectrum(100.0, 100.0, 100.0)));
            lightActor->components.add(std::make_shared<Transform>(Eigen::Vector3d(1, -5, 5)));

            // create a camera
            auto cameraActor = std::make_shared<Actor>("camera");
            auto camera      = std::make_shared<PerspectiveCamera>();
            camera->setLookAt(Eigen::Vector3d(0, 0, 0));
            camera->setPosition(Eigen::Vector3d(0, -2, 2));
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
        }

        /**
         * @brief Advances the simulation one time step forward.
         * @param elapsedTime Elapsed time in milliseconds during the last frame.
         * @param totalTime Total time in milliseconds since the beginning of the first frame.
         * @param timeStep Time step of the simulation. Restarts when resetting the simulation.
         */
        void advance(double elapsedTime, double totalTime, int64_t timeStep) override
        {
            // update ocean mesh
            mOcean->advance(totalTime);

            // notify that the geometry has changed
            mOcean->mesh->markChanged();
        }

        /**
         * @brief Adds graphical user interface elements with imgui.
         */
        void gui() override
        {
        }

    private:
        /**
         * @brief Ocean simualtion.
         */
        std::shared_ptr<Ocean> mOcean;
    };
}

int main()
{
    vislab::Init();

    physsim::PhyssimWindow window(
        800,       // width
        600,       // height
        "8_ocean", // title
        std::make_shared<physsim::OceanSimulation>(),
        false // fullscreen
    );

    return window.run();
}
