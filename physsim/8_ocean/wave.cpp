#include "wave.hpp"

namespace physsim
{
    double Wave::waveNumber() const { return 2 * (float)EIGEN_PI / waveLength; }

    Eigen::Vector2d Wave::direction() const { return Eigen::Vector2d(std::cos(angle), std::sin(angle)); }

    Eigen::Vector2d Wave::waveVector() const { return direction() * waveNumber(); }

    double Wave::angularFrequency() const
    {
        double k = waveNumber();
        return k * dispersion.phaseSpeed(k);
    }

    Eigen::Vector3d Wave::offset(const Eigen::Vector2d& position, const double& t) const
    {
        // TODO: compute the offset for the displaced wave from the given position.
        Eigen::Vector2d ki = waveVector();
        double k = waveNumber();
        Eigen::Vector3d waveHeight;
        waveHeight.z() = amplitude * std::cos(ki.dot(position) - angularFrequency() * t + phase);

        waveHeight.x() = - steepness * ki.x() / k * amplitude * std::sin(ki.dot(position) - angularFrequency() * t + phase);
        waveHeight.y() = - steepness * ki.y() / k * amplitude * std::sin(ki.dot(position) - angularFrequency() * t + phase);
        return waveHeight;
    }
}
