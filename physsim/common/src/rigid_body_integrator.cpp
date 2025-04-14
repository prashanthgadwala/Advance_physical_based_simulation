#include "rigid_body_integrator.hpp"

#include "rigid_body.hpp"

namespace physsim
{
    void explicitEuler(RigidBody& body, double stepSize)
    {
        // get current position and rotation of the body
        Eigen::Vector3d x    = body.position();
        Eigen::Quaterniond q = body.rotation();

        // get current linear and angular velocity of the body
        Eigen::Vector3d v = body.linearVelocity();
        Eigen::Vector3d w = body.angularVelocity();

        // update position of the body using the linear velocity and update body accordingly
        Eigen::Vector3d xnew = x + stepSize * v;
        body.setPosition(xnew);

        // quaternion-based angular velocity update of rotation and update body accordingly
        Eigen::Quaterniond wq(0, w.x(), w.y(), w.z());
        Eigen::Quaterniond qnew = (q + 0.5 * stepSize * wq * q).normalized();
        body.setRotation(qnew);

        // get current linear and angular momentum of body
        Eigen::Vector3d p = body.linearMomentum();
        Eigen::Vector3d l = body.angularMomentum();

        // get force and torque that are currently applied to body
        Eigen::Vector3d f = body.force();
        Eigen::Vector3d t = body.torque();

        // compute new linear momentum
        Eigen::Vector3d pnew = p + stepSize * f;

        // convert from linear momentum to linear velocity and update the body accordingly
        Eigen::Vector3d vnew = body.massInverse() * pnew;
        body.setLinearVelocity(vnew);

        // compute new angular momentum
        Eigen::Matrix3d I    = body.inertiaWorld();
        Eigen::Vector3d lnew = l + stepSize * (t - w.cross(I * w));

        // convert from angular momentum to angular velocity and update the body accordingly
        Eigen::Vector3d wnew = body.inertiaWorldInverse() * lnew;
        body.setAngularVelocity(wnew);
    }

    void symplecticEuler(RigidBody& body, double stepSize)
    {
        // TODO: get current position and rotation of the body

        // TODO: get current linear and angular momentum of body

        // TODO: get force and torque that are currently applied to body

        // TODO: compute new linear momentum

        // TODO: convert from linear momentum to linear velocity and update the body accordingly

        // TODO: compute new angular momentum

        // TODO: convert from angular momentum to angular velocity and update the body accordingly

        // TODO: update position of the body using the linear velocity and update body accordingly

        // TODO: quaternion-based angular velocity update of rotation and update body accordingly
    }

    void implicitEuler(RigidBody& body, double stepSize)
    {
        // See for explanations: https://www.gdcvault.com/play/1022196/Physics-for-Game-Programmers-Numerical

        // TODO: get current position and rotation of the body

        // TODO: get current linear and angular momentum of body

        // TODO: get force and torque that are currently applied to body

        // TODO: compute new linear momentum

        // TODO: convert from linear momentum to linear velocity and update the body accordingly

        // TODO: Convert current angular velocity to body coordinates (initial guess wb0)

        // TODO: Compute residual vector f(wb0) from the from angular velocity in body-coordinates

        // TODO: Compute the Jacobian of f at wb. You can use the function "skew".

        // TODO: Linearly solve for the update step delta_wb, for example using a QR decomposition

        // TODO: Apply the Newton-Raphson iteration by adding delta_wb to the current angular velocity

        // TODO: Transform the angular velocity back to world coordinates

        // TODO: explicitly integrate the torque and update the body accordingly

        // TODO: update position of the body using the linear velocity and update body accordingly

        // TODO: quaternion-based angular velocity update of rotation
    }
}
