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

        Eigen::Vector3d x    = body.position();
        Eigen::Quaterniond q = body.rotation();

        // TODO: get current linear and angular momentum of body
        Eigen::Vector3d p = body.linearMomentum();
        Eigen::Vector3d v = body.linearVelocity();

        // TODO: get force and torque that are currently applied to body
        Eigen::Vector3d f = body.force();
        Eigen::Vector3d t = body.torque();

        // TODO: compute new linear momentum
        Eigen::Vector3d pnew = p + stepSize * f;

        // TODO: convert from linear momentum to linear velocity and update the body accordingly
        Eigen::Vector3d vnew = body.massInverse() * pnew;
        body.setLinearVelocity(vnew);

        // TODO: compute new angular momentum
        Eigen::Matrix3d I    = body.inertiaWorld();
        Eigen::Vector3d l    = body.angularMomentum();
        Eigen::Vector3d lnew = l + stepSize * (t - v.cross(I * v));

        // TODO: convert from angular momentum to angular velocity and update the body accordingly
        Eigen::Vector3d wnew = body.inertiaWorldInverse() * lnew;
        body.setAngularVelocity(wnew);


        // TODO: update position of the body using the linear velocity and update body accordingly
        Eigen::Vector3d xnew = x + stepSize * vnew;
        body.setPosition(xnew);

        // TODO: quaternion-based angular velocity update of rotation and update body accordingly
        Eigen::Quaterniond wq(0, wnew.x(), wnew.y(), wnew.z());
        Eigen::Quaterniond qnew = (q + 0.5 * stepSize * wq * q).normalized();
        body.setRotation(qnew);
    }

    void implicitEuler(RigidBody& body, double stepSize)
    {
        // See for explanations: https://www.gdcvault.com/play/1022196/Physics-for-Game-Programmers-Numerical

        // TODO: get current position and rotation of the body
        Eigen::Vector3d x    = body.position();
        Eigen::Quaterniond q = body.rotation();

        // TODO: get current linear and angular momentum of body
        Eigen::Vector3d p = body.linearMomentum();
        Eigen::Vector3d v = body.linearVelocity();

        // TODO: get force and torque that are currently applied to body
        Eigen::Vector3d f = body.force();
        Eigen::Vector3d t = body.torque();

        // TODO: compute new linear momentum
        Eigen::Vector3d pnew = p + stepSize * f;


        // TODO: convert from linear momentum to linear velocity and update the body accordingly
        Eigen::Vector3d vnew = body.massInverse() * pnew;
        body.setLinearVelocity(vnew);

        // TODO: Convert current angular velocity to body coordinates (initial guess wb0)
        Eigen::Vector3d wb0 = body.inertiaWorldInverse() * body.angularMomentum();
        Eigen::Vector3d wb  = wb0;
        Eigen::Vector3d wbWorld = body.rotation().inverse() * wb;

        // TODO: Compute residual vector f(wb0) from the from angular velocity in body-coordinates
        Eigen::Vector3d fwb = wb.cross(body.inertiaWorld() * wb) - t;
        Eigen::Vector3d residual = fwb;
        Eigen::Matrix3d Ib = body.inertiaWorldInverse();
        Eigen::Matrix3d skew_wb = skew(wbWorld);
        Eigen::Matrix3d skew_Ib_wb = Ib * skew_wb;
        Eigen::Vector3d fwb0 = residual + skew_Ib_wb * wbWorld;
        Eigen::Vector3d delta_wb;
        Eigen::Matrix3d Jwb = Ib * skew_wb;
        Eigen::Matrix3d JwbT = Jwb.transpose();
        Eigen::Matrix3d J = JwbT * Jwb;
        Eigen::Vector3d fwb0T = fwb0.transpose();
        Eigen::Vector3d delta_wb0;
        Eigen::Vector3d delta_wbWorld;
        Eigen::Vector3d delta_wbWorldT;
        Eigen::Vector3d delta_wbWorld0;
        Eigen::Vector3d delta_wb0T;
        Eigen::Vector3d delta_wbWorld0T;
        Eigen::Vector3d delta_wbWorld0T2;
        Eigen::Vector3d delta_wbWorld0T3; 

        // TODO: Compute the Jacobian of f at wb. You can use the function "skew".
        Eigen::Matrix3d J = Ib * skew(wbWorld);
        Eigen::Matrix3d Jt = J.transpose();
        Eigen::Matrix3d JtJ = Jt * J;
        Eigen::Matrix3d JtJInv = JtJ.inverse();

        // TODO: Linearly solve for the update step delta_wb, for example using a QR decomposition
        Eigen::Vector3d delta_wb = J.colPivHouseholderQr().solve(-residual);
        
        // TODO: Apply the Newton-Raphson iteration by adding delta_wb to the current angular velocity
        wb += delta_wb;

        // TODO: Transform the angular velocity back to world coordinates
        Eigen::Vector3d wnew = body.rotation() * wb;
        body.setAngularVelocity(wnew);

        // TODO: explicitly integrate the torque and update the body accordingly
        Eigen::Vector3d lnew = body.inertiaWorld() * wnew;
        body.setAngularMomentum(lnew);

        // TODO: update position of the body using the linear velocity and update body accordingly
        Eigen::Vector3d xnew = x + stepSize * vnew;
        body.setPosition(xnew);

        // TODO: quaternion-based angular velocity update of rotation
        Eigen::Quaterniond wq(0, wnew.x(), wnew.y(), wnew.z());
        Eigen::Quaterniond qnew = (q + 0.5 * stepSize * wq * q).normalized();

    }
}
