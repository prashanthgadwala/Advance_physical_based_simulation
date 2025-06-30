#include "cloth.hpp"

#include <vislab/core/array.hpp>
#include <vislab/geometry/surfaces.hpp>
#include <vislab/geometry/vertex_normals.hpp>
#include <vislab/graphics/trimesh_geometry.hpp>

namespace physsim
{
    Cloth::Cloth(const Eigen::Vector2i& resolution, const ETopology& topology, const Eigen::Vector3f& origin, const Eigen::Vector3f& axis1, const Eigen::Vector3f& axis2)
        : mResolution(resolution)
        , mTopology(topology)
        , mOrigin(origin)
        , mAxis1(axis1)
        , mAxis2(axis2)
    {
        // allocate mesh
        mesh          = std::make_shared<vislab::TrimeshGeometry>();

        // create buffer to hold the spring positions
        mesh->positions = std::make_shared<vislab::Array3f>();
        mesh->positions->setSize(resolution.prod());

        // create buffer to hold the normals
        mesh->normals = std::make_shared<vislab::Array3f>();
        mesh->normals->setSize(resolution.prod());

        // create buffer to hold the spring velocities
        velocities = std::make_shared<vislab::Array3f>();
        velocities->setSize(resolution.prod());

        // create buffer to accumulate accelerations in
        accelerations = std::make_shared<vislab::Array3f>();
        accelerations->setSize(resolution.prod());

        // create index buffer
        mesh->indices = std::make_shared<vislab::Array3u>();
        mesh->indices->setSize((resolution.x() - 1) * (resolution.y() - 1) * 2);
        int index = 0;
        for (int iy = 0; iy < resolution.y() - 1; ++iy)
            for (int ix = 0; ix < resolution.x() - 1; ++ix)
            {
                int i00 = (iy + 0) * resolution.x() + (ix + 0);
                int i01 = (iy + 1) * resolution.x() + (ix + 0);
                int i10 = (iy + 0) * resolution.x() + (ix + 1);
                int i11 = (iy + 1) * resolution.x() + (ix + 1);
                mesh->indices->setValue(index++, Eigen::Vector3u(i00, i01, i10));
                mesh->indices->setValue(index++, Eigen::Vector3u(i01, i11, i10));
            }

        reset();
    }

    void Cloth::advance(double stepSize)
    {
        float Lx  = mAxis1.norm();
        float Ly  = mAxis2.norm();
        float Lxy = (mAxis1 + mAxis2).norm();
        accelerations->setZero();

#ifndef _DEBUG
#pragma omp parallel for
#endif
        for (int64_t i = 0; i < mResolution.prod(); ++i)
        {
            // early out if the start of the spring is fixed
            if (mFixed.find(i) != mFixed.end())
                continue;

            // vertex grid indices
            int ix = i % mResolution.x();
            int iy = i / mResolution.x();

            // add spring forces according to the spring topology
            if (mTopology == ETopology::Structural)
            {
                if (ix > 0)
                    addSpringForce(Eigen::Vector2i(ix - 1, iy), Eigen::Vector2i(ix, iy), Lx);
                if (ix < mResolution.x() - 1)
                    addSpringForce(Eigen::Vector2i(ix + 1, iy), Eigen::Vector2i(ix, iy), Lx);
                if (iy > 0)
                    addSpringForce(Eigen::Vector2i(ix, iy - 1), Eigen::Vector2i(ix, iy), Ly);
                if (iy < mResolution.y() - 1)
                    addSpringForce(Eigen::Vector2i(ix, iy + 1), Eigen::Vector2i(ix, iy), Ly);
            }
            else if (mTopology == ETopology::Diagonal)
            {
                if (ix > 0)
                    addSpringForce(Eigen::Vector2i(ix - 1, iy), Eigen::Vector2i(ix, iy), Lx);
                if (ix < mResolution.x() - 1)
                    addSpringForce(Eigen::Vector2i(ix + 1, iy), Eigen::Vector2i(ix, iy), Lx);
                if (iy > 0)
                    addSpringForce(Eigen::Vector2i(ix, iy - 1), Eigen::Vector2i(ix, iy), Ly);
                if (iy < mResolution.y() - 1)
                    addSpringForce(Eigen::Vector2i(ix, iy + 1), Eigen::Vector2i(ix, iy), Ly);
                if (ix > 0 && iy > 0)
                    addSpringForce(Eigen::Vector2i(ix - 1, iy - 1), Eigen::Vector2i(ix, iy), Lxy);
                if (ix > 0 && iy < mResolution.y() - 1)
                    addSpringForce(Eigen::Vector2i(ix - 1, iy + 1), Eigen::Vector2i(ix, iy), Lxy);
                if (ix < mResolution.x() - 1 && iy > 0)
                    addSpringForce(Eigen::Vector2i(ix + 1, iy - 1), Eigen::Vector2i(ix, iy), Lxy);
                if (ix < mResolution.x() - 1 && iy < mResolution.y() - 1)
                    addSpringForce(Eigen::Vector2i(ix + 1, iy + 1), Eigen::Vector2i(ix, iy), Lxy);
            }
        }

        // symplectic euler update
#ifndef _DEBUG
#pragma omp parallel for
#endif
        for (int64_t i = 0; i < mResolution.prod(); ++i)
        {
            // get position, velocity, and acceleration
            Eigen::Vector3f x = mesh->positions->getValue(i);
            Eigen::Vector3f v = velocities->getValue(i);
            Eigen::Vector3f a = accelerations->getValue(i);

            // TODO: perform symplectic Euler update
            Eigen::Vector3f vnew = v + stepSize * a;
            Eigen::Vector3f xnew = x + stepSize * vnew;

            // TODO: set the new position and the new velocity
            mesh->positions->setValue(i, xnew);
            velocities->setValue(i, vnew);
        }

        // recompute the vertex normals
        vislab::VertexNormals3f::computeNormals(mesh->positions, mesh->indices, mesh->normals);

        // notify that mesh positions have changed
        mesh->markChanged();
    }

    void Cloth::reset()
    {
        for (int iy = 0; iy < mResolution.y(); ++iy)
            for (int ix = 0; ix < mResolution.x(); ++ix)
            {
                mesh->positions->setValue(iy * mResolution.x() + ix, mOrigin + ix * mAxis1 + iy * mAxis2);
            }
        velocities->setZero();
    }

    const Eigen::Vector2i& Cloth::resolution() const
    {
        return mResolution;
    }

    void Cloth::pin(const Eigen::Vector2i& gridIndex)
    {
        mFixed.insert(gridIndex.y() * mResolution.x() + gridIndex.x());
    }

    void Cloth::unpin(const Eigen::Vector2i& gridIndex)
    {
        mFixed.erase(gridIndex.y() * mResolution.x() + gridIndex.x());
    }

    void Cloth::addSpringForce(const Eigen::Vector2i& startPoint, const Eigen::Vector2i& endPoint, const double& L)
    {
        // get start and end position of spring, as well as current velocity
        Eigen::Vector3f startPos = mesh->positions->getValue(startPoint.y() * mResolution.x() + startPoint.x());
        Eigen::Vector3f endPos   = mesh->positions->getValue(endPoint.y() * mResolution.x() + endPoint.x());
        Eigen::Vector3f endVel   = velocities->getValue(endPoint.y() * mResolution.x() + endPoint.x());
        Eigen::Vector3f endAcc   = accelerations->getValue(endPoint.y() * mResolution.x() + endPoint.x());

        // grab parameters
        Eigen::Vector3f spring_dir = (endPos - startPos).normalized();
        double spring_norm         = (endPos - startPos).norm();

        // TODO: compute force

        Eigen::Vector3f force = -stiffness * (spring_norm - L) * spring_dir;

        // TODO: add to acceleration
        Eigen::Vector3f totalForce = force - damping * endVel + 
                                     gravity * mass; // add damping and gravity

        accelerations->setValue(endPoint.y() * mResolution.x() + endPoint.x(), endAcc + totalForce / mass);

        // int startIdx = startPoint.y() * mResolution.x() + startPoint.x();
        // if (mFixed.find(startIdx) == mFixed.end()) {
        //     Eigen::Vector3f startAcc = accelerations->getValue(startIdx);
        //     accelerations->setValue(startIdx, startAcc - force / mass);
        // }
    }
}
