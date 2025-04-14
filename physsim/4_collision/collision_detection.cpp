#include "collision_detection.hpp"

#include "gilbert_johnson_keerthi.hpp"
#include "transformed_mesh.hpp"

#include <vislab/core/array.hpp>
#include <vislab/graphics/actor.hpp>
#include <vislab/graphics/trimesh_geometry.hpp>
#include <vislab/geometry/surfaces.hpp>

namespace physsim
{
    CollisionDetection::CollisionDetection(const std::vector<std::shared_ptr<RigidBody>>& objects)
        : mObjects(objects)
    {
    }

    void CollisionDetection::computeCollisionDetection(EBroadPhaseMethod broadPhaseMethod, double eps, double stepSize)
    {
        clearDataStructures();
        computeBroadPhase(broadPhaseMethod);
        computeNarrowPhase();
        applyImpulse(eps, stepSize);
    }

    void CollisionDetection::computeBroadPhase(EBroadPhaseMethod broadPhaseMethod)
    {
        // compute possible collisions
        mOverlappingBodys.clear();

        switch (broadPhaseMethod)
        {
        case EBroadPhaseMethod::None:
        {
            for (std::size_t i = 0; i < mObjects.size(); i++)
            {
                for (std::size_t j = i + 1; j < mObjects.size(); j++)
                {
                    if (mObjects[i]->type() == RigidBody::EType::Dynamic || mObjects[j]->type() == RigidBody::EType::Dynamic)
                        mOverlappingBodys.push_back(std::make_pair(i, j));
                }
            }
            break;
        }

        case EBroadPhaseMethod::AABB:
        {
            // compute bounding boxes
            std::vector<Eigen::AlignedBox3d> aabbs(mObjects.size());
            for (std::size_t i = 0; i < aabbs.size(); i++)
            {
                aabbs[i] = mObjects[i]->worldBounds();
            }
            for (std::size_t i = 0; i < mObjects.size(); i++)
            {
                for (std::size_t j = i + 1; j < mObjects.size(); j++)
                {
                    // add pair of objects to possible collision if their bounding boxes overlap
                    if (mObjects[i]->type() == RigidBody::EType::Dynamic || mObjects[j]->type() == RigidBody::EType::Dynamic)
                    {
                        if (aabbs[i].intersects(aabbs[j]))
                        {
                            mOverlappingBodys.push_back(std::make_pair(i, j));
                        }
                    }
                }
            }
            break;
        }

        case EBroadPhaseMethod::SweepAndPrune:
        {
            // TODO: compute bounding boxes and create intervals on the 3 main axes

            // TODO: sort intervals in ascending order by beginning of interval

            // TODO: iterate and place overlaps in a set

            // TODO: grab elements that occurred in all containers for the narrow test

            // TODO: pass the intersections on to the narrow phase
            break;
        }
        }
    }

    void CollisionDetection::computeNarrowPhase()
    {
        // iterate through all pairs of possible collisions
        for (auto overlap : mOverlappingBodys)
        {
            mPenetratingVertices.clear();
            mPenetratingEdges.clear();
            std::vector<Contact> temp_contacts[2];
            // compute intersection of a with b first and intersection
            // of b with a and store results in temp_contacts
            for (int switcher = 0; switcher < 2; switcher++)
            {
                RigidBody* a =
                    mObjects[(!switcher) ? overlap.first
                                         : overlap.second]
                        .get();
                RigidBody* b =
                    mObjects[(!switcher) ? overlap.second
                                         : overlap.first]
                        .get();

                auto a_trimesh = a->actor()->components.get<const vislab::TrimeshGeometry>();
                auto b_trimesh = b->actor()->components.get<const vislab::TrimeshGeometry>();

                auto a_transform = a->actor()->components.get<const vislab::Transform>();
                auto b_transform = b->actor()->components.get<const vislab::Transform>();

                TransformedMesh a_tmesh(
                    *a_transform,
                    *a_trimesh->positions.get(),
                    *a_trimesh->indices.get());

                TransformedMesh b_tmesh(
                    *b_transform,
                    *b_trimesh->positions.get(),
                    *b_trimesh->indices.get());

                Contact contact = GilbertJohnsonKeerthi::findCollisionGJK(a_tmesh, b_tmesh);
                if (contact.type != ContactType::None)
                {
                    contact.a = a;
                    contact.b = b;
                    temp_contacts[switcher].push_back(contact);
                }
            }
            // look for vertexFace
            bool found = false;
            for (int i = 0; i < 2; i++)
            {
                for (const auto& cont : temp_contacts[i])
                {
                    if (cont.type == ContactType::Vertex_Face)
                    {
                        mContacts.push_back(cont);
                        found = true;
                        break;
                    }
                }
                if (found)
                {
                    break;
                }
            }
            if (found)
            {
                continue;
            }

            // take single edgeedge if possible
            if (temp_contacts[0].size() > 0 &&
                temp_contacts[0].size() < temp_contacts[1].size())
            {
                for (int i = 0; i < temp_contacts[0].size(); i++)
                {
                    mContacts.push_back(temp_contacts[0][i]);
                }
            }
            else if (temp_contacts[1].size() > 0 &&
                     temp_contacts[0].size() >
                         temp_contacts[1].size())
            {
                for (int i = 0; i < temp_contacts[1].size(); i++)
                {
                    mContacts.push_back(temp_contacts[1][i]);
                }
            }
            else if (temp_contacts[0].size() > 0)
            {
                for (int i = 0; i < temp_contacts[0].size(); i++)
                {
                    mContacts.push_back(temp_contacts[0][i]);
                }
            }
            else if (temp_contacts[1].size() > 0)
            {
                for (int i = 0; i < temp_contacts[1].size(); i++)
                {
                    mContacts.push_back(temp_contacts[1][i]);
                }
            }
        }
    }

    void CollisionDetection::applyImpulse(double eps, double stepSize)
    {
        // see here fore more information: https://en.wikipedia.org/wiki/Collision_response
        if (mContacts.empty())
            return;

        // compute impulse for all contacts
        for (auto contact : mContacts)
        {
            // compute relative velocity and skip if the bodies are already moving apart
            Eigen::Vector3d vrel_vec = contact.b->velocity(contact.p) - contact.a->velocity(contact.p);
            double vrel              = contact.n.dot(vrel_vec);
            if (vrel < 0)
            {
                // bodies are moving apart
                continue;
            }

            // contact points ra and rb translated to local frame
            Eigen::Vector3d ra = contact.p - contact.a->position();
            Eigen::Vector3d rb = contact.p - contact.b->position();

            // TODO: compute impulse response

            // TODO: apply impulse forces to the bodies at the contact point
        }
    }

    void CollisionDetection::clearDataStructures()
    {
        mPenetratingEdges.clear();
        mPenetratingVertices.clear();
        mOverlappingBodys.clear();
        mContacts.clear();
    }
}
