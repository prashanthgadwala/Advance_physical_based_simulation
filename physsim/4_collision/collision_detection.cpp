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

            std::vector<Eigen::AlignedBox3d> aabbs(mObjects.size());
            for (std::size_t i = 0; i < aabbs.size(); i++)
                aabbs[i] = mObjects[i]->worldBounds();

            std::vector<std::pair<std::size_t, std::size_t>> overlaps[3];

            for (int axis = 0; axis < 3; ++axis)
            {
            // TODO: sort intervals in ascending order by beginning of interval
                struct Interval {
                    double min, max;
                    std::size_t idx;
                };

                std::vector<Interval> intervals;

                for (std::size_t i = 0; i < aabbs.size(); i++)
                {
                    intervals.push_back({aabbs[i].min()[axis], aabbs[i].max()[axis], i});
                }
                std::sort(intervals.begin(), intervals.end(),
                        [](const Interval& a, const Interval& b) { return a.min < b.min; });

            // TODO: iterate and place overlaps in a set
                for (std::size_t i = 0; i < intervals.size(); ++i)
                {
                    for (std::size_t j = i + 1; j < intervals.size(); ++j)
                    {
                        if (intervals[j].min > intervals[i].max)
                            break; // No more overlaps possible
                        // Only consider dynamic bodies
                        if (mObjects[intervals[i].idx]->type() == RigidBody::EType::Dynamic ||
                            mObjects[intervals[j].idx]->type() == RigidBody::EType::Dynamic)
                        {
                            // Ensure consistent ordering of pairs (smaller index first)
                            std::size_t idx1 = intervals[i].idx;
                            std::size_t idx2 = intervals[j].idx;
                            if (idx1 > idx2) std::swap(idx1, idx2);
                            overlaps[axis].emplace_back(idx1, idx2);
                        }
                    }
                }
            }

            // TODO: grab elements that occurred in all containers for the narrow test
            std::set<std::pair<std::size_t, std::size_t>> set_x(overlaps[0].begin(), overlaps[0].end());
            std::set<std::pair<std::size_t, std::size_t>> set_y(overlaps[1].begin(), overlaps[1].end());
            std::set<std::pair<std::size_t, std::size_t>> set_z(overlaps[2].begin(), overlaps[2].end());

            std::vector<std::pair<std::size_t, std::size_t>> intersection_xy;
            std::set_intersection(set_x.begin(), set_x.end(),
                                set_y.begin(), set_y.end(),
                                std::back_inserter(intersection_xy));
            std::vector<std::pair<std::size_t, std::size_t>> intersection_xyz;
            std::set_intersection(intersection_xy.begin(), intersection_xy.end(),
                                set_z.begin(), set_z.end(),
                                std::back_inserter(intersection_xyz));

            // TODO: pass the intersections on to the narrow phase
            mOverlappingBodys = intersection_xyz;

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
            if (vrel <= 0)
            {
                // bodies are moving apart or stationary
                continue;
            }

            // contact points ra and rb translated to local frame
            Eigen::Vector3d ra = contact.p - contact.a->position();
            Eigen::Vector3d rb = contact.p - contact.b->position();

            // TODO: compute impulse response
            double invMassA = contact.a->massInverse();
            double invMassB = contact.b->massInverse();
            Eigen::Matrix3d invInertiaA = contact.a->inertiaWorldInverse();
            Eigen::Matrix3d invInertiaB = contact.b->inertiaWorldInverse();

            Eigen::Vector3d n = contact.n;

            // Terms for denominator
            Eigen::Vector3d ra_cross_n = ra.cross(n);
            Eigen::Vector3d rb_cross_n = rb.cross(n);

            double term1 = invMassA;
            double term2 = invMassB;
            double term3 = n.dot((invInertiaA * ra_cross_n).cross(ra));
            double term4 = n.dot((invInertiaB * rb_cross_n).cross(rb));

            double denom = term1 + term2 + term3 + term4;
            if (denom == 0.0)
                continue; // avoid division by zero

            double j = -(1.0 + eps) * vrel / denom;
            Eigen::Vector3d impulse = j * n;

            // TODO: apply impulse forces to the bodies at the contact point
            // For body A: apply -impulse
            contact.a->setLinearVelocity(contact.a->linearVelocity() - impulse * invMassA);
            contact.a->setAngularVelocity(contact.a->angularVelocity() - invInertiaA * ra_cross_n * j);
            
            // For body B: apply +impulse  
            contact.b->setLinearVelocity(contact.b->linearVelocity() + impulse * invMassB);
            contact.b->setAngularVelocity(contact.b->angularVelocity() + invInertiaB * rb_cross_n * j);


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
