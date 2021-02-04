#include "CollisionDetection.h"
#include "lcp_solver.h"
#include "gjk.h"
#include <cmath>

const double THRESHOLD = 0.00001;


void CollisionDetection::computeBroadPhase(int broadPhaseMethod) {
    // compute possible collisions
    m_overlappingBodys.clear();

    switch (broadPhaseMethod) {
    case 0: { // none
        for (size_t i = 0; i < m_objects.size(); i++) {
            for (size_t j = i + 1; j < m_objects.size(); j++) {
                m_overlappingBodys.push_back(std::make_pair(i, j));
            }
        }
        break;
    }

    case 1: {  // AABB
        // compute bounding boxes
        std::vector<AABB> aabbs(m_objects.size());
        for (size_t i = 0; i < aabbs.size(); i++) {
            aabbs[i].computeAABB(m_objects[i]);
        }
        for (size_t i = 0; i < m_objects.size(); i++) {
            for (size_t j = i + 1; j < m_objects.size(); j++) {
                // add pair of objects to possible collision if their
                // bounding boxes overlap
                if (aabbs[i].testCollision(aabbs[j])) {
                    m_overlappingBodys.push_back(std::make_pair(i, j));
                }
            }
        }
        break;
    }

    }
}

void CollisionDetection::computeNarrowPhase(int narrowPhaseMethod) {
    // exhaustive search is from previous exercise, it has been adjusted slighlty (because there were a few bugs, but some still persist)
    // default is our GJK implementation
    switch (narrowPhaseMethod) {
    case 0: {
        // exhaustive
        // iterate through all pairs of possible collisions

        for (auto overlap : m_overlappingBodys) {
            std::vector<Contact> temp_contacts[2];


            // compute intersection of a with b first and intersection
            // of b with a and store results in temp_contacts
            for (int switcher = 0; switcher < 2; switcher++) {
                RigidObject* a =
                    &m_objects[(!switcher) ? overlap.first
                    : overlap.second];
                RigidObject* b =
                    &m_objects[(!switcher) ? overlap.second
                    : overlap.first];

                Eigen::MatrixXd Va, Vb;
                Eigen::MatrixXi Fa, Fb;
                a->getMesh(Va, Fa);
                b->getMesh(Vb, Fb);

                std::set<int> penetratingVertices;
                std::set<std::pair<int, int>> penetratingEdges;

                // iterate through all faces of first object

                for (int face = 0; face < Fa.rows(); face++) {
                    // iterate through all edges of given face                    
                    for (size_t j = 0; j < 3; j++) {

                        int start = Fa(face, j);
                        int end = Fa(face, (j + 1) % 3);

                        // check if there is a collision
                        ContactType ct = isColliding(
                            Va.row(start), Va.row(end), Vb, Fb);


                        // find collision and check for duplicates
                        switch (ct) {
                        case ContactType::VERTEXFACE: {

                            auto ret = penetratingVertices.insert(
                                Fa(face, j));
                            // if not already in set
                            if (ret.second) {
                                Contact temp_col =
                                    findVertexFaceCollision(
                                        Va.row(Fa(face, j)), Vb,
                                        Fb);
                                temp_col.a = a;
                                temp_col.b = b;
                                temp_col.type =
                                    ContactType::VERTEXFACE;

                                temp_contacts[switcher].push_back(
                                    temp_col);
                            }
                            break;
                        }
                        case ContactType::EDGEEDGE: {
                            int orderedStart = std::min(start, end);
                            int orderedEnd = std::max(start, end);
                            auto ret = penetratingEdges.insert(
                                std::make_pair(orderedStart,
                                    orderedEnd));

                            // if not already in set
                            if (ret.second) {
                                Contact temp_col =
                                    findEdgeEdgeCollision(
                                        Va.row(orderedStart),
                                        Va.row(orderedEnd), Vb, Fb);
                                temp_col.a = a;
                                temp_col.b = b;
                                temp_col.type =
                                    ContactType::EDGEEDGE;
                                temp_contacts[switcher].push_back(
                                    temp_col);
                            }
                            break;
                        }
                        case ContactType::NONE: {
                            break;
                        }
                        }
                    }
                }
            }


            // look for vertexFace
            bool found = false;
            for (int i = 0; i < 2; i++) {
                for (auto cont : temp_contacts[i]) {
                    if (cont.type == ContactType::VERTEXFACE) {
                        m_contacts.push_back(cont);
                        found = true;
                    }
                }
                if (found) {
                    break;
                }
            }
            if (found) {
                continue;
            }

            // take single edgeedge if possible

            if (temp_contacts[0].size() > 0 &&
                temp_contacts[0].size() < temp_contacts[1].size()) {
                for (int i = 0; i < temp_contacts[0].size(); i++) {
                    m_contacts.push_back(temp_contacts[0][i]);
                }
            }
            else if (temp_contacts[1].size() > 0 &&
                temp_contacts[0].size() >
                temp_contacts[1].size()) {
                for (int i = 0; i < temp_contacts[1].size(); i++) {
                    m_contacts.push_back(temp_contacts[1][i]);
                }
            }
            else if (temp_contacts[0].size() > 0) {
                for (int i = 0; i < temp_contacts[0].size(); i++) {
                    m_contacts.push_back(temp_contacts[0][i]);
                }
            }
            else if (temp_contacts[1].size() > 0) {
                for (int i = 0; i < temp_contacts[1].size(); i++) {
                    m_contacts.push_back(temp_contacts[1][i]);
                }
            }

        }
        break;
    }

    case 1: {
        // this is used for the simulation
        // GJK

        for (auto overlap : m_overlappingBodys) {
            RigidObject* A =
                &m_objects[overlap.first];
            RigidObject* B =
                &m_objects[overlap.second];


            Eigen::MatrixXd Va, Vb;
            Eigen::MatrixXi Fa, Fb;

            A->getMesh(Va, Fa);

            B->getMesh(Vb, Fb);


            GJK::run(Va, Vb, m_contacts, A, B);

        }
        break;
    }

    }

}

void CollisionDetection::applyImpulse(double eps) {
    // compute impulse for all contacts
    // apply iteratively
    // based on D. Baraff's lecture notes "An Introduction to Physically Based Modeling"

    for (auto contact : m_contacts) {
        Eigen::Vector3d vrel_vec = contact.a->getVelocity(contact.p) -
            contact.b->getVelocity(contact.p);
        double vrel = contact.n.dot(vrel_vec);


        if (vrel > THRESHOLD) {
            // bodies are moving apart
            continue;
        }

        else if (abs(vrel) <= THRESHOLD) {
            // resting contact
            continue;
        }

        else {
            // compute impulse and update momenta

            // compute each of the 4 terms of the denominator (c.f. Baraff lecture notes)
            Eigen::Vector3d numerator = -(1 + eps) * vrel_vec;
            double t1 = contact.a->getMassInv();
            double t2 = contact.b->getMassInv();
            Eigen::Vector3d ra = contact.p - contact.a->getPosition();
            Eigen::Vector3d rb = contact.p - contact.b->getPosition();
            double t3 = contact.n.dot(
                (contact.a->getInertiaInvWorld() * (ra.cross(contact.n)))
                .cross(ra));
            double t4 = contact.n.dot(
                (contact.b->getInertiaInvWorld() * (rb.cross(contact.n)))
                .cross(rb));


            // compute the impulse magnitude
            Eigen::Vector3d j = numerator / (t1 + t2 + t3 + t4);
            Eigen::Vector3d force = j.dot(contact.n) * contact.n;

            // apply the computed impulse to the bodies
            contact.a->setLinearMomentum(contact.a->getLinearMomentum() +
                force);
            contact.b->setLinearMomentum(contact.b->getLinearMomentum() -
                force);
            contact.a->setAngularMomentum(contact.a->getAngularMomentum() +
                ra.cross(force));
            contact.b->setAngularMomentum(contact.b->getAngularMomentum() -
                rb.cross(force));
        }
    }

}


Eigen::Vector3d CollisionDetection::compute_ndot(Contact* c)
{
    // based on D. Baraff's lecture notes "An Introduction to Physically Based Modeling"

    // returns the derivative of the normal vector

    if (c->type == ContactType::VERTEXFACE) // vertex-face contact
    {
        // in this case, the normal vector is attached to B
        return (c->b->getAngularVelocity()).cross(c->n);
    }
    else // edge-edge contact
    {
        // differentiating n with resepect to time
        Eigen::Vector3d eadot = c->a->getAngularVelocity().cross(c->ea);
        Eigen::Vector3d ebdot = c->b->getAngularVelocity().cross(c->eb);

        Eigen::Vector3d n1 = c->ea.cross(c->eb);
        Eigen::Vector3d z = eadot.cross(c->eb) + c->ea.cross(ebdot);
        double l = n1.norm();
        // normalize
        n1 = n1.normalized();
        return (z - ((z.dot(n1)) * (n1))) / l;

    }
}

void CollisionDetection::compute_b_vector(Contact contacts[], int ncontacts, std::vector<double>& b)
{
    // based on D. Baraff's lecture notes "An Introduction to Physically Based Modeling"
    // b vector required for defining LCP problem

    for (int i = 0; i < ncontacts; i++)
    {
        // get information about contact and bodies invovlved
        Contact* c = &contacts[i];
        RigidObject* A = c->a, * B = c->b;
        Eigen::Vector3d n = c->n,
            ra = c->p - A->getPosition(),
            rb = c->p - B->getPosition();

        // get external forces and torques
        Eigen::Vector3d f_ext_a = A->getForce();
        Eigen::Vector3d f_ext_b = B->getForce();
        Eigen::Vector3d	t_ext_a = A->getTorque();
        Eigen::Vector3d	t_ext_b = B->getTorque();
        Eigen::Vector3d a_ext_part, a_vel_part,
            b_ext_part, b_vel_part;

        // compute part of second derivative of p_a (position of body A) due to external force and torque
        // analogously for body B
        a_ext_part = f_ext_a / A->getMass() +
            ((A->getInertiaInvWorld() * t_ext_a).cross(ra)),
            b_ext_part = f_ext_b / B->getMass() +
            ((B->getInertiaInvWorld() * t_ext_b).cross(rb));


        // compute part of second derivative of p_a (position of body A) due to velocity
        // analogously for body B
        a_vel_part = (A->getAngularVelocity().cross((A->getAngularVelocity().cross(ra)))) + ((A->getInertiaInvWorld() * (A->getAngularMomentum().cross(A->getAngularVelocity()))).cross(ra));
        b_vel_part = (B->getAngularVelocity().cross((B->getAngularVelocity().cross(rb)))) + ((B->getInertiaInvWorld() * (B->getAngularMomentum().cross(B->getAngularVelocity()))).cross(rb));


        // combine above results and dot product with normal vector
        double k1 = n.dot((a_ext_part + a_vel_part) - (b_ext_part + b_vel_part));

        Eigen::Vector3d ndot = compute_ndot(c);

        double k2 = 2 * ndot.dot(A->getVelocity(c->p) - B->getVelocity(c->p));

        b[i] = abs(k1 + k2) <= std::numeric_limits<float>::epsilon() ? 0 : k1 + k2;
    }
}


void CollisionDetection::compute_a_matrix(Contact contacts[], int ncontacts, Eigen::MatrixXd& a)
{
    // based on D. Baraff's lecture notes "An Introduction to Physically Based Modeling"

    // A matrix is computed based on resting contacts information
    // required for defining LCP problem

    for (int i = 0; i < ncontacts; i++) {
        for (int j = 0; j < ncontacts; j++) {
            a(i, j) = compute_aij(contacts[i], contacts[j]);
        }
    }
}


double CollisionDetection::compute_aij(Contact ci, Contact cj)
{
    // based on D. Baraff's lecture notes "An Introduction to Physically Based Modeling"
    // there were minor bugs in the original implementation which have been fixed

    // If the bodies involved in the ith and jth contact are distinct, then aij is zero

    if ((ci.a->getID() != cj.a->getID()) && (ci.b->getID() != cj.b->getID()) &&
        (ci.a->getID() != cj.b->getID()) && (ci.b->getID() != cj.a->getID()))
        return 0.0;
    RigidObject* A = ci.a,
        * B = ci.b;
    Eigen::Vector3d ni = ci.n, // contact normal of contact i
        nj = cj.n, // contact normal of contact j
        pi = ci.p, // ith contact point location
        pj = cj.p, // jth contact point location
        ra = pi - A->getPosition(),
        rb = pi - B->getPosition();

    // force and torque that j exerts on body A */
    Eigen::Vector3d force_on_a = Eigen::Vector3d::Zero(),
        torque_on_a = Eigen::Vector3d::Zero();
    if (cj.a->getID() == ci.a->getID())
    {
        // force direction of jth contact force on A
        force_on_a = nj;
        // torque direction
        torque_on_a = (pj - A->getPosition()).cross(force_on_a);
    }
    else if (cj.b->getID() == ci.a->getID())
    {
        force_on_a = -nj;
        torque_on_a = (pj - A->getPosition()).cross(force_on_a);
    }
    // force and torque that j exerts on body B
    Eigen::Vector3d force_on_b = Eigen::Vector3d::Zero(),
        torque_on_b = Eigen::Vector3d::Zero();

    if (cj.a->getID() == ci.b->getID())
    {
        // force direction of jth contact force on B
        force_on_b = nj;
        // torque direction
        torque_on_b = (pj - B->getPosition()).cross(force_on_b);
    }
    else if (cj.b->getID() == ci.b->getID())
    {
        force_on_b = -nj;
        torque_on_b = (pj - B->getPosition()).cross(force_on_b);
    }

    Eigen::Vector3d a_linear;
    Eigen::Vector3d b_linear;
    Eigen::Vector3d a_angular;
    Eigen::Vector3d b_angular;

    // compute how the jth contact force affects the linear and angular acceleration of the contact point on body A
    a_linear = force_on_a / A->getMass();
    a_angular = (A->getInertiaInvWorld() * torque_on_a).cross(ra);

    // analogously for body B
    b_linear = force_on_b / B->getMass();
    b_angular = (B->getInertiaInvWorld() * torque_on_b).cross(rb);

    return ni.dot((a_linear + a_angular) - (b_linear + b_angular));
}


void CollisionDetection::compute_contact_forces(Contact contacts[], int ncontacts) {
    // based on D. Baraff's lecture notes "An Introduction to Physically Based Modeling"

    if (ncontacts > 0) {

        // extract resting contacts from the set of all contacts
        std::vector<Contact> resting_contacts;
        std::vector<Contact> all_contacts;

        for (int i = 0; i < ncontacts; i++) {
            Eigen::Vector3d vrel_vec = contacts[i].a->getVelocity(contacts[i].p) -
                contacts[i].b->getVelocity(contacts[i].p);
            double vrel = contacts[i].n.dot(vrel_vec);

            if (abs(vrel) <= THRESHOLD) {
                // resting contact
                resting_contacts.push_back(contacts[i]);
            }
        }

        if (resting_contacts.size() > 0) {

            int tot_contacts = resting_contacts.size();

            Eigen::MatrixXd amat = Eigen::MatrixXd(tot_contacts, tot_contacts);

            std::vector<double> bvec = std::vector<double>(resting_contacts.size());

            Contact* resting_contacts_arr = &resting_contacts[0];


            // compute A matrix for resting contacts
            compute_a_matrix(resting_contacts_arr, tot_contacts, amat);

            // compute b_vec for resting contacts
            if (!resting_contacts.empty()) {
                int ncontacts_resting_size = resting_contacts.size();
                compute_b_vector(resting_contacts_arr, ncontacts_resting_size, bvec);
            }

            // use Dantzig's algorithm in order to solve LCP problem, i.e. solve for the resting contact forces using A and b which define the LCP
            LCPSolver lcpSolver = LCPSolver(tot_contacts, amat, bvec);

            std::vector<double> fvec = std::vector<double>(tot_contacts);
            fvec = lcpSolver.getForces();


            // add the resting contact forces we just computed into the force and torque field of each rigid body

            for (int i = 0; i < tot_contacts; i++) {
                double f_ = fvec[i]; // retrieve computed contact force

                Eigen::Vector3d n = contacts[i].n;
                RigidObject* A = contacts[i].a,
                    * B = contacts[i].b;

                // apply the force positively to A
                A->setForce(A->getForce() + f_ * n);
                A->setTorque(A->getTorque() + (contacts[i].p - A->getPosition()).cross(f_ * n));

                // apply the force negatively to B
                B->setForce(B->getForce() - f_ * n);
                B->setTorque(B->getTorque() - (contacts[i].p - B->getPosition()).cross(f_ * n));

            }
        }
    }
}


