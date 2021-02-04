#include "RigidBodySim.h"

#include <igl/writeOFF.h>

bool ball_collision = false;

Eigen::Matrix3d skew(const Eigen::Vector3d& a) {
    Eigen::Matrix3d s;
    s << 0, -a.z(), a.y(), a.z(), 0, -a.x(), -a.y(), a.x(), 0;
    return s;
}

bool RigidBodySim::advance() {

    // record: store each mesh in the scene in a separate .obj file in each step of the simulation
    /*
    string path = "../recording/";

    for (int i = 0; i < m_objects.size(); i++) {
        std::string filename =
            path + "obj_" + to_string(i) + "step_" + to_string(m_step) + ".obj";
        auto mesh = m_objects[i].getMesherino();
        bool succ = igl::writeOBJ(filename, mesh.V, mesh.F);
        if (!succ) {
            std::cerr << "Failed to write recording" << std::endl;
        }
    }
    */

    // compute collisions (collision detection)
    m_collisionDetection.computeCollisionDetection(m_broadPhaseMethod, m_narrowPhaseMethod, m_eps);

    // hack: when ball is close enough to wall, turn gravity on for wall fractures
    if (m_objects[m_objects.size() - 1].getPosition()[2] <= 0.65) {
        ball_collision = true;
    }

    // apply gravity to ball
    m_objects[m_objects.size() - 1].applyForceToCOM(m_gravity);

    // apply external force (gravity) to all objects, as soon as ball is about to collide with fall
    if (ball_collision) {
        for (auto& o : m_objects) {
            o.applyForceToCOM(m_gravity);
        }
    }

    // compute resting contact forces
    if (m_collisionDetection.getContacts().size() > 0) {
        m_collisionDetection.compute_contact_forces(&m_collisionDetection.getContacts()[0], m_collisionDetection.getContacts().size());
    }

    for (auto& o : m_objects) {
        // integrate velocities (symplectic euler)
        o.setLinearMomentum(o.getLinearMomentum() + m_dt * o.getForce());
        o.setAngularMomentum(o.getAngularMomentum() + m_dt * o.getTorque());
        o.resetForce();
        o.resetTorque();

        // integrate position (symplectic euler)
        o.setPosition(o.getPosition() + m_dt * o.getLinearVelocity());

        // integrate rotation

        // get angular velocity
        Eigen::Vector3d w = o.getAngularVelocity();

        // update orientation
        switch (m_method) {
        case 0: {
            // matrix-based method was not used for the simulation, default is quaternion-based

            // this code is taken directly from one of the previous exercises 

            // matrix-based angular velocity
            Eigen::Matrix3d r = o.getRotationMatrix();
            Eigen::Matrix3d W;

            // skew-matrix (row-wise)
            W << 0, -w.z(), w.y(),
                w.z(), 0, -w.x(),
                -w.y(), w.x(), 0;

            r = r + m_dt * W * r;

            // orthogonalize rotation matrix to show issue
            // https://math.stackexchange.com/questions/3292034/normalizing-a-rotation-matrix
            // https://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem
            // idea is to find the nearest orthogonral matrix by SVD
            Eigen::JacobiSVD<Eigen::Matrix3d> svd(r, Eigen::ComputeFullU | Eigen::ComputeFullV);
            r = svd.matrixU() * svd.matrixV().transpose();

            o.setRotation(r);
            break;
        }
        default: {
            // quaternion-based
            // solve gyroscopic force

            //Define quaternion, inertia tensor and angular velocity in body coordinates
            Eigen::Quaterniond q = o.getRotation();
            Eigen::Matrix3d I = o.getInertia();
            Eigen::Vector3d wb = q.inverse() * w;

            //Compute f according to (7)
            Eigen::Vector3d f = m_dt * wb.cross(I * wb);

            //Compute J according to (9)
            Eigen::Matrix3d J = I + m_dt * (skew(wb) * I - skew(I * wb));

            //Solve according to (8)
            Eigen::Vector3d delta_wb = J.colPivHouseholderQr().solve(-f);

            //Update and transform back to world coordinates according to (10)
            wb = wb + delta_wb;

            w = q * wb;

            o.setAngularVelocity(w);

            // update orientation

            Eigen::Quaterniond w_prime(0, w.x(), w.y(), w.z());

            Eigen::Quaterniond w_prime_q = w_prime * q;

            Eigen::Quaterniond new_q(q.w() + 0.5 * m_dt * w_prime_q.w(), q.x() + 0.5 * m_dt * w_prime_q.x(), q.y() + 0.5 * m_dt * w_prime_q.y(), q.z() + 0.5 * m_dt * w_prime_q.z());

            o.setRotation(new_q.normalized());

            break;
        }
        }
    }

    // advance time
    m_time += m_dt;
    m_step++;

    return false;
}
