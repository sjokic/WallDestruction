#include "CollisionDetection.h"
#include "Simulation.h"

#include <deque>
#include <sys/stat.h>

using namespace std;

/*
 * Simulation that drops several cubes colliding each other.
 */
class RigidBodySim : public Simulation {
public:
    RigidBodySim() : Simulation(), m_collisionDetection(m_objects) { init(); }

    int NUM_CUBES = 79; // number of fractures (DO NOT MODIFY)

    bool IsPathExist(const std::string& s)
    {
        struct stat buffer;
        return (stat(s.c_str(), &buffer) == 0);
    }

    virtual void init() override {
        string path = "rectangle.off";

        // load first fracture
        m_objects.push_back(RigidObject("Cube_cell.obj"));

        int count_cubes = 0;

        // load remaining fractures
        for (int i = 1; i < NUM_CUBES; i++) {

            if (i < 10) {
                path = "Cube_cell.00" + std::to_string(i) + ".obj";
            }
            else {
                path = "Cube_cell.0" + std::to_string(i) + ".obj";
            }

            if (!IsPathExist("../data/" + path)) {
                std::cout << "COULDN'T FIND CUBE_CELL NUMBER: " << i << " SKIPPING!" << std::endl;
                continue;
            }

            count_cubes++;
            m_objects.push_back(RigidObject(path));

        }

        NUM_CUBES = count_cubes + 1;

        path = "cube_shallow.off";
        m_objects.push_back(RigidObject(path));

        path = "sphere.off";
        m_objects.push_back(RigidObject(path));

        ball = &m_objects[m_objects.size() - 1];

        m_collisionDetection.setObjects(m_objects);

        m_dt = 1e-3 * 3;
        m_gravity << 0, -9.81, 0;
        m_mass = 1.0;
        m_showContacts = false;
        m_broadPhaseMethod = 0;
        m_narrowPhaseMethod = 0;
        m_eps = 1.0;

        reset();
    }

    virtual void resetMembers() override {
        for (auto& o : m_objects) {
            o.reset();
        }

        ball->reset();
        ball->setScale(0.05);
        ball->setMass(15*m_mass);

        ball->setPosition(Eigen::Vector3d(-0.55, 2.5, 8.0));
        ball->setColors(Eigen::RowVector3d(1.0, 0, 0));
        updateVars();

        for (size_t i = 0; i < NUM_CUBES; i++) {
            m_objects[i].setMass(5*m_mass);
            m_objects[i].setPosition(m_objects[i].getPosition() + Eigen::Vector3d(0, 0, 1.0));
        }


        for (size_t i = 0; i < NUM_CUBES; i++) {
            m_objects[i].setMass(m_mass);
        }

        m_objects[NUM_CUBES].setScale(10);
        m_objects[NUM_CUBES].setType(ObjType::STATIC);
        m_objects[NUM_CUBES].setColors(Eigen::RowVector3d(0.5, 0.5, 0.5));
        m_objects[NUM_CUBES].setMass(std::numeric_limits<double>::max());

        m_objects[0 + NUM_CUBES].setPosition(
            Eigen::Vector3d(0, -0.1, 0));

    }

    virtual void updateRenderGeometry() override {
        for (size_t i = 0; i < m_objects.size(); i++) {
            RigidObject& o = m_objects[i];
            if (o.getID() < 0) {
                m_renderVs.emplace_back();
                m_renderFs.emplace_back();
            }

            m_objects[i].getMesh(m_renderVs[i], m_renderFs[i]);
        }
    }

    virtual bool advance() override;

    virtual void renderRenderGeometry(
        igl::opengl::glfw::Viewer& viewer) override {
        for (size_t i = 0; i < m_objects.size(); i++) {
            RigidObject& o = m_objects[i];
            if (o.getID() < 0) {
                int new_id = 0;
                if (i > 0) {
                    new_id = viewer.append_mesh();
                    o.setID(new_id);
                }
                else {
                    o.setID(new_id);
                }

                size_t meshIndex = viewer.mesh_index(o.getID());
                viewer.data_list[meshIndex].show_lines = false;

                viewer.data_list[meshIndex].set_face_based(true);
                viewer.data_list[meshIndex].point_size = 2.0f;
                viewer.data_list[meshIndex].clear();
            }
            size_t meshIndex = viewer.mesh_index(o.getID());

            viewer.data_list[meshIndex].set_mesh(m_renderVs[i], m_renderFs[i]);
            viewer.data_list[meshIndex].compute_normals();

            Eigen::MatrixXd color;
            o.getColors(color);
            viewer.data_list[meshIndex].set_colors(color);
        }

        if (m_showContacts) {
            // number of timesteps to keep showing collision
            int delay = 10;

            // clear old points
            viewer.data_list[1].points = Eigen::MatrixXd(0, 6);
            viewer.data_list[1].point_size = 10.0f;

            // remove expired points
            while (m_contactMemory.size() > 0 &&
                m_contactMemory.front().second + delay < m_step) {
                m_contactMemory.pop_front();
            }

            // get new points and add them to memory
            auto contacts = m_collisionDetection.getContacts();
            for (auto& contact : contacts) {
                m_contactMemory.push_back(std::make_pair(contact, m_step));
            }

            // show points
            for (auto& contact_int_p : m_contactMemory) {
                viewer.data_list[1].add_points(
                    contact_int_p.first.p.transpose(),
                    (contact_int_p.first.type == ContactType::EDGEEDGE)
                    ? Eigen::RowVector3d(0, 1, 0)
                    : Eigen::RowVector3d(0, 0, 1));
            }
        }
    }

#pragma region SettersAndGetters
    void setMethod(int m) { m_method = m; }
    /*
     * Compute magnitude and direction of momentum and apply it to o
     */
    void updateVars() {
        Eigen::Vector3d momentum;
        momentum << std::sin(m_angle), std::cos(m_angle), -6;
        momentum *= 6*m_force;
        ball->setLinearVelocity(momentum / ball->getMass());
    }

    void setAngle(double a) {
        m_angle = a;
        updateVars();
    }

    void setForce(double f) {
        m_force = f;
        updateVars();
    }

    void setMass(double m) { m_mass = m; }

    void showContacts(bool s) {
        if (!s) {
            m_contactMemory.clear();
        }
        m_showContacts = s;
    }

    void setBroadPhaseMethod(int m) { m_broadPhaseMethod = m; }
    void setNarrowPhaseMethod(int m) { m_narrowPhaseMethod = m; }

    void setEps(double eps) { m_eps = eps; }

    Eigen::Vector3d getKineticEnergy() {
        Eigen::Vector3d res;
        res.setZero();
        for (auto o : m_objects) {
            if (o.getType() == ObjType::STATIC) continue;
            Eigen::Vector3d rotE = 0.5 * o.getInertia().diagonal().cwiseProduct(
                o.getAngularVelocity());
            Eigen::Vector3d kinE =
                0.5 * o.getMass() * o.getLinearVelocity().array().square();
            res += rotE + kinE;
        }
        return res;
    }

    Eigen::Vector3d getLinearMomentum() {
        Eigen::Vector3d res;
        res.setZero();
        for (auto o : m_objects) {
            if (o.getType() == ObjType::STATIC) continue;
            res += o.getLinearMomentum();
        }
        return res;
    }

    Eigen::Vector3d getAngularMomentum() {
        Eigen::Vector3d res;
        res.setZero();
        for (auto o : m_objects) {
            if (o.getType() == ObjType::STATIC) continue;
            res += o.getAngularMomentum();
        }
        return res;
    }

#pragma endregion SettersAndGetters

private:
    int m_method;
    double m_angle;
    double m_force;
    double m_mass;

    Eigen::Vector3d m_gravity;

    CollisionDetection m_collisionDetection;
    int m_broadPhaseMethod;
    int m_narrowPhaseMethod;
    double m_eps;

    RigidObject* ball;

    std::vector<Eigen::MatrixXd> m_renderVs;  // vertex positions for rendering
    std::vector<Eigen::MatrixXi> m_renderFs;  // face indices for rendering

    bool m_showContacts;
    std::deque<std::pair<Contact, int>> m_contactMemory;
};