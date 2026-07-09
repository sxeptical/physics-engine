#pragma once

#include "phys/body.hpp"
#include "phys/collision.hpp"
#include "phys/joint.hpp"
#include "phys/raycast.hpp"

#include <functional>
#include <optional>
#include <vector>

namespace phys {

struct WorldSettings {
    Vec3 gravity{0.f, -9.81f, 0.f};
    int velocity_iterations = 10;
    int position_iterations = 4;
    float linear_slop = 0.005f;
    float baumgarte = 0.2f;
    float max_linear_correction = 0.2f;
    float allowed_penetration = 0.01f;
    float restitution_threshold = 1.f;
    bool enable_gyroscopic = false;

    // Sleeping
    bool allow_sleep = true;
    float sleep_linear_threshold = 0.08f;
    float sleep_angular_threshold = 0.1f;
    float sleep_time = 0.5f;

    // Spatial hash broadphase
    float broadphase_cell_size = 2.f;
};

class World {
public:
    explicit World(WorldSettings settings = {});

    BodyId create_body(const BodyDef& def);
    void destroy_body(BodyId id);

    Body* get_body(BodyId id);
    const Body* get_body(BodyId id) const;

    JointId create_distance_joint(const DistanceJointDef& def);
    JointId create_ball_socket_joint(const BallSocketJointDef& def);
    JointId create_mouse_joint(const MouseJointDef& def);
    void destroy_joint(JointId id);
    Joint* get_joint(JointId id);
    void set_mouse_joint_target(JointId id, const Vec3& target);

    [[nodiscard]] const std::vector<Body>& bodies() const { return bodies_; }
    [[nodiscard]] const std::vector<Manifold>& contacts() const { return contacts_; }
    [[nodiscard]] const std::vector<Joint>& joints() const { return joints_; }
    [[nodiscard]] const WorldSettings& settings() const { return settings_; }
    WorldSettings& settings() { return settings_; }

    void step(float dt);
    void clear();

    /// Closest ray hit against all bodies. `skip_static` ignores static bodies.
    [[nodiscard]] std::optional<RaycastHit> raycast(const Ray& ray, float max_distance = 1e6f,
                                                    bool skip_static = false) const;

    using ContactCallback = std::function<void(const Manifold&)>;
    void set_contact_callback(ContactCallback cb) { contact_callback_ = std::move(cb); }

private:
    void integrate_forces(float dt);
    void integrate_velocities(float dt);
    void broadphase_and_narrowphase();
    void solve_velocity_constraints();
    void solve_position_constraints();
    void solve_joints(float dt);
    void solve_mouse_joints(float dt);
    void solve_constraint_joints(float dt);
    void warm_start();
    void update_sleeping(float dt);
    void wake_body(BodyId id);
    [[nodiscard]] bool bodies_jointed(BodyId a, BodyId b) const;

    WorldSettings settings_;
    std::vector<Body> bodies_;
    std::vector<bool> body_alive_;
    std::vector<BodyId> free_ids_;
    std::vector<Manifold> contacts_;
    std::vector<Manifold> previous_contacts_;
    std::vector<Joint> joints_;
    std::vector<bool> joint_alive_;
    std::vector<JointId> free_joint_ids_;
    ContactCallback contact_callback_;
    BodyId next_id_ = 0;
    JointId next_joint_id_ = 0;
};

}  // namespace phys
