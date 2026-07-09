#include "phys/world.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace phys {
namespace {

struct CellKey {
    int x, y, z;
    bool operator==(const CellKey& o) const { return x == o.x && y == o.y && z == o.z; }
};

struct CellKeyHash {
    std::size_t operator()(const CellKey& k) const noexcept {
        // 3D spatial hash
        const std::size_t h1 = static_cast<std::size_t>(k.x) * 73856093u;
        const std::size_t h2 = static_cast<std::size_t>(k.y) * 19349663u;
        const std::size_t h3 = static_cast<std::size_t>(k.z) * 83492791u;
        return h1 ^ h2 ^ h3;
    }
};

float effective_mass(const Body& a, const Body& b, const Vec3& ra, const Vec3& rb, const Vec3& n) {
    const Vec3 ra_n = ra.cross(n);
    const Vec3 rb_n = rb.cross(n);
    return a.inv_mass + b.inv_mass + ra_n.dot(a.inv_inertia_world * ra_n) +
           rb_n.dot(b.inv_inertia_world * rb_n);
}

}  // namespace

World::World(WorldSettings settings) : settings_(settings) {}

BodyId World::create_body(const BodyDef& def) {
    BodyId id;
    if (!free_ids_.empty()) {
        id = free_ids_.back();
        free_ids_.pop_back();
        body_alive_[id] = true;
        bodies_[id] = Body{};
    } else {
        id = next_id_++;
        bodies_.push_back(Body{});
        body_alive_.push_back(true);
    }

    Body& body = bodies_[id];
    body.id = id;
    body.type = def.type;
    body.shape = def.shape;
    body.position = def.position;
    body.orientation = def.orientation.normalized();
    body.linear_velocity = def.linear_velocity;
    body.angular_velocity = def.angular_velocity;
    body.restitution = def.restitution;
    body.friction = def.friction;
    body.fixed_rotation = def.fixed_rotation;
    body.allow_sleep = def.allow_sleep;
    body.awake = true;
    body.sleep_timer = 0.f;
    body.set_mass_from_density(def.density);
    return id;
}

void World::destroy_body(BodyId id) {
    if (id >= bodies_.size() || !body_alive_[id]) {
        return;
    }
    body_alive_[id] = false;
    free_ids_.push_back(id);

    // Remove joints attached to this body
    for (std::size_t j = 0; j < joints_.size(); ++j) {
        if (!joint_alive_[j]) {
            continue;
        }
        if (joints_[j].body_a == id || joints_[j].body_b == id) {
            destroy_joint(static_cast<JointId>(j));
        }
    }
}

Body* World::get_body(BodyId id) {
    if (id >= bodies_.size() || !body_alive_[id]) {
        return nullptr;
    }
    return &bodies_[id];
}

const Body* World::get_body(BodyId id) const {
    if (id >= bodies_.size() || !body_alive_[id]) {
        return nullptr;
    }
    return &bodies_[id];
}

JointId World::create_distance_joint(const DistanceJointDef& def) {
    Body* a = get_body(def.body_a);
    Body* b = get_body(def.body_b);
    if (!a || !b) {
        return kInvalidJoint;
    }

    JointId id;
    if (!free_joint_ids_.empty()) {
        id = free_joint_ids_.back();
        free_joint_ids_.pop_back();
        joint_alive_[id] = true;
        joints_[id] = Joint{};
    } else {
        id = next_joint_id_++;
        joints_.push_back(Joint{});
        joint_alive_.push_back(true);
    }

    Joint& j = joints_[id];
    j.id = id;
    j.type = JointType::Distance;
    j.body_a = def.body_a;
    j.body_b = def.body_b;
    j.local_anchor_a = def.local_anchor_a;
    j.local_anchor_b = def.local_anchor_b;
    j.stiffness = def.stiffness;
    j.damping = def.damping;
    j.max_force = 1e9f;
    j.alive = true;

    const Vec3 wa = a->local_to_world(def.local_anchor_a);
    const Vec3 wb = b->local_to_world(def.local_anchor_b);
    j.length = def.length >= 0.f ? def.length : (wb - wa).length();

    a->wake();
    b->wake();
    return id;
}

JointId World::create_ball_socket_joint(const BallSocketJointDef& def) {
    Body* a = get_body(def.body_a);
    Body* b = get_body(def.body_b);
    if (!a || !b) {
        return kInvalidJoint;
    }

    JointId id;
    if (!free_joint_ids_.empty()) {
        id = free_joint_ids_.back();
        free_joint_ids_.pop_back();
        joint_alive_[id] = true;
        joints_[id] = Joint{};
    } else {
        id = next_joint_id_++;
        joints_.push_back(Joint{});
        joint_alive_.push_back(true);
    }

    Joint& j = joints_[id];
    j.id = id;
    j.type = JointType::BallSocket;
    j.body_a = def.body_a;
    j.body_b = def.body_b;
    j.local_anchor_a = def.local_anchor_a;
    j.local_anchor_b = def.local_anchor_b;
    j.stiffness = def.stiffness;
    j.damping = def.damping;
    j.length = 0.f;
    j.max_force = 1e9f;
    j.alive = true;

    a->wake();
    b->wake();
    return id;
}

JointId World::create_mouse_joint(const MouseJointDef& def) {
    Body* body = get_body(def.body);
    if (!body || body->type != BodyType::Dynamic) {
        return kInvalidJoint;
    }

    JointId id;
    if (!free_joint_ids_.empty()) {
        id = free_joint_ids_.back();
        free_joint_ids_.pop_back();
        joint_alive_[id] = true;
        joints_[id] = Joint{};
    } else {
        id = next_joint_id_++;
        joints_.push_back(Joint{});
        joint_alive_.push_back(true);
    }

    Joint& j = joints_[id];
    j.id = id;
    j.type = JointType::Mouse;
    j.body_a = def.body;
    j.body_b = kInvalidBody;
    j.local_anchor_a = def.local_anchor;
    j.target = def.target;
    j.frequency_hz = std::max(0.1f, def.frequency_hz);
    j.damping_ratio = std::max(0.f, def.damping_ratio);
    j.grab_point_blend = std::clamp(def.grab_point_blend, 0.f, 1.f);
    j.max_force = def.max_force;
    j.length = 0.f;
    j.alive = true;

    body->disturb();
    body->allow_sleep = false;  // keep awake while grabbed
    // Kill residual velocity so the grab doesn't inherit bounce/jitter
    body->linear_velocity = {};
    body->angular_velocity = {};
    return id;
}

void World::destroy_joint(JointId id) {
    if (id >= joints_.size() || !joint_alive_[id]) {
        return;
    }
    // Restore sleep permission on mouse-joint body
    if (joints_[id].type == JointType::Mouse) {
        if (Body* b = get_body(joints_[id].body_a)) {
            b->allow_sleep = true;
            b->wake();
        }
    }
    joint_alive_[id] = false;
    joints_[id].alive = false;
    free_joint_ids_.push_back(id);
}

Joint* World::get_joint(JointId id) {
    if (id >= joints_.size() || !joint_alive_[id]) {
        return nullptr;
    }
    return &joints_[id];
}

void World::set_mouse_joint_target(JointId id, const Vec3& target) {
    if (Joint* j = get_joint(id)) {
        if (j->type == JointType::Mouse) {
            j->target = target;
            if (Body* b = get_body(j->body_a)) {
                b->disturb();
            }
        }
    }
}

void World::clear() {
    bodies_.clear();
    body_alive_.clear();
    free_ids_.clear();
    contacts_.clear();
    previous_contacts_.clear();
    joints_.clear();
    joint_alive_.clear();
    free_joint_ids_.clear();
    next_id_ = 0;
    next_joint_id_ = 0;
}

void World::wake_body(BodyId id) {
    if (Body* b = get_body(id)) {
        b->wake();
    }
}

std::optional<RaycastHit> World::raycast(const Ray& ray, float max_distance,
                                         bool skip_static) const {
    std::optional<RaycastHit> best;
    for (std::size_t i = 0; i < bodies_.size(); ++i) {
        if (!body_alive_[i]) {
            continue;
        }
        const Body& b = bodies_[i];
        if (skip_static && b.type == BodyType::Static) {
            continue;
        }
        if (auto hit = raycast_body(b, ray, max_distance)) {
            if (!best || hit->distance < best->distance) {
                best = hit;
            }
        }
    }
    return best;
}

void World::integrate_forces(float dt) {
    for (std::size_t i = 0; i < bodies_.size(); ++i) {
        if (!body_alive_[i]) {
            continue;
        }
        Body& b = bodies_[i];
        if (b.type != BodyType::Dynamic || !b.awake) {
            continue;
        }

        b.update_world_inertia();
        b.linear_velocity += (settings_.gravity + b.force * b.inv_mass) * dt;

        Vec3 torque = b.torque;
        if (settings_.enable_gyroscopic) {
            const Mat3 I_world = rotate_inertia(b.orientation.to_mat3(), b.inertia_local);
            torque -= b.angular_velocity.cross(I_world * b.angular_velocity);
        }
        b.angular_velocity += b.inv_inertia_world * torque * dt;
    }
}

void World::integrate_velocities(float dt) {
    for (std::size_t i = 0; i < bodies_.size(); ++i) {
        if (!body_alive_[i]) {
            continue;
        }
        Body& b = bodies_[i];
        if (b.type == BodyType::Static || !b.awake) {
            continue;
        }

        b.position += b.linear_velocity * dt;
        if (!b.fixed_rotation) {
            b.orientation.integrate(b.angular_velocity, dt);
            b.update_world_inertia();
        }
        b.clear_forces();
    }
}

void World::broadphase_and_narrowphase() {
    previous_contacts_ = std::move(contacts_);
    contacts_.clear();

    const float cell = std::max(settings_.broadphase_cell_size, 0.5f);
    std::unordered_map<CellKey, std::vector<BodyId>, CellKeyHash> grid;
    std::vector<BodyId> active;
    active.reserve(bodies_.size());

    for (std::size_t i = 0; i < bodies_.size(); ++i) {
        if (!body_alive_[i]) {
            continue;
        }
        active.push_back(static_cast<BodyId>(i));
        const AABB bb = bodies_[i].bounds().expanded(settings_.linear_slop);
        const int x0 = static_cast<int>(std::floor(bb.min.x / cell));
        const int y0 = static_cast<int>(std::floor(bb.min.y / cell));
        const int z0 = static_cast<int>(std::floor(bb.min.z / cell));
        const int x1 = static_cast<int>(std::floor(bb.max.x / cell));
        const int y1 = static_cast<int>(std::floor(bb.max.y / cell));
        const int z1 = static_cast<int>(std::floor(bb.max.z / cell));
        for (int x = x0; x <= x1; ++x) {
            for (int y = y0; y <= y1; ++y) {
                for (int z = z0; z <= z1; ++z) {
                    grid[{x, y, z}].push_back(static_cast<BodyId>(i));
                }
            }
        }
    }

    std::unordered_set<std::uint64_t> tested;
    tested.reserve(active.size() * 4);

    auto pair_key = [](BodyId a, BodyId b) -> std::uint64_t {
        if (a > b) {
            std::swap(a, b);
        }
        return (static_cast<std::uint64_t>(a) << 32) | b;
    };

    for (const auto& entry : grid) {
        const auto& list = entry.second;
        for (std::size_t i = 0; i < list.size(); ++i) {
            for (std::size_t j = i + 1; j < list.size(); ++j) {
                BodyId id_a = list[i];
                BodyId id_b = list[j];
                if (id_a == id_b) {
                    continue;
                }
                const std::uint64_t key = pair_key(id_a, id_b);
                if (!tested.insert(key).second) {
                    continue;
                }

                Body& a = bodies_[id_a];
                Body& b = bodies_[id_b];

                if (a.type == BodyType::Static && b.type == BodyType::Static) {
                    continue;
                }
                // Both sleeping → skip (no relative motion)
                if (!a.awake && !b.awake) {
                    continue;
                }
                // Chain / pendulum links should not collide with each other
                if (bodies_jointed(id_a, id_b)) {
                    continue;
                }

                if (!a.bounds().expanded(settings_.linear_slop)
                         .overlaps(b.bounds().expanded(settings_.linear_slop))) {
                    continue;
                }

                if (auto manifold = collide(a, b)) {
                    // Wake both if one is awake
                    if (a.awake || b.awake) {
                        a.wake();
                        b.wake();
                    }
                    contacts_.push_back(*manifold);
                }
            }
        }
    }

    if (contact_callback_) {
        for (const Manifold& m : contacts_) {
            contact_callback_(m);
        }
    }
}

void World::warm_start() {
    for (Manifold& m : contacts_) {
        for (Manifold& prev : previous_contacts_) {
            if (!((m.a == prev.a && m.b == prev.b) || (m.a == prev.b && m.b == prev.a))) {
                continue;
            }
            const bool swapped = (m.a != prev.a);
            for (int i = 0; i < m.point_count; ++i) {
                for (int j = 0; j < prev.point_count; ++j) {
                    const Vec3 d = m.points[i].position - prev.points[j].position;
                    if (d.length_sq() < 0.01f * 0.01f) {
                        m.points[i].normal_impulse = prev.points[j].normal_impulse;
                        if (swapped) {
                            m.points[i].tangent_impulse[0] = -prev.points[j].tangent_impulse[0];
                            m.points[i].tangent_impulse[1] = -prev.points[j].tangent_impulse[1];
                        } else {
                            m.points[i].tangent_impulse[0] = prev.points[j].tangent_impulse[0];
                            m.points[i].tangent_impulse[1] = prev.points[j].tangent_impulse[1];
                        }
                    }
                }
            }
        }
    }

    for (const Manifold& m : contacts_) {
        Body* a = get_body(m.a);
        Body* b = get_body(m.b);
        if (!a || !b) {
            continue;
        }
        for (int i = 0; i < m.point_count; ++i) {
            const ContactPoint& cp = m.points[i];
            const Vec3 impulse = m.normal * cp.normal_impulse + m.tangent1 * cp.tangent_impulse[0] +
                                 m.tangent2 * cp.tangent_impulse[1];
            a->apply_impulse(-impulse, cp.position);
            b->apply_impulse(impulse, cp.position);
        }
    }
}

void World::solve_velocity_constraints() {
    for (Manifold& m : contacts_) {
        Body* a = get_body(m.a);
        Body* b = get_body(m.b);
        if (!a || !b) {
            continue;
        }
        const float e = std::min(a->restitution, b->restitution);
        for (int i = 0; i < m.point_count; ++i) {
            ContactPoint& cp = m.points[i];
            const Vec3 rv = b->velocity_at_point(cp.position) - a->velocity_at_point(cp.position);
            const float vn = rv.dot(m.normal);
            cp.velocity_bias =
                (vn < -settings_.restitution_threshold) ? -e * vn : 0.f;
        }
    }

    for (int iter = 0; iter < settings_.velocity_iterations; ++iter) {
        for (Manifold& m : contacts_) {
            Body* a = get_body(m.a);
            Body* b = get_body(m.b);
            if (!a || !b) {
                continue;
            }

            const float mu = std::sqrt(std::max(0.f, a->friction * b->friction));

            for (int i = 0; i < m.point_count; ++i) {
                ContactPoint& cp = m.points[i];
                const Vec3 ra = cp.position - a->position;
                const Vec3 rb = cp.position - b->position;

                Vec3 rv = b->velocity_at_point(cp.position) - a->velocity_at_point(cp.position);
                float vn = rv.dot(m.normal);
                float kn = effective_mass(*a, *b, ra, rb, m.normal);
                if (kn > 1e-8f) {
                    float jn = -(vn - cp.velocity_bias) / kn;
                    const float jn_old = cp.normal_impulse;
                    cp.normal_impulse = std::max(jn_old + jn, 0.f);
                    jn = cp.normal_impulse - jn_old;
                    const Vec3 impulse = m.normal * jn;
                    a->apply_impulse(-impulse, cp.position);
                    b->apply_impulse(impulse, cp.position);
                }

                rv = b->velocity_at_point(cp.position) - a->velocity_at_point(cp.position);
                const float max_friction = mu * cp.normal_impulse;

                for (int t = 0; t < 2; ++t) {
                    const Vec3& tangent = (t == 0) ? m.tangent1 : m.tangent2;
                    const float vt = rv.dot(tangent);
                    const float kt = effective_mass(*a, *b, ra, rb, tangent);
                    if (kt < 1e-8f) {
                        continue;
                    }
                    float jt = -vt / kt;
                    const float jt_old = cp.tangent_impulse[t];
                    cp.tangent_impulse[t] =
                        std::clamp(jt_old + jt, -max_friction, max_friction);
                    jt = cp.tangent_impulse[t] - jt_old;
                    const Vec3 impulse = tangent * jt;
                    a->apply_impulse(-impulse, cp.position);
                    b->apply_impulse(impulse, cp.position);
                    rv = b->velocity_at_point(cp.position) - a->velocity_at_point(cp.position);
                }
            }
        }
    }
}

void World::solve_position_constraints() {
    for (int iter = 0; iter < settings_.position_iterations; ++iter) {
        float max_correction = 0.f;

        for (const Manifold& m : contacts_) {
            Body* a = get_body(m.a);
            Body* b = get_body(m.b);
            if (!a || !b) {
                continue;
            }

            auto fresh = collide(*a, *b);
            if (!fresh || fresh->point_count == 0) {
                continue;
            }

            Vec3 normal = fresh->normal;
            if (normal.dot(m.normal) < 0.f) {
                normal = -normal;
            }

            for (int i = 0; i < fresh->point_count; ++i) {
                const float penetration = fresh->points[i].penetration;
                if (penetration <= settings_.allowed_penetration) {
                    continue;
                }

                const Vec3 point = fresh->points[i].position;
                const Vec3 ra = point - a->position;
                const Vec3 rb = point - b->position;

                const float kn = effective_mass(*a, *b, ra, rb, normal);
                if (kn < 1e-8f) {
                    continue;
                }

                float correction =
                    settings_.baumgarte * (penetration - settings_.allowed_penetration) / kn;
                correction = std::min(correction, settings_.max_linear_correction);
                max_correction = std::max(max_correction, correction);

                const Vec3 correction_vec = normal * correction;

                if (a->type == BodyType::Dynamic && a->awake) {
                    a->position -= correction_vec * a->inv_mass;
                    if (!a->fixed_rotation) {
                        const Vec3 dtheta = a->inv_inertia_world * ra.cross(-correction_vec);
                        a->orientation.integrate(dtheta, 1.f);
                        a->update_world_inertia();
                    }
                }
                if (b->type == BodyType::Dynamic && b->awake) {
                    b->position += correction_vec * b->inv_mass;
                    if (!b->fixed_rotation) {
                        const Vec3 dtheta = b->inv_inertia_world * rb.cross(correction_vec);
                        b->orientation.integrate(dtheta, 1.f);
                        b->update_world_inertia();
                    }
                }
            }
        }

        if (max_correction < settings_.linear_slop) {
            break;
        }
    }
}

bool World::bodies_jointed(BodyId a, BodyId b) const {
    if (a > b) {
        std::swap(a, b);
    }
    for (std::size_t ji = 0; ji < joints_.size(); ++ji) {
        if (!joint_alive_[ji]) {
            continue;
        }
        const Joint& j = joints_[ji];
        if (j.type == JointType::Mouse) {
            continue;
        }
        BodyId ja = j.body_a;
        BodyId jb = j.body_b;
        if (ja > jb) {
            std::swap(ja, jb);
        }
        if (ja == a && jb == b) {
            return true;
        }
    }
    return false;
}

void World::solve_joints(float dt) {
    if (dt <= 0.f) {
        return;
    }
    // Mouse joint: budget impulse for velocity solve later
    for (std::size_t ji = 0; ji < joints_.size(); ++ji) {
        if (!joint_alive_[ji]) {
            continue;
        }
        if (joints_[ji].type == JointType::Mouse) {
            joints_[ji].max_impulse_budget = joints_[ji].max_force * dt;
        }
    }
}

void World::solve_constraint_joints(float dt) {
    if (dt <= 0.f) {
        return;
    }

    constexpr int kIters = 8;
    // Map stiffness → soft constraint (critically damped-ish)
    // stiffness acts as spring k; damping as linear damper coefficient

    for (int iter = 0; iter < kIters; ++iter) {
        for (std::size_t ji = 0; ji < joints_.size(); ++ji) {
            if (!joint_alive_[ji]) {
                continue;
            }
            Joint& joint = joints_[ji];
            if (joint.type == JointType::Mouse) {
                continue;
            }

            Body* a = get_body(joint.body_a);
            Body* b = get_body(joint.body_b);
            if (!a || !b) {
                continue;
            }
            if (a->type != BodyType::Dynamic && b->type != BodyType::Dynamic) {
                continue;
            }
            a->wake();
            b->wake();
            a->update_world_inertia();
            b->update_world_inertia();

            const Vec3 wa = a->local_to_world(joint.local_anchor_a);
            const Vec3 wb = b->local_to_world(joint.local_anchor_b);
            const Vec3 ra = wa - a->position;
            const Vec3 rb = wb - b->position;

            auto kn_along = [&](const Vec3& n) -> float {
                const Vec3 ra_n = ra.cross(n);
                const Vec3 rb_n = rb.cross(n);
                return a->inv_mass + b->inv_mass + ra_n.dot(a->inv_inertia_world * ra_n) +
                       rb_n.dot(b->inv_inertia_world * rb_n);
            };

            if (joint.type == JointType::Distance) {
                Vec3 delta = wb - wa;
                float dist = delta.length();
                if (dist < 1e-5f) {
                    continue;
                }
                const Vec3 n = delta / dist;
                const float C = dist - joint.length;  // position error
                const Vec3 va = a->velocity_at_point(wa);
                const Vec3 vb = b->velocity_at_point(wb);
                const float Cdot = (vb - va).dot(n);

                // Soft constraint from k, c (use reduced mass scale)
                const float kn = kn_along(n);
                if (kn < 1e-8f) {
                    continue;
                }
                const float m = 1.f / kn;
                const float k = std::max(joint.stiffness, 1.f);
                const float c = std::max(joint.damping, 0.f);
                float gamma = dt * (c + dt * k);
                gamma = gamma > 1e-8f ? 1.f / gamma : 0.f;
                const float beta = dt * k * gamma;

                const float lambda = -(Cdot + beta * C) / (kn + gamma);
                const Vec3 impulse = n * lambda;
                a->apply_impulse(-impulse, wa);
                b->apply_impulse(impulse, wb);

                // Light position correction on last iters
                if (iter >= kIters - 2 && std::fabs(C) > 0.01f) {
                    const float corr = std::clamp(C * 0.25f, -0.15f, 0.15f);
                    if (a->type == BodyType::Dynamic) {
                        a->position += n * (corr * a->inv_mass / kn);
                    }
                    if (b->type == BodyType::Dynamic) {
                        b->position -= n * (corr * b->inv_mass / kn);
                    }
                    (void)m;
                }
            } else if (joint.type == JointType::BallSocket) {
                // 3D point-to-point soft constraint
                const Vec3 C = wb - wa;
                const Vec3 va = a->velocity_at_point(wa);
                const Vec3 vb = b->velocity_at_point(wb);
                const Vec3 Cdot = vb - va;

                const float k = std::max(joint.stiffness, 1.f);
                const float c = std::max(joint.damping, 0.f);
                // Use average mass for soft params
                const float m_avg =
                    1.f / std::max(a->inv_mass + b->inv_mass, 1e-4f);
                float gamma = dt * (c + dt * k);
                // Scale soft terms to mass so high k works across body sizes
                gamma = dt * (c + dt * k * m_avg);
                gamma = gamma > 1e-8f ? 1.f / gamma : 0.f;
                const float beta = dt * k * m_avg * gamma;

                for (int axis = 0; axis < 3; ++axis) {
                    Vec3 n{};
                    if (axis == 0) {
                        n = {1.f, 0.f, 0.f};
                    } else if (axis == 1) {
                        n = {0.f, 1.f, 0.f};
                    } else {
                        n = {0.f, 0.f, 1.f};
                    }
                    float kn = kn_along(n) + gamma;
                    if (kn < 1e-8f) {
                        continue;
                    }
                    const float lambda = -(Cdot.dot(n) + beta * C.dot(n)) / kn;
                    const Vec3 impulse = n * lambda;
                    a->apply_impulse(-impulse, wa);
                    b->apply_impulse(impulse, wb);
                }

                if (iter >= kIters - 2) {
                    const float err = C.length();
                    if (err > 0.02f) {
                        const Vec3 n = C / err;
                        const float corr = std::min(err * 0.3f, 0.2f);
                        const float kn = kn_along(n);
                        if (kn > 1e-8f) {
                            if (a->type == BodyType::Dynamic) {
                                a->position += n * (corr * a->inv_mass / kn);
                            }
                            if (b->type == BodyType::Dynamic) {
                                b->position -= n * (corr * b->inv_mass / kn);
                            }
                        }
                    }
                }
            }
        }
    }
}

/// Box2D-style soft mouse constraint (velocity impulse). Stable for interactive grab.
void World::solve_mouse_joints(float dt) {
    if (dt <= 0.f) {
        return;
    }

    constexpr float kPi = 3.14159265f;
    constexpr int kIters = 6;

    for (int iter = 0; iter < kIters; ++iter) {
        for (std::size_t ji = 0; ji < joints_.size(); ++ji) {
            if (!joint_alive_[ji]) {
                continue;
            }
            Joint& joint = joints_[ji];
            if (joint.type != JointType::Mouse) {
                continue;
            }

            Body* body = get_body(joint.body_a);
            if (!body || body->type != BodyType::Dynamic) {
                continue;
            }
            if (!body->awake) {
                body->disturb();
            }
            body->update_world_inertia();

            const Vec3 anchor = body->local_to_world(joint.local_anchor_a);
            const Vec3 r = anchor - body->position;
            const Vec3 C = anchor - joint.target;  // position error
            const Vec3 v = body->velocity_at_point(anchor);

            // Soft constraint params from frequency + damping ratio
            const float omega = 2.f * kPi * joint.frequency_hz;
            const float zeta = joint.damping_ratio;
            // Use body mass for soft params (stable scale)
            const float m = std::max(body->mass, 1e-4f);
            const float k = m * omega * omega;
            const float d = 2.f * m * zeta * omega;
            float gamma = dt * (d + dt * k);
            gamma = gamma > 1e-8f ? 1.f / gamma : 0.f;
            const float beta = dt * k * gamma;

            // Effective mass at grab point (approx, isotropic for softness)
            // For each direction n: K = inv_m + n·(I⁻¹(r×n)×r) = inv_m + (r×n)·I⁻¹(r×n)
            // Soft: solve (K + gamma) * lambda = -(v + beta*C)
            // We solve full 3D with diagonal approximation using average K
            auto effective_mass_dir = [&](const Vec3& n) -> float {
                const Vec3 rn = r.cross(n);
                return body->inv_mass + rn.dot(body->inv_inertia_world * rn);
            };

            // Desired: correct velocity of anchor so it tracks target (stationary this substep)
            Vec3 impulse{};
            for (int axis = 0; axis < 3; ++axis) {
                Vec3 n{};
                if (axis == 0) {
                    n = {1.f, 0.f, 0.f};
                } else if (axis == 1) {
                    n = {0.f, 1.f, 0.f};
                } else {
                    n = {0.f, 0.f, 1.f};
                }

                float kn = effective_mass_dir(n) + gamma;
                if (kn < 1e-8f) {
                    continue;
                }
                const float cdot = v.dot(n) + beta * C.dot(n);
                const float lambda = -cdot / kn;
                impulse += n * lambda;
            }

            // Clamp to remaining force budget for this step
            const float mag = impulse.length();
            if (mag > joint.max_impulse_budget && mag > 1e-8f) {
                impulse *= joint.max_impulse_budget / mag;
                joint.max_impulse_budget = 0.f;
            } else {
                joint.max_impulse_budget = std::max(0.f, joint.max_impulse_budget - mag);
            }

            // Blend COM vs grab-point application — less torque thrash when off-center
            const float blend = joint.grab_point_blend;
            const Vec3 at_point = impulse * blend;
            const Vec3 at_com = impulse * (1.f - blend);
            body->linear_velocity += at_com * body->inv_mass;
            body->apply_impulse(at_point, anchor);

            // Mild angular damping while held keeps grab from spinning spastically
            body->angular_velocity *= 0.92f;
        }
    }
}

void World::update_sleeping(float dt) {
    if (!settings_.allow_sleep) {
        return;
    }

    for (std::size_t i = 0; i < bodies_.size(); ++i) {
        if (!body_alive_[i]) {
            continue;
        }
        Body& b = bodies_[i];
        if (b.type != BodyType::Dynamic || !b.allow_sleep) {
            continue;
        }
        if (!b.awake) {
            continue;
        }

        // Snap micro-velocities so gravity/contact noise doesn't prevent sleep
        const float lin = b.linear_velocity.length();
        const float ang = b.angular_velocity.length();
        if (lin < settings_.sleep_linear_threshold * 0.5f) {
            b.linear_velocity *= 0.8f;
        }
        if (ang < settings_.sleep_angular_threshold * 0.5f) {
            b.angular_velocity *= 0.8f;
        }

        const float lin2 = b.linear_velocity.length();
        const float ang2 = b.angular_velocity.length();
        if (lin2 < settings_.sleep_linear_threshold &&
            ang2 < settings_.sleep_angular_threshold) {
            b.sleep_timer += dt;
            if (b.sleep_timer >= settings_.sleep_time) {
                b.awake = false;
                b.linear_velocity = {};
                b.angular_velocity = {};
                b.clear_forces();
                b.sleep_timer = 0.f;
            }
        } else {
            b.sleep_timer = 0.f;
        }
    }
}

void World::step(float dt) {
    if (dt <= 0.f) {
        return;
    }

    solve_joints(dt);  // init mouse impulse budgets
    integrate_forces(dt);
    broadphase_and_narrowphase();
    warm_start();
    solve_velocity_constraints();
    solve_constraint_joints(dt);  // distance + ball-socket (stable impulse form)
    solve_mouse_joints(dt);       // grab after contacts
    integrate_velocities(dt);
    solve_position_constraints();
    update_sleeping(dt);
}

}  // namespace phys
