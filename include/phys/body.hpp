#pragma once

#include "phys/mat3.hpp"
#include "phys/quat.hpp"
#include "phys/shapes.hpp"
#include "phys/vec3.hpp"

#include <algorithm>
#include <cstdint>

namespace phys {

using BodyId = std::uint32_t;
constexpr BodyId kInvalidBody = static_cast<BodyId>(-1);

enum class BodyType {
    Dynamic,
    Static,
    Kinematic,
};

struct BodyDef {
    BodyType type = BodyType::Dynamic;
    Shape shape = Sphere{0.5f};
    Vec3 position{0.f, 0.f, 0.f};
    Quat orientation = Quat::identity();
    Vec3 linear_velocity{0.f, 0.f, 0.f};
    Vec3 angular_velocity{0.f, 0.f, 0.f};
    float density = 1.f;
    float restitution = 0.2f;
    float friction = 0.3f;
    bool fixed_rotation = false;
    bool allow_sleep = true;
};

struct Body {
    BodyId id = kInvalidBody;
    BodyType type = BodyType::Dynamic;
    Shape shape = Sphere{0.5f};

    Vec3 position{};
    Quat orientation = Quat::identity();

    Vec3 linear_velocity{};
    Vec3 angular_velocity{};

    Vec3 force{};
    Vec3 torque{};

    float mass = 1.f;
    float inv_mass = 1.f;
    Mat3 inertia_local = Mat3::identity();
    Mat3 inv_inertia_local = Mat3::identity();
    Mat3 inv_inertia_world = Mat3::identity();

    float restitution = 0.2f;
    float friction = 0.3f;
    bool fixed_rotation = false;
    bool allow_sleep = true;
    bool awake = true;
    float sleep_timer = 0.f;

    void set_mass_from_density(float density) {
        if (type != BodyType::Dynamic) {
            mass = 0.f;
            inv_mass = 0.f;
            inertia_local = Mat3::zero();
            inv_inertia_local = Mat3::zero();
            inv_inertia_world = Mat3::zero();
            return;
        }

        const float vol = std::max(shape_volume(shape), 1e-6f);
        mass = std::max(vol * density, 1e-6f);
        inv_mass = 1.f / mass;

        if (fixed_rotation) {
            inertia_local = Mat3::zero();
            inv_inertia_local = Mat3::zero();
            inv_inertia_world = Mat3::zero();
        } else {
            inertia_local = shape_inertia(shape, mass);
            inv_inertia_local = inertia_local.inverted();
            update_world_inertia();
        }
    }

    void update_world_inertia() {
        if (type != BodyType::Dynamic || fixed_rotation) {
            inv_inertia_world = Mat3::zero();
            return;
        }
        const Mat3 r = orientation.to_mat3();
        inv_inertia_world = rotate_inertia(r, inv_inertia_local);
    }

    [[nodiscard]] AABB bounds() const {
        return shape_bounds(shape, position, orientation);
    }

    [[nodiscard]] Vec3 velocity_at_point(const Vec3& world_point) const {
        const Vec3 r = world_point - position;
        return linear_velocity + angular_velocity.cross(r);
    }

    [[nodiscard]] Vec3 world_to_local(const Vec3& world_point) const {
        return orientation.conjugate().rotate(world_point - position);
    }

    [[nodiscard]] Vec3 local_to_world(const Vec3& local_point) const {
        return position + orientation.rotate(local_point);
    }

    /// Wake a sleeping body. Safe to call every frame on already-awake bodies
    /// (does not reset the sleep timer unless the body was actually asleep).
    void wake() {
        if (type == BodyType::Static) {
            return;
        }
        if (!awake) {
            awake = true;
            sleep_timer = 0.f;
        }
    }

    /// Force-reset activity (e.g. user grab, large impulse) — restarts sleep timer.
    void disturb() {
        if (type == BodyType::Static) {
            return;
        }
        awake = true;
        sleep_timer = 0.f;
    }

    void apply_force(const Vec3& f) {
        if (type != BodyType::Dynamic) {
            return;
        }
        if (f.length_sq() > 1e-12f) {
            disturb();
        }
        force += f;
    }

    void apply_force_at_point(const Vec3& f, const Vec3& world_point) {
        if (type != BodyType::Dynamic) {
            return;
        }
        if (f.length_sq() > 1e-12f) {
            disturb();
        }
        force += f;
        torque += (world_point - position).cross(f);
    }

    void apply_impulse(const Vec3& impulse, const Vec3& world_point) {
        if (type != BodyType::Dynamic) {
            return;
        }
        // Solver impulses must not reset sleep_timer every frame — only un-sleep.
        if (!awake) {
            awake = true;
            sleep_timer = 0.f;
        }
        linear_velocity += impulse * inv_mass;
        const Vec3 r = world_point - position;
        angular_velocity += inv_inertia_world * r.cross(impulse);
    }

    void clear_forces() {
        force = {};
        torque = {};
    }
};

}  // namespace phys
