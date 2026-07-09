#pragma once

#include "phys/aabb.hpp"
#include "phys/mat3.hpp"
#include "phys/quat.hpp"
#include "phys/vec3.hpp"

#include <algorithm>
#include <cmath>
#include <type_traits>
#include <variant>

namespace phys {

struct Sphere {
    float radius = 0.5f;

    [[nodiscard]] float volume() const {
        return (4.f / 3.f) * 3.14159265f * radius * radius * radius;
    }

    [[nodiscard]] Mat3 inertia_tensor(float mass) const {
        const float i = 0.4f * mass * radius * radius;
        return Mat3::diagonal(i, i, i);
    }

    [[nodiscard]] AABB bounds(const Vec3& position) const {
        const Vec3 r{radius, radius, radius};
        return {position - r, position + r};
    }
};

struct Box {
    Vec3 half_extents{0.5f, 0.5f, 0.5f};

    [[nodiscard]] float volume() const {
        return 8.f * half_extents.x * half_extents.y * half_extents.z;
    }

    [[nodiscard]] Mat3 inertia_tensor(float mass) const {
        const float x2 = half_extents.x * half_extents.x;
        const float y2 = half_extents.y * half_extents.y;
        const float z2 = half_extents.z * half_extents.z;
        return Mat3::diagonal(
            mass * (y2 + z2) / 3.f,
            mass * (x2 + z2) / 3.f,
            mass * (x2 + y2) / 3.f);
    }

    [[nodiscard]] Vec3 corner(int index, const Vec3& position, const Quat& orientation) const {
        const float sx = (index & 1) ? half_extents.x : -half_extents.x;
        const float sy = (index & 2) ? half_extents.y : -half_extents.y;
        const float sz = (index & 4) ? half_extents.z : -half_extents.z;
        return position + orientation.rotate({sx, sy, sz});
    }

    [[nodiscard]] AABB bounds(const Vec3& position, const Quat& orientation) const {
        Vec3 min_p{1e30f, 1e30f, 1e30f};
        Vec3 max_p{-1e30f, -1e30f, -1e30f};
        for (int i = 0; i < 8; ++i) {
            const Vec3 p = corner(i, position, orientation);
            min_p.x = std::min(min_p.x, p.x);
            min_p.y = std::min(min_p.y, p.y);
            min_p.z = std::min(min_p.z, p.z);
            max_p.x = std::max(max_p.x, p.x);
            max_p.y = std::max(max_p.y, p.y);
            max_p.z = std::max(max_p.z, p.z);
        }
        return {min_p, max_p};
    }

    [[nodiscard]] Mat3 axes(const Quat& orientation) const { return orientation.to_mat3(); }
};

/// Capsule aligned with local Y axis: hemisphere caps + cylinder.
struct Capsule {
    float radius = 0.25f;
    float half_height = 0.5f;  // distance from center to start of hemisphere along Y

    [[nodiscard]] float volume() const {
        const float cyl = 3.14159265f * radius * radius * (2.f * half_height);
        const float sph = (4.f / 3.f) * 3.14159265f * radius * radius * radius;
        return cyl + sph;
    }

    [[nodiscard]] Mat3 inertia_tensor(float mass) const {
        // Approximate as cylinder of total length 2*(half_height+radius)
        const float h = 2.f * (half_height + radius);
        const float r2 = radius * radius;
        const float ix = mass * (3.f * r2 + h * h) / 12.f;
        const float iy = 0.5f * mass * r2;
        return Mat3::diagonal(ix, iy, ix);
    }

    [[nodiscard]] Vec3 axis_start(const Vec3& position, const Quat& orientation) const {
        return position + orientation.rotate({0.f, -half_height, 0.f});
    }

    [[nodiscard]] Vec3 axis_end(const Vec3& position, const Quat& orientation) const {
        return position + orientation.rotate({0.f, half_height, 0.f});
    }

    [[nodiscard]] AABB bounds(const Vec3& position, const Quat& orientation) const {
        const Vec3 a = axis_start(position, orientation);
        const Vec3 b = axis_end(position, orientation);
        const Vec3 r{radius, radius, radius};
        return {
            {std::min(a.x, b.x) - radius, std::min(a.y, b.y) - radius, std::min(a.z, b.z) - radius},
            {std::max(a.x, b.x) + radius, std::max(a.y, b.y) + radius, std::max(a.z, b.z) + radius},
        };
    }
};

using Shape = std::variant<Sphere, Box, Capsule>;

inline float shape_volume(const Shape& shape) {
    return std::visit([](const auto& s) { return s.volume(); }, shape);
}

inline Mat3 shape_inertia(const Shape& shape, float mass) {
    return std::visit([mass](const auto& s) { return s.inertia_tensor(mass); }, shape);
}

inline AABB shape_bounds(const Shape& shape, const Vec3& position, const Quat& orientation) {
    return std::visit(
        [&](const auto& s) -> AABB {
            using T = std::decay_t<decltype(s)>;
            if constexpr (std::is_same_v<T, Sphere>) {
                return s.bounds(position);
            } else {
                return s.bounds(position, orientation);
            }
        },
        shape);
}

}  // namespace phys
