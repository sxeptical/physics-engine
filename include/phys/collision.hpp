#pragma once

#include "phys/body.hpp"
#include "phys/vec3.hpp"

#include <array>
#include <optional>

namespace phys {

struct ContactPoint {
    Vec3 position{};
    float penetration = 0.f;
    float normal_impulse = 0.f;
    float tangent_impulse[2] = {0.f, 0.f};
    float velocity_bias = 0.f;
};

struct Manifold {
    BodyId a = kInvalidBody;
    BodyId b = kInvalidBody;
    Vec3 normal{};  // from A toward B
    int point_count = 0;
    std::array<ContactPoint, 4> points{};  // boxes can generate up to 4 contacts
    Vec3 tangent1{};
    Vec3 tangent2{};
};

[[nodiscard]] std::optional<Manifold> collide(const Body& a, const Body& b);

}  // namespace phys
