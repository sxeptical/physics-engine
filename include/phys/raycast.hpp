#pragma once

#include "phys/body.hpp"
#include "phys/vec3.hpp"

#include <optional>

namespace phys {

struct Ray {
    Vec3 origin{};
    Vec3 direction{};  // should be normalized
};

struct RaycastHit {
    BodyId body = kInvalidBody;
    Vec3 point{};
    Vec3 normal{};
    float distance = 0.f;
};

/// Ray vs single body. Returns closest hit along the ray (t >= 0).
[[nodiscard]] std::optional<RaycastHit> raycast_body(const Body& body, const Ray& ray,
                                                     float max_distance = 1e6f);

}  // namespace phys
