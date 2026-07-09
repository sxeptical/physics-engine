#pragma once

#include "phys/vec3.hpp"

#include <algorithm>

namespace phys {

struct AABB {
    Vec3 min;
    Vec3 max;

    [[nodiscard]] constexpr Vec3 center() const { return (min + max) * 0.5f; }
    [[nodiscard]] constexpr Vec3 size() const { return max - min; }
    [[nodiscard]] constexpr Vec3 half_extents() const { return (max - min) * 0.5f; }

    [[nodiscard]] constexpr bool overlaps(const AABB& o) const {
        return min.x <= o.max.x && max.x >= o.min.x && min.y <= o.max.y && max.y >= o.min.y &&
               min.z <= o.max.z && max.z >= o.min.z;
    }

    [[nodiscard]] static AABB from_center_size(const Vec3& center, const Vec3& size) {
        const Vec3 half = size * 0.5f;
        return {center - half, center + half};
    }

    [[nodiscard]] AABB expanded(float margin) const {
        const Vec3 m{margin, margin, margin};
        return {min - m, max + m};
    }
};

}  // namespace phys
