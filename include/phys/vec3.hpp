#pragma once

#include <cmath>
#include <ostream>

namespace phys {

struct Vec3 {
    float x = 0.f;
    float y = 0.f;
    float z = 0.f;

    constexpr Vec3() = default;
    constexpr Vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}

    constexpr Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    constexpr Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    constexpr Vec3 operator*(float s) const { return {x * s, y * s, z * s}; }
    constexpr Vec3 operator/(float s) const { return {x / s, y / s, z / s}; }
    constexpr Vec3 operator-() const { return {-x, -y, -z}; }

    Vec3& operator+=(const Vec3& o) {
        x += o.x;
        y += o.y;
        z += o.z;
        return *this;
    }
    Vec3& operator-=(const Vec3& o) {
        x -= o.x;
        y -= o.y;
        z -= o.z;
        return *this;
    }
    Vec3& operator*=(float s) {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }

    [[nodiscard]] constexpr float length_sq() const { return x * x + y * y + z * z; }
    [[nodiscard]] float length() const { return std::sqrt(length_sq()); }

    [[nodiscard]] Vec3 normalized() const {
        const float len = length();
        if (len < 1e-8f) {
            return {};
        }
        return *this / len;
    }

    [[nodiscard]] constexpr float dot(const Vec3& o) const {
        return x * o.x + y * o.y + z * o.z;
    }

    [[nodiscard]] constexpr Vec3 cross(const Vec3& o) const {
        return {
            y * o.z - z * o.y,
            z * o.x - x * o.z,
            x * o.y - y * o.x,
        };
    }

    [[nodiscard]] constexpr Vec3 cwise_mul(const Vec3& o) const {
        return {x * o.x, y * o.y, z * o.z};
    }

    [[nodiscard]] constexpr float max_abs_component() const {
        const float ax = x < 0.f ? -x : x;
        const float ay = y < 0.f ? -y : y;
        const float az = z < 0.f ? -z : z;
        float m = ax;
        if (ay > m) {
            m = ay;
        }
        if (az > m) {
            m = az;
        }
        return m;
    }
};

inline constexpr Vec3 operator*(float s, const Vec3& v) { return v * s; }

inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    return os << '(' << v.x << ", " << v.y << ", " << v.z << ')';
}

/// Build an orthonormal basis with `n` as the first axis.
inline void build_orthonormal_basis(const Vec3& n, Vec3& t1, Vec3& t2) {
    if (std::fabs(n.x) >= 0.57735f) {
        t1 = Vec3{n.y, -n.x, 0.f}.normalized();
    } else {
        t1 = Vec3{0.f, n.z, -n.y}.normalized();
    }
    t2 = n.cross(t1);
}

}  // namespace phys
