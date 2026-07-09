#pragma once

#include "phys/mat3.hpp"
#include "phys/vec3.hpp"

#include <cmath>

namespace phys {

struct Quat {
    float w = 1.f;
    float x = 0.f;
    float y = 0.f;
    float z = 0.f;

    constexpr Quat() = default;
    constexpr Quat(float w_, float x_, float y_, float z_) : w(w_), x(x_), y(y_), z(z_) {}

    static constexpr Quat identity() { return {1.f, 0.f, 0.f, 0.f}; }

    static Quat from_axis_angle(const Vec3& axis, float radians) {
        const Vec3 n = axis.normalized();
        const float half = radians * 0.5f;
        const float s = std::sin(half);
        return {std::cos(half), n.x * s, n.y * s, n.z * s};
    }

    static Quat from_euler_xyz(float rx, float ry, float rz) {
        const Quat qx = from_axis_angle({1.f, 0.f, 0.f}, rx);
        const Quat qy = from_axis_angle({0.f, 1.f, 0.f}, ry);
        const Quat qz = from_axis_angle({0.f, 0.f, 1.f}, rz);
        return qz * qy * qx;
    }

    [[nodiscard]] Quat operator*(const Quat& o) const {
        return {
            w * o.w - x * o.x - y * o.y - z * o.z,
            w * o.x + x * o.w + y * o.z - z * o.y,
            w * o.y - x * o.z + y * o.w + z * o.x,
            w * o.z + x * o.y - y * o.x + z * o.w,
        };
    }

    [[nodiscard]] Quat conjugate() const { return {w, -x, -y, -z}; }

    [[nodiscard]] float length_sq() const { return w * w + x * x + y * y + z * z; }

    [[nodiscard]] Quat normalized() const {
        const float len = std::sqrt(length_sq());
        if (len < 1e-8f) {
            return identity();
        }
        const float inv = 1.f / len;
        return {w * inv, x * inv, y * inv, z * inv};
    }

    /// Rotate vector by this quaternion (unit quat assumed).
    [[nodiscard]] Vec3 rotate(const Vec3& v) const {
        // v' = q * (0,v) * q^-1
        const Vec3 qv{x, y, z};
        const Vec3 t = 2.f * qv.cross(v);
        return v + w * t + qv.cross(t);
    }

    [[nodiscard]] Mat3 to_mat3() const {
        const float xx = x * x, yy = y * y, zz = z * z;
        const float xy = x * y, xz = x * z, yz = y * z;
        const float wx = w * x, wy = w * y, wz = w * z;

        Mat3 m;
        m.cols[0] = {1.f - 2.f * (yy + zz), 2.f * (xy + wz), 2.f * (xz - wy)};
        m.cols[1] = {2.f * (xy - wz), 1.f - 2.f * (xx + zz), 2.f * (yz + wx)};
        m.cols[2] = {2.f * (xz + wy), 2.f * (yz - wx), 1.f - 2.f * (xx + yy)};
        return m;
    }

    /// Integrate orientation with angular velocity ω (rad/s): q += 0.5 * Ω(ω) * q * dt
    void integrate(const Vec3& omega, float dt) {
        const Quat omega_q{0.f, omega.x, omega.y, omega.z};
        const Quat dq = omega_q * (*this);
        w += 0.5f * dq.w * dt;
        x += 0.5f * dq.x * dt;
        y += 0.5f * dq.y * dt;
        z += 0.5f * dq.z * dt;
        *this = normalized();
    }
};

}  // namespace phys
