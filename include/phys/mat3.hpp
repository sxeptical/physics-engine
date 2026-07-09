#pragma once

#include "phys/vec3.hpp"

#include <cmath>

namespace phys {

namespace detail {
inline constexpr float& vec_at(Vec3& v, int i) {
    return i == 0 ? v.x : (i == 1 ? v.y : v.z);
}
inline constexpr float vec_at(const Vec3& v, int i) {
    return i == 0 ? v.x : (i == 1 ? v.y : v.z);
}
}  // namespace detail

/// Column-major 3×3 matrix.
struct Mat3 {
    Vec3 cols[3]{};

    constexpr Mat3() = default;
    constexpr Mat3(const Vec3& c0, const Vec3& c1, const Vec3& c2) : cols{c0, c1, c2} {}

    static constexpr Mat3 identity() {
        return {
            {1.f, 0.f, 0.f},
            {0.f, 1.f, 0.f},
            {0.f, 0.f, 1.f},
        };
    }

    static constexpr Mat3 zero() { return {}; }

    static constexpr Mat3 diagonal(float xx, float yy, float zz) {
        return {
            {xx, 0.f, 0.f},
            {0.f, yy, 0.f},
            {0.f, 0.f, zz},
        };
    }

    [[nodiscard]] constexpr float operator()(int r, int c) const {
        return detail::vec_at(cols[c], r);
    }
    float& operator()(int r, int c) { return detail::vec_at(cols[c], r); }

    [[nodiscard]] constexpr Vec3 operator*(const Vec3& v) const {
        return cols[0] * v.x + cols[1] * v.y + cols[2] * v.z;
    }

    [[nodiscard]] Mat3 operator*(const Mat3& o) const {
        return {
            (*this) * o.cols[0],
            (*this) * o.cols[1],
            (*this) * o.cols[2],
        };
    }

    [[nodiscard]] Mat3 operator*(float s) const {
        return {cols[0] * s, cols[1] * s, cols[2] * s};
    }

    [[nodiscard]] Mat3 transposed() const {
        return {
            {(*this)(0, 0), (*this)(0, 1), (*this)(0, 2)},
            {(*this)(1, 0), (*this)(1, 1), (*this)(1, 2)},
            {(*this)(2, 0), (*this)(2, 1), (*this)(2, 2)},
        };
    }

    [[nodiscard]] Mat3 inverted() const {
        const float a = (*this)(0, 0), b = (*this)(0, 1), c = (*this)(0, 2);
        const float d = (*this)(1, 0), e = (*this)(1, 1), f = (*this)(1, 2);
        const float g = (*this)(2, 0), h = (*this)(2, 1), i = (*this)(2, 2);

        const float A = e * i - f * h;
        const float B = f * g - d * i;
        const float C = d * h - e * g;
        const float D = c * h - b * i;
        const float E = a * i - c * g;
        const float F = b * g - a * h;
        const float G = b * f - c * e;
        const float H = c * d - a * f;
        const float I = a * e - b * d;

        const float det = a * A + b * B + c * C;
        if (std::fabs(det) < 1e-12f) {
            return zero();
        }
        const float inv = 1.f / det;
        Mat3 out;
        out(0, 0) = A * inv;
        out(1, 0) = B * inv;
        out(2, 0) = C * inv;
        out(0, 1) = D * inv;
        out(1, 1) = E * inv;
        out(2, 1) = F * inv;
        out(0, 2) = G * inv;
        out(1, 2) = H * inv;
        out(2, 2) = I * inv;
        return out;
    }
};

/// Rotate a local-space diagonal inertia into world space: R * I * R^T
inline Mat3 rotate_inertia(const Mat3& rotation, const Mat3& local_inertia) {
    return rotation * local_inertia * rotation.transposed();
}

}  // namespace phys
