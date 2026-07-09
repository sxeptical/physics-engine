#include "phys/raycast.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

namespace phys {
namespace {

constexpr float kEpsilon = 1e-6f;

std::optional<RaycastHit> ray_sphere(const Body& body, const Sphere& s, const Ray& ray,
                                     float max_distance) {
    const Vec3 oc = ray.origin - body.position;
    const float b = oc.dot(ray.direction);
    const float c = oc.dot(oc) - s.radius * s.radius;
    const float disc = b * b - c;
    if (disc < 0.f) {
        return std::nullopt;
    }
    const float sqrt_d = std::sqrt(disc);
    float t = -b - sqrt_d;
    if (t < 0.f) {
        t = -b + sqrt_d;
    }
    if (t < 0.f || t > max_distance) {
        return std::nullopt;
    }
    RaycastHit hit;
    hit.body = body.id;
    hit.distance = t;
    hit.point = ray.origin + ray.direction * t;
    hit.normal = (hit.point - body.position).normalized();
    return hit;
}

std::optional<RaycastHit> ray_box(const Body& body, const Box& box, const Ray& ray,
                                  float max_distance) {
    // Transform ray to local box space
    const Mat3 rot = body.orientation.to_mat3();
    const Mat3 rot_t = rot.transposed();
    const Vec3 local_origin = rot_t * (ray.origin - body.position);
    const Vec3 local_dir = rot_t * ray.direction;

    float tmin = 0.f;
    float tmax = max_distance;
    const float he[3] = {box.half_extents.x, box.half_extents.y, box.half_extents.z};
    const float o[3] = {local_origin.x, local_origin.y, local_origin.z};
    const float d[3] = {local_dir.x, local_dir.y, local_dir.z};
    int hit_axis = 0;
    float hit_sign = 1.f;

    for (int i = 0; i < 3; ++i) {
        if (std::fabs(d[i]) < kEpsilon) {
            if (o[i] < -he[i] || o[i] > he[i]) {
                return std::nullopt;
            }
            continue;
        }
        float t1 = (-he[i] - o[i]) / d[i];
        float t2 = (he[i] - o[i]) / d[i];
        float sign1 = -1.f;
        float sign2 = 1.f;
        if (t1 > t2) {
            std::swap(t1, t2);
            std::swap(sign1, sign2);
        }
        if (t1 > tmin) {
            tmin = t1;
            hit_axis = i;
            hit_sign = sign1;
        }
        tmax = std::min(tmax, t2);
        if (tmin > tmax) {
            return std::nullopt;
        }
    }

    if (tmin < 0.f || tmin > max_distance) {
        return std::nullopt;
    }

    RaycastHit hit;
    hit.body = body.id;
    hit.distance = tmin;
    hit.point = ray.origin + ray.direction * tmin;
    Vec3 local_n{};
    if (hit_axis == 0) {
        local_n = {hit_sign, 0.f, 0.f};
    } else if (hit_axis == 1) {
        local_n = {0.f, hit_sign, 0.f};
    } else {
        local_n = {0.f, 0.f, hit_sign};
    }
    hit.normal = rot * local_n;
    return hit;
}

std::optional<RaycastHit> ray_capsule(const Body& body, const Capsule& cap, const Ray& ray,
                                      float max_distance) {
    // Approximate: ray vs inflated segment (two spheres + cylinder)
    // Sample by testing spheres at ends and closest approach to axis
    const Vec3 a = cap.axis_start(body.position, body.orientation);
    const Vec3 b = cap.axis_end(body.position, body.orientation);

    std::optional<RaycastHit> best;

    auto consider_sphere = [&](const Vec3& center) {
        Body tmp = body;
        tmp.position = center;
        if (auto h = ray_sphere(tmp, Sphere{cap.radius}, ray, max_distance)) {
            h->body = body.id;
            if (!best || h->distance < best->distance) {
                best = h;
            }
        }
    };
    consider_sphere(a);
    consider_sphere(b);

    // Cylinder (infinite, then clamp) — quadratic against axis
    const Vec3 ba = b - a;
    const float ba_len_sq = ba.length_sq();
    if (ba_len_sq > kEpsilon) {
        const Vec3 ro = ray.origin - a;
        const Vec3 ba_n = ba / std::sqrt(ba_len_sq);
        // Project ray onto plane perpendicular to axis
        const Vec3 d = ray.direction - ba_n * ray.direction.dot(ba_n);
        const Vec3 o = ro - ba_n * ro.dot(ba_n);
        const float A = d.dot(d);
        const float B = 2.f * o.dot(d);
        const float C = o.dot(o) - cap.radius * cap.radius;
        if (A > kEpsilon) {
            const float disc = B * B - 4.f * A * C;
            if (disc >= 0.f) {
                const float t = (-B - std::sqrt(disc)) / (2.f * A);
                if (t >= 0.f && t <= max_distance) {
                    const Vec3 p = ray.origin + ray.direction * t;
                    const float along = (p - a).dot(ba_n);
                    if (along >= 0.f && along * along <= ba_len_sq) {
                        const Vec3 axis_pt = a + ba_n * along;
                        RaycastHit hit;
                        hit.body = body.id;
                        hit.distance = t;
                        hit.point = p;
                        hit.normal = (p - axis_pt).normalized();
                        if (!best || hit.distance < best->distance) {
                            best = hit;
                        }
                    }
                }
            }
        }
    }
    return best;
}

}  // namespace

std::optional<RaycastHit> raycast_body(const Body& body, const Ray& ray, float max_distance) {
    Ray r = ray;
    const float dir_len = r.direction.length();
    if (dir_len < kEpsilon) {
        return std::nullopt;
    }
    r.direction = r.direction / dir_len;

    return std::visit(
        [&](const auto& shape) -> std::optional<RaycastHit> {
            using T = std::decay_t<decltype(shape)>;
            if constexpr (std::is_same_v<T, Sphere>) {
                return ray_sphere(body, shape, r, max_distance);
            } else if constexpr (std::is_same_v<T, Box>) {
                return ray_box(body, shape, r, max_distance);
            } else if constexpr (std::is_same_v<T, Capsule>) {
                return ray_capsule(body, shape, r, max_distance);
            }
            return std::nullopt;
        },
        body.shape);
}

}  // namespace phys
