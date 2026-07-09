#include "phys/collision.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

namespace phys {
namespace {

constexpr float kEpsilon = 1e-6f;

void finalize_manifold(Manifold& m) {
    build_orthonormal_basis(m.normal, m.tangent1, m.tangent2);
}

std::optional<Manifold> sphere_sphere(const Body& a, const Body& b, const Sphere& sa,
                                      const Sphere& sb) {
    const Vec3 delta = b.position - a.position;
    const float dist_sq = delta.length_sq();
    const float radius_sum = sa.radius + sb.radius;

    if (dist_sq > radius_sum * radius_sum) {
        return std::nullopt;
    }

    Manifold m;
    m.a = a.id;
    m.b = b.id;
    m.point_count = 1;

    if (dist_sq < kEpsilon) {
        m.normal = {0.f, 1.f, 0.f};
        m.points[0].position = a.position;
        m.points[0].penetration = radius_sum;
        finalize_manifold(m);
        return m;
    }

    const float dist = std::sqrt(dist_sq);
    m.normal = delta / dist;
    m.points[0].penetration = radius_sum - dist;
    m.points[0].position = a.position + m.normal * sa.radius;
    finalize_manifold(m);
    return m;
}

std::optional<Manifold> sphere_box(const Body& sphere_body, const Body& box_body,
                                   const Sphere& sphere, const Box& box, bool swapped) {
    const Mat3 rot = box_body.orientation.to_mat3();
    const Mat3 rot_t = rot.transposed();
    const Vec3 local = rot_t * (sphere_body.position - box_body.position);

    const Vec3 closest{
        std::clamp(local.x, -box.half_extents.x, box.half_extents.x),
        std::clamp(local.y, -box.half_extents.y, box.half_extents.y),
        std::clamp(local.z, -box.half_extents.z, box.half_extents.z),
    };

    Vec3 local_delta = local - closest;
    float dist_sq = local_delta.length_sq();

    Manifold m;
    m.point_count = 1;

    if (dist_sq < kEpsilon) {
        // Sphere center inside box — push out along nearest face
        const float dx = box.half_extents.x - std::fabs(local.x);
        const float dy = box.half_extents.y - std::fabs(local.y);
        const float dz = box.half_extents.z - std::fabs(local.z);

        Vec3 local_n;
        float pen;
        Vec3 local_contact;

        if (dx <= dy && dx <= dz) {
            const float sign = local.x >= 0.f ? 1.f : -1.f;
            local_n = {sign, 0.f, 0.f};
            pen = dx + sphere.radius;
            local_contact = {sign * box.half_extents.x, local.y, local.z};
        } else if (dy <= dx && dy <= dz) {
            const float sign = local.y >= 0.f ? 1.f : -1.f;
            local_n = {0.f, sign, 0.f};
            pen = dy + sphere.radius;
            local_contact = {local.x, sign * box.half_extents.y, local.z};
        } else {
            const float sign = local.z >= 0.f ? 1.f : -1.f;
            local_n = {0.f, 0.f, sign};
            pen = dz + sphere.radius;
            local_contact = {local.x, local.y, sign * box.half_extents.z};
        }

        const Vec3 world_n = rot * local_n;
        m.points[0].penetration = pen;
        m.points[0].position = box_body.position + rot * local_contact;

        if (swapped) {
            m.a = box_body.id;
            m.b = sphere_body.id;
            m.normal = world_n;
        } else {
            m.a = sphere_body.id;
            m.b = box_body.id;
            m.normal = -world_n;
        }
        finalize_manifold(m);
        return m;
    }

    if (dist_sq > sphere.radius * sphere.radius) {
        return std::nullopt;
    }

    const float dist = std::sqrt(dist_sq);
    const Vec3 local_n = local_delta / dist;
    const Vec3 world_n = rot * local_n;

    m.points[0].penetration = sphere.radius - dist;
    m.points[0].position = box_body.position + rot * closest;

    if (swapped) {
        m.a = box_body.id;
        m.b = sphere_body.id;
        m.normal = world_n;
    } else {
        m.a = sphere_body.id;
        m.b = box_body.id;
        m.normal = -world_n;
    }
    finalize_manifold(m);
    return m;
}

// ---------------------------------------------------------------------------
// Box–box SAT
// ---------------------------------------------------------------------------

struct FaceAxisResult {
    float separation = -std::numeric_limits<float>::infinity();
    int axis_index = 0;  // 0,1,2 for local axes
    Vec3 axis{};
};

float project_box(const Box& box, const Mat3& axes, const Vec3& axis) {
    return box.half_extents.x * std::fabs(axis.dot(axes.cols[0])) +
           box.half_extents.y * std::fabs(axis.dot(axes.cols[1])) +
           box.half_extents.z * std::fabs(axis.dot(axes.cols[2]));
}

FaceAxisResult query_face_axes(const Body& a, const Box& ba, const Mat3& axes_a, const Body& b,
                               const Box& bb, const Mat3& axes_b) {
    FaceAxisResult best;
    best.separation = -std::numeric_limits<float>::infinity();
    const Vec3 d = b.position - a.position;

    for (int i = 0; i < 3; ++i) {
        Vec3 axis = axes_a.cols[i];
        const float dist = d.dot(axis);
        // Ensure axis points from A toward B
        if (dist < 0.f) {
            axis = -axis;
        }

        const float ra = project_box(ba, axes_a, axis);
        const float rb = project_box(bb, axes_b, axis);
        const float separation = std::fabs(d.dot(axis)) - ra - rb;

        if (separation > best.separation) {
            best.separation = separation;
            best.axis_index = i;
            best.axis = axis;
        }
    }
    return best;
}

struct EdgeAxisResult {
    float separation = -std::numeric_limits<float>::infinity();
    int edge_a = 0;
    int edge_b = 0;
    Vec3 axis{};
};

EdgeAxisResult query_edge_axes(const Body& a, const Box& ba, const Mat3& axes_a, const Body& b,
                               const Box& bb, const Mat3& axes_b) {
    EdgeAxisResult best;
    best.separation = -std::numeric_limits<float>::infinity();
    const Vec3 d = b.position - a.position;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Vec3 axis = axes_a.cols[i].cross(axes_b.cols[j]);
            const float len_sq = axis.length_sq();
            if (len_sq < 1e-8f) {
                continue;  // parallel edges
            }
            axis = axis / std::sqrt(len_sq);
            if (d.dot(axis) < 0.f) {
                axis = -axis;
            }

            const float ra = project_box(ba, axes_a, axis);
            const float rb = project_box(bb, axes_b, axis);
            const float separation = std::fabs(d.dot(axis)) - ra - rb;

            if (separation > best.separation) {
                best.separation = separation;
                best.edge_a = i;
                best.edge_b = j;
                best.axis = axis;
            }
        }
    }
    return best;
}

// Support point of a box in a world-space direction
Vec3 box_support(const Body& body, const Box& box, const Mat3& axes, const Vec3& dir) {
    return body.position + axes.cols[0] * (dir.dot(axes.cols[0]) >= 0.f ? box.half_extents.x
                                                                        : -box.half_extents.x) +
           axes.cols[1] * (dir.dot(axes.cols[1]) >= 0.f ? box.half_extents.y
                                                        : -box.half_extents.y) +
           axes.cols[2] * (dir.dot(axes.cols[2]) >= 0.f ? box.half_extents.z
                                                        : -box.half_extents.z);
}

// Closest points between two line segments
void closest_points_segments(const Vec3& p1, const Vec3& q1, const Vec3& p2, const Vec3& q2,
                             Vec3& c1, Vec3& c2) {
    const Vec3 d1 = q1 - p1;
    const Vec3 d2 = q2 - p2;
    const Vec3 r = p1 - p2;
    const float a = d1.dot(d1);
    const float e = d2.dot(d2);
    const float f = d2.dot(r);

    float s, t;
    if (a <= kEpsilon && e <= kEpsilon) {
        c1 = p1;
        c2 = p2;
        return;
    }
    if (a <= kEpsilon) {
        s = 0.f;
        t = std::clamp(f / e, 0.f, 1.f);
    } else {
        const float c = d1.dot(r);
        if (e <= kEpsilon) {
            t = 0.f;
            s = std::clamp(-c / a, 0.f, 1.f);
        } else {
            const float b = d1.dot(d2);
            const float denom = a * e - b * b;
            s = denom != 0.f ? std::clamp((b * f - c * e) / denom, 0.f, 1.f) : 0.f;
            t = (b * s + f) / e;
            if (t < 0.f) {
                t = 0.f;
                s = std::clamp(-c / a, 0.f, 1.f);
            } else if (t > 1.f) {
                t = 1.f;
                s = std::clamp((b - c) / a, 0.f, 1.f);
            }
        }
    }
    c1 = p1 + d1 * s;
    c2 = p2 + d2 * t;
}

// Get the 4 vertices of the face whose outward normal best matches `axis`
void get_face_vertices(const Body& body, const Box& box, const Mat3& axes, const Vec3& axis,
                       Vec3 out[4], Vec3& face_normal) {
    int best = 0;
    float best_dot = -std::numeric_limits<float>::infinity();
    float best_sign = 1.f;
    for (int i = 0; i < 3; ++i) {
        const float d = axes.cols[i].dot(axis);
        if (d > best_dot) {
            best_dot = d;
            best = i;
            best_sign = 1.f;
        }
        if (-d > best_dot) {
            best_dot = -d;
            best = i;
            best_sign = -1.f;
        }
    }

    const int u = (best + 1) % 3;
    const int v = (best + 2) % 3;
    const float he[3] = {box.half_extents.x, box.half_extents.y, box.half_extents.z};
    const Vec3 center_pt = body.position + axes.cols[best] * (best_sign * he[best]);
    face_normal = axes.cols[best] * best_sign;

    out[0] = center_pt + axes.cols[u] * he[u] + axes.cols[v] * he[v];
    out[1] = center_pt - axes.cols[u] * he[u] + axes.cols[v] * he[v];
    out[2] = center_pt - axes.cols[u] * he[u] - axes.cols[v] * he[v];
    out[3] = center_pt + axes.cols[u] * he[u] - axes.cols[v] * he[v];
}

// Clip polygon against plane (keep points with n·x <= d)
int clip_poly(const Vec3* in, int in_count, const Vec3& n, float d, Vec3* out) {
    int out_count = 0;
    for (int i = 0; i < in_count; ++i) {
        const Vec3& a = in[i];
        const Vec3& b = in[(i + 1) % in_count];
        const float da = a.dot(n) - d;
        const float db = b.dot(n) - d;

        if (da <= 0.f && db <= 0.f) {
            out[out_count++] = b;
        } else if (da <= 0.f && db > 0.f) {
            const float t = da / (da - db);
            out[out_count++] = a + (b - a) * t;
        } else if (da > 0.f && db <= 0.f) {
            const float t = da / (da - db);
            out[out_count++] = a + (b - a) * t;
            out[out_count++] = b;
        }
    }
    return out_count;
}

std::optional<Manifold> box_box(const Body& a, const Body& b, const Box& ba, const Box& bb) {
    const Mat3 axes_a = a.orientation.to_mat3();
    const Mat3 axes_b = b.orientation.to_mat3();

    const FaceAxisResult face_a = query_face_axes(a, ba, axes_a, b, bb, axes_b);
    if (face_a.separation > 0.f) {
        return std::nullopt;
    }
    const FaceAxisResult face_b = query_face_axes(b, bb, axes_b, a, ba, axes_a);
    if (face_b.separation > 0.f) {
        return std::nullopt;
    }
    const EdgeAxisResult edges = query_edge_axes(a, ba, axes_a, b, bb, axes_b);
    if (edges.separation > 0.f) {
        return std::nullopt;
    }

    constexpr float kRelTol = 0.95f;
    constexpr float kAbsTol = 0.01f;

    enum class AxisType { FaceA, FaceB, Edge };
    AxisType type = AxisType::FaceA;
    float best_sep = face_a.separation;
    Vec3 normal = face_a.axis;

    if (face_b.separation > kRelTol * best_sep + kAbsTol * std::min(1.f, -best_sep)) {
        type = AxisType::FaceB;
        best_sep = face_b.separation;
        normal = -face_b.axis;  // from A to B
    }
    // Prefer face over edge slightly for stability
    if (edges.separation > kRelTol * best_sep + kAbsTol * std::min(1.f, -best_sep)) {
        type = AxisType::Edge;
        best_sep = edges.separation;
        normal = edges.axis;
    }

    Manifold m;
    m.a = a.id;
    m.b = b.id;
    m.normal = normal;

    if (type == AxisType::Edge) {
        // Find supporting edges and closest points
        const Vec3 axis_a = axes_a.cols[edges.edge_a];
        const Vec3 axis_b = axes_b.cols[edges.edge_b];

        // Edge midpoints via support in directions orthogonal to edge toward other
        const Vec3 mid_a = box_support(a, ba, axes_a, normal);
        const Vec3 mid_b = box_support(b, bb, axes_b, -normal);

        // Build edge segments through support points along edge directions
        const float he_a[3] = {ba.half_extents.x, ba.half_extents.y, ba.half_extents.z};
        const float he_b[3] = {bb.half_extents.x, bb.half_extents.y, bb.half_extents.z};
        const Vec3 p1 = mid_a - axis_a * he_a[edges.edge_a];
        const Vec3 q1 = mid_a + axis_a * he_a[edges.edge_a];
        const Vec3 p2 = mid_b - axis_b * he_b[edges.edge_b];
        const Vec3 q2 = mid_b + axis_b * he_b[edges.edge_b];

        Vec3 c1, c2;
        closest_points_segments(p1, q1, p2, q2, c1, c2);

        m.point_count = 1;
        m.points[0].position = (c1 + c2) * 0.5f;
        m.points[0].penetration = -best_sep;
        // Ensure normal points A → B
        if ((b.position - a.position).dot(m.normal) < 0.f) {
            m.normal = -m.normal;
        }
        finalize_manifold(m);
        return m;
    }

    // Face contact: clip incident face against reference face side planes
    const bool a_is_ref = (type == AxisType::FaceA);
    const Body& ref_body = a_is_ref ? a : b;
    const Body& inc_body = a_is_ref ? b : a;
    const Box& ref_box = a_is_ref ? ba : bb;
    const Box& inc_box = a_is_ref ? bb : ba;
    const Mat3& ref_axes = a_is_ref ? axes_a : axes_b;
    const Mat3& inc_axes = a_is_ref ? axes_b : axes_a;
    Vec3 ref_normal = a_is_ref ? normal : -normal;

    Vec3 ref_face[4];
    Vec3 ref_n;
    get_face_vertices(ref_body, ref_box, ref_axes, ref_normal, ref_face, ref_n);
    ref_normal = ref_n;

    Vec3 inc_face[4];
    Vec3 inc_n;
    get_face_vertices(inc_body, inc_box, inc_axes, -ref_normal, inc_face, inc_n);

    // Side planes of reference face
    Vec3 clipped[16];
    int count = 4;
    Vec3 poly[16];
    for (int i = 0; i < 4; ++i) {
        poly[i] = inc_face[i];
    }

    for (int i = 0; i < 4; ++i) {
        const Vec3& v0 = ref_face[i];
        const Vec3& v1 = ref_face[(i + 1) % 4];
        const Vec3 edge = v1 - v0;
        const Vec3 plane_n = edge.cross(ref_normal).normalized();
        // Plane faces inward: keep points on the inside of the face polygon
        // plane_n should point inward toward face center
        const Vec3 face_center =
            (ref_face[0] + ref_face[1] + ref_face[2] + ref_face[3]) * 0.25f;
        Vec3 inward = plane_n;
        if (inward.dot(face_center - v0) < 0.f) {
            inward = -inward;
        }
        // Clip keeping n·x <= d where d = n·v0 for outward... use inward normal
        // Keep points with inward·(x - v0) >= 0  ⇒  -inward · x <= -inward·v0
        count = clip_poly(poly, count, -inward, -inward.dot(v0), clipped);
        if (count < 1) {
            return std::nullopt;
        }
        for (int k = 0; k < count; ++k) {
            poly[k] = clipped[k];
        }
    }

    // Keep points below reference plane (penetrating)
    const float plane_d = ref_normal.dot(ref_face[0]);
    m.point_count = 0;
    m.normal = a_is_ref ? ref_normal : -ref_normal;

    for (int i = 0; i < count && m.point_count < 4; ++i) {
        const float sep = ref_normal.dot(poly[i]) - plane_d;
        if (sep <= 0.f) {
            m.points[m.point_count].position = poly[i];
            m.points[m.point_count].penetration = -sep;
            ++m.point_count;
        }
    }

    if (m.point_count == 0) {
        // Fallback single contact
        m.point_count = 1;
        m.points[0].position = (a.position + b.position) * 0.5f;
        m.points[0].penetration = -best_sep;
    }

    finalize_manifold(m);
    return m;
}

// Closest point on segment ab to point p
Vec3 closest_on_segment(const Vec3& a, const Vec3& b, const Vec3& p) {
    const Vec3 ab = b - a;
    const float len_sq = ab.length_sq();
    if (len_sq < kEpsilon) {
        return a;
    }
    const float t = std::clamp((p - a).dot(ab) / len_sq, 0.f, 1.f);
    return a + ab * t;
}

std::optional<Manifold> sphere_capsule(const Body& sphere_body, const Body& cap_body,
                                       const Sphere& sphere, const Capsule& cap, bool swapped) {
    const Vec3 a = cap.axis_start(cap_body.position, cap_body.orientation);
    const Vec3 b = cap.axis_end(cap_body.position, cap_body.orientation);
    const Vec3 closest = closest_on_segment(a, b, sphere_body.position);
    const Vec3 delta = sphere_body.position - closest;
    const float dist_sq = delta.length_sq();
    const float radius_sum = sphere.radius + cap.radius;

    if (dist_sq > radius_sum * radius_sum) {
        return std::nullopt;
    }

    Manifold m;
    m.point_count = 1;

    if (dist_sq < kEpsilon) {
        // Degenerate — use capsule axis as normal fallback
        Vec3 axis = (b - a).normalized();
        if (axis.length_sq() < kEpsilon) {
            axis = {0.f, 1.f, 0.f};
        }
        const Vec3 n = axis;  // arbitrary
        m.points[0].position = closest + n * cap.radius;
        m.points[0].penetration = radius_sum;
        if (swapped) {
            m.a = cap_body.id;
            m.b = sphere_body.id;
            m.normal = n;
        } else {
            m.a = sphere_body.id;
            m.b = cap_body.id;
            m.normal = -n;
        }
        finalize_manifold(m);
        return m;
    }

    const float dist = std::sqrt(dist_sq);
    const Vec3 n = delta / dist;  // from capsule axis toward sphere

    m.points[0].penetration = radius_sum - dist;
    m.points[0].position = closest + n * cap.radius;

    if (swapped) {
        m.a = cap_body.id;
        m.b = sphere_body.id;
        m.normal = n;
    } else {
        m.a = sphere_body.id;
        m.b = cap_body.id;
        m.normal = -n;
    }
    finalize_manifold(m);
    return m;
}

std::optional<Manifold> capsule_capsule(const Body& a, const Body& b, const Capsule& ca,
                                        const Capsule& cb) {
    const Vec3 a0 = ca.axis_start(a.position, a.orientation);
    const Vec3 a1 = ca.axis_end(a.position, a.orientation);
    const Vec3 b0 = cb.axis_start(b.position, b.orientation);
    const Vec3 b1 = cb.axis_end(b.position, b.orientation);

    // Closest points between two segments
    const Vec3 d1 = a1 - a0;
    const Vec3 d2 = b1 - b0;
    const Vec3 r = a0 - b0;
    const float A = d1.dot(d1);
    const float E = d2.dot(d2);
    const float F = d2.dot(r);
    float s, t;
    if (A <= kEpsilon && E <= kEpsilon) {
        s = t = 0.f;
    } else if (A <= kEpsilon) {
        s = 0.f;
        t = std::clamp(F / E, 0.f, 1.f);
    } else {
        const float C = d1.dot(r);
        if (E <= kEpsilon) {
            t = 0.f;
            s = std::clamp(-C / A, 0.f, 1.f);
        } else {
            const float B = d1.dot(d2);
            const float denom = A * E - B * B;
            s = denom != 0.f ? std::clamp((B * F - C * E) / denom, 0.f, 1.f) : 0.f;
            t = (B * s + F) / E;
            if (t < 0.f) {
                t = 0.f;
                s = std::clamp(-C / A, 0.f, 1.f);
            } else if (t > 1.f) {
                t = 1.f;
                s = std::clamp((B - C) / A, 0.f, 1.f);
            }
        }
    }
    const Vec3 p1 = a0 + d1 * s;
    const Vec3 p2 = b0 + d2 * t;
    const Vec3 delta = p2 - p1;
    const float dist_sq = delta.length_sq();
    const float radius_sum = ca.radius + cb.radius;
    if (dist_sq > radius_sum * radius_sum) {
        return std::nullopt;
    }

    Manifold m;
    m.a = a.id;
    m.b = b.id;
    m.point_count = 1;
    if (dist_sq < kEpsilon) {
        m.normal = {0.f, 1.f, 0.f};
        m.points[0].penetration = radius_sum;
        m.points[0].position = p1;
    } else {
        const float dist = std::sqrt(dist_sq);
        m.normal = delta / dist;
        m.points[0].penetration = radius_sum - dist;
        m.points[0].position = p1 + m.normal * ca.radius;
    }
    finalize_manifold(m);
    return m;
}

std::optional<Manifold> capsule_box(const Body& cap_body, const Body& box_body,
                                    const Capsule& cap, const Box& box, bool swapped) {
    // Sample capsule as sphere at closest point on axis to box center (good approx)
    // Better: find closest point on capsule segment to box (iterate)
    const Vec3 a = cap.axis_start(cap_body.position, cap_body.orientation);
    const Vec3 b = cap.axis_end(cap_body.position, cap_body.orientation);

    const Mat3 rot = box_body.orientation.to_mat3();
    const Mat3 rot_t = rot.transposed();

    // Closest point on segment to box (approximate: closest to box center, then refine)
    Vec3 best_on_seg = closest_on_segment(a, b, box_body.position);
    // Project box center-relative and clamp for better contact
    for (int iter = 0; iter < 3; ++iter) {
        const Vec3 local = rot_t * (best_on_seg - box_body.position);
        const Vec3 closest_local{
            std::clamp(local.x, -box.half_extents.x, box.half_extents.x),
            std::clamp(local.y, -box.half_extents.y, box.half_extents.y),
            std::clamp(local.z, -box.half_extents.z, box.half_extents.z),
        };
        const Vec3 closest_world = box_body.position + rot * closest_local;
        best_on_seg = closest_on_segment(a, b, closest_world);
    }

    // Treat as sphere at best_on_seg
    Body pseudo = cap_body;
    pseudo.position = best_on_seg;
    Sphere sph{cap.radius};
    auto m = sphere_box(pseudo, box_body, sph, box, false);
    if (!m) {
        return std::nullopt;
    }
    // Fix body ids / normal for swap
    if (swapped) {
        // Want A=box, B=capsule — sphere_box(false) gave A=sphere(pseudo), B=box
        std::swap(m->a, m->b);
        m->normal = -m->normal;
        m->a = box_body.id;
        m->b = cap_body.id;
    } else {
        m->a = cap_body.id;
        m->b = box_body.id;
    }
    return m;
}

}  // namespace

std::optional<Manifold> collide(const Body& a, const Body& b) {
    return std::visit(
        [&](const auto& sa, const auto& sb) -> std::optional<Manifold> {
            using A = std::decay_t<decltype(sa)>;
            using B = std::decay_t<decltype(sb)>;

            if constexpr (std::is_same_v<A, Sphere> && std::is_same_v<B, Sphere>) {
                return sphere_sphere(a, b, sa, sb);
            } else if constexpr (std::is_same_v<A, Sphere> && std::is_same_v<B, Box>) {
                return sphere_box(a, b, sa, sb, false);
            } else if constexpr (std::is_same_v<A, Box> && std::is_same_v<B, Sphere>) {
                return sphere_box(b, a, sb, sa, true);
            } else if constexpr (std::is_same_v<A, Box> && std::is_same_v<B, Box>) {
                return box_box(a, b, sa, sb);
            } else if constexpr (std::is_same_v<A, Sphere> && std::is_same_v<B, Capsule>) {
                return sphere_capsule(a, b, sa, sb, false);
            } else if constexpr (std::is_same_v<A, Capsule> && std::is_same_v<B, Sphere>) {
                return sphere_capsule(b, a, sb, sa, true);
            } else if constexpr (std::is_same_v<A, Capsule> && std::is_same_v<B, Capsule>) {
                return capsule_capsule(a, b, sa, sb);
            } else if constexpr (std::is_same_v<A, Capsule> && std::is_same_v<B, Box>) {
                return capsule_box(a, b, sa, sb, false);
            } else if constexpr (std::is_same_v<A, Box> && std::is_same_v<B, Capsule>) {
                return capsule_box(b, a, sb, sa, true);
            }
            return std::nullopt;
        },
        a.shape, b.shape);
}

}  // namespace phys
