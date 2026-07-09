// Interactive 3D physics sandbox (raylib).
//
// Controls:
//   Left click + drag     — grab and move a body (throw on release)
//   Right click           — spawn sphere
//   Shift + right         — spawn box
//   Ctrl  + right         — spawn capsule
//   Middle mouse drag     — orbit camera
//   Alt + left drag       — orbit camera (alternate)
//   Scroll / +/-          — zoom in and out
//   Space                 — pause / resume
//   R / 1                 — default scene
//   2                     — box pyramid
//   3                     — sphere pile
//   4                     — hanging chain (distance joints)
//   5                     — pendulum (ball-socket)
//   B                     — place / move black hole under cursor
//   H                     — dismiss black hole
//   G                     — toggle gravity
//   Esc                   — quit

#include "phys/phys.hpp"

#include "raylib.h"
#include "raymath.h"
#include "rcamera.h"
#include "rlgl.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
#include <variant>

namespace {

constexpr int kScreenW = 1280;
constexpr int kScreenH = 720;
constexpr float kFixedDt = 1.f / 60.f;

// Playable arena (half-extents on X/Z). Walls keep pieces from falling off.
constexpr float kMapHalf = 12.f;
constexpr float kWallHeight = 2.5f;
constexpr float kWallThickness = 0.35f;
constexpr float kFloorThickness = 0.25f;
constexpr float kZoomMin = 2.5f;
constexpr float kZoomMax = 55.f;

Color body_color(const phys::Body& b, int index) {
    if (b.type == phys::BodyType::Static) {
        return Color{70, 72, 80, 255};
    }
    if (!b.awake) {
        return Color{90, 90, 110, 255};  // sleeping = dimmer
    }
    static const Color kPalette[] = {
        Color{80, 180, 255, 255},  Color{255, 120, 100, 255}, Color{120, 220, 140, 255},
        Color{255, 200, 80, 255},  Color{200, 140, 255, 255}, Color{100, 230, 220, 255},
    };
    return kPalette[index % 6];
}

Vector3 to_rl(phys::Vec3 v) { return {v.x, v.y, v.z}; }
phys::Vec3 to_phys(Vector3 v) { return {v.x, v.y, v.z}; }

Matrix quat_to_matrix(const phys::Quat& q) {
    const phys::Mat3 m = q.to_mat3();
    Matrix out = MatrixIdentity();
    out.m0 = m.cols[0].x;
    out.m1 = m.cols[0].y;
    out.m2 = m.cols[0].z;
    out.m4 = m.cols[1].x;
    out.m5 = m.cols[1].y;
    out.m6 = m.cols[1].z;
    out.m8 = m.cols[2].x;
    out.m9 = m.cols[2].y;
    out.m10 = m.cols[2].z;
    return out;
}

void draw_body(const phys::Body& b, int index, bool highlight) {
    Color fill = body_color(b, index);
    if (highlight) {
        fill = Color{255, 255, 120, 255};
    }
    const Color wire = Color{20, 20, 28, 255};

    if (const auto* sphere = std::get_if<phys::Sphere>(&b.shape)) {
        DrawSphere(to_rl(b.position), sphere->radius, fill);
        DrawSphereWires(to_rl(b.position), sphere->radius, 12, 12, wire);
    } else if (const auto* box = std::get_if<phys::Box>(&b.shape)) {
        const Vector3 size = {box->half_extents.x * 2.f, box->half_extents.y * 2.f,
                              box->half_extents.z * 2.f};
        rlPushMatrix();
        rlTranslatef(b.position.x, b.position.y, b.position.z);
        rlMultMatrixf(MatrixToFloat(quat_to_matrix(b.orientation)));
        DrawCube({0.f, 0.f, 0.f}, size.x, size.y, size.z, fill);
        DrawCubeWires({0.f, 0.f, 0.f}, size.x, size.y, size.z, wire);
        rlPopMatrix();
    } else if (const auto* cap = std::get_if<phys::Capsule>(&b.shape)) {
        rlPushMatrix();
        rlTranslatef(b.position.x, b.position.y, b.position.z);
        rlMultMatrixf(MatrixToFloat(quat_to_matrix(b.orientation)));
        // Cylinder along Y + end spheres
        DrawCylinderEx({0.f, -cap->half_height, 0.f}, {0.f, cap->half_height, 0.f}, cap->radius,
                       cap->radius, 12, fill);
        DrawSphere({0.f, -cap->half_height, 0.f}, cap->radius, fill);
        DrawSphere({0.f, cap->half_height, 0.f}, cap->radius, fill);
        DrawCylinderWiresEx({0.f, -cap->half_height, 0.f}, {0.f, cap->half_height, 0.f},
                            cap->radius, cap->radius, 8, wire);
        rlPopMatrix();
    }
}

void draw_contacts(const phys::World& world) {
    for (const phys::Manifold& m : world.contacts()) {
        for (int i = 0; i < m.point_count; ++i) {
            DrawSphere(to_rl(m.points[i].position), 0.05f, RED);
            DrawLine3D(to_rl(m.points[i].position),
                       to_rl(m.points[i].position + m.normal * 0.3f), ORANGE);
        }
    }
}

void draw_joints(const phys::World& world) {
    for (const phys::Joint& j : world.joints()) {
        if (!j.alive) {
            continue;
        }
        if (j.type == phys::JointType::Mouse) {
            const phys::Body* b = world.get_body(j.body_a);
            if (!b) {
                continue;
            }
            const phys::Vec3 anchor = b->local_to_world(j.local_anchor_a);
            DrawLine3D(to_rl(anchor), to_rl(j.target), YELLOW);
            DrawSphere(to_rl(j.target), 0.08f, YELLOW);
            DrawSphere(to_rl(anchor), 0.06f, GOLD);
            continue;
        }
        const phys::Body* a = world.get_body(j.body_a);
        const phys::Body* b = world.get_body(j.body_b);
        if (!a || !b) {
            continue;
        }
        const Vector3 wa = to_rl(a->local_to_world(j.local_anchor_a));
        const Vector3 wb = to_rl(b->local_to_world(j.local_anchor_b));
        DrawLine3D(wa, wb, Color{180, 220, 255, 255});
        DrawSphere(wa, 0.05f, SKYBLUE);
        DrawSphere(wb, 0.05f, SKYBLUE);
    }
}

void add_arena(phys::World& world) {
    // Floor
    phys::BodyDef ground;
    ground.type = phys::BodyType::Static;
    ground.shape = phys::Box{{kMapHalf, kFloorThickness, kMapHalf}};
    ground.position = {0.f, -kFloorThickness, 0.f};
    ground.friction = 0.6f;
    ground.restitution = 0.1f;
    world.create_body(ground);

    // Four walls forming a border around the map
    const float wall_y = kWallHeight * 0.5f;
    const float inset = kMapHalf + kWallThickness;  // outer face sits at map edge+

    auto add_wall = [&](phys::Vec3 half, phys::Vec3 pos) {
        phys::BodyDef wall;
        wall.type = phys::BodyType::Static;
        wall.shape = phys::Box{half};
        wall.position = pos;
        wall.friction = 0.45f;
        wall.restitution = 0.15f;
        world.create_body(wall);
    };

    // North / South (along X)
    add_wall({kMapHalf + kWallThickness * 2.f, kWallHeight * 0.5f, kWallThickness},
             {0.f, wall_y, inset});
    add_wall({kMapHalf + kWallThickness * 2.f, kWallHeight * 0.5f, kWallThickness},
             {0.f, wall_y, -inset});
    // East / West (along Z)
    add_wall({kWallThickness, kWallHeight * 0.5f, kMapHalf}, {inset, wall_y, 0.f});
    add_wall({kWallThickness, kWallHeight * 0.5f, kMapHalf}, {-inset, wall_y, 0.f});
}

void draw_map_border() {
    // Bright outline on the floor so the play area is obvious
    const float y = 0.02f;
    const float h = kMapHalf;
    const Color edge = Color{100, 200, 255, 220};
    const Color rim = Color{60, 90, 120, 180};

    DrawLine3D({-h, y, -h}, {h, y, -h}, edge);
    DrawLine3D({h, y, -h}, {h, y, h}, edge);
    DrawLine3D({h, y, h}, {-h, y, h}, edge);
    DrawLine3D({-h, y, h}, {-h, y, -h}, edge);

    // Top rim of walls
    const float top = kWallHeight;
    const float o = kMapHalf + kWallThickness;
    DrawLine3D({-o, top, -o}, {o, top, -o}, rim);
    DrawLine3D({o, top, -o}, {o, top, o}, rim);
    DrawLine3D({o, top, o}, {-o, top, o}, rim);
    DrawLine3D({-o, top, o}, {-o, top, -o}, rim);
}

void build_default_scene(phys::World& world) {
    world.clear();
    add_arena(world);

    for (int i = 0; i < 6; ++i) {
        phys::BodyDef def;
        def.type = phys::BodyType::Dynamic;
        def.position = {-3.f + i * 1.2f, 2.f + i * 0.3f, 0.f};
        def.restitution = 0.25f;
        def.friction = 0.4f;
        def.density = 1.f;
        def.orientation = phys::Quat::from_euler_xyz(0.1f * i, 0.15f * i, 0.05f * i);
        if (i % 3 == 0) {
            def.shape = phys::Sphere{0.35f + 0.05f * (i % 3)};
        } else if (i % 3 == 1) {
            def.shape = phys::Box{{0.35f, 0.35f, 0.35f}};
        } else {
            def.shape = phys::Capsule{0.2f, 0.35f};
        }
        world.create_body(def);
    }
}

void build_pyramid(phys::World& world) {
    world.clear();
    add_arena(world);

    const float size = 0.4f;
    const int rows = 6;
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < rows - row; ++col) {
            phys::BodyDef def;
            def.type = phys::BodyType::Dynamic;
            def.shape = phys::Box{{size, size, size}};
            def.position = {
                (col - (rows - row - 1) * 0.5f) * (size * 2.1f),
                size + row * (size * 2.1f),
                0.f,
            };
            def.friction = 0.5f;
            def.restitution = 0.05f;
            def.density = 1.f;
            world.create_body(def);
        }
    }
}

void build_sphere_pile(phys::World& world) {
    world.clear();
    add_arena(world);

    for (int y = 0; y < 4; ++y) {
        for (int x = 0; x < 4; ++x) {
            for (int z = 0; z < 4; ++z) {
                phys::BodyDef def;
                def.type = phys::BodyType::Dynamic;
                def.shape = phys::Sphere{0.3f};
                def.position = {
                    (x - 1.5f) * 0.65f,
                    0.35f + y * 0.65f,
                    (z - 1.5f) * 0.65f,
                };
                def.restitution = 0.2f;
                def.friction = 0.3f;
                def.density = 1.f;
                world.create_body(def);
            }
        }
    }
}

void build_chain(phys::World& world) {
    world.clear();
    add_arena(world);

    // Ceiling anchor (static)
    phys::BodyDef anchor;
    anchor.type = phys::BodyType::Static;
    anchor.shape = phys::Box{{0.35f, 0.2f, 0.35f}};
    anchor.position = {0.f, 7.5f, 0.f};
    const phys::BodyId anchor_id = world.create_body(anchor);

    constexpr float kLinkR = 0.22f;
    constexpr float kSpacing = 0.7f;  // center-to-center > 2*r so links don't stack-fight
    phys::BodyId prev = anchor_id;
    for (int i = 0; i < 7; ++i) {
        phys::BodyDef link;
        link.type = phys::BodyType::Dynamic;
        link.shape = phys::Sphere{kLinkR};
        link.position = {0.15f * static_cast<float>(i % 2), 7.0f - i * kSpacing, 0.f};
        link.density = 1.f;
        link.friction = 0.4f;
        link.restitution = 0.05f;
        link.allow_sleep = false;
        const phys::BodyId id = world.create_body(link);

        phys::DistanceJointDef jd;
        jd.body_a = prev;
        jd.body_b = id;
        jd.local_anchor_a = (i == 0) ? phys::Vec3{0.f, -0.2f, 0.f} : phys::Vec3{0.f, -kLinkR, 0.f};
        jd.local_anchor_b = {0.f, kLinkR, 0.f};
        jd.length = kSpacing;
        // Soft impulse joints: high k + damping holds under gravity
        jd.stiffness = 800.f;
        jd.damping = 40.f;
        world.create_distance_joint(jd);
        prev = id;
    }
}

void build_pendulum(phys::World& world) {
    world.clear();
    add_arena(world);

    phys::BodyDef pivot;
    pivot.type = phys::BodyType::Static;
    pivot.shape = phys::Box{{0.25f, 0.25f, 0.25f}};
    pivot.position = {0.f, 6.5f, 0.f};
    const phys::BodyId pivot_id = world.create_body(pivot);

    phys::BodyDef bob;
    bob.type = phys::BodyType::Dynamic;
    bob.shape = phys::Sphere{0.45f};
    bob.position = {2.5f, 6.5f, 0.f};  // offset → swings under gravity
    bob.density = 1.f;
    bob.friction = 0.2f;
    bob.restitution = 0.15f;
    bob.allow_sleep = false;
    const phys::BodyId bob_id = world.create_body(bob);

    phys::BallSocketJointDef jd;
    jd.body_a = pivot_id;
    jd.body_b = bob_id;
    jd.local_anchor_a = {0.f, 0.f, 0.f};
    if (phys::Body* b = world.get_body(bob_id)) {
        // Attach at bob surface toward pivot so rod length ≈ 2.5
        jd.local_anchor_b = b->world_to_local(pivot.position);
    }
    jd.stiffness = 1200.f;
    jd.damping = 50.f;
    world.create_ball_socket_joint(jd);
}

void spawn_at(phys::World& world, phys::Vec3 pos, int kind) {
    // Keep spawns inside the bordered map
    const float margin = 1.f;
    pos.x = std::clamp(pos.x, -kMapHalf + margin, kMapHalf - margin);
    pos.z = std::clamp(pos.z, -kMapHalf + margin, kMapHalf - margin);

    phys::BodyDef def;
    def.type = phys::BodyType::Dynamic;
    def.position = {pos.x, std::max(pos.y, 0.f) + 1.2f, pos.z};
    def.restitution = 0.25f;
    def.friction = 0.4f;
    def.density = 1.f;
    def.orientation = phys::Quat::from_euler_xyz(
        static_cast<float>(GetRandomValue(0, 62)) / 20.f,
        static_cast<float>(GetRandomValue(0, 62)) / 20.f,
        static_cast<float>(GetRandomValue(0, 62)) / 20.f);
    if (kind == 1) {
        const float h = 0.25f + static_cast<float>(GetRandomValue(0, 20)) / 50.f;
        def.shape = phys::Box{{h, h, h}};
    } else if (kind == 2) {
        def.shape = phys::Capsule{0.2f, 0.35f};
    } else {
        def.shape = phys::Sphere{0.25f + static_cast<float>(GetRandomValue(0, 20)) / 50.f};
    }
    world.create_body(def);
}

phys::Ray screen_ray(Camera3D camera) {
    const Ray r = GetScreenToWorldRay(GetMousePosition(), camera);
    return {to_phys(r.position), to_phys(r.direction).normalized()};
}

/// Target point at fixed depth along the mouse ray (stable while dragging).
phys::Vec3 mouse_at_depth(Camera3D camera, float depth) {
    const phys::Ray ray = screen_ray(camera);
    const float d = std::max(depth, 0.5f);
    return ray.origin + ray.direction * d;
}

phys::Vec3 ray_ground_hit(Camera3D camera) {
    const phys::Ray ray = screen_ray(camera);
    if (std::fabs(ray.direction.y) < 1e-6f) {
        return {0.f, 2.f, 0.f};
    }
    const float t = -ray.origin.y / ray.direction.y;
    if (t < 0.f) {
        return ray.origin + ray.direction * 5.f;
    }
    return ray.origin + ray.direction * t;
}

int count_alive(const phys::World& world) {
    int n = 0;
    for (const phys::Body& b : world.bodies()) {
        if (world.get_body(b.id)) {
            ++n;
        }
    }
    return n;
}

int count_sleeping(const phys::World& world) {
    int n = 0;
    for (const phys::Body& b : world.bodies()) {
        if (world.get_body(b.id) && b.type == phys::BodyType::Dynamic && !b.awake) {
            ++n;
        }
    }
    return n;
}

struct GrabState {
    phys::JointId joint = phys::kInvalidJoint;
    phys::BodyId body = phys::kInvalidBody;
    float grab_depth = 0.f;  // fixed distance along camera ray at grab time
    phys::Vec3 last_target{};
    phys::Vec3 prev_target{};
    phys::Vec3 smooth_target{};
    bool active = false;
};

struct BlackHole {
    bool active = false;
    phys::Vec3 position{0.f, 1.2f, 0.f};
    float base_horizon = 0.55f;
    float horizon = 0.55f;       // destroy radius (grows with meals)
    float influence = 14.f;      // max pull range (scales with size)
    float strength = 55.f;       // gravitational strength
    float swirl = 8.f;           // tangential spin force
    float visual_pulse = 0.f;
    int eaten = 0;
    float mass_eaten = 0.f;

    void recompute_size() {
        // Keep growing forever with diminishing returns (no hard stop around ~40 meals).
        const float meals = static_cast<float>(std::max(eaten, 0));
        const float from_count = 0.55f * std::log1p(meals * 0.35f);
        const float from_mass = 0.25f * std::log1p(mass_eaten * 0.15f);
        horizon = base_horizon + from_count + from_mass;

        constexpr float kSoftMax = 8.f;
        if (horizon > kSoftMax) {
            horizon = kSoftMax;
        }

        position.y = std::max(1.2f, horizon * 1.05f);

        // Gravity / suction scale strongly with size (≈ mass ~ r³ feel, soft-capped)
        const float s = horizon / base_horizon;  // 1 at birth, ~3–4 when large
        const float s2 = s * s;
        const float s3 = s2 * s;

        // Pull range grows with size so a bigger hole dominates more of the arena
        influence = 9.f * s + horizon * 5.f;
        influence = std::min(influence, kMapHalf * 2.2f);

        // Strength grows super-linearly — bigger hole = much stronger suction
        strength = 35.f * s2 + 20.f * s3 + mass_eaten * 5.f + meals * 1.5f;
        swirl = 5.f * s2 + 3.f * s3;
    }
};

void draw_black_hole(const BlackHole& bh, float time) {
    if (!bh.active) {
        return;
    }

    const Vector3 p = to_rl(bh.position);
    const float pulse = 1.f + 0.06f * std::sin(time * 4.f + bh.visual_pulse);
    const float core_r = bh.horizon * 0.9f * pulse;
    const float disc_r = bh.horizon * 2.6f;

    // Accretion disc (flat rings) — scales with horizon
    for (int i = 0; i < 5; ++i) {
        const float t = static_cast<float>(i) / 4.f;
        const float r = disc_r * (0.55f + 0.45f * t);
        const unsigned char a = static_cast<unsigned char>(200 - i * 35);
        Color c = Color{
            static_cast<unsigned char>(255),
            static_cast<unsigned char>(120 + i * 20),
            static_cast<unsigned char>(40 + i * 10),
            a,
        };
        const float spin = time * (1.5f + i * 0.3f);
        DrawCircle3D(p, r, {0.f, 1.f, 0.f}, 90.f, c);
        for (int k = 0; k < 3; ++k) {
            const float ang = spin + k * 2.094f + i * 0.4f;
            const Vector3 spot = {
                p.x + r * std::cos(ang),
                p.y,
                p.z + r * std::sin(ang),
            };
            // Keep hotspots small so they don't look like uneaten physics balls
            DrawSphere(spot, 0.06f * (1.f - t * 0.4f), Color{255, 200, 80, 160});
        }
    }

    // Event horizon (near-black core)
    DrawSphere(p, core_r, Color{5, 5, 8, 255});
    DrawSphereWires(p, core_r * 1.05f, 12, 12, Color{40, 20, 60, 200});

    // Photon-sphere glow
    DrawSphereWires(p, core_r * 1.35f, 10, 10, Color{120, 60, 180, 90});
    DrawSphereWires(p, std::min(bh.influence * 0.15f, bh.horizon * 2.5f), 8, 8,
                    Color{80, 40, 120, 40});

    // Pull indicator rays (inward)
    for (int i = 0; i < 12; ++i) {
        const float ang = time * 0.8f + i * (6.283185f / 12.f);
        const float outer = bh.horizon * 2.4f;
        const Vector3 a = {
            p.x + outer * std::cos(ang),
            p.y + 0.1f * bh.horizon * std::sin(time * 3.f + i),
            p.z + outer * std::sin(ang),
        };
        DrawLine3D(a, p, Color{160, 80, 220, 70});
    }
}

/// Apply attraction and swallow dynamic bodies. Returns number eaten this tick.
int update_black_hole(BlackHole& bh, phys::World& world, GrabState& grab, float dt) {
    if (!bh.active) {
        return 0;
    }

    bh.visual_pulse += dt;
    int eaten_now = 0;

    bh.recompute_size();

    // Consume zone larger than the visual core so objects actually vanish when "eaten"
    // (previously large holes left balls orbiting inside the disc, outside a small horizon)
    const float consume_r = std::max(bh.horizon * 1.4f, bh.horizon + 0.35f);
    const float disc_r = bh.horizon * 2.2f;
    const float size_scale = std::max(1.f, bh.horizon / bh.base_horizon);

    std::vector<phys::BodyId> to_destroy;
    to_destroy.reserve(16);

    for (const phys::Body& body_ref : world.bodies()) {
        phys::Body* body = world.get_body(body_ref.id);
        if (!body || body->type != phys::BodyType::Dynamic) {
            continue;
        }

        // Use body id from slot, not stale copy
        const phys::BodyId id = body->id;
        phys::Vec3 delta = bh.position - body->position;
        float dist = delta.length();
        if (dist < 1e-4f) {
            dist = 1e-4f;
            delta = {0.f, 1.f, 0.f};
        }

        // Swallowed — use expanded consume radius so large holes still eat cleanly
        if (dist < consume_r) {
            to_destroy.push_back(id);
            continue;
        }

        if (dist > bh.influence) {
            continue;
        }

        const phys::Vec3 dir = delta / dist;

        // Softening stays modest so pull keeps rising near the hole (not flat for large r)
        const float soft_eps = 0.12f + 0.06f * bh.horizon * bh.horizon;
        const float soft = dist * dist + soft_eps;
        float pull = bh.strength * body->mass / soft;

        // Extra vacuum once inside the accretion disc
        if (dist < disc_r) {
            pull *= 2.5f;
        }

        const float edge = 1.f - (dist / bh.influence);
        const float falloff = edge * edge;

        body->apply_force(dir * (pull * falloff));

        // Swirl only in the outer region — near the hole it caused permanent orbits
        // so balls never crossed the event horizon once the hole got big.
        const float swirl_fade =
            std::clamp((dist - consume_r) / std::max(bh.horizon, 0.3f), 0.f, 1.f);
        if (swirl_fade > 0.05f) {
            phys::Vec3 up{0.f, 1.f, 0.f};
            phys::Vec3 tangential = up.cross(dir);
            if (tangential.length_sq() < 1e-6f) {
                tangential = phys::Vec3{1.f, 0.f, 0.f}.cross(dir);
            }
            tangential = tangential.normalized();
            body->apply_force(tangential * (bh.swirl * 0.65f * body->mass * falloff *
                                            swirl_fade / (dist + 0.5f)));
        }

        // Kill tangential velocity near the hole so objects fall in instead of orbiting
        if (dist < disc_r) {
            const float radial_speed = body->linear_velocity.dot(dir);
            phys::Vec3 radial_v = dir * radial_speed;
            phys::Vec3 tang_v = body->linear_velocity - radial_v;
            const float kill = std::clamp(1.5f * dt * size_scale, 0.f, 0.85f);
            body->linear_velocity = radial_v + tang_v * (1.f - kill);
            // Ensure net inward motion
            if (radial_speed < 2.f * size_scale) {
                body->linear_velocity += dir * ((4.f + 3.f * size_scale) * dt);
            }
        }

        // General radial vacuum
        const float vacuum = (6.f + 5.f * size_scale) * falloff;
        body->linear_velocity += dir * (vacuum * dt);
        body->angular_velocity *= (1.f - std::min(0.9f, 0.4f * size_scale * falloff * dt));
        body->disturb();
    }

    for (phys::BodyId id : to_destroy) {
        if (grab.active && grab.body == id) {
            world.destroy_joint(grab.joint);
            grab = {};
        }
        if (const phys::Body* b = world.get_body(id)) {
            bh.mass_eaten += b->mass;
            world.destroy_body(id);
            ++bh.eaten;
            ++eaten_now;
        }
    }

    return eaten_now;
}

}  // namespace

int main() {
    SetConfigFlags(FLAG_MSAA_4X_HINT | FLAG_VSYNC_HINT);
    InitWindow(kScreenW, kScreenH, "phys — 3D Physics Sandbox");
    SetTargetFPS(60);

    Camera3D camera = {};
    camera.position = {8.f, 6.f, 10.f};
    camera.target = {0.f, 1.f, 0.f};
    camera.up = {0.f, 1.f, 0.f};
    camera.fovy = 55.f;
    camera.projection = CAMERA_PERSPECTIVE;

    phys::World world;
    world.settings().gravity = {0.f, -9.81f, 0.f};
    world.settings().velocity_iterations = 12;
    world.settings().position_iterations = 4;
    world.settings().allow_sleep = true;
    build_default_scene(world);

    bool paused = false;
    float accumulator = 0.f;
    bool gravity_on = true;
    GrabState grab;
    BlackHole black_hole;
    float anim_time = 0.f;

    // Fixed starting view — camera only moves while orbit bind is held
    float yaw = 0.7f;
    float pitch = 0.5f;
    float dist = 14.f;
    const auto apply_camera = [&]() {
        camera.position = {
            camera.target.x + dist * std::cos(pitch) * std::sin(yaw),
            camera.target.y + dist * std::sin(pitch),
            camera.target.z + dist * std::cos(pitch) * std::cos(yaw),
        };
    };
    apply_camera();

    while (!WindowShouldClose()) {
        const Vector2 mouse_delta = GetMouseDelta();
        const bool alt_down = IsKeyDown(KEY_LEFT_ALT) || IsKeyDown(KEY_RIGHT_ALT);
        const bool orbiting =
            IsMouseButtonDown(MOUSE_BUTTON_MIDDLE) ||
            (alt_down && IsMouseButtonDown(MOUSE_BUTTON_LEFT));

        // --- Camera orbit (MMB drag, or Alt+LMB) — stays fixed otherwise ---
        if (orbiting) {
            yaw += mouse_delta.x * 0.005f;
            pitch += mouse_delta.y * 0.005f;
            pitch = Clamp(pitch, 0.08f, 1.45f);
        }

        // Zoom: mouse wheel and +/- keys
        const float wheel = GetMouseWheelMove();
        if (wheel != 0.f) {
            dist *= std::pow(0.88f, wheel);
        }
        if (IsKeyDown(KEY_EQUAL) || IsKeyDown(KEY_KP_ADD) || IsKeyDown(KEY_RIGHT_BRACKET)) {
            dist *= 0.97f;
        }
        if (IsKeyDown(KEY_MINUS) || IsKeyDown(KEY_KP_SUBTRACT) || IsKeyDown(KEY_LEFT_BRACKET)) {
            dist *= 1.03f;
        }
        dist = Clamp(dist, kZoomMin, kZoomMax);
        apply_camera();

        // --- Grab / drag (disabled while orbiting with Alt+LMB) ---
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && !grab.active && !alt_down && !orbiting) {
            const phys::Ray ray = screen_ray(camera);
            if (auto hit = world.raycast(ray, 200.f, true)) {
                phys::Body* body = world.get_body(hit->body);
                if (body && body->type == phys::BodyType::Dynamic) {
                    phys::MouseJointDef md;
                    md.body = hit->body;
                    md.local_anchor = body->world_to_local(hit->point);
                    md.target = hit->point;
                    // Generous force budget; soft constraint prevents explosion
                    md.max_force = std::max(2000.f, body->mass * 800.f);
                    md.frequency_hz = 7.f;
                    md.damping_ratio = 1.0f;   // critically damped — no bounce
                    md.grab_point_blend = 0.3f; // mostly translate COM, less spin
                    grab.joint = world.create_mouse_joint(md);
                    grab.body = hit->body;
                    grab.grab_depth = hit->distance;
                    grab.last_target = hit->point;
                    grab.prev_target = hit->point;
                    grab.smooth_target = hit->point;
                    grab.active = grab.joint != phys::kInvalidJoint;
                }
            }
        }

        if (grab.active && IsMouseButtonDown(MOUSE_BUTTON_LEFT) && !alt_down) {
            // Fixed-depth ray targeting (stable) + light smoothing (kills pixel jitter)
            const phys::Vec3 raw = mouse_at_depth(camera, grab.grab_depth);
            constexpr float kSmooth = 0.4f;
            grab.smooth_target = grab.smooth_target + (raw - grab.smooth_target) * kSmooth;
            grab.prev_target = grab.last_target;
            grab.last_target = grab.smooth_target;
            world.set_mouse_joint_target(grab.joint, grab.smooth_target);
        }

        if (grab.active && (IsMouseButtonReleased(MOUSE_BUTTON_LEFT) || alt_down)) {
            if (phys::Body* body = world.get_body(grab.body)) {
                // Throw from smoothed motion; clamp so release doesn't fling wildly
                phys::Vec3 throw_v = (grab.last_target - grab.prev_target) * (1.f / std::max(GetFrameTime(), 1e-4f));
                const float speed = throw_v.length();
                constexpr float kMaxThrow = 18.f;
                if (speed > kMaxThrow) {
                    throw_v = throw_v * (kMaxThrow / speed);
                }
                body->linear_velocity += throw_v * 0.45f;
                body->angular_velocity *= 0.5f;
                body->wake();
            }
            world.destroy_joint(grab.joint);
            grab = {};
        }

        // --- Keyboard ---
        if (IsKeyPressed(KEY_SPACE)) {
            paused = !paused;
        }
        if (IsKeyPressed(KEY_R) || IsKeyPressed(KEY_ONE)) {
            if (grab.active) {
                world.destroy_joint(grab.joint);
                grab = {};
            }
            build_default_scene(world);
        }
        if (IsKeyPressed(KEY_TWO)) {
            if (grab.active) {
                world.destroy_joint(grab.joint);
                grab = {};
            }
            build_pyramid(world);
        }
        if (IsKeyPressed(KEY_THREE)) {
            if (grab.active) {
                world.destroy_joint(grab.joint);
                grab = {};
            }
            build_sphere_pile(world);
        }
        if (IsKeyPressed(KEY_FOUR)) {
            if (grab.active) {
                world.destroy_joint(grab.joint);
                grab = {};
            }
            build_chain(world);
        }
        if (IsKeyPressed(KEY_FIVE)) {
            if (grab.active) {
                world.destroy_joint(grab.joint);
                grab = {};
            }
            build_pendulum(world);
        }
        if (IsKeyPressed(KEY_G)) {
            gravity_on = !gravity_on;
            world.settings().gravity =
                gravity_on ? phys::Vec3{0.f, -9.81f, 0.f} : phys::Vec3{0.f, 0.f, 0.f};
            // Wake everything when gravity toggles
            for (const phys::Body& b : world.bodies()) {
                if (phys::Body* p = world.get_body(b.id)) {
                    p->wake();
                }
            }
        }

        // Black hole: B places/moves under cursor, H dismisses
        if (IsKeyPressed(KEY_B)) {
            phys::Vec3 pos = ray_ground_hit(camera);
            pos.x = std::clamp(pos.x, -kMapHalf + 1.5f, kMapHalf - 1.5f);
            pos.z = std::clamp(pos.z, -kMapHalf + 1.5f, kMapHalf - 1.5f);
            pos.y = 1.2f;
            black_hole.active = true;
            black_hole.position = pos;
            black_hole.recompute_size();
            // Wake everything so sleeping piles get pulled
            for (const phys::Body& b : world.bodies()) {
                if (phys::Body* p = world.get_body(b.id)) {
                    p->disturb();
                }
            }
        }
        if (IsKeyPressed(KEY_H)) {
            black_hole = {};
        }

        // Reset black hole meal stats when loading a new scene (keep active if on)
        // Spawn on right click (not while grabbing)
        if (!grab.active && IsMouseButtonPressed(MOUSE_BUTTON_RIGHT)) {
            int kind = 0;
            if (IsKeyDown(KEY_LEFT_SHIFT) || IsKeyDown(KEY_RIGHT_SHIFT)) {
                kind = 1;
            } else if (IsKeyDown(KEY_LEFT_CONTROL) || IsKeyDown(KEY_RIGHT_CONTROL)) {
                kind = 2;
            }
            spawn_at(world, ray_ground_hit(camera), kind);
        }

        const float frame_dt = GetFrameTime();
        anim_time += frame_dt;

        if (!paused) {
            accumulator += frame_dt;
            if (accumulator > 0.25f) {
                accumulator = 0.25f;
            }
            while (accumulator >= kFixedDt) {
                update_black_hole(black_hole, world, grab, kFixedDt);
                world.step(kFixedDt);
                accumulator -= kFixedDt;
            }
        }

        BeginDrawing();
        ClearBackground(Color{16, 16, 22, 255});

        BeginMode3D(camera);
        DrawGrid(static_cast<int>(kMapHalf * 2.f), 1.f);
        draw_map_border();

        int index = 0;
        for (const phys::Body& body : world.bodies()) {
            if (!world.get_body(body.id)) {
                continue;
            }
            const bool hl = grab.active && body.id == grab.body;
            draw_body(body, index++, hl);
        }
        draw_contacts(world);
        draw_joints(world);
        draw_black_hole(black_hole, anim_time);
        EndMode3D();

        DrawRectangle(0, 0, kScreenW, 86, Color{10, 10, 14, 220});
        DrawText("phys 3D sandbox", 16, 10, 22, RAYWHITE);
        DrawText(
            "LMB grab | MMB/Alt+LMB orbit | Scroll zoom | RMB spawn | B black hole | H dismiss | 1-5 | Space | G",
            16, 38, 14, Color{170, 170, 190, 255});
        DrawText("Black hole sucks dynamic objects and deletes them at the event horizon.", 16, 58,
                 14, Color{140, 140, 160, 255});

        char stats[256];
        if (black_hole.active) {
            std::snprintf(stats, sizeof(stats),
                          "bodies: %d  ate: %d  size: %.1f  pull: %.0f  zoom: %.1f  %s%s",
                          count_alive(world), black_hole.eaten, black_hole.horizon,
                          black_hole.strength, dist, paused ? "PAUSED" : "RUNNING",
                          grab.active ? "  GRAB" : "");
        } else {
            std::snprintf(stats, sizeof(stats),
                          "bodies: %d  sleep: %d  zoom: %.1f  %s  g=%s%s",
                          count_alive(world), count_sleeping(world), dist,
                          paused ? "PAUSED" : "RUNNING", gravity_on ? "on" : "off",
                          grab.active ? "  GRAB" : "");
        }
        DrawText(stats, kScreenW - 480, 28, 16, Color{150, 220, 150, 255});

        // Crosshair
        DrawCircleLines(kScreenW / 2, kScreenH / 2, 6, Color{255, 255, 255, 80});

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
