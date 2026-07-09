#include "phys/phys.hpp"

#include <cmath>
#include <iostream>

namespace {

int g_failed = 0;
int g_passed = 0;

#define EXPECT_TRUE(cond)                                                                              \
    do {                                                                                               \
        if (cond) {                                                                                    \
            ++g_passed;                                                                                \
        } else {                                                                                       \
            ++g_failed;                                                                                \
            std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << "  " << #cond << '\n';              \
        }                                                                                              \
    } while (0)

#define EXPECT_NEAR(a, b, tol)                                                                         \
    do {                                                                                               \
        const auto _a = (a);                                                                           \
        const auto _b = (b);                                                                           \
        if (std::fabs(_a - _b) <= (tol)) {                                                             \
            ++g_passed;                                                                                \
        } else {                                                                                       \
            ++g_failed;                                                                                \
            std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << "  " << #a << " ≈ " << #b            \
                      << "  got " << _a << " vs " << _b << '\n';                                       \
        }                                                                                              \
    } while (0)

void test_vec3() {
    using phys::Vec3;
    Vec3 a{1.f, 2.f, 2.f};
    EXPECT_NEAR(a.length(), 3.f, 1e-5f);
    EXPECT_NEAR(a.normalized().length(), 1.f, 1e-5f);
    const float cross_z = Vec3(1.f, 0.f, 0.f).cross(Vec3(0.f, 1.f, 0.f)).z;
    EXPECT_NEAR(cross_z, 1.f, 1e-5f);
}

void test_quat_rotate() {
    using phys::Quat;
    using phys::Vec3;
    const Quat q = Quat::from_axis_angle({0.f, 1.f, 0.f}, 3.14159265f * 0.5f);
    const Vec3 r = q.rotate({1.f, 0.f, 0.f});
    EXPECT_NEAR(r.x, 0.f, 1e-4f);
    EXPECT_NEAR(r.y, 0.f, 1e-4f);
    EXPECT_NEAR(r.z, -1.f, 1e-4f);
}

void test_sphere_sphere() {
    phys::Body a;
    a.id = 0;
    a.shape = phys::Sphere{1.f};
    a.position = {0.f, 0.f, 0.f};

    phys::Body b;
    b.id = 1;
    b.shape = phys::Sphere{1.f};
    b.position = {1.5f, 0.f, 0.f};

    auto m = phys::collide(a, b);
    EXPECT_TRUE(m.has_value());
    EXPECT_NEAR(m->points[0].penetration, 0.5f, 1e-4f);
}

void test_capsule_sphere() {
    phys::Body cap;
    cap.id = 0;
    cap.shape = phys::Capsule{0.5f, 1.f};
    cap.position = {0.f, 0.f, 0.f};

    phys::Body sph;
    sph.id = 1;
    sph.shape = phys::Sphere{0.5f};
    sph.position = {0.8f, 0.f, 0.f};  // radii 0.5+0.5 → overlap of 0.2

    auto m = phys::collide(cap, sph);
    EXPECT_TRUE(m.has_value());
    EXPECT_TRUE(m->points[0].penetration > 0.1f);
}

void test_ball_bounce() {
    phys::World world;
    world.settings().gravity = {0.f, -10.f, 0.f};
    world.settings().restitution_threshold = 0.2f;
    world.settings().allow_sleep = false;

    phys::BodyDef ground;
    ground.type = phys::BodyType::Static;
    ground.shape = phys::Box{{20.f, 0.5f, 20.f}};
    ground.position = {0.f, -0.5f, 0.f};
    ground.restitution = 0.8f;
    world.create_body(ground);

    phys::BodyDef ball;
    ball.shape = phys::Sphere{0.5f};
    ball.position = {0.f, 5.f, 0.f};
    ball.restitution = 0.8f;
    ball.density = 1.f;
    const phys::BodyId ball_id = world.create_body(ball);

    const float dt = 1.f / 60.f;
    float max_h = -1e9f;
    bool hit = false;
    for (int i = 0; i < 180; ++i) {
        world.step(dt);
        const phys::Body* b = world.get_body(ball_id);
        if (b->position.y < 1.f) {
            hit = true;
        }
        if (hit && b->linear_velocity.y > 0.f) {
            max_h = std::max(max_h, b->position.y);
        }
    }
    EXPECT_TRUE(hit);
    EXPECT_TRUE(max_h > 1.5f);
}

void test_sleeping() {
    phys::World world;
    world.settings().gravity = {0.f, -10.f, 0.f};
    world.settings().allow_sleep = true;
    world.settings().sleep_time = 0.25f;
    world.settings().sleep_linear_threshold = 0.25f;
    world.settings().sleep_angular_threshold = 0.35f;
    world.settings().velocity_iterations = 16;
    world.settings().position_iterations = 6;

    phys::BodyDef ground;
    ground.type = phys::BodyType::Static;
    ground.shape = phys::Box{{20.f, 0.5f, 20.f}};
    ground.position = {0.f, -0.5f, 0.f};
    ground.restitution = 0.f;
    ground.friction = 0.9f;
    world.create_body(ground);

    phys::BodyDef ball;
    ball.shape = phys::Sphere{0.5f};  // spheres rest more stably than boxes here
    ball.position = {0.f, 0.5f, 0.f};
    ball.restitution = 0.f;
    ball.friction = 0.9f;
    ball.density = 1.f;
    const phys::BodyId id = world.create_body(ball);

    const float dt = 1.f / 60.f;
    for (int i = 0; i < 480; ++i) {
        world.step(dt);
    }

    const phys::Body* b = world.get_body(id);
    EXPECT_TRUE(!b->awake);
    EXPECT_NEAR(b->linear_velocity.length(), 0.f, 1e-5f);
}

void test_distance_joint() {
    phys::World world;
    world.settings().gravity = {0.f, -10.f, 0.f};
    world.settings().allow_sleep = false;

    phys::BodyDef anchor;
    anchor.type = phys::BodyType::Static;
    anchor.shape = phys::Box{{0.2f, 0.2f, 0.2f}};
    anchor.position = {0.f, 5.f, 0.f};
    const phys::BodyId a = world.create_body(anchor);

    phys::BodyDef bob;
    bob.shape = phys::Sphere{0.3f};
    bob.position = {0.f, 3.f, 0.f};
    bob.density = 1.f;
    bob.allow_sleep = false;
    const phys::BodyId b = world.create_body(bob);

    phys::DistanceJointDef jd;
    jd.body_a = a;
    jd.body_b = b;
    jd.local_anchor_a = {};
    jd.local_anchor_b = {};
    jd.length = 2.f;
    jd.stiffness = 800.f;
    jd.damping = 40.f;
    world.create_distance_joint(jd);

    const float dt = 1.f / 60.f;
    for (int i = 0; i < 240; ++i) {
        world.step(dt);
    }

    const phys::Body* bob_b = world.get_body(b);
    const float dist = (bob_b->position - phys::Vec3{0.f, 5.f, 0.f}).length();
    EXPECT_NEAR(dist, 2.f, 0.45f);
    EXPECT_TRUE(bob_b->position.y < 5.f);
    EXPECT_TRUE(bob_b->position.y > 2.2f);  // not collapsed to floor
}

void test_mouse_joint_moves_body() {
    phys::World world;
    world.settings().gravity = {0.f, 0.f, 0.f};
    world.settings().allow_sleep = false;

    phys::BodyDef def;
    def.shape = phys::Sphere{0.5f};
    def.position = {0.f, 0.f, 0.f};
    def.density = 1.f;
    const phys::BodyId id = world.create_body(def);

    phys::MouseJointDef md;
    md.body = id;
    md.local_anchor = {};
    md.target = {2.f, 0.f, 0.f};
    md.max_force = 5000.f;
    md.frequency_hz = 8.f;
    md.damping_ratio = 1.f;
    world.create_mouse_joint(md);

    const float dt = 1.f / 60.f;
    for (int i = 0; i < 180; ++i) {
        world.step(dt);
    }

    const phys::Body* b = world.get_body(id);
    EXPECT_TRUE(b->position.x > 1.0f);
}

void test_raycast_sphere() {
    phys::World world;
    phys::BodyDef def;
    def.shape = phys::Sphere{1.f};
    def.position = {0.f, 0.f, 5.f};
    def.type = phys::BodyType::Static;
    world.create_body(def);

    phys::Ray ray{{0.f, 0.f, 0.f}, {0.f, 0.f, 1.f}};
    auto hit = world.raycast(ray, 100.f);
    EXPECT_TRUE(hit.has_value());
    EXPECT_NEAR(hit->distance, 4.f, 0.05f);
}

void test_momentum_transfer() {
    phys::World world;
    world.settings().gravity = {0.f, 0.f, 0.f};
    world.settings().allow_sleep = false;

    phys::BodyDef a_def;
    a_def.shape = phys::Sphere{0.5f};
    a_def.position = {0.f, 0.f, 0.f};
    a_def.linear_velocity = {5.f, 0.f, 0.f};
    a_def.restitution = 1.f;
    a_def.friction = 0.f;
    a_def.density = 1.f;
    const phys::BodyId a = world.create_body(a_def);

    phys::BodyDef b_def;
    b_def.shape = phys::Sphere{0.5f};
    b_def.position = {2.f, 0.f, 0.f};
    b_def.restitution = 1.f;
    b_def.friction = 0.f;
    b_def.density = 1.f;
    const phys::BodyId b = world.create_body(b_def);

    const float dt = 1.f / 120.f;
    for (int i = 0; i < 120; ++i) {
        world.step(dt);
    }

    EXPECT_NEAR(world.get_body(a)->linear_velocity.x, 0.f, 1.5f);
    EXPECT_NEAR(world.get_body(b)->linear_velocity.x, 5.f, 1.5f);
}

}  // namespace

int main() {
    test_vec3();
    test_quat_rotate();
    test_sphere_sphere();
    test_capsule_sphere();
    test_ball_bounce();
    test_sleeping();
    test_distance_joint();
    test_mouse_joint_moves_body();
    test_raycast_sphere();
    test_momentum_transfer();

    std::cout << "Passed: " << g_passed << "  Failed: " << g_failed << '\n';
    return g_failed == 0 ? 0 : 1;
}
