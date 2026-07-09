#pragma once

#include "phys/body.hpp"
#include "phys/vec3.hpp"

#include <cstdint>
#include <vector>

namespace phys {

using JointId = std::uint32_t;
constexpr JointId kInvalidJoint = static_cast<JointId>(-1);

enum class JointType {
    Distance,
    BallSocket,
    Mouse,
};

struct DistanceJointDef {
    BodyId body_a = kInvalidBody;
    BodyId body_b = kInvalidBody;
    Vec3 local_anchor_a{};  // in body A local space
    Vec3 local_anchor_b{};
    float length = -1.f;  // < 0 → use current distance at creation
    float stiffness = 50.f;
    float damping = 5.f;
};

struct BallSocketJointDef {
    BodyId body_a = kInvalidBody;
    BodyId body_b = kInvalidBody;
    Vec3 local_anchor_a{};
    Vec3 local_anchor_b{};
    float stiffness = 200.f;
    float damping = 10.f;
};

struct MouseJointDef {
    BodyId body = kInvalidBody;
    Vec3 local_anchor{};   // grab point in body local space
    Vec3 target{};         // world-space target
    float max_force = 1000.f;
    /// Response speed in Hz (higher = snappier). Soft constraint, not spring k.
    float frequency_hz = 6.f;
    /// 1 = critically damped (stable). <1 oscillates, >1 overdamps.
    float damping_ratio = 1.0f;
    /// Blend of impulse at grab point vs center of mass (1 = all at point).
    /// Lower values reduce wild spinning when grabbing off-center.
    float grab_point_blend = 0.35f;
};

struct Joint {
    JointId id = kInvalidJoint;
    JointType type = JointType::Distance;
    BodyId body_a = kInvalidBody;
    BodyId body_b = kInvalidBody;  // unused for Mouse
    Vec3 local_anchor_a{};
    Vec3 local_anchor_b{};
    Vec3 target{};  // Mouse joint world target
    float length = 0.f;
    float stiffness = 50.f;       // distance / ball-socket spring k
    float damping = 5.f;          // distance / ball-socket damper
    float frequency_hz = 6.f;     // mouse joint
    float damping_ratio = 1.f;    // mouse joint
    float grab_point_blend = 0.35f;
    float max_force = 1e9f;
    float max_impulse_budget = 0.f;  // remaining impulse this step (mouse)
    bool alive = true;
};

}  // namespace phys
