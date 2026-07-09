# phys — 3D Rigid Body Physics Engine

A from-scratch **3D rigid-body physics engine** in modern C++17, with an interactive sandbox.

## Features

| Area | Implementation |
|------|----------------|
| **Math** | `Vec3`, `Mat3`, `Quat` |
| **Shapes** | Sphere, oriented box (OBB), **capsule** |
| **Bodies** | Dynamic / static / kinematic, density → mass & inertia |
| **Collision** | Sphere/box/capsule pairs, SAT box–box, manifolds |
| **Solver** | Sequential impulses, warm start, dual-tangent friction |
| **Joints** | Distance, ball-socket, **mouse** (for grabbing) |
| **Sleeping** | Bodies rest and sleep; wake on contact / grab / force |
| **Broadphase** | Spatial hash grid |
| **Queries** | World raycast against all shapes |

## Build

```bash
cd Documents/code/physics-engine
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Tests only (no raylib):

```bash
cmake -B build -DPHYS_BUILD_DEMO=OFF && cmake --build build -j
```

## Run

```bash
./build/phys_tests
./build/phys_sandbox
```

### Sandbox controls

| Input | Action |
|-------|--------|
| **LMB drag** | Grab, move, and throw a body |
| **RMB** | Spawn sphere |
| **Shift+RMB** | Spawn box |
| **Ctrl+RMB** | Spawn capsule |
| **Mouse move** | Orbit camera (when not grabbing) |
| **Scroll** | Zoom |
| **Space** | Pause |
| **1** | Default mixed scene |
| **2** | Box pyramid |
| **3** | Sphere pile |
| **4** | Hanging chain (distance joints) |
| **5** | Pendulum (ball-socket) |
| **G** | Toggle gravity |

Dim bodies are **sleeping**. Grab them to wake and move.

## Quick API

```cpp
#include "phys/phys.hpp"

phys::World world;
world.settings().gravity = {0.f, -9.81f, 0.f};

// Grab with a mouse joint
phys::MouseJointDef grab;
grab.body = id;
grab.local_anchor = body->world_to_local(hit_point);
grab.target = hit_point;
grab.max_force = 1000.f;
auto jid = world.create_mouse_joint(grab);
world.set_mouse_joint_target(jid, new_target);
world.destroy_joint(jid);

// Raycast
auto hit = world.raycast({{0,5,0}, {0,-1,0}}, 100.f);

// Distance joint
phys::DistanceJointDef d;
d.body_a = a; d.body_b = b;
d.length = 1.5f;
world.create_distance_joint(d);
```

## Layout

```
include/phys/   public headers
src/            collision, raycast, world
tests/          unit tests
demos/          3D sandbox (raylib)
```

This project was made with ai for fun and learning purposes
