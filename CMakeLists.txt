cmake_minimum_required(VERSION 3.16)
project(phys VERSION 0.2.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(PHYS_BUILD_TESTS "Build unit tests" ON)
option(PHYS_BUILD_DEMO "Build interactive sandbox demo (downloads raylib)" ON)

# ---------------------------------------------------------------------------
# Core library
# ---------------------------------------------------------------------------
add_library(phys
    src/collision.cpp
    src/raycast.cpp
    src/world.cpp
)
add_library(phys::phys ALIAS phys)

target_include_directories(phys
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

if(MSVC)
    target_compile_options(phys PRIVATE /W4)
else()
    target_compile_options(phys PRIVATE -Wall -Wextra -Wpedantic)
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(PHYS_BUILD_TESTS)
    add_executable(phys_tests tests/test_physics.cpp)
    target_link_libraries(phys_tests PRIVATE phys::phys)
    enable_testing()
    add_test(NAME phys_tests COMMAND phys_tests)
endif()

# ---------------------------------------------------------------------------
# Interactive demo (raylib via FetchContent)
# ---------------------------------------------------------------------------
if(PHYS_BUILD_DEMO)
    include(FetchContent)
    set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(BUILD_GAMES OFF CACHE BOOL "" FORCE)

    FetchContent_Declare(
        raylib
        GIT_REPOSITORY https://github.com/raysan5/raylib.git
        GIT_TAG 5.5
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(raylib)

    add_executable(phys_sandbox demos/sandbox.cpp)
    target_link_libraries(phys_sandbox PRIVATE phys::phys raylib)

    if(APPLE)
        target_link_libraries(phys_sandbox PRIVATE
            "-framework IOKit"
            "-framework Cocoa"
            "-framework OpenGL"
        )
    endif()
endif()

message(STATUS "phys ${PROJECT_VERSION}: tests=${PHYS_BUILD_TESTS} demo=${PHYS_BUILD_DEMO}")
