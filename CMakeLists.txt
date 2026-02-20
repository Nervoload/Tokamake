cmake_minimum_required(VERSION 3.20)

project(tokamakfusion LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(tokamak_core STATIC
    src/artifact_export.cpp
    src/collision.cpp
    src/electrostatic_field.cpp
    src/electrostatic_models.cpp
    src/engine.cpp
    src/magnetic_field.cpp
    src/particle_system.cpp
    src/particle_push.cpp
    src/spatial_grid.cpp
    src/telemetry.cpp)

target_include_directories(tokamak_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(tokamakfusion src/main.cpp)
target_link_libraries(tokamakfusion PRIVATE tokamak_core)

include(CTest)
if(BUILD_TESTING)
    add_subdirectory(third_party/googletest)

    add_executable(tokamak_tests
        tests/engine_tests.cpp
        tests/integration_tests.cpp
        tests/milestone2_magnetic_field_tests.cpp
        tests/milestone3_electrostatic_field_tests.cpp
        tests/milestone6_validation_tests.cpp
        tests/particle_system_tests.cpp)
    target_link_libraries(tokamak_tests PRIVATE tokamak_core gtest_main)
    target_compile_definitions(tokamak_tests PRIVATE TOKAMAKFUSION_PATH="$<TARGET_FILE:tokamakfusion>")

    add_test(NAME tokamak_tests COMMAND tokamak_tests)
endif()
