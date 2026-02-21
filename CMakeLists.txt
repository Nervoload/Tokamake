cmake_minimum_required(VERSION 3.20)

project(tokamakfusion LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_VIEWER "Build tokamak_viewer OpenGL replay executable" OFF)
option(VIEWER_STRICT_VENDOR_CHECK "Require vendored ImGui/GLAD files when BUILD_VIEWER is ON" ON)

add_library(tokamak_core STATIC
    src/artifact_export.cpp
    src/collision.cpp
    src/electrostatic_field.cpp
    src/electrostatic_models.cpp
    src/engine.cpp
    src/magnetic_field.cpp
    src/particle_system.cpp
    src/particle_push.cpp
    src/reactivity.cpp
    src/spatial_grid.cpp
    src/telemetry.cpp
    src/viewer/replay_loader.cpp
    src/viewer/replay_manifest.cpp
    src/viewer/replay_snapshot.cpp)

target_include_directories(tokamak_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(tokamakfusion src/main.cpp)
target_link_libraries(tokamakfusion PRIVATE tokamak_core)

if(BUILD_VIEWER)
    find_package(glfw3 REQUIRED)

    set(VIEWER_VENDOR_REQUIRED_FILES
        third_party/glad/src/glad.c
        third_party/glad/include/glad/glad.h
        third_party/glad/include/KHR/khrplatform.h
        third_party/imgui/imgui.cpp
        third_party/imgui/imgui_draw.cpp
        third_party/imgui/imgui_tables.cpp
        third_party/imgui/imgui_widgets.cpp
        third_party/imgui/imgui.h
        third_party/imgui/backends/imgui_impl_glfw.cpp
        third_party/imgui/backends/imgui_impl_glfw.h
        third_party/imgui/backends/imgui_impl_opengl3.cpp
        third_party/imgui/backends/imgui_impl_opengl3.h)

    if(VIEWER_STRICT_VENDOR_CHECK)
        set(VIEWER_VENDOR_MISSING_FILES "")
        foreach(vendor_file IN LISTS VIEWER_VENDOR_REQUIRED_FILES)
            if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${vendor_file}")
                list(APPEND VIEWER_VENDOR_MISSING_FILES "${vendor_file}")
            endif()
        endforeach()

        if(VIEWER_VENDOR_MISSING_FILES)
            string(JOIN "\n  " VIEWER_VENDOR_MISSING_TEXT ${VIEWER_VENDOR_MISSING_FILES})
            message(FATAL_ERROR
                "BUILD_VIEWER=ON requires vendored ImGui/GLAD files in third_party. Missing:\n"
                "  ${VIEWER_VENDOR_MISSING_TEXT}\n"
                "See docs/VIEWER_SETUP.md for setup instructions.")
        endif()
    endif()

    add_executable(tokamak_viewer
        src/viewer/main.cpp
        src/viewer/viewer_app.cpp
        src/viewer/gl_renderer.cpp
        src/viewer/camera.cpp
        third_party/glad/src/glad.c
        third_party/imgui/imgui.cpp
        third_party/imgui/imgui_draw.cpp
        third_party/imgui/imgui_tables.cpp
        third_party/imgui/imgui_widgets.cpp
        third_party/imgui/backends/imgui_impl_glfw.cpp
        third_party/imgui/backends/imgui_impl_opengl3.cpp)

    target_link_libraries(tokamak_viewer PRIVATE tokamak_core glfw)
    target_include_directories(tokamak_viewer PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/glad/include
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/imgui
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/imgui/backends)

    if(APPLE)
        target_link_libraries(tokamak_viewer PRIVATE "-framework OpenGL")
    else()
        find_package(OpenGL REQUIRED)
        target_link_libraries(tokamak_viewer PRIVATE OpenGL::GL)
    endif()
endif()

include(CTest)
if(BUILD_TESTING)
    add_subdirectory(third_party/googletest)

    add_executable(tokamak_tests
        tests/engine_tests.cpp
        tests/integration_tests.cpp
        tests/milestone2_magnetic_field_tests.cpp
        tests/milestone3_electrostatic_field_tests.cpp
        tests/milestone6_validation_tests.cpp
        tests/particle_system_tests.cpp
        tests/reactivity_tests.cpp
        tests/viewer_replay_tests.cpp)
    target_link_libraries(tokamak_tests PRIVATE tokamak_core gtest_main)
    target_compile_definitions(tokamak_tests PRIVATE
        TOKAMAKFUSION_PATH="$<TARGET_FILE:tokamakfusion>"
        TOKAMAK_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

    add_test(NAME tokamak_tests COMMAND tokamak_tests)
endif()
