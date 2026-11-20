# src/CMakeLists.txt — the BladeGen core library.
# Layers (each depends only on the ones above it):
#   math     B-splines, interpolation, quadrature, roots, least squares (Eigen)
#   profile  2-D C2 blade profile, its parameters and quality checks
#   blade    meridional channel, conformal m' mapping, 3-D sections and blades
#   fit      inverse design: profile parameters from published coordinates
#   occt     OpenCASCADE: NURBS surfaces, solids, hub fillet, STEP/IGES/STL
#   io       JSON parameter files, published data readers
#   pipeline one-call blade generation (CLI and Python)

set(BLADEGEN_SOURCES
    profile/ProfileParams.cpp
    profile/Profile2D.cpp
    profile/Quality.cpp
    blade/Meridional.cpp
    blade/Blade.cpp
    occt/OcctBlade.cpp
    occt/Fillet.cpp
)

add_library(bladegen_lib STATIC ${BLADEGEN_SOURCES})
add_library(bladegen::lib ALIAS bladegen_lib)
target_include_directories(bladegen_lib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(bladegen_lib
    PUBLIC  Eigen3::Eigen bladegen::occt bladegen::json
    PRIVATE bladegen::compiler_flags)
set_target_properties(bladegen_lib PROPERTIES OUTPUT_NAME bladegen)
