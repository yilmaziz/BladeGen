# src/CMakeLists.txt — the BladeGen core library.
# Layers (each depends only on the ones above it):
#   math     B-splines, interpolation, quadrature, roots, least squares (Eigen)
#   profile  2-D C2 blade profile, its parameters and quality checks
#   blade    meridional channel, conformal m' mapping, 3-D sections and blades
#   fit      inverse design: profile parameters from published coordinates
#   occt     OpenCASCADE: NURBS surfaces, solids, hub fillet, STEP/IGES/STL
#   io       JSON parameter files, published data readers
#   design   the designer's blade: meridional channel, five layers, mesh and fields
#   pipeline one-call blade generation (CLI and Python)

set(BLADEGEN_SOURCES
    profile/ProfileParams.cpp
    profile/Profile2D.cpp
    profile/Quality.cpp
    profile/Camber.cpp
    blade/Meridional.cpp
    blade/Blade.cpp
    blade/FlowPath.cpp
    occt/OcctBlade.cpp
    occt/Fillet.cpp
    fit/InverseFit.cpp
    fit/SectionFit.cpp
    fit/NacaFit.cpp
    io/Rotor67.cpp
    io/BladeJson.cpp
    naca/Naca65.cpp
    design/BladeDesign.cpp
    design/DesignSolid.cpp
)

add_library(bladegen_lib STATIC ${BLADEGEN_SOURCES})
add_library(bladegen::lib ALIAS bladegen_lib)
target_include_directories(bladegen_lib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(bladegen_lib
    PUBLIC  Eigen3::Eigen bladegen::occt bladegen::json
    PRIVATE bladegen::compiler_flags)
set_target_properties(bladegen_lib PROPERTIES OUTPUT_NAME bladegen_core)
