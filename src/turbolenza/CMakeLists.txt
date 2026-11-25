# src/turbolenza — TURBOlenza, the blade-to-blade compressible potential-flow solver after
# Aungier (2003) Ch. 5: the C++ port of the verified TURBOlenza v1 Python package
# (turbolenza.b2b).  A library of its own, without OpenCASCADE or Qt:
#   turbolenza_core      cascade, grid, metrics .. solver, postprocess (Eigen; OpenMP, MKL PARDISO optional)
#   turbolenza_bladegen  the BladeGen adapters: a design layer's profile -> the cascade's surface curves,
#                        and a design -> the through-flow's quasi-normals and blade row (ThroughflowBladeGen.cpp)
# Every .cpp in this folder except the adapters belongs to the core (globbed, so modules can be
# added without editing this file).

file(GLOB TURBOLENZA_CORE_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
list(FILTER TURBOLENZA_CORE_SOURCES EXCLUDE REGEX "(BladeGenSource|ThroughflowBladeGen)\.cpp$")
if(NOT TURBOLENZA_CORE_SOURCES)
    set(TURBOLENZA_CORE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Version.cpp)
endif()

add_library(turbolenza_core STATIC ${TURBOLENZA_CORE_SOURCES})
add_library(turbolenza::core ALIAS turbolenza_core)
target_include_directories(turbolenza_core PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)
target_link_libraries(turbolenza_core PUBLIC Eigen3::Eigen bladegen::json PRIVATE bladegen::compiler_flags)

# OpenMP: the node loops of a pass (coefficients, flow, gas state) run on the chosen cores.
find_package(OpenMP QUIET COMPONENTS CXX)
if(OpenMP_CXX_FOUND)
    target_link_libraries(turbolenza_core PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(turbolenza_core PUBLIC TURBOLENZA_HAVE_OPENMP)
    message(STATUS "TURBOlenza: OpenMP ${OpenMP_CXX_VERSION}")
else()
    message(STATUS "TURBOlenza: OpenMP not found, the solver runs on one core")
endif()

# Intel MKL PARDISO: the parallel sparse direct solver, selectable beside Eigen's SparseLU.
set(TURBOLENZA_WITH_PARDISO AUTO CACHE STRING "Use Intel MKL PARDISO when found (AUTO, ON, OFF)")
set_property(CACHE TURBOLENZA_WITH_PARDISO PROPERTY STRINGS AUTO ON OFF)
if(NOT TURBOLENZA_WITH_PARDISO STREQUAL "OFF")
    set(MKL_INTERFACE lp64)
    # MKL's TBB threading layer: the OpenMP loops of a pass keep the compiler's OpenMP runtime, and
    # PARDISO runs on TBB (it runs between those loops, never beside them).  MKL's own OpenMP layer
    # (libiomp5md) cannot share MSVC's /openmp runtime.
    set(MKL_THREADING tbb_thread)
    set(MKL_LINK dynamic)
    if(TURBOLENZA_WITH_PARDISO STREQUAL "ON")
        find_package(MKL CONFIG REQUIRED)
    else()
        find_package(MKL CONFIG QUIET)
    endif()
    if(MKL_FOUND)
        target_link_libraries(turbolenza_core PUBLIC MKL::MKL)
        target_compile_definitions(turbolenza_core PUBLIC TURBOLENZA_HAVE_PARDISO)
        message(STATUS "TURBOlenza: Intel MKL PARDISO enabled (${MKL_THREADING})")
    else()
        message(STATUS "TURBOlenza: Intel MKL not found, PARDISO disabled (Eigen SparseLU only)")
    endif()
endif()

# The BladeGen adapter (links the BladeGen core library).
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/BladeGenSource.cpp)
    set(TURBOLENZA_BLADEGEN_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/BladeGenSource.cpp)
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ThroughflowBladeGen.cpp)
        list(APPEND TURBOLENZA_BLADEGEN_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/ThroughflowBladeGen.cpp)
    endif()
    add_library(turbolenza_bladegen STATIC ${TURBOLENZA_BLADEGEN_SOURCES})
    add_library(turbolenza::bladegen ALIAS turbolenza_bladegen)
    target_link_libraries(turbolenza_bladegen PUBLIC turbolenza::core bladegen::lib PRIVATE bladegen::compiler_flags)
endif()
