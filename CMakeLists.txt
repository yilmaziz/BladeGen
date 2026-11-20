# ------------------------------------------------------------------------------
# BladeGen — parametric turbomachinery blade generator (C++23, OpenCASCADE)
# ------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.25)

project(BladeGen
    VERSION 1.0.0
    DESCRIPTION "Parametric C2 compressor/turbine blade generator on OpenCASCADE"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)   # the static library also goes into the Python module

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option(BLADEGEN_BUILD_TESTS   "Build the Catch2 test suite"                       ON)
option(BLADEGEN_BUILD_PYTHON  "Build the pybind11 Python module (bladegen)"       OFF)
option(BLADEGEN_WARNINGS_AS_ERRORS "Treat compiler warnings as errors"            OFF)
option(BLADEGEN_ENABLE_ASAN   "Enable AddressSanitizer + UBSan (GCC/Clang)"       OFF)

include(CompilerFlags)
include(Dependencies)

add_subdirectory(src)
add_subdirectory(apps)

if(BLADEGEN_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(BLADEGEN_BUILD_PYTHON)
    add_subdirectory(python)
endif()
