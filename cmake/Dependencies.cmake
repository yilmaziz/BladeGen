# cmake/Dependencies.cmake
# Third-party libraries.  Each one is looked up as an installed package first
# (conda-forge on Windows, apt on Linux, or anything on CMAKE_PREFIX_PATH);
# the header-only ones fall back to FetchContent so a bare machine still builds.
#
# Exposed targets:
#   bladegen::occt   OpenCASCADE geometry kernel + data exchange
#   bladegen::occt_viewer  OCCT visualization (V3d, AIS, OpenGL) when installed (else empty)
#   Eigen3::Eigen    dense / sparse linear algebra (header only)
#   bladegen::json   nlohmann/json (header only)
#   Catch2::Catch2WithMain   (tests only)

include(FetchContent)

# ── OpenCASCADE ───────────────────────────────────────────────────────────────
find_package(OpenCASCADE CONFIG REQUIRED)
message(STATUS "OpenCASCADE ${OpenCASCADE_VERSION} from ${OpenCASCADE_DIR}")

add_library(bladegen_occt INTERFACE)
add_library(bladegen::occt ALIAS bladegen_occt)
target_include_directories(bladegen_occt SYSTEM INTERFACE ${OpenCASCADE_INCLUDE_DIR})

set(_occt_libs
    TKernel TKMath TKG2d TKG3d TKGeomBase TKGeomAlgo TKBRep TKTopAlgo TKPrim
    TKBO TKBool TKFillet TKOffset TKShHealing TKMesh TKXSBase)
# OCCT 7.8 renamed the data-exchange toolkits (TKSTEP -> TKDESTEP, ...).
if(TARGET TKDESTEP)
    list(APPEND _occt_libs TKDE TKDESTEP TKDEIGES TKDESTL)
else()
    list(APPEND _occt_libs TKSTEP TKSTEPBase TKSTEPAttr TKSTEP209 TKIGES TKSTL)
endif()
foreach(_lib IN LISTS _occt_libs)
    if(NOT TARGET ${_lib})
        message(FATAL_ERROR "OpenCASCADE toolkit ${_lib} not found in ${OpenCASCADE_DIR}")
    endif()
endforeach()
target_link_libraries(bladegen_occt INTERFACE ${_occt_libs})

# Visualization toolkits (V3d / AIS / OpenGL), optional: the designer's 3-D solid view.
add_library(bladegen_occt_viewer INTERFACE)
add_library(bladegen::occt_viewer ALIAS bladegen_occt_viewer)
set(BLADEGEN_HAVE_OCCT_VIEWER OFF)
if(TARGET TKService AND TARGET TKV3d AND TARGET TKOpenGl)
    set(BLADEGEN_HAVE_OCCT_VIEWER ON)
    target_link_libraries(bladegen_occt_viewer INTERFACE bladegen::occt TKService TKV3d TKOpenGl)
    target_compile_definitions(bladegen_occt_viewer INTERFACE BLADEGEN_HAVE_OCCT_VIEWER)
    message(STATUS "OpenCASCADE visualization toolkits found: the designer gets the 3-D solid view")
else()
    message(STATUS "OpenCASCADE visualization toolkits not found: the designer's 3-D solid view is disabled")
endif()

# Directory holding the OCCT shared libraries (Windows: needed on PATH for tests).
get_target_property(_tkernel_loc TKernel IMPORTED_LOCATION_RELEASE)
if(NOT _tkernel_loc)
    get_target_property(_tkernel_loc TKernel IMPORTED_LOCATION)
endif()
if(_tkernel_loc)
    get_filename_component(BLADEGEN_OCCT_BIN_DIR "${_tkernel_loc}" DIRECTORY)
endif()
set(BLADEGEN_OCCT_BIN_DIR "${BLADEGEN_OCCT_BIN_DIR}" CACHE INTERNAL "OCCT runtime directory")

# ── Eigen ─────────────────────────────────────────────────────────────────────
find_package(Eigen3 3.3 CONFIG QUIET)
if(NOT TARGET Eigen3::Eigen)
    message(STATUS "Eigen3 not installed: fetching 3.4.0")
    FetchContent_Declare(eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG        3.4.0
        GIT_SHALLOW    TRUE)
    FetchContent_GetProperties(eigen)
    if(NOT eigen_POPULATED)
        FetchContent_Populate(eigen)
    endif()
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${eigen_SOURCE_DIR}")
endif()

# ── nlohmann/json ─────────────────────────────────────────────────────────────
find_package(nlohmann_json 3.10 CONFIG QUIET)
if(NOT TARGET nlohmann_json::nlohmann_json)
    message(STATUS "nlohmann_json not installed: fetching v3.11.3")
    FetchContent_Declare(nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG        v3.11.3
        GIT_SHALLOW    TRUE)
    set(JSON_BuildTests OFF CACHE BOOL "" FORCE)
    set(JSON_Install    OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(nlohmann_json)
endif()
add_library(bladegen_json INTERFACE)
add_library(bladegen::json ALIAS bladegen_json)
target_link_libraries(bladegen_json INTERFACE nlohmann_json::nlohmann_json)

# ── Catch2 ────────────────────────────────────────────────────────────────────
if(BLADEGEN_BUILD_TESTS)
    find_package(Catch2 3 CONFIG QUIET)
    if(NOT TARGET Catch2::Catch2WithMain)
        message(STATUS "Catch2 3 not installed: fetching v3.6.0")
        FetchContent_Declare(Catch2
            GIT_REPOSITORY https://github.com/catchorg/Catch2.git
            GIT_TAG        v3.6.0
            GIT_SHALLOW    TRUE)
        FetchContent_MakeAvailable(Catch2)
        list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")
    endif()
endif()
