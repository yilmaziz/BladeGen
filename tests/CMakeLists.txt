# tests/CMakeLists.txt — one Catch2 executable per layer.
include(Catch)
# BLADEGEN_TEST_SPEC: a Catch2 test spec, e.g. "~[slow]" (the ASan preset skips the slow fits)
set(BLADEGEN_TEST_SPEC "" CACHE STRING "Catch2 test spec for test discovery")

function(bladegen_add_test name)
    add_executable(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE bladegen::lib bladegen::compiler_flags Catch2::Catch2WithMain)
    target_compile_definitions(${name} PRIVATE BLADEGEN_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
                                               BLADEGEN_TEST_OUTPUT="${CMAKE_BINARY_DIR}/test_output")
    set(_env "")
    if(WIN32 AND BLADEGEN_OCCT_BIN_DIR)
        set(_env "PATH=${BLADEGEN_OCCT_BIN_DIR}\;$ENV{PATH}")
        string(REPLACE ";" "\;" _env "${_env}")
    endif()
    catch_discover_tests(${name} DISCOVERY_MODE PRE_TEST TEST_SPEC "${BLADEGEN_TEST_SPEC}" PROPERTIES ENVIRONMENT "${_env}")
endfunction()

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/test_output)
bladegen_add_test(test_math    math/test_bspline.cpp math/test_numerics.cpp)
bladegen_add_test(test_profile profile/test_profile.cpp profile/test_quality.cpp profile/test_blade2d_reference.cpp profile/test_camber.cpp)
bladegen_add_test(test_blade   blade/test_blade.cpp blade/test_flowpath.cpp)
bladegen_add_test(test_occt    occt/test_occt.cpp)
bladegen_add_test(test_fit     fit/test_fit.cpp)
bladegen_add_test(test_io      io/test_io.cpp io/test_section_export.cpp)
bladegen_add_test(test_naca    naca/test_naca65.cpp)
bladegen_add_test(test_design  design/test_design.cpp design/test_restagger.cpp design/test_scene.cpp design/test_cg.cpp design/test_stacking.cpp)

# TURBOlenza (src/turbolenza): one executable, every tests/turbolenza/test_*.cpp (globbed).
# The v1 reference comparisons SKIP when tests/turbolenza/data has not been generated.
file(GLOB TURBOLENZA_TESTS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/turbolenza/test_*.cpp)
if(TURBOLENZA_TESTS)
    add_executable(test_turbolenza ${TURBOLENZA_TESTS})
    target_link_libraries(test_turbolenza PRIVATE turbolenza::core bladegen::compiler_flags Catch2::Catch2WithMain)
    if(TARGET turbolenza::bladegen)
        target_link_libraries(test_turbolenza PRIVATE turbolenza::bladegen)
    endif()
    target_include_directories(test_turbolenza PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/turbolenza)
    target_compile_definitions(test_turbolenza PRIVATE BLADEGEN_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
                                                       BLADEGEN_TEST_OUTPUT="${CMAKE_BINARY_DIR}/test_output")
    set(_env "")
    if(WIN32 AND BLADEGEN_OCCT_BIN_DIR)
        set(_env "PATH=${BLADEGEN_OCCT_BIN_DIR}\;$ENV{PATH}")
        string(REPLACE ";" "\;" _env "${_env}")
    endif()
    # a v1-comparison test without its reference data SKIPs (Catch2 exit code 4): ctest reports it as skipped
    catch_discover_tests(test_turbolenza DISCOVERY_MODE PRE_TEST TEST_SPEC "${BLADEGEN_TEST_SPEC}"
                         PROPERTIES ENVIRONMENT "${_env}" SKIP_RETURN_CODE 4)
endif()
