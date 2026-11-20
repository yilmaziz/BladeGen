# tests/CMakeLists.txt — one Catch2 executable per layer.
include(Catch)

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
    catch_discover_tests(${name} DISCOVERY_MODE PRE_TEST PROPERTIES ENVIRONMENT "${_env}")
endfunction()

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/test_output)
bladegen_add_test(test_math    math/test_bspline.cpp math/test_numerics.cpp)
bladegen_add_test(test_profile profile/test_profile.cpp profile/test_quality.cpp profile/test_blade2d_reference.cpp)
