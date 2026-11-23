# apps/CMakeLists.txt — command-line tools
#   bladegen      profile check, NACA 65 sections, blade build (STEP/IGES/STL), inverse fit,
#                 section export (export-sections, JSON), demo
#   rotor67       NASA Rotor 67 reconstruction (fit of the 14 published sections, solid, fillet)
#   profile_dump  development: side samples of a built-in profile as CSV
#   fillet_lab    development: fillet resolution study
foreach(_app bladegen_cli rotor67 profile_dump fillet_lab)
    add_executable(${_app} ${_app}.cpp)
    target_link_libraries(${_app} PRIVATE bladegen::lib bladegen::compiler_flags)
endforeach()
set_target_properties(bladegen_cli PROPERTIES OUTPUT_NAME bladegen)
add_subdirectory(designer)
target_compile_definitions(bladegen_cli PRIVATE BLADEGEN_VERSION="${PROJECT_VERSION}")

if(BLADEGEN_BUILD_TESTS)
    set(_env "")
    if(WIN32 AND BLADEGEN_OCCT_BIN_DIR)
        set(_env "PATH=${BLADEGEN_OCCT_BIN_DIR}\;$ENV{PATH}")
        string(REPLACE ";" "\;" _env "${_env}")
    endif()
    set(_out ${CMAKE_BINARY_DIR}/test_output/cli)
    add_test(NAME cli_version COMMAND bladegen_cli version)
    set_tests_properties(cli_version PROPERTIES PASS_REGULAR_EXPRESSION "BladeGen ${PROJECT_VERSION}")
    add_test(NAME cli_profile COMMAND bladegen_cli profile ${PROJECT_SOURCE_DIR}/tests/data/naca65_1210_blade_params.json --csv ${_out}_curvature.csv)
    # the blade2d NACA fit has an S-shaped TE connector: built, exit code 3 (geometry errors)
    set_tests_properties(cli_profile PROPERTIES WILL_FAIL FALSE PASS_REGULAR_EXPRESSION "S-shaped")
    add_test(NAME cli_demo COMMAND bladegen_cli demo ${_out})
    add_test(NAME cli_build COMMAND bladegen_cli build ${_out}/demo_rotor.json ${_out}/rebuilt.stp --unit m)
    set_tests_properties(cli_build PROPERTIES DEPENDS cli_demo)
    add_test(NAME cli_naca65 COMMAND bladegen_cli naca65 1.2 0.10 --variant emery_base --csv ${_out}_naca65.csv)
    set_tests_properties(cli_naca65 PROPERTIES PASS_REGULAR_EXPRESSION "NACA 65-.12A10.10 .emery_base.*[^1-9]0 errors")
    add_test(NAME cli_naca65_bad COMMAND bladegen_cli naca65 1.2 0.10 --variant dca)
    set_tests_properties(cli_naca65_bad PROPERTIES WILL_FAIL TRUE)
    # export-sections: a design file at the hub, a NACA 65 cascade, and refused input
    add_test(NAME cli_export_sections_design COMMAND bladegen_cli export-sections
             ${PROJECT_SOURCE_DIR}/data/rotor67/rotor67_bladegen.json --span 0 --out ${_out}_rotor67_hub_sections.json)
    set_tests_properties(cli_export_sections_design PROPERTIES
                         PASS_REGULAR_EXPRESSION "NASA_Rotor_67: 1 section.s., 238 points per side, 22 blades")
    add_test(NAME cli_export_sections_naca65 COMMAND bladegen_cli export-sections naca65 1.2 0.10
             --variant kovach_sandercock --chord 0.5 --stagger 33 --solidity 1.0 --stations 101
             --out ${_out}_naca65_cascade.json)
    set_tests_properties(cli_export_sections_naca65 PROPERTIES
                         PASS_REGULAR_EXPRESSION "NACA 65-.12A10.10 kovach_sandercock, 179 SS . 179 PS points, pitch 0.5 m")
    add_test(NAME cli_export_sections_stdout COMMAND bladegen_cli export-sections naca65 0.8 0.08
             --variant emery_base --chord 1 --stagger 25 --solidity 1.25)
    set_tests_properties(cli_export_sections_stdout PROPERTIES
                         PASS_REGULAR_EXPRESSION "\"format\": \"bladegen.sections\",[^{}]*\"version\": 1,")
    add_test(NAME cli_export_sections_bad_span COMMAND bladegen_cli export-sections
             ${PROJECT_SOURCE_DIR}/data/rotor67/rotor67_bladegen.json --span 1.5)
    add_test(NAME cli_export_sections_bad_number COMMAND bladegen_cli export-sections naca65 1.2 0.10
             --variant kovach_sandercock --chord 0.5x --stagger 33 --solidity 1)
    add_test(NAME cli_export_sections_missing COMMAND bladegen_cli export-sections naca65 1.2 0.10
             --variant kovach_sandercock --chord 0.5 --stagger 33)
    set_tests_properties(cli_export_sections_bad_span cli_export_sections_bad_number cli_export_sections_missing
                         PROPERTIES WILL_FAIL TRUE)
    set_tests_properties(cli_export_sections_design cli_export_sections_naca65 cli_export_sections_stdout
                         cli_export_sections_bad_span cli_export_sections_bad_number cli_export_sections_missing
                         PROPERTIES ENVIRONMENT "${_env}")
    add_test(NAME cli_bad_usage COMMAND bladegen_cli frobnicate)
    set_tests_properties(cli_bad_usage PROPERTIES WILL_FAIL TRUE)
    set_tests_properties(cli_version cli_profile cli_naca65 cli_naca65_bad cli_demo cli_build cli_bad_usage PROPERTIES ENVIRONMENT "${_env}")
endif()
