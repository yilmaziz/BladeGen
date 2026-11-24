# apps/designer — bladegen_designer, the Qt Widgets MDI designer (C++ views, no Python).
# Built when Qt 6 (Widgets, Concurrent) is found: conda-forge qt6-main, apt qt6-base-dev,
# Homebrew qt.  BLADEGEN_BUILD_DESIGNER=OFF skips it; =ON makes Qt required.
set(BLADEGEN_BUILD_DESIGNER AUTO CACHE STRING "Build the Qt designer (AUTO, ON, OFF)")
set_property(CACHE BLADEGEN_BUILD_DESIGNER PROPERTY STRINGS AUTO ON OFF)
if(BLADEGEN_BUILD_DESIGNER STREQUAL "OFF")
    return()
endif()
if(BLADEGEN_BUILD_DESIGNER STREQUAL "ON")
    find_package(Qt6 6.4 REQUIRED COMPONENTS Widgets Concurrent)
else()
    find_package(Qt6 6.4 QUIET COMPONENTS Widgets Concurrent)
    if(NOT Qt6_FOUND)
        message(STATUS "Qt 6 not found: the designer (bladegen_designer) is not built")
        return()
    endif()
endif()
message(STATUS "Qt ${Qt6_VERSION}: building bladegen_designer")

add_executable(bladegen_designer
    main.cpp DesignerWindow.cpp DesignerWindow.hpp ContentView.cpp ContentView.hpp
    ParameterPanel.cpp ParameterPanel.hpp PlotWidget.cpp PlotWidget.hpp View3D.cpp View3D.hpp
    OcctView.cpp OcctView.hpp Engine.cpp Engine.hpp)
# The TURBOlenza window (the blade-to-blade solver of src/turbolenza on a copy of one layer).
target_sources(bladegen_designer PRIVATE
    turbolenza/TurbolenzaWindow.cpp turbolenza/TurbolenzaWindow.hpp turbolenza/TurbolenzaPages.cpp
    turbolenza/TurbolenzaPlots.cpp turbolenza/TurbolenzaPlots.hpp)
target_include_directories(bladegen_designer PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(bladegen_designer PROPERTIES AUTOMOC ON WIN32_EXECUTABLE OFF)
target_link_libraries(bladegen_designer PRIVATE bladegen::lib bladegen::occt_viewer bladegen::compiler_flags Qt6::Widgets Qt6::Concurrent
                                                turbolenza::core turbolenza::bladegen)
target_compile_definitions(bladegen_designer PRIVATE BLADEGEN_VERSION="${PROJECT_VERSION}")

if(BLADEGEN_BUILD_TESTS)
    set(_env "QT_QPA_PLATFORM=offscreen")
    if(WIN32 AND BLADEGEN_OCCT_BIN_DIR)
        list(APPEND _env "PATH=${BLADEGEN_OCCT_BIN_DIR}\;$ENV{PATH}")
        string(REPLACE ";" "\;" _p "${BLADEGEN_OCCT_BIN_DIR};$ENV{PATH}")
        set(_env "QT_QPA_PLATFORM=offscreen;PATH=${_p}")
    endif()
    add_test(NAME designer_self_test
             COMMAND bladegen_designer --self-test ${CMAKE_BINARY_DIR}/test_output/designer)
    set_tests_properties(designer_self_test PROPERTIES ENVIRONMENT "${_env}" PASS_REGULAR_EXPRESSION "view renders OK" TIMEOUT 300)
    add_test(NAME designer_turbolenza_self_test
             COMMAND bladegen_designer --turbolenza-self-test ${CMAKE_BINARY_DIR}/test_output/turbolenza_ui)
    set_tests_properties(designer_turbolenza_self_test PROPERTIES ENVIRONMENT "${_env}" PASS_REGULAR_EXPRESSION "turbolenza self-test OK" TIMEOUT 900)
    add_test(NAME designer_bad_file COMMAND bladegen_designer ${PROJECT_SOURCE_DIR}/README.md --screenshot ${CMAKE_BINARY_DIR}/test_output/bad.png)
    set_tests_properties(designer_bad_file PROPERTIES ENVIRONMENT "${_env}" WILL_FAIL TRUE TIMEOUT 120)
endif()
