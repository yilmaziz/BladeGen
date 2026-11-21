# python/CMakeLists.txt — the `bladegen` Python package (pybind11 extension + helpers).
#
# Layout written to the build tree (put <build>/python on PYTHONPATH):
#   <build>/python/bladegen/__init__.py   re-exports the extension, finds OCCT DLLs on Windows
#   <build>/python/bladegen/_paths.py     build-time paths (OCCT runtime directory)
#   <build>/python/bladegen/_core.*.pyd   the extension module (C++ engine)
#   <build>/python/bladegen/notebook.py   plotting / widget helpers for the notebooks

find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module)
execute_process(COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
                OUTPUT_VARIABLE _pybind11_dir OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
find_package(pybind11 CONFIG REQUIRED HINTS "${_pybind11_dir}")

set(_pkg ${CMAKE_BINARY_DIR}/python/bladegen)
pybind11_add_module(bladegen_py MODULE bindings.cpp bind_profile.cpp bind_naca.cpp bind_blade.cpp bind_fit.cpp bind_occt.cpp bind_design.cpp)
target_compile_definitions(bladegen_py PRIVATE BLADEGEN_VERSION="${PROJECT_VERSION}")
target_link_libraries(bladegen_py PRIVATE bladegen::lib bladegen::compiler_flags)
set_target_properties(bladegen_py PROPERTIES
    OUTPUT_NAME _core
    LIBRARY_OUTPUT_DIRECTORY ${_pkg}
    LIBRARY_OUTPUT_DIRECTORY_RELEASE ${_pkg}
    LIBRARY_OUTPUT_DIRECTORY_DEBUG ${_pkg}
    LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO ${_pkg})

file(MAKE_DIRECTORY ${_pkg})
configure_file(bladegen/_paths.py.in ${_pkg}/_paths.py @ONLY)
foreach(_f __init__.py notebook.py)
    configure_file(bladegen/${_f} ${_pkg}/${_f} COPYONLY)
endforeach()

if(BLADEGEN_BUILD_TESTS)
    add_test(NAME python_bindings
             COMMAND ${Python_EXECUTABLE} -m pytest -q ${CMAKE_CURRENT_SOURCE_DIR}/tests -p no:cacheprovider)
    set_tests_properties(python_bindings PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python;MPLBACKEND=Agg;BLADEGEN_SOURCE_DIR=${PROJECT_SOURCE_DIR};BLADEGEN_TEST_OUTPUT=${CMAKE_BINARY_DIR}/test_output")
    execute_process(COMMAND ${Python_EXECUTABLE} -c "import nbclient, ipykernel, matplotlib"
                    RESULT_VARIABLE _nb_rc OUTPUT_QUIET ERROR_QUIET)
    if(_nb_rc EQUAL 0)
        add_test(NAME python_notebooks
                 COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/run_notebooks.py)
        set_tests_properties(python_notebooks PROPERTIES
            ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python;MPLBACKEND=Agg;BLADEGEN_NB_STATIC=1;BLADEGEN_SOURCE_DIR=${PROJECT_SOURCE_DIR};BLADEGEN_TEST_OUTPUT=${CMAKE_BINARY_DIR}/test_output"
            TIMEOUT 1800)
    endif()
endif()
