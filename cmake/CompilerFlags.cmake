# cmake/CompilerFlags.cmake
# Strict warnings for every BladeGen target (not for third-party headers, which
# are included as SYSTEM).  Targets link bladegen::compiler_flags.

add_library(bladegen_compiler_flags INTERFACE)
add_library(bladegen::compiler_flags ALIAS bladegen_compiler_flags)

set(_msvc  $<CXX_COMPILER_ID:MSVC>)
set(_gnu   $<CXX_COMPILER_ID:GNU>)
set(_clang $<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>)
set(_gcclike $<OR:${_gnu},${_clang}>)

target_compile_options(bladegen_compiler_flags INTERFACE
    $<${_msvc}:
        /W4 /permissive- /Zc:__cplusplus /utf-8 /EHsc /bigobj
        /external:anglebrackets /external:W0
        /wd4251 /wd4275          # dll-interface warnings from OCCT export macros
    >
    $<${_gcclike}:
        -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
        -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused
        -Woverloaded-virtual -Wdouble-promotion -Wformat=2
    >
    $<${_gnu}:
        -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op
    >
)

if(BLADEGEN_WARNINGS_AS_ERRORS)
    target_compile_options(bladegen_compiler_flags INTERFACE
        $<${_msvc}:/WX>
        $<${_gcclike}:-Werror>
    )
endif()

target_compile_definitions(bladegen_compiler_flags INTERFACE
    $<${_msvc}:NOMINMAX WIN32_LEAN_AND_MEAN _USE_MATH_DEFINES>
)

if(BLADEGEN_ENABLE_ASAN)
    if(MSVC)
        message(FATAL_ERROR "BLADEGEN_ENABLE_ASAN is for GCC/Clang (use a linux-asan preset)")
    endif()
    target_compile_options(bladegen_compiler_flags INTERFACE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(bladegen_compiler_flags INTERFACE -fsanitize=address,undefined)
endif()
