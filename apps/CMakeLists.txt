# apps/CMakeLists.txt — command-line tools
add_executable(profile_dump profile_dump.cpp)
target_link_libraries(profile_dump PRIVATE bladegen::lib bladegen::compiler_flags)
