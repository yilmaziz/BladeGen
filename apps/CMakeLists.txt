# apps/CMakeLists.txt — command-line tools
add_executable(profile_dump profile_dump.cpp)
target_link_libraries(profile_dump PRIVATE bladegen::lib bladegen::compiler_flags)
add_executable(fillet_lab fillet_lab.cpp)
target_link_libraries(fillet_lab PRIVATE bladegen::lib bladegen::compiler_flags)
add_executable(rotor67 rotor67.cpp)
target_link_libraries(rotor67 PRIVATE bladegen::lib bladegen::compiler_flags)
