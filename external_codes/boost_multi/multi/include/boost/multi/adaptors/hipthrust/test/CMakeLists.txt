cmake_minimum_required(VERSION 3.24)

project(boost-multi-adaptor-thrust-hip LANGUAGES CXX HIP)

find_package(Boost     REQUIRED COMPONENTS headers)
find_package(rocthrust)
if(NOT rocthrust_FOUND)
	message(FATAL_ERROR "Could not find rocThrust. On Ubuntu/Debian with the ROCm apt repo configured, install the development package with:\n  sudo apt-get install rocthrust-dev\nor point CMake at an existing install via -Drocthrust_DIR=<prefix>/lib/cmake/rocthrust .")
endif()

include(CMakePrintHelpers)

cmake_print_properties(TARGETS roc::rocthrust PROPERTIES INTERFACE_LINK_LIBRARIES INTERFACE_INCLUDE_DIRECTORIES)

enable_testing()

include(CTest)

add_executable(vector.hip.amd.x vector.hip)
#set_source_files_properties(vector.hip PROPERTIES LANGUAGE HIP)
target_link_libraries(vector.hip.amd.x PRIVATE multi roc::rocthrust)

add_test(NAME vector.hip.amd.x COMMAND $<TARGET_FILE:vector.hip.amd.x>)

add_executable(array.hip.amd.x array.hip)
#set_source_files_properties(array.hip PROPERTIES LANGUAGE HIP)
target_link_libraries(array.hip.amd.x PRIVATE multi roc::rocthrust)

add_test(NAME array.hip.amd.x COMMAND $<TARGET_FILE:array.hip.amd.x>)

add_executable(speed.hip.amd.x speed.hip)
#set_source_files_properties(speed.hip PROPERTIES LANGUAGE HIP)
target_link_libraries(speed.hip.amd.x PRIVATE multi roc::rocthrust)

add_test(NAME speed.hip.amd.x COMMAND $<TARGET_FILE:speed.hip.amd.x>)

set_tests_properties(speed.hip.amd.x PROPERTIES RUN_SERIAL TRUE)

