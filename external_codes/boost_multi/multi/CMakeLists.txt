cmake_minimum_required(VERSION 3.11)
project(boost-multi CXX)
set(CMAKE_CXX_STANDARD 14)

add_library(${PROJECT_NAME} INTERFACE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake") 

#configure_file(config.h.in config.h)

set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_RULE_MESSAGES OFF)

enable_testing()
find_program(MEMORYCHECK_COMMAND valgrind)
set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --error-exitcode=1")
include(CTest)

add_subdirectory(test)
add_subdirectory(adaptors/blas)

