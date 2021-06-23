# Test that the compiler is configured with a C++ 14 standard library

set(TEST_CXX14_SOURCE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/try_cxx14_library.cpp)
file(
  WRITE ${TEST_CXX14_SOURCE}
  "// Test for C++14 standard library support
#include <tuple>

int main(int argc, char **argv)
{
    std::tuple<int,double> t(1, 3.5);
    // Accessing a tuple member by type is a C++14 feature
    int j = std::get<int>(t);
    return 0;
}
")

try_compile(
  CXX14_LIBRARY_OKAY
  ${CMAKE_BINARY_DIR}
  ${TEST_CXX14_SOURCE}
  CXX_STANDARD
  14
  CXX_STANDARD_REQUIRED
  ON
  OUTPUT_VARIABLE COMPILE_OUTPUT)

if(NOT CXX14_LIBRARY_OKAY)
  set(COMPILE_FAIL_OUTPUT cpp14_compile_fail.txt)
  file(WRITE "${CMAKE_BINARY_DIR}/${COMPILE_FAIL_OUTPUT}" "${COMPILE_OUTPUT}")

  message(STATUS "C++14 standard library support not found")
  message("compiler is ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message("Compiler detected is g++.\n  Use version 5.0 or newer for a C++14 standard library")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(
      "Compiler detected is clang++.\n  If not using libcxx, ensure a g++ version greater than 5.0 is on the path")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    message(
      "Compiler detected is icpc.\n  Ensure a gcc version greater than 5.0 is on the path.  Or use the -cxxlib switch to point to a newer gcc install."
    )
  endif()
  message("  Output of test compile is in ${COMPILE_FAIL_OUTPUT}")
  message(FATAL_ERROR "stopping")
else()
  message(STATUS "C++14 standard library supported")
endif()
