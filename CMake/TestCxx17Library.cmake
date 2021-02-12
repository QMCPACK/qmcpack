
# Test that the compiler is configured with a C++17 standard library

set(TEST_CXX17_SOURCE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/try_cxx17_library.cpp)
file(WRITE ${TEST_CXX17_SOURCE}
"// Test for C++17 standard library support
#include <variant>
#include <string>

int main(int argc, char **argv)
{
    std::variant<int, float, std::string> intFloatString;
    return 0;
}
")


try_compile(CXX17_LIBRARY_OKAY ${CMAKE_BINARY_DIR}
            ${TEST_CXX17_SOURCE}
            CXX_STANDARD 17
            CXX_STANDARD_REQUIRED ON
            OUTPUT_VARIABLE COMPILE_OUTPUT)


if (NOT CXX17_LIBRARY_OKAY)
  set(COMPILE_FAIL_OUTPUT cpp17_compile_fail.txt)
  file(WRITE "${CMAKE_BINARY_DIR}/${COMPILE_FAIL_OUTPUT}" "${COMPILE_OUTPUT}")

  message(STATUS "C++17 standard library support not found")
  message("compiler is ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU") 
    message("Compiler detected is g++.\n  Use version 7.0 or newer for C++17 standard library support.")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message("Compiler detected is clang++.\n  If not using libcxx, ensure a g++ version greater than 7.0 is also on the path so that its C++17 library can be used.")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    message("Compiler detected is icpc.\n  Ensure a gcc version greater than 7.0 is also on the path so that its C++17 library can be used.  Or use the -cxxlib switch to point to a newer gcc install.")
  endif()
  message("  Output of test compile is in ${COMPILE_FAIL_OUTPUT}")
  message(FATAL_ERROR "stopping")
else()
  message(STATUS "C++17 standard library supported")
endif()
