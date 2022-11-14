# Test that the compiler is configured with a C++17 standard library

set(TEST_CXX17_SOURCE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/try_cxx17_library.cpp)
file(
  WRITE ${TEST_CXX17_SOURCE}
  "// Test for C++17 standard library support
#include <array>
#include <cstddef>
#include <memory_resource>

int main(int argc, char **argv)
{
    // allocate memory on the stack
    std::array<std::byte, 20000> buf;

    // without fallback memory allocation on heap
    std::pmr::monotonic_buffer_resource pool{ buf.data(), buf.size(),
                                              std::pmr::null_memory_resource() };
    return 0;
}
")

try_compile(
  CXX17_LIBRARY_OKAY
  ${CMAKE_BINARY_DIR}
  ${TEST_CXX17_SOURCE}
  CXX_STANDARD
  17
  CXX_STANDARD_REQUIRED
  ON
  OUTPUT_VARIABLE COMPILE_OUTPUT)

if(NOT CXX17_LIBRARY_OKAY)
  set(COMPILE_FAIL_OUTPUT cpp17_compile_fail.txt)
  file(WRITE "${CMAKE_BINARY_DIR}/${COMPILE_FAIL_OUTPUT}" "${COMPILE_OUTPUT}")

  message(STATUS "C++17 standard library support not found or incomplete")
  message("compiler is ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message("Compiler detected is g++.\n  Use version 9.0 or newer for complete C++17 standard library support.")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
    message(
      "Compiler detected is <Clang> namely clang++ or a vendor variant (icpx, amdclang++, armclang++).\n  If not using libcxx, ensure a GCC toolchain version equal or greater "
      "than 9.0 gets picked up. Check with '<Clang> -v'. Or use the --gcc-toolchain compiler option "
      "(added to both CMAKE_C_FLAGS and CMAKE_CXX_FLAGS) to point to a newer GCC installation."
    )
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    message(
      "Compiler detected is icpc.\n  Ensure a GCC version equal or greater than 9.0 is also on the PATH "
       "such that its C++17 library can be used. Check with 'icpc -v'. Or use the -cxxlib compiler option "
       "(added to CMAKE_CXX_FLAGS) to point to a newer GCC installation."
    )
  endif()
  message("  Output of test compile is in ${COMPILE_FAIL_OUTPUT}")
  message(FATAL_ERROR "stopping")
else()
  message(STATUS "C++17 standard library supported")
endif()
