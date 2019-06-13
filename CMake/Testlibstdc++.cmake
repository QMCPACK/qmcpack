
# Test that if a C++ compiler is compatiable with the libstdc++ in use

set(TEST_LIBSTDCXX_SOURCE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/try_libstdcxx.cpp)
file(WRITE ${TEST_LIBSTDCXX_SOURCE}
"// Test the compatibility between the compiler and the libstdc++ from GNU
#include <cstdio>

int main(int argc, char **argv)
{
#if ( ( __INTEL_COMPILER == 1800 ) && (  _GLIBCXX_RELEASE > 7 ) )
#error too new libstdc++ from GNU for Intel 18, use GNU version <= 7
#endif
    return 0;
}
")


if (${CMAKE_VERSION} VERSION_LESS 3.8.0)
  try_compile(LIBSTDCXX_OKAY ${CMAKE_BINARY_DIR}
              ${TEST_LIBSTDCXX_SOURCE}
              COMPILE_DEFINITIONS ${CMAKE_CXX14_STANDARD_COMPILE_OPTION}
              OUTPUT_VARIABLE COMPILE_OUTPUT)
else()
  try_compile(LIBSTDCXX_OKAY ${CMAKE_BINARY_DIR}
              ${TEST_LIBSTDCXX_SOURCE}
              CXX_STANDARD 14
              CXX_STANDARD_REQUIRED ON
              OUTPUT_VARIABLE COMPILE_OUTPUT)
endif()


if (NOT LIBSTDCXX_OKAY)
  set(COMPILE_FAIL_OUTPUT libstdc++_compile_fail.txt)
  file(WRITE "${CMAKE_BINARY_DIR}/${COMPILE_FAIL_OUTPUT}" "${COMPILE_OUTPUT}")
  message(STATUS "libstdc++/C++ compiler version compatibility check failed")
  message(FATAL_ERROR "${COMPILE_OUTPUT}")
else()
  message(STATUS "libstdc++/C++ compiler version compatibility check pass")
endif()
