
# Test that if a C++ compiler is compatiable with the libstdc++ in use

set(TEST_LIBSTDCXX_SOURCE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/try_libstdcxx.cpp)
file(WRITE ${TEST_LIBSTDCXX_SOURCE}
"// Test the compatibility between the compiler and the libstdc++ from GNU
#include <cstdio>

int main(int argc, char **argv)
{
// Unfortunately this check doesn't work because _GLIBCXX_RELEASE appeared in the GCC 7.1 release.
// It is kept here as an exmple for the future.
#if ( defined(__INTEL_COMPILER) && (  _GLIBCXX_RELEASE < 7 ) )
#error too old libstdc++ from GNU for Intel compilers, use GNU version >= 7
#endif
// libstdc++ from GCC 8 is bad for Intel 19 in both C++14 and C++17
#if ( ( __INTEL_COMPILER == 1900 ) && (  _GLIBCXX_RELEASE > 7 ) )
#error too new libstdc++ from GNU for Intel 19, use GNU version <= 7
#endif
#if ( ( __INTEL_COMPILER == 1910 ) && (  _GLIBCXX_RELEASE > 9 ) )
#error too new libstdc++ from GNU for Intel 20, use GNU version <= 9
#endif
    return 0;
}
")


try_compile(LIBSTDCXX_OKAY ${CMAKE_BINARY_DIR}
            ${TEST_LIBSTDCXX_SOURCE}
            CXX_STANDARD 14
            CXX_STANDARD_REQUIRED ON
            OUTPUT_VARIABLE COMPILE_OUTPUT)


if (NOT LIBSTDCXX_OKAY)
  set(COMPILE_FAIL_OUTPUT libstdc++_compile_fail.txt)
  file(WRITE "${CMAKE_BINARY_DIR}/${COMPILE_FAIL_OUTPUT}" "${COMPILE_OUTPUT}")
  message(STATUS "libstdc++/C++ compiler version compatibility check failed")
  message(FATAL_ERROR "${COMPILE_OUTPUT}")
else()
  message(STATUS "libstdc++/C++ compiler version compatibility check pass")
endif()
