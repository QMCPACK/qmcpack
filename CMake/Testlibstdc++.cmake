# Test that if a C++ compiler is compatiable with the libstdc++ in use

# Test "#include <cstdio>" before version compatibility checks.
include(CheckIncludeFileCXX)
check_include_file_cxx(cstdio INCLUDE_CSTDIO_WORKS)
if(NOT INCLUDE_CSTDIO_WORKS)
  unset(INCLUDE_CSTDIO_WORKS CACHE)
  message(FATAL_ERROR "`#include <cstdio>` test failed! Please provide a working C++ compiler.")
endif()

try_compile(
  LIBSTDCXX_OKAY
  ${CMAKE_BINARY_DIR}
  ${PROJECT_CMAKE}/try_compile_sources/check_libstdcxx.cpp
  CXX_STANDARD
  ${QMC_CXX_STANDARD}
  CXX_STANDARD_REQUIRED
  ON
  OUTPUT_VARIABLE COMPILE_OUTPUT)

if(NOT LIBSTDCXX_OKAY)
  set(COMPILE_FAIL_OUTPUT libstdc++_compile_fail.txt)
  file(WRITE "${CMAKE_BINARY_DIR}/${COMPILE_FAIL_OUTPUT}" "${COMPILE_OUTPUT}")
  message(STATUS "libstdc++/C++ compiler version compatibility check failed")
  message(FATAL_ERROR "${COMPILE_OUTPUT}")
else()
  message(STATUS "libstdc++/C++ compiler version compatibility check pass")
endif()
