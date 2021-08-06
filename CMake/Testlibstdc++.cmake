# Test that if a C++ compiler is compatiable with the libstdc++ in use

try_compile(
  LIBSTDCXX_OKAY
  ${CMAKE_BINARY_DIR}
  ${PROJECT_CMAKE}/try_compile_sources/check_libstdcxx.cpp
  CXX_STANDARD
  14
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

find_package(Filesystem
  COMPONENTS Experimental Final
  REQUIRED)
