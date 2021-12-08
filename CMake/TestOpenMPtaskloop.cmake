# Test if C++ compiler supports OpenMP taskloop construct

try_compile(OMP_TASKLOOP_OKAY ${CMAKE_BINARY_DIR} ${PROJECT_CMAKE}/try_compile_sources/check_openmp_taskloop.cpp
            OUTPUT_VARIABLE COMPILE_OUTPUT)

if(NOT OMP_TASKLOOP_OKAY)
  set(COMPILE_FAIL_OUTPUT omp_taskloop_compile_fail.txt)
  file(WRITE "${CMAKE_BINARY_DIR}/${COMPILE_FAIL_OUTPUT}" "${COMPILE_OUTPUT}")
  message(STATUS "OpenMP taskloop functionality check failed!" "See compiler output at ${COMPILE_FAIL_OUTPUT}")
else()
  message(STATUS "OpenMP taskloop functionality check pass")
endif()
