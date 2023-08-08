#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2023 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(OMP_RT_SRCS OMPallocator.cpp OMPDeviceManager.cpp)
set(OMP_LA_SRCS ompBLAS.cpp)

add_library(platform_omptarget_runtime ${OMP_RT_SRCS})
if(USE_OBJECT_TARGET)
  add_library(platform_omptarget_LA OBJECT ${OMP_LA_SRCS})
else()
  add_library(platform_omptarget_LA ${OMP_LA_SRCS})
endif()
target_link_libraries(platform_omptarget_runtime PRIVATE platform_host_runtime)
target_link_libraries(platform_omptarget_LA PUBLIC platform_omptarget_runtime)

target_compile_options(platform_omptarget_LA
                       PRIVATE "$<$<BOOL:${ENABLE_OFFLOAD_CLANG_DEBUG_O3}>:$<$<CONFIG:DEBUG>:-O3>>")

if(ENABLE_CUDA AND QMC_OFFLOAD_MEM_ASSOCIATED)
  target_link_libraries(platform_omptarget_runtime PUBLIC platform_cuda_runtime)
endif()

if(NOT QMC_COMPLEX)
  target_compile_definitions(platform_omptarget_LA PUBLIC OPENMP_NO_COMPLEX)
else()
  # Test if C++ compiler supports OpenMP user defined reduction on complex type
  try_compile(
    OMP_UDR_COMPLEX_OKAY ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/test_user_defined_reduction_complex.cpp
    CMAKE_FLAGS "${CMAKE_CXX_FLAGS}"
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  if(NOT OMP_UDR_COMPLEX_OKAY)
    set(COMPILE_FAIL_OUTPUT omp_UDR_complex_fail.txt)
    file(WRITE "${CMAKE_BINARY_DIR}/${COMPILE_FAIL_OUTPUT}" "${COMPILE_OUTPUT}")
    message(STATUS "OpenMP user-defined reduction functionality compilation check failed!"
                   "See compiler output at ${COMPILE_FAIL_OUTPUT}")

    # Test if C++ compiler supports OpenMP implicit reduction on complex type
    try_compile(
      OMP_IMPLICIT_REDUCTION_COMPLEX_OKAY ${CMAKE_CURRENT_BINARY_DIR}
      ${CMAKE_CURRENT_SOURCE_DIR}/test_user_defined_reduction_complex.cpp
      CMAKE_FLAGS "${CMAKE_CXX_FLAGS}"
      COMPILE_DEFINITIONS "-DOPENMP_NO_UDR"
      OUTPUT_VARIABLE COMPILE_OUTPUT)

    if(NOT OMP_IMPLICIT_REDUCTION_COMPLEX_OKAY)
      set(COMPILE_FAIL_OUTPUT omp_implicit_reduction_complex_fail.txt)
      file(WRITE "${CMAKE_BINARY_DIR}/${COMPILE_FAIL_OUTPUT}" "${COMPILE_OUTPUT}")
      message(FATAL_ERROR "OpenMP implicit reduction on complex functionality compilation check failed!"
                          "See compiler output at ${COMPILE_FAIL_OUTPUT}")
    else()
      target_compile_definitions(platform_omptarget_LA PUBLIC OPENMP_NO_UDR)
      message(STATUS "OpenMP implicit reduction functionality compilation check pass")
    endif()
  else()
    message(STATUS "OpenMP user-defined reduction functionality compilation check pass")
  endif()
endif()
