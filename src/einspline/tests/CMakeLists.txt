#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////

set(SRC_DIR einspline)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

set(SRCS test_one.cpp test_3d.cpp)
if(QMC_CUDA)
  set(SRCS ${SRCS} test_cuda.cu)
  if(NOT QMC_CUDA2HIP)
    add_library(cudatests test_cuda.cu)
  else()
    hip_add_library(cudatests test_cuda.cu)
  endif(NOT QMC_CUDA2HIP)
  set(CUDATESTLIB cudatests)
endif()

add_executable(${UTEST_EXE} ${SRCS})
target_link_libraries(${UTEST_EXE} catch_main ${CUDATESTLIB} einspline qmcutil)

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
