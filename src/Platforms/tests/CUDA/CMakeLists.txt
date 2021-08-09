#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2019 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

set(SRC_DIR CUDA)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

cuda_add_library(cuda_device_value_test_kernels test_device_value_kernels.cu)
target_link_libraries(cuda_device_value_test_kernels PUBLIC platform_device)

add_executable(${UTEST_EXE} test_CUDAallocator.cpp)
target_link_libraries(${UTEST_EXE} platform_device containers catch_main)
if(USE_OBJECT_TARGET)
  target_link_libraries(${UTEST_EXE} platform_omptarget)
endif()

add_unit_test(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
