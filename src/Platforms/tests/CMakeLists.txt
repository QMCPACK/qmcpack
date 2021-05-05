#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////


IF (BUILD_UNIT_TESTS)
  IF(ENABLE_CUDA)
    SUBDIRS(CUDA)
  ENDIF()
  SUBDIRS(OMPTarget)
  SUBDIRS(CPU)
ENDIF()

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

if (ENABLE_CUDA)
  CUDA_ADD_LIBRARY(cuda_test_kernels test_kernels.cu)
  TARGET_LINK_LIBRARIES(cuda_test_kernels PUBLIC platform_device)
  SET(UTEST_EXE test_DeviceValue)
  SET(UTEST_NAME deterministic-unit_${UTEST_EXE})
  ADD_EXECUTABLE(${UTEST_EXE} test_DeviceValue.cpp )
  TARGET_LINK_LIBRARIES(${UTEST_EXE} platform_device cuda_test_kernels catch_main qmcutil)
  IF(USE_OBJECT_TARGET)
    TARGET_LINK_LIBRARIES(${UTEST_EXE} platform_omptarget)
  ENDIF()
  TARGET_INCLUDE_DIRECTORIES(${UTEST_EXE} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/..")
  ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
endif()
