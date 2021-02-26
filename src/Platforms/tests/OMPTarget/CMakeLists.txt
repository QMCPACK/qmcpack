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


IF(ENABLE_OFFLOAD)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

SET(SRC_DIR omptarget)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME deterministic-unit_test_${SRC_DIR})

ADD_EXECUTABLE(${UTEST_EXE} test_vector.cpp test_math.cpp test_deep_copy.cpp test_class_member.cpp test_runtime_mem.cpp)
TARGET_LINK_LIBRARIES(${UTEST_EXE} platform_device catch_main)
IF(USE_OBJECT_TARGET)
  TARGET_LINK_LIBRARIES(${UTEST_EXE}  platform_omptarget)
ENDIF()

ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")

ENDIF(ENABLE_OFFLOAD)
