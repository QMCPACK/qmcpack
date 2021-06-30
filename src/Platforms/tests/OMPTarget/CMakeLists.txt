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

if(ENABLE_OFFLOAD)

  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

  set(SRC_DIR omptarget)
  set(UTEST_EXE test_${SRC_DIR})
  set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

  add_executable(${UTEST_EXE} test_vector.cpp test_math.cpp test_deep_copy.cpp test_class_member.cpp
                              test_runtime_mem.cpp)
  target_link_libraries(${UTEST_EXE} platform_device catch_main)
  if(USE_OBJECT_TARGET)
    target_link_libraries(${UTEST_EXE} platform_omptarget)
  endif()

  add_unit_test(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")

endif(ENABLE_OFFLOAD)
