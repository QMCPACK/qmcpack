#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////
    
INCLUDE("${qmcpack_SOURCE_DIR}/CMake/macros.cmake")
    

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

SET(SRC_DIR type_traits)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME deterministic-unit_test_${SRC_DIR})
SET(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/type_traits)

EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")

#this is the only unit test program for type_traits
#I'm not going to figure out how to deal with that fact.
#If you write another test program in here you can make
#test_cuda_type_aliases.cpp an optional src for the add_executable command
IF(QMC_CUDA)
  ADD_EXECUTABLE(${UTEST_EXE} test_cuda_types.cpp)
  ADD_UNIT_TEST(${UTEST_NAME} "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
  SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS "unit" WORKING_DIRECTORY ${UTEST_DIR})
ENDIF(QMC_CUDA)
