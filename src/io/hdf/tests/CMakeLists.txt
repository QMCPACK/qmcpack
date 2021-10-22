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


set(SRC_DIR io_hdf5)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

# Directory where input file is copied to, and working directory for unit test
set(UTEST_DIR ${CMAKE_CURRENT_BINARY_DIR})

set(UTEST_SRC test_hdf_archive.cpp test_hdf_parallel.cpp test_hdf_reshape.cpp test_hdf_hyperslab.cpp)
if(BUILD_AFQMC)
  set(UTEST_SRC ${UTEST_SRC} test_hdf_multiarray.cpp)
endif(BUILD_AFQMC)

add_executable(${UTEST_EXE} ${UTEST_SRC})
target_link_libraries(${UTEST_EXE} catch_main qmcio)

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
