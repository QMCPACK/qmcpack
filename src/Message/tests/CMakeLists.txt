#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Mark Dewing, mdewing@ganl.gov Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(SRC_DIR message)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

add_executable(${UTEST_EXE} test_communciate.cpp)
target_link_libraries(${UTEST_EXE} PUBLIC message catch_main)

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)

if(HAVE_MPI)
  set(UTEST_EXE test_${SRC_DIR}_mpi)
  set(UTEST_NAME deterministic-unit_test_${SRC_DIR}_mpi)
  #this is dependent on the directory creation and sym linking of earlier driver tests
  set(MPI_UTILITY_TEST_SRC test_mpi_exception_wrapper.cpp)
  add_executable(${UTEST_EXE} ${MPI_UTILITY_TEST_SRC})
  #Way too many depenedencies make for very slow test linking
  target_link_libraries(${UTEST_EXE} PUBLIC message catch_main)
  # Right now the unified driver mpi tests are hard coded for 3 MPI ranks
  add_unit_test(${UTEST_NAME} 3 1 $<TARGET_FILE:${UTEST_EXE}>)
endif(HAVE_MPI)
