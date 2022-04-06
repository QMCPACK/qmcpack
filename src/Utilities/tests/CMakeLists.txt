#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////

set(SRC_DIR utilities)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

add_executable(
  ${UTEST_EXE}
  test_rng.cpp
  test_parser.cpp
  test_timer.cpp
  test_runtime_manager.cpp
  test_prime_set.cpp
  test_partition.cpp
  test_ResourceCollection.cpp
  test_infostream.cpp
  test_project_data.cpp
  test_rng_control.cpp
  test_output_manager.cpp
  test_ModernStringUtils.cpp
  test_StlPrettyPrint.cpp
  test_StdRandom.cpp)
target_link_libraries(${UTEST_EXE} catch_main qmcutil)

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)

subdirs(for_testing)
