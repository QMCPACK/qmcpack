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

set(SRC_DIR io_ohmmsdata)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

# Directory where input file is copied to, and working directory for unit test
set(UTEST_DIR ${CMAKE_CURRENT_BINARY_DIR})

add_executable(${UTEST_EXE} test_xml.cpp)
target_include_directories(${UTEST_EXE} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../")
target_link_libraries(${UTEST_EXE} catch_main qmcutil qmcio_xml)

execute_process(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/simple.xml" ${UTEST_DIR})

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
