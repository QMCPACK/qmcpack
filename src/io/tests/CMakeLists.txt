# //////////////////////////////////////////////////////////////////////////////////////
# // This file is distributed under the University of Illinois/NCSA Open Source License.
# // See LICENSE file in top directory for details.
# //
# // Copyright (c) 2022 QMCPACK developers
# //
# // File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
# //
# // File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
# //////////////////////////////////////////////////////////////////////////////////////

set(SRC_DIR io)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

# Directory where input file is copied to, and working directory for unit test
set(UTEST_DIR ${CMAKE_CURRENT_BINARY_DIR})

add_executable(${UTEST_EXE} test_InputNode.cpp)
target_include_directories(${UTEST_EXE} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../")
target_link_libraries(${UTEST_EXE} catch_main qmcutil)

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
