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

set(SRC_DIR numerics)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

set(UTEST_SRCS
    test_grid_functor.cpp
    test_stdlib.cpp
    test_bessel.cpp
    test_ylm.cpp
    test_gaussian_basis.cpp
    test_cartesian_tensor.cpp
    test_soa_cartesian_tensor.cpp
    test_transform.cpp
    test_min_oned.cpp
    test_OneDimCubicSplineLinearGrid.cpp
    test_one_dim_cubic_spline.cpp)

# Run gen_gto.py to create these files.  They may take a long time to compile.
#SET(UTEST_SRCS ${UTEST_SRCS} test_full_cartesian_tensor.cpp test_full_soa_cartesian_tensor.cpp)

add_executable(${UTEST_EXE} ${UTEST_SRCS})

target_link_libraries(${UTEST_EXE} catch_main containers qmcio)

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
