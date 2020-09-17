#////////////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2020 QMCPACK developers.
#//
#// File developed by: Fionn Malone malone14@llnl.gov, Lawrence Livermore National Laboratory
#//
#// File created by: Fionn Malone malone14@llnl.gov, Lawrence Livermore National Laboratory
#////////////////////////////////////////////////////////////////////////////////////////////

MESSAGE("Building AFQMC performance executable ")

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${qmcpack_BINARY_DIR}/bin)

ADD_EXECUTABLE(afqmc-qmc-performance performance.cpp)
TARGET_LINK_LIBRARIES(afqmc-qmc-performance afqmc Math::BLAS_LAPACK)
