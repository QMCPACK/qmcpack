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

message("Building AFQMC performance executable ")

add_executable(qmc-afqmc-performance performance.cpp)
target_link_libraries(qmc-afqmc-performance afqmc platform_LA)
