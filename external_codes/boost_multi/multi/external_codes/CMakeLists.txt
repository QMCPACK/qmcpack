#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////


add_library(qmc_external_mpi_wrapper INTERFACE)
target_include_directories(qmc_external_mpi_wrapper INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/mpi_wrapper/mpi3/include)
target_link_libraries(qmc_external_mpi_wrapper INTERFACE MPI::MPI_CXX Boost::boost)

add_library(qmc_external_multi INTERFACE)
target_include_directories(qmc_external_multi INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/boost_multi/multi/include)

add_library(qmc_external_thrust INTERFACE)
target_include_directories(qmc_external_thrust INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/thrust)

add_library(qmc_external_catch2 INTERFACE)
target_include_directories(qmc_external_catch2 INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/catch)
