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

# Hide potentially lengthy initial Catch2 output such as
# -- Performing Test HAVE_FLAG__ffile_prefix_map__workspaces_qmcpack_external_codes_Catch2__
# Avoid changing CMAKE_REQUIRED_QUIET defined status or value
if(DEFINED CMAKE_REQUIRED_QUIET)
  set(QMC_CMAKE_REQUIRED_QUIET_SAVED "${CMAKE_REQUIRED_QUIET}")
  set(QMC_CMAKE_REQUIRED_QUIET_WAS_DEFINED TRUE)
endif()
set(CMAKE_REQUIRED_QUIET TRUE)
add_subdirectory(Catch2 EXCLUDE_FROM_ALL)
# Catch2's floating-point matchers use infinities, which -ffast-math otherwise assumes cannot occur.
target_compile_options(Catch2 PRIVATE $<$<CXX_COMPILER_ID:Clang,GNU>:-fno-finite-math-only>)
if(QMC_CMAKE_REQUIRED_QUIET_WAS_DEFINED)
  set(CMAKE_REQUIRED_QUIET "${QMC_CMAKE_REQUIRED_QUIET_SAVED}")
else()
  unset(CMAKE_REQUIRED_QUIET)
endif()
unset(QMC_CMAKE_REQUIRED_QUIET_SAVED)
unset(QMC_CMAKE_REQUIRED_QUIET_WAS_DEFINED)

add_library(qmc_external_catch2 INTERFACE)
target_link_libraries(qmc_external_catch2 INTERFACE Catch2::Catch2)
