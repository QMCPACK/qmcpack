#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(COMM_SRCS Communicate.cpp AppAbort.cpp MPIObjectBase.cpp)

add_library(message ${COMM_SRCS})
target_link_libraries(message PUBLIC platform_host_runtime)

if(HAVE_MPI)
  target_link_libraries(message PUBLIC MPI::MPI_CXX Boost::boost)
  target_include_directories(message PUBLIC ${qmcpack_SOURCE_DIR}/external_codes/mpi_wrapper)
endif()

add_library(catch_main catch_main.cpp)
target_include_directories(catch_main PUBLIC "${PROJECT_SOURCE_DIR}/external_codes/catch")
target_compile_definitions(catch_main PUBLIC "CATCH_CONFIG_ENABLE_BENCHMARKING")
target_link_libraries(catch_main PRIVATE message platform_runtime)
if(HAVE_MPI)
  target_compile_definitions(catch_main PRIVATE "CATCH_MAIN_HAVE_MPI")
endif()

add_library(catch_main_no_mpi catch_main.cpp)
target_include_directories(catch_main_no_mpi PUBLIC "${PROJECT_SOURCE_DIR}/external_codes/catch")
target_link_libraries(catch_main_no_mpi PRIVATE platform_runtime)

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
