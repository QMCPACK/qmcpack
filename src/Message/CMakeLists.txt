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
target_link_libraries(message PUBLIC platform_host)

if(HAVE_MPI)
  target_link_libraries(message PUBLIC MPI::MPI_CXX Boost::boost)
  set_property(
    TARGET message
    APPEND
    PROPERTY INCLUDE_DIRECTORIES "${qmcpack_SOURCE_DIR}/external_codes/mpi_wrapper")
  set_property(TARGET message APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                              "${qmcpack_SOURCE_DIR}/external_codes/mpi_wrapper")
endif()

add_library(catch_main catch_main.cpp)
target_include_directories(catch_main PUBLIC "${PROJECT_SOURCE_DIR}/external_codes/catch")
# This plus using the definition of the catch_benchmark_main as INTERFACE target
# would allow just building catch main once for benchmarking and mpi enabled testing
# But it causes a segv on power9 with clang12 and on alcf-gnu-mkl18
# TARGET_COMPILE_DEFINITIONS(catch_main PRIVATE "CATCH_CONFIG_ENABLE_BENCHMARKING")
target_link_libraries(catch_main message)
if(HAVE_MPI)
  set_property(TARGET catch_main APPEND PROPERTY COMPILE_DEFINITIONS "CATCH_MAIN_HAVE_MPI")
endif()

add_library(catch_main_no_mpi catch_main.cpp)
target_include_directories(catch_main_no_mpi PUBLIC "${PROJECT_SOURCE_DIR}/external_codes/catch")

# ADD_LIBRARY(catch_benchmark_main INTERFACE)
# SET_TARGET_PROPERTIES(catch_benchmark_main PROPERTIES INTERFACE_COMPILE_DEFINITIONS "CATCH_CONFIG_ENABLE_BENCHMARKING")
# TARGET_LINK_LIBRARIES(catch_benchmark_main INTERFACE catch_main)

add_library(catch_benchmark_main catch_main.cpp)
target_compile_definitions(catch_benchmark_main PUBLIC "CATCH_CONFIG_ENABLE_BENCHMARKING")
target_link_libraries(catch_benchmark_main message)
if(HAVE_MPI)
  set_property(TARGET catch_benchmark_main APPEND PROPERTY COMPILE_DEFINITIONS "CATCH_MAIN_HAVE_MPI")
endif()
target_include_directories(catch_benchmark_main PUBLIC "${PROJECT_SOURCE_DIR}/external_codes/catch")
