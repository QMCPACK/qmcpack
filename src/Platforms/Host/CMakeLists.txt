#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2022 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////


# platform_host_runtime is the target for host runtime system which includes
# interaction with OS libraries
set(HOST_SRCS sysutil.cpp InfoStream.cpp OutputManager.cpp)
add_library(platform_host_runtime ${HOST_SRCS})
target_include_directories(platform_host_runtime PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
