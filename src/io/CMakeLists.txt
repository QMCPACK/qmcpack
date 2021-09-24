#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2020 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

add_subdirectory(hdf)
add_subdirectory(OhmmsData)

add_library(qmcio INTERFACE)
target_include_directories(qmcio INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(qmcio INTERFACE qmcio_hdf qmcio_xml)
