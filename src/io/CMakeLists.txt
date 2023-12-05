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

# qmcio is intended to be an interface library servicing codes outside src/io.
# header files inclusion from outside src/io needs to treat src/io being the root for
# the header file inclusion.
# qmcio_hdf qmcio_xml are targets for internal use only.
add_library(qmcio INTERFACE)
target_include_directories(qmcio INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(qmcio INTERFACE qmcio_hdf qmcio_xml)
