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


SET(IO_HDF5_SRC hdf_archive.cpp)

ADD_LIBRARY(qmcio_hdf ${IO_HDF5_SRC})

#TARGET_INCLUDE_DIRECTORIES(qmcio_hdf PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
TARGET_LINK_LIBRARIES(qmcio_hdf PUBLIC containers IO::HDF5 message)
