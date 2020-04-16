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


SET(HOST_SRCS
    sysutil.cpp
)

SET(DEVICE_SRCS
    accelerators.cpp
)

ADD_LIBRARY(platform_host ${HOST_SRCS})

ADD_LIBRARY(platform_device ${DEVICE_SRCS})
TARGET_LINK_LIBRARIES(platform_device PUBLIC message)
