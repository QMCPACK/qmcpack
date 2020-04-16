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


SUBDIRS(OpenMP)

SET(PLATFORMS_SRCS
    accelerators.cpp
)

ADD_LIBRARY(qmcplatforms ${PLATFORMS_SRCS})
TARGET_LINK_LIBRARIES(qmcplatforms PUBLIC qmcomp qmcutil)
