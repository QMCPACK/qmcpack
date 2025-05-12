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


add_library(spline2_omptarget OBJECT MultiBsplineOffload.cpp)
target_link_libraries(spline2_omptarget PUBLIC platform_runtime)

if(USE_OBJECT_TARGET)
  add_library(spline2 OBJECT MultiBspline.cpp)
else()
  add_library(spline2 MultiBspline.cpp)
endif()
target_link_libraries(spline2 PUBLIC spline2_omptarget)

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
