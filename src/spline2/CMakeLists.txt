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
target_link_libraries(spline2_omptarget PUBLIC einspline platform_runtime)

set(SPLINE2_SRC MultiBspline.cpp MultiBspline1D.cpp SplineUtils.cpp)
if(HAVE_MPI)
  list(APPEND SPLINE2_SRC MultiBsplineMPIShared.cpp)
endif()

if(USE_OBJECT_TARGET)
  add_library(spline2 OBJECT ${SPLINE2_SRC})
else()
  add_library(spline2 ${SPLINE2_SRC})
endif()
target_link_libraries(spline2 PUBLIC spline2_omptarget message qmcio)

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
