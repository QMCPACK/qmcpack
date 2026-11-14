#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
#//		       Anouar Benali, benali@anl.gov, Argonne National Laboratory
#//		       Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////

set(SRCS bspline_create.cpp bspline_data.cpp multi_bspline_create.cpp multi_bspline_copy.cpp bspline_eval_d_std.cpp)

add_library(einspline ${SRCS})
target_compile_definitions(einspline INTERFACE "HAVE_EINSPLINE")

#ADD_EXECUTABLE(time_multi time_multi_new.c)
#target_link_libraries(time_multi einspline)
#add_dependencies(time_multi ${PROJECT_BINARY_HOME}/include/einspline/bspline.h)
#ADD_EXECUTABLE(test_multi_single test_multi_single.cpp)
#target_link_libraries(test_multi_single einspline)

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
