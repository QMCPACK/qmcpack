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

set(SRCS )

SET(SRCS ${SRCS}
  bspline_create.c
  bspline_data.c
  multi_bspline_create.c
  multi_nubspline_create.c
  nubspline_create.c
  nubasis.c
  nugrid.c
  multi_bspline_copy.c
)

SET(SRCS ${SRCS}
  bspline_eval_d_std.cpp
  nubspline_eval_d_std.cpp
  multi_bspline_eval_s_std3.cpp
  multi_bspline_eval_d_std3.cpp
  multi_bspline_eval_z_std3.cpp
  multi_nubspline_eval_z_std.cpp
)


if(HAVE_CUDA)
  SET(SRCS  ${SRCS}
    multi_bspline_create_cuda.cu
    bspline_create_cuda.cu
    )
  CUDA_ADD_LIBRARY(einspline ${SRCS})
else()
  ADD_LIBRARY(einspline ${SRCS})
endif()

#ADD_EXECUTABLE(time_multi time_multi_new.c)
#target_link_libraries(time_multi einspline)
#add_dependencies(time_multi ${PROJECT_BINARY_HOME}/include/einspline/bspline.h)
#ADD_EXECUTABLE(test_multi_single test_multi_single.cpp)
#target_link_libraries(test_multi_single einspline)
