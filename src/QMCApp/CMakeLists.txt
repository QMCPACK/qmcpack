#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
#//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
#//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
#//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
#//                    Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
#//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//
#// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////

####################################
# create libqmc
####################################
set(QMCAPPDIR QMCAppBase.cpp QMCMainState.cpp QMCMain.cpp CoupledMC.cpp)
if(USE_OBJECT_TARGET)
  add_library(qmc OBJECT ${QMCAPPDIR})
else()
  add_library(qmc ${QMCAPPDIR})
endif()
target_link_libraries(qmc PUBLIC qmcdriver)
if(BUILD_AFQMC)
  target_link_libraries(qmc PUBLIC afqmc)
endif(BUILD_AFQMC)

####################################
# QMCPACK application
####################################
add_executable(qmcpack qmcapp.cpp)
set_target_properties(qmcpack PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${qmcpack_BINARY_DIR}/bin)
if(QMC_COMPLEX)
  set_target_properties(qmcpack PROPERTIES OUTPUT_NAME "qmcpack_complex")
endif(QMC_COMPLEX)

# Clean up git-rev-tmp.h after the build is finished
add_custom_command(
  TARGET qmcpack
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E remove ${GITREV_TMP})

target_link_libraries(qmcpack qmc qmcdriver)
if(USE_OBJECT_TARGET)
  target_link_libraries(
    qmcpack
    qmcestimators
    qmcham
    qmcwfs
    qmcparticle
    qmcutil
    platform_omptarget_LA)
endif()

add_test_target_in_output_location(qmcpack bin)
install(TARGETS qmcpack RUNTIME DESTINATION bin)
