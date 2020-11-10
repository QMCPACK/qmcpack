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
SET(QMCAPPDIR
    QMCAppBase.cpp
    QMCMainState.cpp
    QMCMain.cpp
    CoupledMC.cpp
    accelerators.cpp
)
IF(USE_OBJECT_TARGET)
  ADD_LIBRARY(qmc OBJECT ${QMCAPPDIR})
ELSE()
  ADD_LIBRARY(qmc ${QMCAPPDIR})
ENDIF()
TARGET_LINK_LIBRARIES(qmc PUBLIC qmcdriver)

####################################
# QMCPACK application
####################################
ADD_EXECUTABLE(qmcpack qmcapp.cpp)

# Clean up git-rev-tmp.h after the build is finished
ADD_CUSTOM_COMMAND(TARGET qmcpack
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E remove ${GITREV_TMP})

TARGET_LINK_LIBRARIES(qmcpack qmc qmcdriver)
IF(USE_OBJECT_TARGET)
  TARGET_LINK_LIBRARIES(qmcpack qmcestimators qmcham qmcwfs qmcparticle qmcutil containers platform_device platform_omptarget)
ENDIF()

IF(BUILD_AFQMC)
  TARGET_LINK_LIBRARIES(qmcpack afqmc)
ENDIF(BUILD_AFQMC)

IF(BUILD_FCIQMC)
  TARGET_LINK_LIBRARIES(qmcpack fciqmc)
ENDIF(BUILD_FCIQMC)

INSTALL(TARGETS qmcpack
        RUNTIME DESTINATION bin)
