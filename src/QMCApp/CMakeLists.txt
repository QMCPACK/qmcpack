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
    
    
PROJECT(qmcapp)
# formic dependencies
#FIND_PATH(LMYENGINE_INCLUDE_DIR formic/fqmc/engine/var_dependencies.h)
#MESSAGE("-- LMYENGINE_INCLUDE_DIR=${LMYENGINE_INCLUDE_DIR}")
#INCLUDE_DIRECTORIES(${LMYENGINE_INCLUDE_DIR})

####################################
# create libqmc 
####################################
SET(QMCAPPDIR 
    ParticleSetPool.cpp
    WaveFunctionPool.cpp
    HamiltonianPool.cpp
    InitMolecularSystem.cpp
    QMCAppBase.cpp
    QMCDriverFactory.cpp
    QMCMain.cpp
    CoupledMC.cpp
)
ADD_LIBRARY(qmc ${QMCAPPDIR})

####################################
# list of QMC applications
####################################
SET(QMCSRCS qmcpack)

FOREACH(p ${QMCSRCS})

# QMCPACK historically called qmcapp, and still built from qmcapp.cpp
# Override source name for executable in this case, otherwise
# executablename.cpp used for main()
 IF ( ${p} STREQUAL "qmcpack" )
  ADD_EXECUTABLE(${p} qmcapp.cpp) 
ELSE()
  ADD_EXECUTABLE(${p} ${p}.cpp) 
ENDIF()

  TARGET_LINK_LIBRARIES(${p} qmc qmcdriver qmcham qmcwfs qmcbase qmcutil adios_config)

  FOREACH(l ${QMC_UTIL_LIBS})
    TARGET_LINK_LIBRARIES(${p} ${l})
  ENDFOREACH(l ${QMC_UTIL_LIBS})

  if(ENABLE_TAU_PROFILE)
    TARGET_LINK_LIBRARIES(${p} tau)
  endif(ENABLE_TAU_PROFILE)

  IF(MPI_LIBRARY)
    TARGET_LINK_LIBRARIES(${p} ${MPI_LIBRARY})
  ENDIF(MPI_LIBRARY)

  IF(BUILD_LMYENGINE_INTERFACE)
    TARGET_LINK_LIBRARIES(${p} formic_utils)
  ENDIF(BUILD_LMYENGINE_INTERFACE)

ENDFOREACH(p ${QMCSRCS})

