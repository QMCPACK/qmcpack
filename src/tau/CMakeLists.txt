#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
#//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
#//
#// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////
    
    

include_directories(${TAUDIR}/include)

ADD_DEFINITIONS (-DTAU_LIBRARY_SOURCE)

#-DTAU_MPIFILE
SET(tausrcs
  ${TAU_HOME}/src/Profile/Profiler.cpp
  ${TAU_HOME}/src/Profile/TulipTimers.cpp
  ${TAU_HOME}/src/Profile/UserEvent.cpp
  ${TAU_HOME}/src/Profile/FunctionInfo.cpp
  ${TAU_HOME}/src/Profile/RtsLayer.cpp
  ${TAU_HOME}/src/Profile/RtsThread.cpp
  ${TAU_HOME}/src/Profile/TauCAPI.cpp
  ${TAU_HOME}/src/Profile/TauFAPI.cpp
  ${TAU_HOME}/src/Profile/TauMapping.cpp
  ${TAU_HOME}/src/Profile/TauHooks.cpp
  ${TAU_HOME}/src/Profile/TauHandler.cpp
  ${TAU_HOME}/src/Profile/TauMemory.cpp
  ${TAU_HOME}/src/Profile/MetaData.cpp
  ${TAU_HOME}/src/Profile/TauEnv.cpp
  ${TAU_HOME}/src/Profile/TauInit.cpp
  ${TAU_HOME}/src/Profile/TauCallPath.cpp
  ${TAU_HOME}/src/Profile/OpenMPLayer.cpp
  ${TAU_HOME}/src/Profile/TauKojakOpariDisable.cpp
  ${TAU_HOME}/src/Profile/Comp_gnu.cpp
 )
#  ${TAU_HOME}/src/Profile/Comp_pgi.cpp
#  ${TAU_HOME}/src/Profile/Comp_xl.cpp
#  ${TAU_HOME}/src/Profile/Comp_intel9.cpp

ADD_LIBRARY(tau ${tausrcs})

# TauMpi.c
# ClockSync.c
# ProfileMerge.c
# TauFMpi.c
# TauMpiExtensions.c
