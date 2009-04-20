include_directories(${TAUDIR}/include)

ADD_DEFINITIONS (-DPROFILING_ON
  -DTAU_STDCXXLIB
  -DTAU_DOT_H_LESS_HEADERS
  -DTAU_MPI
  -DTAU_MPI_THREADED
  -DTAU_CALLPATH
  -DTAU_MPIGREQUEST
  -DTAU_MPIDATAREP
  -DTAU_MPIERRHANDLER
  -DTAU_MPIATTRFUNCTION
  -DTAU_MPITYPEEX
  -DTAU_MPIADDERROR
  -DTAU_LARGEFILE64_SOURCE
  -DTAU_BFD
  )

IF(ENABLE_OPENMP)
  ADD_DEFINITIONS (-DTAU_OPENMP)
IF(ENABLE_OPENMP)

ADD_DEFINITIONS (-DTAU_LIBRARY_SOURCE)

#-DTAU_MPIFILE
SET(tausrcs
  ${TAUDIR}/src/Profile/Profiler.cpp
  ${TAUDIR}/src/Profile/TulipTimers.cpp
  ${TAUDIR}/src/Profile/UserEvent.cpp
  ${TAUDIR}/src/Profile/FunctionInfo.cpp
  ${TAUDIR}/src/Profile/RtsLayer.cpp
  ${TAUDIR}/src/Profile/RtsThread.cpp
  ${TAUDIR}/src/Profile/TauCAPI.cpp
  ${TAUDIR}/src/Profile/TauFAPI.cpp
  ${TAUDIR}/src/Profile/TauMapping.cpp
  ${TAUDIR}/src/Profile/TauHooks.cpp
  ${TAUDIR}/src/Profile/TauHandler.cpp
  ${TAUDIR}/src/Profile/TauMemory.cpp
  ${TAUDIR}/src/Profile/MetaData.cpp
  ${TAUDIR}/src/Profile/TauEnv.cpp
  ${TAUDIR}/src/Profile/TauInit.cpp
  ${TAUDIR}/src/Profile/TauCallPath.cpp
  ${TAUDIR}/src/Profile/OpenMPLayer.cpp
  ${TAUDIR}/src/Profile/TauKojakOpariDisable.cpp
  ${TAUDIR}/src/Profile/Comp_intel9.cpp
 )
# /home/jnkim/src/tau-2.18.1p1/src/Profile/Tracer.cpp

# /home/jnkim/src/tau-2.18.1p1/src/Comp_pgi.cpp
# /home/jnkim/src/tau-2.18.1p1/src/Comp_xl.cpp
# /home/jnkim/src/tau-2.18.1p1/src/Comp_gnu.cpp
ADD_LIBRARY(tau ${tausrcs})

# TauMpi.c
# ClockSync.c
# ProfileMerge.c
# TauFMpi.c
# TauMpiExtensions.c
