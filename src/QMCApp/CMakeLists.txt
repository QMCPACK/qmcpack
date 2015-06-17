PROJECT(qmcapp)

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
# qmcapp is the only one
####################################
SET(QMCSRCS qmcapp)

FOREACH(p ${QMCSRCS})

  ADD_EXECUTABLE(${p} ${p}.cpp)
  TARGET_LINK_LIBRARIES(${p} qmc qmcdriver qmcham qmcwfs qmcbase qmcutil adios_config)

  FOREACH(l ${QMC_UTIL_LIBS})
    TARGET_LINK_LIBRARIES(${p} ${l})
  ENDFOREACH(l ${QMC_UTIL_LIBS})

  TARGET_LINK_LIBRARIES(${p} ${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${FORTRAN_LIBRARIES})

  if(ENABLE_TAU_PROFILE)
    TARGET_LINK_LIBRARIES(${p} tau)
  endif(ENABLE_TAU_PROFILE)

  IF(MPI_LIBRARY)
    TARGET_LINK_LIBRARIES(${p} ${MPI_LIBRARY})
  ENDIF(MPI_LIBRARY)

ENDFOREACH(p ${QMCSRCS})

####################################
# Add test subdirectory
####################################
ADD_SUBDIRECTORY(test)

