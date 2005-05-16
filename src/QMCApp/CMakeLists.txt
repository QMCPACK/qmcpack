PROJECT(newqmc)

####################################
# common wavefunctions and hamiltonians
####################################
SET(WFSRCS
  ../QMCWaveFunctions/OrbitalBuilderBase.cpp
  ../QMCWaveFunctions/AtomicOrbitals/HFAtomicSTOSet.cpp 
  ../QMCWaveFunctions/AtomicOrbitals/HFAtomicSTOSetBuilder.cpp 
  ../QMCWaveFunctions/AtomicOrbitals/HeSTOClementiRottie.cpp
  ../QMCWaveFunctions/MolecularOrbitals/STOMolecularOrbitals.cpp
  ../QMCWaveFunctions/MolecularOrbitals/RGFBuilderBase.cpp
  ../QMCWaveFunctions/MolecularOrbitals/STO2GridBuilder.cpp
  ../QMCWaveFunctions/MolecularOrbitals/GTO2GridBuilder.cpp
  ../QMCWaveFunctions/MolecularOrbitals/Any2GridBuilder.cpp
  ../QMCWaveFunctions/MolecularOrbitals/NumericalRGFBuilder.cpp
  ../QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.cpp
  ../QMCWaveFunctions/MolecularOrbitals/MolecularOrbitalBuilder.cpp
  ../QMCWaveFunctions/MolecularOrbitals/NumericalMolecularOrbitals.cpp
  ../QMCHamiltonians/LocalPPotential.cpp
  ../QMCHamiltonians/NonLocalPPotential.cpp
  ../QMCHamiltonians/GeCorePolPotential.cpp
)

####################################
# QMC Drivers which implements various QMC algorithms
####################################
SET(QMCDRIVERS
  ../QMCDrivers/QMCDriver.cpp
  ../QMCDrivers/VMC.cpp
  ../QMCDrivers/VMCParticleByParticle.cpp
  ../QMCDrivers/VMCMultiple.cpp
  ../QMCDrivers/VMC_OPT.cpp
  ../QMCDrivers/MolecuDMC.cpp
  ../QMCDrivers/DMCParticleByParticle.cpp
  ../QMCDrivers/ReptationMC.cpp
  ../Estimators/ScalarEstimatorManager.cpp
)

####################################
#sources for new qmcapp, general QMC code to handle multiple
#configurations, hamiltonians and wavefunctions.
####################################
SET(MAINSRCS 
  ${WFSRCS}
  ${QMCDRIVERS}
  ParticleSetPool.cpp
  WaveFunctionPool.cpp
  HamiltonianPool.cpp
  QMCAppBase.cpp
  QMCMain.cpp
  qmcapp.cpp
)

####################################
# sources for mo2grid
####################################
SET(MOSRCS 
  ${WFSRCS}
  ${QMCDRIVERS}
  QMCAppBase.cpp
  MO2Grid3D.cpp 
  mo2grid3d_main.cpp
)

####################################
# link lapack/blas and fortran library 
####################################
LINK_LIBRARIES(
  ${LAPACK_LIBRARY} 
  ${BLAS_LIBRARY} 
  ${FORTRAN_LIBS} 
)

####################################
# adding thread library: e.g., -lguide -lpthread
####################################
IF(${BLAS_EXTRA_LIBRARY} MATCHES "guide")
  LINK_LIBRARIES(${BLAS_EXTRA_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
ELSE(${BLAS_EXTRA_LIBRARY} MATCHES "guide")
  LINK_LIBRARIES(${CMAKE_THREAD_LIBS_INIT})
ENDIF(${BLAS_EXTRA_LIBRARY} MATCHES "guide")

####################################
# build mo2grid3d executable
####################################
ADD_EXECUTABLE(mo2grid3d ${MOSRCS})
TARGET_LINK_LIBRARIES(mo2grid3d qmc qmcbase)

####################################
#general qmcapp for multiple configurations, hamiltonians and wavefunctions.
####################################
ADD_EXECUTABLE(oops ${MAINSRCS})
TARGET_LINK_LIBRARIES(oops qmc qmcbase)
