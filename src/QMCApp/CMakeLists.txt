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
  ../QMCWaveFunctions/MolecularOrbitals/GTOMolecularOrbitals.cpp
  ../QMCWaveFunctions/MolecularOrbitals/RGFBuilderBase.cpp
  ../QMCWaveFunctions/MolecularOrbitals/STO2GridBuilder.cpp
  ../QMCWaveFunctions/MolecularOrbitals/GTO2GridBuilder.cpp
  ../QMCWaveFunctions/MolecularOrbitals/Any2GridBuilder.cpp
  ../QMCWaveFunctions/MolecularOrbitals/NumericalRGFBuilder.cpp
  ../QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.cpp
  ../QMCWaveFunctions/MolecularOrbitals/MolecularOrbitalBuilder.cpp
  ../QMCWaveFunctions/NumericalOrbitalSetBuilder.cpp
  ../QMCWaveFunctions/ElectronGasOrbitalBuilder.cpp
  ../QMCHamiltonians/LocalPPotential.cpp
  ../QMCHamiltonians/NonLocalPPotential.cpp
  ../QMCHamiltonians/LocalCorePolPotential.cpp
  ../LongRange/LPQHIBasis.cpp
)

####################################
# QMC Drivers which implements various QMC algorithms
####################################
SET(QMCDRIVERS
  ../QMCDrivers/SimpleFixedNodeBranch.cpp
  ../QMCDrivers/QMCDriver.cpp
  ../QMCDrivers/DummyQMC.cpp
  ../QMCDrivers/VMC.cpp
  ../QMCDrivers/VMCParticleByParticle.cpp
  ../QMCDrivers/MolecuDMC.cpp
  ../QMCDrivers/ReptationMC.cpp
  ../QMCDrivers/PolymerEstimator.cpp
  ../QMCDrivers/VMCMultiple.cpp
  ../QMCDrivers/VMCPbyPMultiple.cpp
  ../QMCDrivers/RQMCMultiple.cpp
  ../QMCDrivers/MultiChain.cpp
  ../QMCDrivers/DMCParticleByParticle.cpp
  ../QMCDrivers/WalkerControlBase.cpp
  ../Estimators/ScalarEstimatorManager.cpp
  ../Estimators/MultipleEnergyEstimator.cpp
  ../QMCDrivers/QMCOptimize.cpp
  ../QMCDrivers/QMCCostFunction.cpp
)

#include GSLOptimize.cpp if libgsl is enabled
#IF(HAVE_LIBGSL)
#  SET(QMCDRIVERS ${QMCDRIVERS} ../QMCDrivers/GSLOptimize.cpp)
#ENDIF(HAVE_LIBGSL)

IF(HAVE_MPI) 
SET(QMCDRIVERS ${QMCDRIVERS} 
    ../QMCDrivers/GlobalWalkerControl.cpp
    )
ENDIF(HAVE_MPI) 


####################################
# Test code to debug DMCParticleByParticle
#  ../QMCDrivers/DMCParticleByParticleDebug.cpp
####################################

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
  InitMolecularSystem.cpp
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
  ParticleSetPool.cpp
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
IF(INTEL_MKL)
  LINK_LIBRARIES(${INTEL_GUIDE_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
ELSE(INTEL_MKL)
  LINK_LIBRARIES(${CMAKE_THREAD_LIBS_INIT})
ENDIF(INTEL_MKL)

####################################
# build mo2grid3d executable
####################################
IF(NOT HAVE_MPI)
ADD_EXECUTABLE(mo2grid3d ${MOSRCS})
TARGET_LINK_LIBRARIES(mo2grid3d qmc qmcbase)
ENDIF(NOT HAVE_MPI)

####################################
#general qmcapp for multiple configurations, hamiltonians and wavefunctions.
####################################
ADD_EXECUTABLE(qmcapp ${MAINSRCS})
TARGET_LINK_LIBRARIES(qmcapp qmc qmcbase)

