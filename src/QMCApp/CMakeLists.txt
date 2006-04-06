PROJECT(newqmc)

####################################
# common wavefunctions and hamiltonians
####################################
SET(WFSRCS
  ../QMCWaveFunctions/OrbitalBuilderBase.cpp
  ../QMCWaveFunctions/JAABuilder.cpp
  ../QMCWaveFunctions/NJAABuilder.cpp
  ../QMCWaveFunctions/JABBuilder.cpp
  ../QMCWaveFunctions/NJABBuilder.cpp
  ../QMCWaveFunctions/AtomicOrbitals/HFAtomicSTOSet.cpp 
  ../QMCWaveFunctions/AtomicOrbitals/HFAtomicSTOSetBuilder.cpp 
  ../QMCWaveFunctions/AtomicOrbitals/HeSTOClementiRottie.cpp
  ../QMCWaveFunctions/MolecularOrbitals/STOMolecularOrbitals.cpp
  ../QMCWaveFunctions/MolecularOrbitals/GTOMolecularOrbitals.cpp
  ../QMCWaveFunctions/MolecularOrbitals/RGFBuilderBase.cpp
  ../QMCWaveFunctions/MolecularOrbitals/NumericalRGFBuilder.cpp
  ../QMCWaveFunctions/MolecularOrbitals/STO2GridBuilder.cpp
  ../QMCWaveFunctions/MolecularOrbitals/GTO2GridBuilder.cpp
  ../QMCWaveFunctions/MolecularOrbitals/Any2GridBuilder.cpp
  ../QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.cpp
  ../QMCWaveFunctions/MolecularOrbitals/MolecularOrbitalBuilder.cpp
  ../QMCWaveFunctions/ElectronGasOrbitalBuilder.cpp
  ../QMCWaveFunctions/ThreeBodyGeminal.cpp
  ../QMCWaveFunctions/ThreeBodyGeminalBuilder.cpp
  ../QMCWaveFunctions/WaveFunctionFactory.cpp
  ../QMCHamiltonians/LocalPPotential.cpp
  ../QMCHamiltonians/NonLocalPPotential.cpp
  ../QMCHamiltonians/LocalCorePolPotential.cpp
  ../QMCHamiltonians/HamiltonianFactory.cpp
)

IF(NOT QMCPLUSPLUS_RELEASE)
SET(WFSRCS ${WFSRCS}
  ../QMCWaveFunctions/NumericalOrbitalSetBuilder.cpp
  ../LongRange/LPQHIBasis.cpp
  )

ENDIF(NOT QMCPLUSPLUS_RELEASE)

####################################
# QMC Drivers which implements various QMC algorithms
####################################
SET(DMCDRIVERS
  ../QMCDrivers/WalkerControlBase.cpp
  ../QMCDrivers/DMC/DMCMoveAll.cpp
  ../QMCDrivers/DMC/DMCPbyP.cpp
  ../QMCDrivers/DMC/DMCFactory.cpp
  ../QMCDrivers/DMC/DMCUpdateBase.cpp
  ../QMCDrivers/DMC/DMCUpdateAll.cpp
  ../QMCDrivers/DMC/DMCUpdatePbyP.cpp
  ../QMCDrivers/DMC/WalkerControlFactory.cpp
  ../QMCDrivers/DMC/WalkerReconfiguration.cpp
)

IF(ENABLE_OPENMP) 
  SET(DMCDRIVERS ${DMCDRIVERS} ../QMCDrivers/DMC/DMCPbyPOMP.cpp)
ENDIF(ENABLE_OPENMP) 

IF(HAVE_MPI) 
  SET(DMCDRIVERS ${DMCDRIVERS} 
  ../QMCDrivers/DMC/WalkerControlMPI.cpp
  ../QMCDrivers/DMC/WalkerReconfigurationMPI.cpp)
ENDIF(HAVE_MPI) 


SET(QMCDRIVERS 
  ../QMCDrivers/SimpleFixedNodeBranch.cpp
  ../QMCDrivers/QMCDriver.cpp
  ../QMCDrivers/DummyQMC.cpp
  ../QMCDrivers/VMC.cpp
  ../QMCDrivers/VMCParticleByParticle.cpp
  ../Estimators/ScalarEstimatorManager.cpp
  ../QMCDrivers/QMCOptimize.cpp
  ../QMCDrivers/QMCCostFunction.cpp
  ../QMCDrivers/WaveFunctionTester.cpp
  ${DMCDRIVERS}
)

IF(NOT QMCPLUSPLUS_RELEASE)
  SET(QMCDRIVERS ${QMCDRIVERS}
    ../QMCDrivers/VMCMultiple.cpp
    ../QMCDrivers/VMCPbyPMultiple.cpp
    ../QMCDrivers/PolymerEstimator.cpp
    ../QMCDrivers/RQMCMultiple.cpp
    ../QMCDrivers/MultiChain.cpp
    ../Estimators/MultipleEnergyEstimator.cpp
    ../QMCWaveFunctions/AGPDeterminant.cpp
    ../QMCWaveFunctions/AGPDeterminantBuilder.cpp
    ../QMCDrivers/SpaceWarp.cpp
    ../QMCDrivers/VMCMultipleWarp.cpp
    ../QMCDrivers/VMCPbyPMultiWarp.cpp
    ../QMCDrivers/RQMCMultiWarp.cpp
  )

#include GSLOptimize.cpp if libgsl is enabled
  IF(HAVE_LIBGSL)
    SET(QMCDRIVERS ${QMCDRIVERS} ../QMCDrivers/GSLOptimize.cpp)
  ENDIF(HAVE_LIBGSL)
ENDIF(NOT QMCPLUSPLUS_RELEASE)



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
  QMCDriverFactory.cpp
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
#IF(NOT HAVE_MPI)
#ADD_EXECUTABLE(mo2grid3d ${MOSRCS})
#TARGET_LINK_LIBRARIES(mo2grid3d qmc qmcbase)
#ENDIF(NOT HAVE_MPI)

####################################
#general qmcapp for multiple configurations, hamiltonians and wavefunctions.
####################################
ADD_EXECUTABLE(qmcapp ${MAINSRCS})
TARGET_LINK_LIBRARIES(qmcapp qmc qmcbase)

