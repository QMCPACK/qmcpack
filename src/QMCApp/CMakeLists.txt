PROJECT(newqmc)

####################################
# common wavefunctions and hamiltonians
####################################
SET(WFSRCS
  ../QMCWaveFunctions/OrbitalBuilderBase.cpp
  ../QMCWaveFunctions/JAABuilder.cpp
  ../QMCWaveFunctions/JABBuilder.cpp
  ../QMCWaveFunctions/NJAABuilder.cpp
  ../QMCWaveFunctions/NJABBuilder.cpp
  ../QMCWaveFunctions/LRTwoBodyJastrow.cpp
  ../QMCWaveFunctions/JAAPBCBuilder.cpp
  ../QMCWaveFunctions/BasisSetBase.cpp
  ../QMCWaveFunctions/ComboOrbital.cpp
  ../QMCWaveFunctions/OrbitalConstraintsBase.cpp
  ../QMCWaveFunctions/PadeConstraints.cpp
  ../QMCWaveFunctions/TwoBodyJastrowBuilder.cpp
  ../QMCWaveFunctions/RPAConstraints.cpp
  ../QMCWaveFunctions/WMConstraints.cpp
  ../QMCWaveFunctions/WaveFunctionFactory.cpp
  ../QMCWaveFunctions/SPOSetBase.cpp
  ../QMCWaveFunctions/BasisSetFactory.cpp
  ../QMCWaveFunctions/Fermion/DiracDeterminantBase.cpp
  ../QMCWaveFunctions/Fermion/SlaterDet.cpp
  ../QMCWaveFunctions/Fermion/SlaterDetBuilder.cpp
  ../QMCWaveFunctions/MolecularOrbitals/STOBuilder.cpp
  ../QMCWaveFunctions/MolecularOrbitals/GTOBuilder.cpp
  ../QMCWaveFunctions/MolecularOrbitals/NGOBuilder.cpp
  ../LongRange/LPQHIBasis.cpp
  ../QMCHamiltonians/CoulombPBCAATemp.cpp
  ../QMCHamiltonians/CoulombPBCABTemp.cpp
  ../QMCHamiltonians/LocalCorePolPotential.cpp
  ../QMCHamiltonians/HamiltonianFactory.cpp
  ../QMCHamiltonians/Ylm.cpp
  ../Estimators/LocalEnergyEstimator.cpp
  ../Estimators/ScalarEstimatorManager.cpp
)

IF(QMC_COMPLEX)
  SET(WFSRCS ${WFSRCS}
    ../QMCWaveFunctions/ElectronGasComplexOrbitalBuilder.cpp
    ../QMCWaveFunctions/PlaneWaveOrbitalBuilder.cpp
    ../QMCWaveFunctions/PlaneWave/PWOrbitalBuilder.cpp
    ../QMCWaveFunctions/PlaneWave/PWOrbitalSet.cpp
    ../QMCHamiltonians/LocalECPotential.cpp
    ../QMCHamiltonians/NonLocalECPComponent.cpp
    ../QMCHamiltonians/NonLocalECPotential.cpp
    ../QMCHamiltonians/ECPotentialBuilder.cpp
    ../QMCHamiltonians/ECPComponentBuilder.cpp
  )
ELSE(QMC_COMPLEX)
  SET(WFSRCS ${WFSRCS}
    ../QMCWaveFunctions/NumericalOrbitalSetBuilder.cpp
    ../QMCWaveFunctions/SplineSetBuilder.cpp
    ../QMCWaveFunctions/ElectronGasOrbitalBuilder.cpp
    ../QMCWaveFunctions/PlaneWaveOrbitalBuilder.cpp
    ../QMCWaveFunctions/ThreeBodyGeminal.cpp
    ../QMCWaveFunctions/ThreeBodyGeminalBuilder.cpp
    ../QMCWaveFunctions/ThreeBodyPadeBuilder.cpp
    ../QMCWaveFunctions/MolecularOrbitals/STOMolecularOrbitals.cpp
    ../QMCWaveFunctions/MolecularOrbitals/GTOMolecularOrbitals.cpp
    ../QMCWaveFunctions/MolecularOrbitals/RGFBuilderBase.cpp
    ../QMCWaveFunctions/MolecularOrbitals/NumericalRGFBuilder.cpp
    ../QMCWaveFunctions/MolecularOrbitals/STO2GridBuilder.cpp
    ../QMCWaveFunctions/MolecularOrbitals/GTO2GridBuilder.cpp
    ../QMCWaveFunctions/MolecularOrbitals/Any2GridBuilder.cpp
    ../QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.cpp
    ../QMCWaveFunctions/MolecularOrbitals/MolecularOrbitalBuilder.cpp
    ../QMCWaveFunctions/AtomicOrbitals/HFAtomicSTOSet.cpp 
    ../QMCWaveFunctions/AtomicOrbitals/HFAtomicSTOSetBuilder.cpp 
    ../QMCWaveFunctions/AtomicOrbitals/HeSTOClementiRottie.cpp
    ../QMCHamiltonians/NonLocalPPotential.cpp
    ../QMCHamiltonians/LocalECPotential.cpp
    ../QMCHamiltonians/NonLocalECPComponent.cpp
    ../QMCHamiltonians/NonLocalECPotential.cpp
    ../QMCHamiltonians/ECPotentialBuilder.cpp
    ../QMCHamiltonians/ECPComponentBuilder.cpp
  )
ENDIF(QMC_COMPLEX)

####################################
# QMC Drivers which implements various QMC algorithms
####################################
SET(DMCDRIVERS
  ../QMCHamiltonians/NonLocalTOperator.cpp
  ../QMCDrivers/WalkerControlBase.cpp
  ../QMCDrivers/DMC/DMCMoveAll.cpp
  ../QMCDrivers/DMC/DMCPbyP.cpp
  ../QMCDrivers/DMC/DMCFactory.cpp
  ../QMCDrivers/DMC/DMCUpdateBase.cpp
  ../QMCDrivers/DMC/DMCUpdateAll.cpp
  ../QMCDrivers/DMC/DMCUpdatePbyP.cpp
  ../QMCDrivers/DMC/DMCNonLocalUpdate.cpp
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
  ../QMCDrivers/ReptationMC.cpp
  ../QMCDrivers/QMCOptimize.cpp
  ../QMCDrivers/QMCCostFunction.cpp
  ../QMCDrivers/WaveFunctionTester.cpp
  ${DMCDRIVERS}
)

IF(NOT QMCPLUSPLUS_RELEASE)

  IF(NOT QMC_COMPLEX)
    SET(QMCDRIVERS ${QMCDRIVERS}
      ../QMCWaveFunctions/AGPDeterminant.cpp
      ../QMCWaveFunctions/AGPDeterminantBuilder.cpp
    )
  ENDIF(NOT QMC_COMPLEX)

  SET(QMCDRIVERS ${QMCDRIVERS}
    ../QMCDrivers/VMCMultiple.cpp
    ../QMCDrivers/VMCPbyPMultiple.cpp
    ../QMCDrivers/PolymerEstimator.cpp
    ../QMCDrivers/RQMCMultiple.cpp
    ../QMCDrivers/MultiChain.cpp
    ../QMCDrivers/SpaceWarp.cpp
    ../QMCDrivers/VMCMultipleWarp.cpp
    ../QMCDrivers/VMCPbyPMultiWarp.cpp
    ../QMCDrivers/RQMCMultiWarp.cpp
    ../Estimators/MultipleEnergyEstimator.cpp
  )

#include GSLOptimize.cpp if libgsl is enabled
  IF(HAVE_LIBGSL)
    SET(QMCDRIVERS ${QMCDRIVERS} ../QMCDrivers/GSLOptimize.cpp)
  ENDIF(HAVE_LIBGSL)
ENDIF(NOT QMCPLUSPLUS_RELEASE)

SET(QMCAPPDIR 
    ${WFSRCS} 
    ${QMCDRIVERS}
    ParticleSetPool.cpp
    WaveFunctionPool.cpp
    HamiltonianPool.cpp
    InitMolecularSystem.cpp
    QMCAppBase.cpp
    QMCDriverFactory.cpp
    QMCMain.cpp
		QMCInterface.cpp
)

####################################
# create libqmc 
####################################
ADD_LIBRARY(qmc ${QMCAPPDIR})

####################################
#sources for new qmcapp, general QMC code to handle multiple
#configurations, hamiltonians and wavefunctions.
####################################
SET(QMCMAIN 
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
#LINK_LIBRARIES(
#  ${LAPACK_LIBRARY} 
#  ${BLAS_LIBRARY} 
#  ${FORTRAN_LIBS} 
#)

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
ADD_EXECUTABLE(qmcapp ${QMCMAIN})

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  TARGET_LINK_LIBRARIES(qmcapp qmc qmcbase qmcutil)
ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  TARGET_LINK_LIBRARIES(qmcapp qmc qmcbase qmcutil 
  ${LAPACK_LIBRARY}
  ${BLAS_LIBRARY} ${FORTRAN_LIBS})
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

