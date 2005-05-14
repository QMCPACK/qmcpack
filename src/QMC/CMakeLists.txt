PROJECT(qmc)

CONFIGURE_FILE(${qmc_SOURCE_DIR}/main.cpp.in
              ${qmc_SOURCE_DIR}/main.cpp)

#choose classes for MolecuApps
IF(QMCAPPS MATCHES "MolecuApps")
  SET(WAVEFUNCTIONS
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
    ../QMCHamiltonians/WOS/Device.cpp
    ../QMCHamiltonians/WOS/WOSPotential.cpp
  )
  SET(APPSRC MolecuApps.cpp)
ENDIF(QMCAPPS MATCHES "MolecuApps")

#choose classes for QDApps
IF(QMCAPPS MATCHES "QDApps")
  SET(WAVEFUNCTIONS
    ../QMCWaveFunctions/OrbitalBuilderBase.cpp
    ../Numerics/Spline3D/Grid1D.cpp
    ../Numerics/Spline3D/Grid3D.cpp
    ../Numerics/Spline3D/CubicSpline.cpp
    ../Numerics/Spline3D/TriCubicSpline.cpp
    ../Numerics/Spline3D/TriCubicSplineSet.cpp
    ../QMCWaveFunctions/TriCubicSplineBuilder.cpp 
    ../QMCWaveFunctions/QDwfBuilder.cpp 
    ../QMCWaveFunctions/QDwf.cpp 
  )
  SET(APPSRC QDApps.cpp)
ENDIF(QMCAPPS MATCHES "QDApps")

SET(QMCSRCS 
  ../QMCDrivers/QMCDriver.cpp
  ../QMCDrivers/VMC.cpp
  ../QMCDrivers/VMCParticleByParticle.cpp
  ../QMCDrivers/DMCParticleByParticle.cpp
  ../QMCDrivers/VMC_OPT.cpp
  ../QMCDrivers/MolecuDMC.cpp
  ../QMCDrivers/ReptationMC.cpp
  ../Estimators/ScalarEstimatorManager.cpp
  QMCApps.cpp
  ${APPSRC}
)

LINK_LIBRARIES(
${LAPACK_LIBRARY} 
${BLAS_LIBRARY} 
${FORTRAN_LIBS} 
)

IF(${BLAS_EXTRA_LIBRARY} MATCHES "guide")
  LINK_LIBRARIES(${BLAS_EXTRA_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
ELSE(${BLAS_EXTRA_LIBRARY} MATCHES "guide")
  LINK_LIBRARIES(${CMAKE_THREAD_LIBS_INIT})
ENDIF(${BLAS_EXTRA_LIBRARY} MATCHES "guide")

ADD_EXECUTABLE(qmcapp ${WAVEFUNCTIONS} ${QMCSRCS} main.cpp)
TARGET_LINK_LIBRARIES(qmcapp qmc qmcbase)
