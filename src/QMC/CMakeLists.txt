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
    ../QMCWaveFunctions/MolecularOrbitals/RadialGridFunctorBuilder.cpp
    ../QMCWaveFunctions/MolecularOrbitals/STOMolecularOrbitals.cpp
    ../QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.cpp
    ../QMCWaveFunctions/MolecularOrbitals/MolecularOrbitalBuilder.cpp
    ../QMCHamiltonians/CorePolPotential.cpp
    ../QMCHamiltonians/LocalPPotential.cpp
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
  QMCDriver.cpp
  VMC.cpp
  VMCParticleByParticle.cpp
  VMC_OPT.cpp
  MolecuDMC.cpp
  WaveFunctionTester.cpp
  ../Estimators/ScalarEstimatorManager.cpp
  QMCApps.cpp
  ${APPSRC}
  main.cpp
)

LINK_LIBRARIES(${FORTRAN_LIBS} 
${LAPACK_LIBRARY} 
${BLAS_LIBRARY} 
${BLAS_EXTRA_LIBRARY}
${CMAKE_THREAD_LIBS_INIT}
)

ADD_EXECUTABLE(qmcapp ${WAVEFUNCTIONS} ${QMCSRCS})
TARGET_LINK_LIBRARIES(qmcapp qmcbase qmc)
