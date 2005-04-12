PROJECT(qmctools)

SET(QMCTOOLS
  ../QMCWaveFunctions/OrbitalBuilderBase.cpp
  ../QMCWaveFunctions/MolecularOrbitals/RGFBuilderBase.cpp
  ../QMCWaveFunctions/MolecularOrbitals/GTO2GridBuilder.cpp
  gto2grid.cpp
)

SET(CASINO
  ../QMCWaveFunctions/OrbitalBuilderBase.cpp
  ../QMCWaveFunctions/MolecularOrbitals/RGFBuilderBase.cpp
  ../QMCWaveFunctions/MolecularOrbitals/GTO2GridBuilder.cpp
  QMCGaussianParserBase.cpp
  CasinoParser.cpp
  parse.cpp
)

LINK_LIBRARIES(${FORTRAN_LIBS} 
${LAPACK_LIBRARY} 
${BLAS_LIBRARY} 
${BLAS_EXTRA_LIBRARY}
${CMAKE_THREAD_LIBS_INIT}
)

#ADD_EXECUTABLE(gto2grid ${QMCTOOLS})
#TARGET_LINK_LIBRARIES(gto2grid qmc qmcbase)
#ADD_EXECUTABLE(casino ${CASINO})
#TARGET_LINK_LIBRARIES(casino qmc qmcbase)
#ADD_EXECUTABLE(molinit ../QMCWaveFunctions/OrbitalBuilderBase.cpp Assign.cpp)
#TARGET_LINK_LIBRARIES(molinit qmc qmcbase)
ADD_EXECUTABLE(spline spline.cpp)
TARGET_LINK_LIBRARIES(spline qmc qmcbase)
ADD_EXECUTABLE(spline3d spline3.cpp)
TARGET_LINK_LIBRARIES(spline3d qmc qmcbase)
ADD_EXECUTABLE(spline3d_io spline3_io.cpp)
TARGET_LINK_LIBRARIES(spline3d_io qmc qmcbase)
