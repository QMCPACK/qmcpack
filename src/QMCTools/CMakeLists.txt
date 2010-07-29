PROJECT(qmctools)

SET(QMCTOOLS
  ../QMCWaveFunctions/OrbitalBuilderBase.cpp
  RGFBuilderBase.cpp
  GTO2GridBuilder.cpp
)

SET(CASINO
  ../QMCWaveFunctions/OrbitalBuilderBase.cpp
  RGFBuilderBase.cpp
  GTO2GridBuilder.cpp
  ../QMCApp/InitMolecularSystem.cpp
  QMCGaussianParserBase.cpp
  GaussianFCHKParser.cpp
  GamesXmlParser.cpp
  GamesAsciiParser.cpp
  CasinoParser.cpp
  BParser.cpp
  BMakeFunc.cpp
  convert4qmc.cpp
)

LINK_LIBRARIES(${FORTRAN_LIBS} 
${LAPACK_LIBRARY} 
${BLAS_LIBRARY} 
)

IF(INTEL_MKL)
  LINK_LIBRARIES(${INTEL_GUIDE_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
ELSE(INTEL_MKL)
  LINK_LIBRARIES(${CMAKE_THREAD_LIBS_INIT})
ENDIF(INTEL_MKL)

#ADD_EXECUTABLE(gto2grid gto2grid.cpp ${QMCTOOLS})
#TARGET_LINK_LIBRARIES(gto2grid qmc qmcbase)
#ADD_EXECUTABLE(gto2slater gto2slater.cpp ${QMCTOOLS})
#TARGET_LINK_LIBRARIES(gto2slater qmcbase)
#ADD_EXECUTABLE(boostrandom random.cpp)
#TARGET_LINK_LIBRARIES(boostrandom qmcbase)
ADD_EXECUTABLE(convert4qmc ${CASINO})
TARGET_LINK_LIBRARIES(convert4qmc qmc qmcbase qmcutil
${LAPACK_LIBRARY}
${BLAS_LIBRARY} ${FORTRAN_LIBS})

#ADD_EXECUTABLE(density2D densityplot.cpp)
#TARGET_LINK_LIBRARIES(density2D qmc qmcbase qmcutil)
#
#ADD_EXECUTABLE(gofr gofr.cpp)
#TARGET_LINK_LIBRARIES(gofr qmc qmcbase qmcutil)

#ADD_EXECUTABLE(spline spline.cpp)
#TARGET_LINK_LIBRARIES(spline qmc qmcbase)
#ADD_EXECUTABLE(spline3d spline3.cpp)
#TARGET_LINK_LIBRARIES(spline3d qmc qmcbase)
#ADD_EXECUTABLE(spline3d_io spline3_io.cpp)
#TARGET_LINK_LIBRARIES(spline3d_io qmc qmcbase)
#ADD_EXECUTABLE(molinit ../QMCWaveFunctions/OrbitalBuilderBase.cpp Assign.cpp)
#TARGET_LINK_LIBRARIES(molinit qmc qmcbase)
#ADD_EXECUTABLE(h5merge HDFWalkerMerger h5merge.cpp)
#TARGET_LINK_LIBRARIES(h5merge qmc qmcbase qmcutil)
