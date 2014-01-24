PROJECT(qmctools)

SET(QMCTOOLS
  ../QMCWaveFunctions/OrbitalBuilderBase.cpp
  RGFBuilderBase.cpp
  GTO2GridBuilder.cpp
)

SET(MOSRCS
  ../QMCWaveFunctions/OrbitalBuilderBase.cpp
  ../QMCApp/InitMolecularSystem.cpp
  ../QMCApp/ParticleSetPool.cpp
  RGFBuilderBase.cpp
  GTO2GridBuilder.cpp
  QMCGaussianParserBase.cpp
  GaussianFCHKParser.cpp
  GamesXmlParser.cpp
  GamesAsciiParser.cpp
  CasinoParser.cpp
  BParser.cpp
  BMakeFunc.cpp
)

# create libmocommon 
ADD_LIBRARY(mocommon ${MOSRCS})

set(QTOOLS convert4qmc MSDgenerator extract-eshdf-kvectors trace_density)
ADD_EXECUTABLE(getSupercell getSupercell.cpp)

FOREACH(p ${QTOOLS})

  ADD_EXECUTABLE(${p} ${p}.cpp)
  TARGET_LINK_LIBRARIES(${p} mocommon qmcbase qmcutil)

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

ENDFOREACH(p ${QTOOLS})

