PROJECT(external)

SET(PC
  ../QMCWaveFunctions/OrbitalBuilderBase.cpp 
  PC.cpp
)

SET(ASSIGN 
  ../QMCWaveFunctions/OrbitalBuilderBase.cpp
  Assign.cpp
)

ADD_EXECUTABLE(pc ${PC})
TARGET_LINK_LIBRARIES(pc qmc qmcbase)
ADD_EXECUTABLE(assign ${ASSIGN})
TARGET_LINK_LIBRARIES(assign qmc qmcbase)
LINK_LIBRARIES(${FORTRAN_LIBS} 
${LAPACK_LIBRARY} 
${BLAS_LIBRARY} 
${BLAS_EXTRA_LIBRARY}
${CMAKE_THREAD_LIBS_INIT}
)


