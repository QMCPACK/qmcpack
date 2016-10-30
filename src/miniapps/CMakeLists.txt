PROJECT(miniapps)

# add apps XYZ.cpp, e.g., qmc_particles.cpp
#SET(ESTEST einspline_smp einspline_spo qmc_particles moveonsphere twobody ptclset)
SET(ESTEST diff_distancetables diff_j2 diff_j1 einspline_spo einspline_spo_nested determinant)

FOREACH(p ${ESTEST})
  ADD_EXECUTABLE( ${p}  ${p}.cpp)
  TARGET_LINK_LIBRARIES(${p} qmcbase qmcutil ${QMC_UTIL_LIBS} ${MPI_LIBRARY})
ENDFOREACH(p ${ESTEST})

ADD_LIBRARY(miniwfs AoSWaveFunction.cpp SoAWaveFunction.cpp)

ADD_EXECUTABLE(miniqmc miniqmc.cpp)
TARGET_LINK_LIBRARIES(miniqmc miniwfs qmcbase qmcutil ${QMC_UTIL_LIBS} ${MPI_LIBRARY})

ADD_EXECUTABLE(diff_wfs diff_wfs.cpp)
TARGET_LINK_LIBRARIES(diff_wfs miniwfs qmcbase qmcutil ${QMC_UTIL_LIBS} ${MPI_LIBRARY})


