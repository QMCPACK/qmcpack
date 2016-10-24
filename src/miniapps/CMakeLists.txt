PROJECT(miniapps)

# add apps XYZ.cpp, e.g., qmc_particles.cpp
#SET(ESTEST einspline_smp einspline_spo qmc_particles moveonsphere twobody ptclset)
SET(ESTEST distancetable_soa)

FOREACH(p ${ESTEST})
  ADD_EXECUTABLE( ${p}  ${p}.cpp)
  TARGET_LINK_LIBRARIES(${p} ${QMC_UTIL_LIBS}  qmcutil qmcbase ${MPI_LIBRARY})
ENDFOREACH(p ${ESTEST})


