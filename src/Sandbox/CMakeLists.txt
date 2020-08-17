PROJECT(Sandbox)

# add apps XYZ.cpp, e.g., qmc_particles.cpp
SET(ESTEST einspline_spo einspline_spo_nested determinant restart determinant_delayed_update)

FOREACH(p ${ESTEST})
  ADD_EXECUTABLE( ${p}  ${p}.cpp)
  TARGET_LINK_LIBRARIES(${p} einspline qmcparticle)
ENDFOREACH(p ${ESTEST})
