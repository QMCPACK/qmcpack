project(Sandbox)

# add apps XYZ.cpp, e.g., qmc_particles.cpp
set(ESTEST diff_distancetables einspline_spo einspline_spo_nested determinant restart determinant_delayed_update)

foreach(p ${ESTEST})
  add_executable(${p} ${p}.cpp)
  target_link_libraries(${p} einspline qmcparticle)
  add_unit_test(sandbox_${p} 1 3 $<TARGET_FILE:${p}>)

endforeach(p ${ESTEST})
