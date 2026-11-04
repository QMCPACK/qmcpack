project(Sandbox)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${qmcpack_BINARY_DIR}/bin)
# add apps XYZ.cpp, e.g., qmc_particles.cpp
set(ESTEST diff_distancetables einspline_spo einspline_spo_nested determinant restart determinant_delayed_update)

add_library(sandbox_helper ParticleIOUtility.cpp)
target_link_libraries(sandbox_helper PUBLIC qmcparticle)

foreach(p ${ESTEST})
  add_executable(${p} ${p}.cpp)
  target_link_libraries(${p} spline2 qmcparticle sandbox_helper)
  add_unit_test(BASENAME sandbox_${p} THREADS 3 TEST_BINARY $<TARGET_FILE:${p}>)
endforeach(p ${ESTEST})

if(HAVE_MPI)
  add_unit_test(BASENAME sandbox_restart PROCS 2 THREADS 1 TEST_BINARY $<TARGET_FILE:restart>)
  add_test_check_file_existence(${ADDED_UNIT_TEST_NAME} restart2.config.h5 TRUE)
endif()

get_property(CURRENT_TESTS DIRECTORY PROPERTY TESTS)
add_labels_for_tests("${CURRENT_TESTS}")
