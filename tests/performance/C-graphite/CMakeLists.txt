# C-graphite benchmark tests

# Input variables
# QMC_DATA - C-graphite subdirectory should contain the *.h5 files

function(ADD_C_GRAPHITE_TEST TEST_NAME TEST_DIR INPUT_FILE H5_FILE)
  message(VERBOSE "Adding test ${TEST_NAME}")
  set(WDIR "${qmcpack_BINARY_DIR}/tests/performance/C-graphite/${TEST_DIR}")
  set(SDIR "${qmcpack_SOURCE_DIR}/tests/performance/C-graphite/sample/${TEST_DIR}")
  set(QMC_APP $<TARGET_FILE:qmcpack>)
  if(ENABLE_TIMERS)
    set(PERF_ARGS "--enable-timers=fine")
  endif()

  set(PPDIR "${qmcpack_SOURCE_DIR}/tests/performance/C-graphite")

  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${SDIR} ${WDIR})

  file(GLOB PP "${PPDIR}/*.xml")
  foreach(F IN LISTS PP)
    file(RELATIVE_PATH R ${PPDIR} ${F})
    maybe_symlink(${F} ${WDIR}/../${R})
  endforeach()

  maybe_symlink("${QMC_DATA}/C-graphite/${H5_FILE}" "${WDIR}/../${H5_FILE}")

  set(PROCS 1)
  set(THREADS 16)
  math(EXPR TOT_PROCS "${PROCS} * ${THREADS}")
  if(HAVE_MPI)
    add_test(NAME ${TEST_NAME} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
                                       ${QMC_APP} ${PERF_ARGS} ${INPUT_FILE})
  else()
    add_test(NAME ${TEST_NAME} COMMAND ${QMC_APP} ${PERF_ARGS} ${INPUT_FILE})
  endif()

  set_tests_properties(${TEST_NAME} PROPERTIES LABELS "performance")
  set_tests_properties(${TEST_NAME} PROPERTIES WORKING_DIRECTORY "${WDIR}")
  set_tests_properties(${TEST_NAME} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=${THREADS})
  set_tests_properties(${TEST_NAME} PROPERTIES PROCESSORS ${TOT_PROCS} PROCESSOR_AFFINITY TRUE)

  if(ENABLE_TIMERS)
    add_test(NAME "${TEST_NAME}-time" COMMAND ${Python3_EXECUTABLE} ../process_perf.py ${INPUT_FILE})
    set_tests_properties("${TEST_NAME}-time" PROPERTIES LABELS "performance")
    set_tests_properties("${TEST_NAME}-time" PROPERTIES WORKING_DIRECTORY "${WDIR}")
    set_tests_properties("${TEST_NAME}-time" PROPERTIES DEPENDS ${TEST_NAME})
  endif()
endfunction()

if(NOT QMC_DATA)
  message(VERBOSE "QMC_DATA not set. Performance tests not added.")
elseif(NOT EXISTS ${QMC_DATA}/C-graphite)
  message("C-graphite directory under QMC_DATA does not exist. C-graphite performance tests not added.")
else()

  # *.h5 files and md5 sums
  # 30323351a2636ba4987e216be5c383a1  lda.pwscf.h5

  # Size is number of electrons
  set(C_SIZES 256)

  if(QMC_C_MAX_SIZE)
    foreach(SIZE IN LISTS C_SIZES)
      if(SIZE GREATER QMC_C_MAX_SIZE)
        list(REMOVE_ITEM C_SIZES ${SIZE})
      endif()
    endforeach()
    message("C-graphite sizes to benchmark: ${C_SIZES}")
  endif()

  if(QMC_CUDA)
    set(COMPUTE_TYPE gpu)
  else()
    set(COMPUTE_TYPE cpu)
  endif()

  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../NiO/process_perf.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  foreach(SIZE IN LISTS C_SIZES)
    math(EXPR ATOM_COUNT "${SIZE} / 4")
    set(PERF_TEST performance-C-graphite-${COMPUTE_TYPE}-a${ATOM_COUNT}-e${SIZE}-1-16)
    set(TEST_DIR dmc-a${ATOM_COUNT}-e${SIZE}-${COMPUTE_TYPE})
    set(INPUT_FILE C-graphite-S${SIZE}-dmc.xml)
    # H5 filename is fixed consistent with CORAL name, but should be templated if a range of carbon sizes is added
    set(H5_FILE lda.pwscf.h5)
    #    SET(H5_FILE C-graphite-supertwist111-supershift000-S${SIZE}.h5)
    set(H5_FULL_PATH "${QMC_DATA}/C-graphite/${H5_FILE}")
    if(EXISTS ${H5_FULL_PATH})
      add_c_graphite_test(${PERF_TEST} ${TEST_DIR} ${INPUT_FILE} ${H5_FILE})
    else()
      message("${PERF_TEST} not added because the corresponding h5 file not found: ${H5_FULL_PATH}")
    endif()
  endforeach()

endif()
