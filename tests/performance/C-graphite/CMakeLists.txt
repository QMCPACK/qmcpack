# C-graphite benchmark tests

# Input variables
# QMC_DATA - C-graphite subdirectory should contain the *.h5 files

function(
  ADD_C_GRAPHITE_TEST
  TEST_NAME
  PROCS
  THREADS
  TEST_DIR
  TEST_SOURCE_DIR
  INPUT_FILE
  H5_FILE
  ADJUST_INPUT)
  message(VERBOSE "Adding test ${TEST_NAME}")
  set(WDIR "${qmcpack_BINARY_DIR}/tests/performance/C-graphite/${TEST_DIR}")
  set(SDIR "${qmcpack_SOURCE_DIR}/tests/performance/C-graphite/sample/${TEST_SOURCE_DIR}")
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

  set(TEST_NAME ${TEST_NAME}-${PROCS}-${THREADS})
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
  if(ENABLE_CUDA OR ENABLE_OFFLOAD)
    set_tests_properties(${TEST_NAME} PROPERTIES RESOURCE_LOCK exclusively_owned_gpus)
  endif()

  if(ENABLE_TIMERS)
    add_test(NAME "${TEST_NAME}-time" COMMAND ${Python3_EXECUTABLE} ../process_perf.py ${INPUT_FILE})
    set_tests_properties("${TEST_NAME}-time" PROPERTIES LABELS "performance")
    set_tests_properties("${TEST_NAME}-time" PROPERTIES WORKING_DIRECTORY "${WDIR}")
    set_tests_properties("${TEST_NAME}-time" PROPERTIES DEPENDS ${TEST_NAME})
  endif()

  separate_arguments(ADJUST_INPUT)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} ${qmcpack_SOURCE_DIR}/tests/performance/adjust_qmcpack_input.py ${ADJUST_INPUT}
            ${TEST_DIR}/${INPUT_FILE} WORKING_DIRECTORY "${qmcpack_BINARY_DIR}/tests/performance/C-graphite")

endfunction()

if(NOT QMC_DATA)
  message(VERBOSE "QMC_DATA not set. Performance tests not added.")
elseif(NOT EXISTS ${QMC_DATA}/C-graphite)
  message("C-graphite directory under QMC_DATA does not exist. C-graphite performance tests not added.")
else()

  # *.h5 files and md5 sums
  # 30323351a2636ba4987e216be5c383a1  lda.pwscf.h5

  # Size is number of atoms
  set(C_SIZES 4 16 36 64)
  set(WALKER_COUNTS 128 64 32 16)

  if(QMC_PERFORMANCE_C_GRAPHITE_MAX_ATOMS)
    foreach(SIZE IN LISTS C_SIZES)
      if(SIZE GREATER QMC_PERFORMANCE_C_GRAPHITE_MAX_ATOMS)
        list(REMOVE_ITEM C_SIZES ${SIZE})
      endif()
    endforeach()
    message("C-graphite sizes to benchmark: ${C_SIZES}")
  endif()

  set(COMPUTE_TYPE cpu)

  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../NiO/process_perf.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  list(LENGTH C_SIZES LENGTH_MAX)
  math(EXPR LENGTH_MAX "${LENGTH_MAX} - 1")
  foreach(INDEX RANGE ${LENGTH_MAX})
    list(GET C_SIZES ${INDEX} ATOM_COUNT)
    math(EXPR ELECTRON_COUNT "${ATOM_COUNT} * 4")
    math(EXPR SIZE "${ATOM_COUNT} / 4")
    set(TEST_SOURCE_DIR dmc-a${ATOM_COUNT}-e${ELECTRON_COUNT}-${COMPUTE_TYPE})
    set(INPUT_FILE C-graphite-S${SIZE}-dmc.xml)
    # H5 filename is fixed consistent with CORAL name, but should be templated if a range of carbon sizes is added
    set(H5_FILE lda.pwscf.h5)
    #    SET(H5_FILE C-graphite-supertwist111-supershift000-S${SIZE}.h5)
    set(H5_FULL_PATH "${QMC_DATA}/C-graphite/${H5_FILE}")
    if(EXISTS ${H5_FULL_PATH})
      list(GET WALKER_COUNTS ${INDEX} WALKER_COUNT)
      set(ADJUST_INPUT "-i")
      #Legacy driver
      set(DRIVER_TYPE cpu_driver)
      set(PERF_TEST performance-C-graphite-a${ATOM_COUNT}-e${ELECTRON_COUNT}-${DRIVER_TYPE}-w${WALKER_COUNT})
      set(TEST_DIR dmc-a${ATOM_COUNT}-e${ELECTRON_COUNT}-${DRIVER_TYPE})
      add_c_graphite_test(
        ${PERF_TEST}
        1
        16
        ${TEST_DIR}
        ${TEST_SOURCE_DIR}
        ${INPUT_FILE}
        ${H5_FILE}
        "${ADJUST_INPUT} -w ${WALKER_COUNT}")
      #Batched driver
      set(DRIVER_TYPE batched_driver)
      set(PERF_TEST performance-C-graphite-a${ATOM_COUNT}-e${ELECTRON_COUNT}-${DRIVER_TYPE}-w${WALKER_COUNT})
      set(TEST_DIR dmc-a${ATOM_COUNT}-e${ELECTRON_COUNT}-${DRIVER_TYPE})
      add_c_graphite_test(
        ${PERF_TEST}
        1
        4
        ${TEST_DIR}
        ${TEST_SOURCE_DIR}
        ${INPUT_FILE}
        ${H5_FILE}
        "${ADJUST_INPUT} -w ${WALKER_COUNT} -u --detbatched")

    else()
      message("${PERF_TEST} not added because the corresponding h5 file not found: ${H5_FULL_PATH}")
    endif()
  endforeach()

endif()
