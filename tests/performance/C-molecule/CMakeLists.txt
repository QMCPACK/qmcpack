# C-molecule benchmark tests

# Input variables
# QMC_DATA - C-molecule subdirectory should contain the *.h5 file

if(NOT QMC_COMPLEX)
  function(
    ADD_C_MOLECULE_TEST
    TEST_NAME
    PROCS
    THREADS
    TEST_DIR
    TEST_SOURCE_DIR
    INPUT_FILE
    H5_FILE
    ADJUST_INPUT)
    message(VERBOSE "Adding test ${TEST_NAME}")
    set(WDIR "${qmcpack_BINARY_DIR}/tests/performance/C-molecule/${TEST_DIR}")
    set(SDIR "${qmcpack_SOURCE_DIR}/tests/performance/C-molecule/sample/${TEST_SOURCE_DIR}")

    set(QMC_APP $<TARGET_FILE:qmcpack>)
    if(ENABLE_TIMERS)
      set(PERF_ARGS "--enable-timers=fine")
    endif()

    set(PPDIR "${qmcpack_SOURCE_DIR}/tests/performance/C-molecule")

    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${SDIR} ${WDIR})

    file(GLOB PP "${PPDIR}/*.xml")
    foreach(F IN LISTS PP)
      file(RELATIVE_PATH R ${PPDIR} ${F})
      maybe_symlink(${F} ${WDIR}/../${R})
    endforeach()

    maybe_symlink("${QMC_DATA}/C-molecule/${H5_FILE}" "${WDIR}/../${H5_FILE}")

    separate_arguments(ADJUST_INPUT)
    execute_process(
      COMMAND ${Python3_EXECUTABLE} ${qmcpack_SOURCE_DIR}/tests/performance/adjust_qmcpack_input.py ${ADJUST_INPUT}
              ${TEST_DIR}/${INPUT_FILE} WORKING_DIRECTORY "${qmcpack_BINARY_DIR}/tests/performance/C-molecule")

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
  endfunction()

  if(NOT QMC_DATA)
    message(VERBOSE "QMC_DATA not set. Performance tests not added.")
  elseif(NOT EXISTS ${QMC_DATA}/C-molecule)
    message("C-molecule directory under QMC_DATA does not exist. C-molecule performance tests not added.")
  else()

    set(C_SIZES 12 18 24 30 60)
    set(WALKER_COUNTS 256 128 64 32 16)

    if(QMC_PERFORMANCE_C_MOLECULE_MAX_ATOMS)
      foreach(SIZE IN LISTS C_SIZES)
        if(SIZE GREATER QMC_PERFORMANCE_C_MOLECULE_MAX_ATOMS)
          list(REMOVE_ITEM C_SIZES ${SIZE})
        endif()
      endforeach()
      message("C-molecule sizes to benchmark: ${C_SIZES}")
    endif()

    set(ADJUST_INPUT "-i")

    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../NiO/process_perf.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    list(LENGTH C_SIZES LENGTH_MAX)
    math(EXPR LENGTH_MAX "${LENGTH_MAX} - 1")
    foreach(INDEX RANGE ${LENGTH_MAX})
      foreach(POTENTIAL_TYPE "ae" "pp")
        list(GET C_SIZES ${INDEX} SIZE)
        list(GET WALKER_COUNTS ${INDEX} WALKER_COUNT)
        math(EXPR ATOM_COUNT "${SIZE}")
        if(POTENTIAL_TYPE MATCHES "ae")
          math(EXPR ELECTRON_COUNT "${SIZE} * 6")
        else()
          math(EXPR ELECTRON_COUNT "${SIZE} * 4")
        endif()
        set(TEST_SOURCE_DIR dmc-C${ATOM_COUNT}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE})
        set(INPUT_FILE C${SIZE}-dmc.xml)
        set(H5_FILE C${SIZE}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE}.h5)
        set(H5_FULL_PATH "${QMC_DATA}/C-molecule/${H5_FILE}")
        if(EXISTS ${H5_FULL_PATH})
          set(DRIVER_TYPE cpu_driver)
          set(PERF_TEST
              performance-C-molecule-C${ATOM_COUNT}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE}-${DRIVER_TYPE}-w${WALKER_COUNT}
          )
          set(TEST_DIR dmc-C${ATOM_COUNT}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE}-${DRIVER_TYPE})
          add_c_molecule_test(
            ${PERF_TEST}
            1
            16
            ${TEST_DIR}
            ${TEST_SOURCE_DIR}
            ${INPUT_FILE}
            ${H5_FILE}
            "${ADJUST_INPUT} -w ${WALKER_COUNT}")
          # J3
          set(PERF_TEST
              performance-C-molecule-C${ATOM_COUNT}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE}-J3-${DRIVER_TYPE}-w${WALKER_COUNT}
          )
          set(TEST_DIR dmc-C${ATOM_COUNT}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE}-J3-${DRIVER_TYPE})
          add_c_molecule_test(
            ${PERF_TEST}
            1
            16
            ${TEST_DIR}
            ${TEST_SOURCE_DIR}
            ${INPUT_FILE}
            ${H5_FILE}
            "${ADJUST_INPUT} -w ${WALKER_COUNT} -j ${TEST_DIR}/J123.xml")
          # Batched-driver
          set(DRIVER_TYPE batched_driver)
          set(PERF_TEST
              performance-C-molecule-C${ATOM_COUNT}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE}-${DRIVER_TYPE}-w${WALKER_COUNT}
          )
          set(TEST_DIR dmc-C${ATOM_COUNT}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE}-${DRIVER_TYPE})
          add_c_molecule_test(
            ${PERF_TEST}
            1
            4
            ${TEST_DIR}
            ${TEST_SOURCE_DIR}
            ${INPUT_FILE}
            ${H5_FILE}
            "${ADJUST_INPUT} -w ${WALKER_COUNT} -u --detbatched")
          # J3
          set(PERF_TEST
              performance-C-molecule-C${ATOM_COUNT}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE}-J3-${DRIVER_TYPE}-w${WALKER_COUNT}
          )
          set(TEST_DIR dmc-C${ATOM_COUNT}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE}-J3-${DRIVER_TYPE})
          add_c_molecule_test(
            ${PERF_TEST}
            1
            4
            ${TEST_DIR}
            ${TEST_SOURCE_DIR}
            ${INPUT_FILE}
            ${H5_FILE}
            "${ADJUST_INPUT} -w ${WALKER_COUNT} -u --detbatched -j ${TEST_DIR}/J123.xml")

        elseif(NOT (SIZE MATCHES "60" AND POTENTIAL_TYPE MATCHES "ae"))
          message(
            VERBOSE
            "performance-C-molecule-C${ATOM_COUNT}-e${ELECTRON_COUNT} not added because the corresponding h5 file not found: ${H5_FULL_PATH}"
          )
        endif()
      endforeach()
    endforeach()

  endif()
else()
  message(
    VERBOSE
    "Skipping C-molecule performance tests because the basis sets are not supported by complex build (QMC_COMPLEX=1)")
endif()
