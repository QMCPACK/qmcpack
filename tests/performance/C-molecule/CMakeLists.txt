# C-molecule benchmark tests

# Input variables
# QMC_DATA - C-molecule subdirectory should contain the *.h5 files

if(NOT QMC_CUDA)
  if(NOT QMC_COMPLEX)
    function(
      ADD_C_MOLECULE_TEST
      TEST_NAME
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
    elseif(NOT EXISTS ${QMC_DATA}/C-molecule)
      message("C-molecule directory under QMC_DATA does not exist. C-molecule performance tests not added.")
    else()

      set(C_SIZES 12 18 24 30 60)

      if(QMC_C_MAX_SIZE)
        foreach(SIZE IN LISTS C_SIZES)
          if(SIZE GREATER QMC_C_MAX_SIZE)
            list(REMOVE_ITEM C_SIZES ${SIZE})
          endif()
        endforeach()
        message(VERBOSE "C-molecule sizes to benchmark: ${C_SIZES}")
      endif()

      set(ADJUST_INPUT "-i")

      file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../NiO/process_perf.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

      list(LENGTH C_SIZES LENGTH_MAX)
      math(EXPR LENGTH_MAX "${LENGTH_MAX} - 1")
      foreach(INDEX RANGE ${LENGTH_MAX})
        foreach(POTENTIAL_TYPE "ae" "pp")
          list(GET C_SIZES ${INDEX} SIZE)
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
            set(PERF_TEST performance-C-molecule-C${ATOM_COUNT}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE}-1-16)
            set(TEST_DIR dmc-C${ATOM_COUNT}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE})
            add_c_molecule_test(${PERF_TEST} ${TEST_DIR} ${TEST_SOURCE_DIR} ${INPUT_FILE} ${H5_FILE} "${ADJUST_INPUT}")
            # J3
            set(PERF_TEST performance-C-molecule-C${ATOM_COUNT}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE}-J3-1-16)
            set(TEST_DIR dmc-C${ATOM_COUNT}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE}-J3)
            add_c_molecule_test(${PERF_TEST} ${TEST_DIR} ${TEST_SOURCE_DIR} ${INPUT_FILE} ${H5_FILE}
                                "${ADJUST_INPUT} -j ${TEST_DIR}/J123.xml")
          elseif(NOT (SIZE MATCHES "60" AND POTENTIAL_TYPE MATCHES "ae"))
            message(VERBOSE
              "performance-C-molecule-C${ATOM_COUNT}-e${ELECTRON_COUNT} not added because the corresponding h5 file not found: ${H5_FULL_PATH}"
            )
          endif()
        endforeach()
      endforeach()

    endif()
  else()
    message(VERBOSE
      "Skipping C-molecule performance tests because the basis sets are not supported by complex build (QMC_COMPLEX=1)")
  endif()
else()
  message(VERBOSE
    "Skipping C-molecule performance tests because the basis sets are not supported by CUDA build (QMC_CUDA=1)")
endif()
