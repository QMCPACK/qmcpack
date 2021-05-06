
# C-molecule benchmark tests

# Input variables
# QMC_DATA - C-molecule subdirectory should contain the *.h5 files

IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)
  FUNCTION(ADD_C_MOLECULE_TEST TEST_NAME TEST_DIR TEST_SOURCE_DIR INPUT_FILE H5_FILE ADJUST_INPUT)
    MESSAGE_VERBOSE("Adding test ${TEST_NAME}")
    SET (WDIR "${qmcpack_BINARY_DIR}/tests/performance/C-molecule/${TEST_DIR}")
    SET (SDIR "${qmcpack_SOURCE_DIR}/tests/performance/C-molecule/sample/${TEST_SOURCE_DIR}")

    GET_TARGET_PROPERTY( QMC_APP_NAME qmcpack OUTPUT_NAME )
    SET ( QMC_APP "${qmcpack_BINARY_DIR}/bin/${QMC_APP_NAME}" )
    IF (ENABLE_TIMERS)
      SET(PERF_ARGS "--enable-timers=fine")
    ENDIF()

    SET (PPDIR "${qmcpack_SOURCE_DIR}/tests/performance/C-molecule")

    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_directory ${SDIR} ${WDIR})

    FILE(GLOB PP "${PPDIR}/*.xml")
    FOREACH(F IN LISTS PP)
      FILE(RELATIVE_PATH R ${PPDIR} ${F})
      MAYBE_SYMLINK(${F} ${WDIR}/../${R})
    ENDFOREACH()

    MAYBE_SYMLINK("${QMC_DATA}/C-molecule/${H5_FILE}" "${WDIR}/../${H5_FILE}")

    SEPARATE_ARGUMENTS(ADJUST_INPUT)
    EXECUTE_PROCESS(COMMAND ${qmcpack_SOURCE_DIR}/tests/performance/adjust_qmcpack_input.py ${ADJUST_INPUT} ${TEST_DIR}/${INPUT_FILE}
                   WORKING_DIRECTORY "${qmcpack_BINARY_DIR}/tests/performance/C-molecule")

    SET(PROCS 1)
    SET(THREADS 16)
    MATH(EXPR TOT_PROCS "${PROCS} * ${THREADS}")
    IF (HAVE_MPI)
      ADD_TEST(NAME ${TEST_NAME} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS} ${QMC_APP} ${PERF_ARGS} ${INPUT_FILE})
    ELSE()
      ADD_TEST(NAME ${TEST_NAME} COMMAND ${QMC_APP} ${PERF_ARGS} ${INPUT_FILE})
    ENDIF()

    SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES LABELS "performance")
    SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES WORKING_DIRECTORY "${WDIR}")
    SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=${THREADS})
    SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES PROCESSORS ${TOT_PROCS} PROCESSOR_AFFINITY TRUE)

    IF (ENABLE_TIMERS)
      ADD_TEST(NAME "${TEST_NAME}-time" COMMAND ../process_perf.py ${INPUT_FILE})
      SET_TESTS_PROPERTIES("${TEST_NAME}-time" PROPERTIES LABELS "performance")
      SET_TESTS_PROPERTIES("${TEST_NAME}-time" PROPERTIES WORKING_DIRECTORY "${WDIR}")
      SET_TESTS_PROPERTIES("${TEST_NAME}-time" PROPERTIES DEPENDS ${TEST_NAME})
    ENDIF()
  ENDFUNCTION()

  IF (NOT QMC_DATA)
    MESSAGE_VERBOSE("QMC_DATA not set. Performance tests not added.")
  ELSEIF(NOT EXISTS ${QMC_DATA}/C-molecule)
    MESSAGE("C-molecule directory under QMC_DATA does not exist. C-molecule performance tests not added.")
  ELSE()

    SET(C_SIZES 12 18 24 30 60)

    IF (QMC_C_MAX_SIZE)
      FOREACH(SIZE IN LISTS C_SIZES)
        IF (SIZE GREATER QMC_C_MAX_SIZE)
          LIST(REMOVE_ITEM C_SIZES ${SIZE})
        ENDIF()
      ENDFOREACH()
      MESSAGE_VERBOSE("C-molecule sizes to benchmark: ${C_SIZES}")
    ENDIF()

    SET(ADJUST_INPUT "-i")

    FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../NiO/process_perf.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    LIST(LENGTH C_SIZES LENGTH_MAX)
    MATH(EXPR LENGTH_MAX "${LENGTH_MAX} - 1")
    FOREACH(INDEX RANGE ${LENGTH_MAX})
     FOREACH(POTENTIAL_TYPE "ae" "pp")
      LIST(GET C_SIZES ${INDEX} SIZE)
      MATH(EXPR ATOM_COUNT "${SIZE}")
      IF (POTENTIAL_TYPE MATCHES "ae")
       MATH(EXPR ELECTRON_COUNT "${SIZE} * 6")
      ELSE()
       MATH(EXPR ELECTRON_COUNT "${SIZE} * 4")
      ENDIF()   
      SET(TEST_SOURCE_DIR dmc-C${ATOM_COUNT}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE})
      SET(INPUT_FILE C${SIZE}-dmc.xml)
      SET(H5_FILE C${SIZE}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE}.h5)
      SET(H5_FULL_PATH "${QMC_DATA}/C-molecule/${H5_FILE}")
      IF (EXISTS ${H5_FULL_PATH})
        SET(PERF_TEST performance-C-molecule-C${ATOM_COUNT}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE}-1-16)
        SET(TEST_DIR dmc-C${ATOM_COUNT}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE})
        ADD_C_MOLECULE_TEST(${PERF_TEST} ${TEST_DIR} ${TEST_SOURCE_DIR} ${INPUT_FILE} ${H5_FILE} "${ADJUST_INPUT}")
        # J3
        SET(PERF_TEST performance-C-molecule-C${ATOM_COUNT}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE}-J3-1-16)
        SET(TEST_DIR dmc-C${ATOM_COUNT}-e${ELECTRON_COUNT}-${POTENTIAL_TYPE}-J3)
        ADD_C_MOLECULE_TEST(${PERF_TEST} ${TEST_DIR} ${TEST_SOURCE_DIR} ${INPUT_FILE} ${H5_FILE} "${ADJUST_INPUT} -j ${TEST_DIR}/J123.xml")
      ELSEIF(NOT (SIZE MATCHES "60" AND POTENTIAL_TYPE MATCHES "ae"))   
          MESSAGE_VERBOSE("performance-C-molecule-C${ATOM_COUNT}-e${ELECTRON_COUNT} not added because the corresponding h5 file not found: ${H5_FULL_PATH}")
      ENDIF()
     ENDFOREACH()
    ENDFOREACH()

  ENDIF()
 ELSE()
    MESSAGE_VERBOSE("Skipping C-molecule performance tests because the basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
  MESSAGE_VERBOSE("Skipping C-molecule performance tests because the basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

