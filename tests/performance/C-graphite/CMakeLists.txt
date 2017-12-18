
# C-graphite benchmark tests

# Input variables
# QMC_DATA - C-graphite subdirectory should contain the *.h5 files

FUNCTION(ADD_C_GRAPHITE_TEST TEST_NAME TEST_DIR INPUT_FILE H5_FILE)
  MESSAGE("Adding test ${TEST_NAME}")
  SET (WDIR "${qmcpack_BINARY_DIR}/tests/performance/C-graphite/${TEST_DIR}")
  SET (SDIR "${CMAKE_SOURCE_DIR}/tests/performance/C-graphite/sample/${TEST_DIR}")

  SET(QMC_APP "${qmcpack_BINARY_DIR}/bin/qmcpack")
  IF (ENABLE_TIMERS)
    SET(PERF_ARGS "--enable-timers=fine")
  ENDIF()

  SET (PPDIR "${CMAKE_SOURCE_DIR}/tests/performance/C-graphite")

  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_directory ${SDIR} ${WDIR})

  FILE(GLOB PP "${PPDIR}/*.xml")
  FOREACH(F IN LISTS PP)
    FILE(RELATIVE_PATH R ${PPDIR} ${F})
    MAYBE_SYMLINK(${F} ${WDIR}/../${R})
  ENDFOREACH()

  MAYBE_SYMLINK("${QMC_DATA}/C-graphite/${H5_FILE}" "${WDIR}/../${H5_FILE}")

  SET(PROCS 1)
  SET(THREADS 16)
  MATH(EXPR TOT_PROCS "${PROCS} * ${THREADS}")
  IF (USE_MPI)
    ADD_TEST(NAME ${TEST_NAME} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${QMC_APP} ${PERF_ARGS} ${INPUT_FILE})
  ELSE()
    ADD_TEST(NAME ${TEST_NAME} COMMAND ${QMC_APP} ${PERF_ARGS} ${INPUT_FILE})
  ENDIF()

  SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES LABELS "performance")
  SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES WORKING_DIRECTORY "${WDIR}")
  SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=${THREADS})
  SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES PROCESSORS ${TOT_PROCS})

  IF (ENABLE_TIMERS)
    ADD_TEST(NAME "${TEST_NAME}-time" COMMAND python ../process_perf.py ${INPUT_FILE})
    SET_TESTS_PROPERTIES("${TEST_NAME}-time" PROPERTIES LABELS "performance")
    SET_TESTS_PROPERTIES("${TEST_NAME}-time" PROPERTIES WORKING_DIRECTORY "${WDIR}")
    SET_TESTS_PROPERTIES("${TEST_NAME}-time" PROPERTIES DEPENDS ${TEST_NAME})
  ENDIF()
ENDFUNCTION()

IF (NOT QMC_DATA)
  MESSAGE("QMC_DATA not set. Performance tests not added.")
ELSEIF(NOT EXISTS ${QMC_DATA}/C-graphite)
  MESSAGE("C-graphite directory under QMC_DATA does not exist. C-graphite performance tests not added.")
ELSE()

  # *.h5 files and md5 sums
  # 30323351a2636ba4987e216be5c383a1  lda.pwscf.h5

  # Size is number of electrons
  SET(C_SIZES 256)

  IF (QMC_C_MAX_SIZE)
    FOREACH(SIZE IN LISTS C_SIZES)
      IF (SIZE GREATER QMC_C_MAX_SIZE)
        LIST(REMOVE_ITEM C_SIZES ${SIZE})
      ENDIF()
    ENDFOREACH()
    MESSAGE("C-graphite sizes to benchmark: ${C_SIZES}")
  ENDIF()

  IF (QMC_CUDA)
    SET(COMPUTE_TYPE gpu)
  ELSE()
    SET(COMPUTE_TYPE cpu)
  ENDIF()

  FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../NiO/process_perf.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  FOREACH(SIZE IN LISTS C_SIZES)
    MATH(EXPR ATOM_COUNT "${SIZE} / 4")
    SET(PERF_TEST performance-C-graphite-${COMPUTE_TYPE}-a${ATOM_COUNT}-e${SIZE}-1-16)
    SET(TEST_DIR dmc-a${ATOM_COUNT}-e${SIZE}-${COMPUTE_TYPE})
    SET(INPUT_FILE C-graphite-S${SIZE}-dmc.xml)
# H5 filename is fixed consistent with CORAL name, but should be templated if a range of carbon sizes is added
    SET(H5_FILE lda.pwscf.h5)
#    SET(H5_FILE C-graphite-supertwist111-supershift000-S${SIZE}.h5)
    SET(H5_FULL_PATH "${QMC_DATA}/C-graphite/${H5_FILE}")
    IF (EXISTS ${H5_FULL_PATH})
      ADD_C_GRAPHITE_TEST(${PERF_TEST} ${TEST_DIR} ${INPUT_FILE} ${H5_FILE})
    ELSE()
      MESSAGE("${PERF_TEST} not added because the corresponding h5 file not found: ${H5_FULL_PATH}")
    ENDIF()
  ENDFOREACH()

ENDIF()
