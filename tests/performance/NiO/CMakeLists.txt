
# NiO benchmark tests

# Input variables
# QMC_DATA_NIO_H5 - location of the *.h5 files

FUNCTION(ADD_NIO_TEST TEST_NAME TEST_DIR INPUT_FILE H5_FILE)
  SET (WDIR "${qmcpack_BINARY_DIR}/tests/performance/NiO/${TEST_DIR}")
  SET (SDIR "${CMAKE_SOURCE_DIR}/tests/performance/NiO/${TEST_DIR}")

  SET(QMC_APP "${qmcpack_BINARY_DIR}/bin/qmcpack")
  SET(PERF_ARGS "--enable-timers=fine")

  SET (PPDIR "${CMAKE_SOURCE_DIR}/tests/performance/NiO")

  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_directory ${SDIR} ${WDIR})

  FILE(GLOB PP "${PPDIR}/*.xml")
  FOREACH(F IN LISTS PP)
    FILE(RELATIVE_PATH R ${PPDIR} ${F})
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink "${F}" ${WDIR}/${R})
  ENDFOREACH()

  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink "${QMC_DATA_NIO_H5}/${H5_FILE}" "${WDIR}/${H5_FILE}")

  ADD_TEST(NAME ${TEST_NAME} COMMAND ${QMC_APP} ${PERF_ARGS} ${INPUT_FILE})
  SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES LABELS "performance")
  SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES WORKING_DIRECTORY "${WDIR}")
  SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
ENDFUNCTION()

IF (NOT QMC_DATA_NIO_H5)
  MESSAGE("QMC_DATA_NIO_h5 not set. NiO performance tests not added")
ELSE()

  # Need some control to trim this list for time/memory for particular platforms
  SET(NIO_SIZES 8 16 32 64 128 256)

  IF (QMC_CUDA)
    SET(COMPUTE_TYPE gpu)
  ELSE()
    SET(COMPUTE_TYPE cpu)
  ENDIF()

  FOREACH(SIZE IN LISTS NIO_SIZES)
    SET(PERF_TEST NiO-cpu-S${SIZE}-${COMPUTE_TYPE})
    SET(TEST_DIR dmc-S${SIZE}-${COMPUTE_TYPE})
    SET(INPUT_FILE NiO-fcc-S${SIZE}-dmc.xml)
    SET(H5_FILE NiO-fcc-supertwist111-supershift000-S${SIZE}.h5)
    ADD_NIO_TEST(${PERF_TEST} ${TEST_DIR} ${INPUT_FILE} ${H5_FILE})
  ENDFOREACH()

  # Example with SIZE and COMPUTE_TYPE expanded
  #SET(PERF_TEST NiO-cpu-S8-dmc)
  #SET(TEST_DIR dmc-S8-cpu)
  #SET(INPUT_FILE NiO-fcc-S8-dmc.xml)
  #SET(H5_FILE NiO-fcc-supertwist111-supershift000-S8.h5)
  #ADD_NIO_TEST(${PERF_TEST} ${TEST_DIR} ${INPUT_FILE} ${H5_FILE})

ENDIF()
