
# NiO benchmark tests

# Input variables
# QMC_DATA - NiO subdirectory should contain the *.h5 files

FUNCTION(ADD_NIO_TEST TEST_NAME TEST_DIR INPUT_FILE H5_FILE)
  MESSAGE("Adding test ${TEST_NAME}")
  SET (WDIR "${qmcpack_BINARY_DIR}/tests/performance/NiO/${TEST_DIR}")
  SET (SDIR "${CMAKE_SOURCE_DIR}/tests/performance/NiO/sample/${TEST_DIR}")

  SET(QMC_APP "${qmcpack_BINARY_DIR}/bin/qmcpack")
  IF (ENABLE_TIMERS)
    SET(PERF_ARGS "--enable-timers=fine")
  ENDIF()

  SET (PPDIR "${CMAKE_SOURCE_DIR}/tests/performance/NiO")

  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_directory ${SDIR} ${WDIR})

  FILE(GLOB PP "${PPDIR}/*.xml")
  FOREACH(F IN LISTS PP)
    FILE(RELATIVE_PATH R ${PPDIR} ${F})
    MAYBE_SYMLINK(${F} ${WDIR}/../${R})
  ENDFOREACH()

  MAYBE_SYMLINK("${QMC_DATA}/NiO/${H5_FILE}" "${WDIR}/../${H5_FILE}")

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
ELSEIF(NOT EXISTS ${QMC_DATA}/NiO)
  MESSAGE("NiO directory under QMC_DATA does not exist. NiO performance tests not added.")
ELSE()

  # *.h5 files and md5 sums
  # 6476972b54b58c89d15c478ed4e10317  NiO-fcc-supertwist111-supershift000-S8.h5
  # b47f4be12f98f8a3d4b65d0ae048b837  NiO-fcc-supertwist111-supershift000-S16.h5
  # ee1f6c6699a24e30d7e6c122cde55ac1  NiO-fcc-supertwist111-supershift000-S32.h5
  # 40ecaf05177aa4bbba7d3bf757994548  NiO-fcc-supertwist111-supershift000-S64.h5
  # 0a530594a3c7eec4f0155b5b2ca92eb0  NiO-fcc-supertwist111-supershift000-S128.h5
  # cff0101debb11c8c215e9138658fbd21  NiO-fcc-supertwist111-supershift000-S256.h5


  # Rough guide to time/memory
  # Intel Haswell Xeon, 1 thread
  # S8   - 30 s     1.3 GB (spline) + 18 MB * nwalker
  # S16  - 60 s     2.1 GB (spline) + 71 MB * nwalker
  # S32  - 130 s    3.6 GB (spline) + 280 MB * nwalker
  # S64  - 400 s    7.2 GB (spline) + 1.1 GB * nwalker
  # S128 - 1900 s  14.3 GB (spline) + 4.5 GB * nwalker
  # S256 - 14800 s 27.1 GB (spline) + 18.1 GB * nwalker


  # Size is number of supercells. Multiply by 4 to get the number of atoms.
  SET(NIO_SIZES 8 16 32 64 128 256)

  IF (QMC_NIO_MAX_SIZE)
    FOREACH(SIZE IN LISTS NIO_SIZES)
      IF (SIZE GREATER QMC_NIO_MAX_SIZE)
        LIST(REMOVE_ITEM NIO_SIZES ${SIZE})
      ENDIF()
    ENDFOREACH()
    MESSAGE("NiO sizes to benchmark: ${NIO_SIZES}")
  ENDIF()

  IF (QMC_CUDA)
    SET(COMPUTE_TYPE gpu)
  ELSE()
    SET(COMPUTE_TYPE cpu)
  ENDIF()

  FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/process_perf.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  FOREACH(SIZE IN LISTS NIO_SIZES)
    MATH(EXPR ATOM_COUNT "${SIZE} * 4")
    MATH(EXPR ELECTRON_COUNT "${SIZE} * 48")    
    SET(PERF_TEST performance-NiO-${COMPUTE_TYPE}-a${ATOM_COUNT}-e${ELECTRON_COUNT}-1-16)
    SET(TEST_DIR dmc-a${ATOM_COUNT}-e${ELECTRON_COUNT}-${COMPUTE_TYPE})
    SET(INPUT_FILE NiO-fcc-S${SIZE}-dmc.xml)
    SET(H5_FILE NiO-fcc-supertwist111-supershift000-S${SIZE}.h5)
    SET(H5_FULL_PATH "${QMC_DATA}/NiO/${H5_FILE}")
    IF (EXISTS ${H5_FULL_PATH})
      ADD_NIO_TEST(${PERF_TEST} ${TEST_DIR} ${INPUT_FILE} ${H5_FILE})
    ELSE()
      MESSAGE("${PERF_TEST} not added because the corresponding h5 file not found: ${H5_FULL_PATH}")
    ENDIF()
  ENDFOREACH()

  IF (NOT QMC_CUDA)
    FOREACH(SIZE IN LISTS NIO_SIZES)
      MATH(EXPR ATOM_COUNT "${SIZE} * 4")
      MATH(EXPR ELECTRON_COUNT "${SIZE} * 48")    
      SET(PERF_TEST performance-NiO-${COMPUTE_TYPE}-a${ATOM_COUNT}-e${ELECTRON_COUNT}-J3-1-16)
      SET(TEST_DIR dmc-a${ATOM_COUNT}-e${ELECTRON_COUNT}-${COMPUTE_TYPE}-J3)
      SET(INPUT_FILE NiO-fcc-S${SIZE}-dmc.xml)
      SET(H5_FILE NiO-fcc-supertwist111-supershift000-S${SIZE}.h5)
      SET(H5_FULL_PATH "${QMC_DATA}/NiO/${H5_FILE}")
      IF (EXISTS ${H5_FULL_PATH})
        ADD_NIO_TEST(${PERF_TEST} ${TEST_DIR} ${INPUT_FILE} ${H5_FILE})
      ELSE()
        MESSAGE("${PERF_TEST} not added because the corresponding h5 file not found: ${H5_FULL_PATH}")
      ENDIF()
    ENDFOREACH()
  ENDIF()

  # Example with SIZE and COMPUTE_TYPE expanded
  #SET(PERF_TEST NiO-cpu-a32-e384-dmc)
  #SET(TEST_DIR dmc-a32-e384-cpu)
  #SET(INPUT_FILE NiO-fcc-S8-dmc.xml)
  #SET(H5_FILE NiO-fcc-supertwist111-supershift000-S8.h5)
  #ADD_NIO_TEST(${PERF_TEST} ${TEST_DIR} ${INPUT_FILE} ${H5_FILE})
ENDIF()
