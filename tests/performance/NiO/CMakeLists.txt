
# NiO benchmark tests

# Input variables
# QMC_DATA - NiO subdirectory should contain the *.h5 files

FUNCTION(ADD_NIO_TEST TEST_NAME PROCS THREADS TEST_DIR TEST_SOURCE_DIR INPUT_FILE H5_FILE ADJUST_INPUT)
  MESSAGE_VERBOSE("Adding test ${TEST_NAME}")
  SET (WDIR "${qmcpack_BINARY_DIR}/tests/performance/NiO/${TEST_DIR}")
  SET (SDIR "${qmcpack_SOURCE_DIR}/tests/performance/NiO/sample/${TEST_SOURCE_DIR}")

  SET(QMC_APP "${qmcpack_BINARY_DIR}/bin/qmcpack")
  IF (ENABLE_TIMERS)
    SET(PERF_ARGS "--enable-timers=fine")
  ENDIF()

  SET (PPDIR "${qmcpack_SOURCE_DIR}/tests/performance/NiO")

  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_directory ${SDIR} ${WDIR})

  FILE(GLOB PP "${PPDIR}/*.xml")
  FOREACH(F IN LISTS PP)
    FILE(RELATIVE_PATH R ${PPDIR} ${F})
    MAYBE_SYMLINK(${F} ${WDIR}/../${R})
  ENDFOREACH()

  MAYBE_SYMLINK("${QMC_DATA}/NiO/${H5_FILE}" "${WDIR}/../${H5_FILE}")

  SEPARATE_ARGUMENTS(ADJUST_INPUT)
  EXECUTE_PROCESS(COMMAND ${qmcpack_SOURCE_DIR}/tests/performance/adjust_qmcpack_input.py ${ADJUST_INPUT} ${TEST_DIR}/${INPUT_FILE}
                  WORKING_DIRECTORY "${qmcpack_BINARY_DIR}/tests/performance/NiO")

  SET(TEST_NAME ${TEST_NAME}-${PROCS}-${THREADS})
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
  IF (QMC_CUDA OR ENABLE_CUDA OR ENABLE_OFFLOAD)
    SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES RESOURCE_LOCK exclusively_owned_gpus)
  ENDIF()

  IF (ENABLE_TIMERS)
    ADD_TEST(NAME "${TEST_NAME}-time" COMMAND ../process_perf.py ${INPUT_FILE})
    SET_TESTS_PROPERTIES("${TEST_NAME}-time" PROPERTIES LABELS "performance")
    SET_TESTS_PROPERTIES("${TEST_NAME}-time" PROPERTIES WORKING_DIRECTORY "${WDIR}")
    SET_TESTS_PROPERTIES("${TEST_NAME}-time" PROPERTIES DEPENDS ${TEST_NAME})
  ENDIF()
ENDFUNCTION()

IF (NOT QMC_DATA)
  MESSAGE_VERBOSE("QMC_DATA not set. Performance tests not added.")
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
  SET(NIO_SIZES 1 2 4 8 16 24 32 48 64 128 256)
  SET(DELAY_RANKS 8 16 16 32 32 32 32 32 64 64 64)

  IF (QMC_NIO_MAX_SIZE)
    FOREACH(SIZE IN LISTS NIO_SIZES)
      IF (SIZE GREATER QMC_NIO_MAX_SIZE)
        LIST(REMOVE_ITEM NIO_SIZES ${SIZE})
      ENDIF()
    ENDFOREACH()
    MESSAGE("NiO sizes to benchmark: ${NIO_SIZES}")
  ENDIF()

  FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/process_perf.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  LIST(LENGTH NIO_SIZES LENGTH_MAX)
  MATH(EXPR LENGTH_MAX "${LENGTH_MAX} - 1")
  FOREACH(INDEX RANGE ${LENGTH_MAX})
    LIST(GET NIO_SIZES ${INDEX} SIZE)
    MATH(EXPR ATOM_COUNT "${SIZE} * 4")
    MATH(EXPR ELECTRON_COUNT "${SIZE} * 48")    
    SET(TEST_SOURCE_DIR dmc-a${ATOM_COUNT}-e${ELECTRON_COUNT}-cpu)
    SET(INPUT_FILE NiO-fcc-S${SIZE}-dmc.xml)
    SET(H5_FILE NiO-fcc-supertwist111-supershift000-S${SIZE}.h5)
    SET(H5_FULL_PATH "${QMC_DATA}/NiO/${H5_FILE}")
    IF (EXISTS ${H5_FULL_PATH})
      SET(ADJUST_INPUT "-i")
      IF (NOT QMC_CUDA)
        # CPU driver
        SET(DRIVER_TYPE cpu_driver)
        SET(PERF_TEST performance-NiO-a${ATOM_COUNT}-e${ELECTRON_COUNT}-${DRIVER_TYPE}-w16-SM1)
        SET(TEST_DIR dmc-a${ATOM_COUNT}-e${ELECTRON_COUNT}-${DRIVER_TYPE})
        ADD_NIO_TEST(${PERF_TEST} 1 16 ${TEST_DIR} ${TEST_SOURCE_DIR} ${INPUT_FILE} ${H5_FILE} "${ADJUST_INPUT} -w 16 -d 1")
        ## delayed update
        LIST(GET DELAY_RANKS ${INDEX} NDELAY)
        SET(PERF_TEST performance-NiO-a${ATOM_COUNT}-e${ELECTRON_COUNT}-${DRIVER_TYPE}-w16-DU${NDELAY})
        SET(TEST_DIR dmc-a${ATOM_COUNT}-e${ELECTRON_COUNT}-${DRIVER_TYPE}-DU${NDELAY})
        ADD_NIO_TEST(${PERF_TEST} 1 16 ${TEST_DIR} ${TEST_SOURCE_DIR} ${INPUT_FILE} ${H5_FILE} "${ADJUST_INPUT} -w 16 -d ${NDELAY}")
        ## J3
        SET(PERF_TEST performance-NiO-a${ATOM_COUNT}-e${ELECTRON_COUNT}-${DRIVER_TYPE}-w16-J3)
        SET(TEST_DIR dmc-a${ATOM_COUNT}-e${ELECTRON_COUNT}-${DRIVER_TYPE}-J3)
        ADD_NIO_TEST(${PERF_TEST} 1 16 ${TEST_DIR} ${TEST_SOURCE_DIR} ${INPUT_FILE} ${H5_FILE} "${ADJUST_INPUT} -w 16 -j ${TEST_DIR}/J123.xml")

        # Batched driver
        SET(DRIVER_TYPE batched_driver)
        ## delayed update
        SET(PERF_TEST performance-NiO-a${ATOM_COUNT}-e${ELECTRON_COUNT}-${DRIVER_TYPE}-w16-DU${NDELAY})
        SET(TEST_DIR dmc-a${ATOM_COUNT}-e${ELECTRON_COUNT}-${DRIVER_TYPE}-DU${NDELAY})
        ADD_NIO_TEST(${PERF_TEST} 1 4 ${TEST_DIR} ${TEST_SOURCE_DIR} ${INPUT_FILE} ${H5_FILE} "${ADJUST_INPUT} -w 16 -d ${NDELAY} -u --detbatched")
      ELSE (NOT QMC_CUDA)
        # legacy CUDA driver
        SET(DRIVER_TYPE legacy_cuda)
        SET(PERF_TEST performance-NiO-a${ATOM_COUNT}-e${ELECTRON_COUNT}-${DRIVER_TYPE}-w32)
        SET(TEST_DIR dmc-a${ATOM_COUNT}-e${ELECTRON_COUNT}-${DRIVER_TYPE})
        ADD_NIO_TEST(${PERF_TEST} 1 16 ${TEST_DIR} ${TEST_SOURCE_DIR} ${INPUT_FILE} ${H5_FILE} "${ADJUST_INPUT} -w 32")
      ENDIF (NOT QMC_CUDA)
    ELSE()
      MESSAGE("NiO-a${ATOM_COUNT}-e${ELECTRON_COUNT} performance test not added because the corresponding h5 file not found: ${H5_FULL_PATH}")
    ENDIF()
  ENDFOREACH()

  # Example with SIZE and DRIVER_TYPE expanded
  #SET(PERF_TEST NiO-cpu-a32-e384)
  #SET(TEST_DIR dmc-a32-e384-cpu)
  #SET(INPUT_FILE NiO-fcc-S8-dmc.xml)
  #SET(H5_FILE NiO-fcc-supertwist111-supershift000-S8.h5)
  #ADD_NIO_TEST(${PERF_TEST} ${TEST_DIR} ${INPUT_FILE} ${H5_FILE})
ENDIF()
