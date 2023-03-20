# NiO benchmark tests

# Input variables
# QMC_DATA - NiO subdirectory should contain the *.h5 files

function(
  ADD_NIO_TEST
  TEST_NAME
  PROCS
  THREADS
  TEST_DIR
  TEST_SOURCE_DIR
  INPUT_FILE
  H5_FILE
  ADJUST_INPUT)
  message(VERBOSE "Adding test ${TEST_NAME}")
  set(WDIR "${qmcpack_BINARY_DIR}/tests/performance/NiO/${TEST_DIR}")
  set(SDIR "${qmcpack_SOURCE_DIR}/tests/performance/NiO/sample/${TEST_SOURCE_DIR}")
  set(QMC_APP $<TARGET_FILE:qmcpack>)
  if(ENABLE_TIMERS)
    set(PERF_ARGS "--enable-timers=fine")
  endif()

  file(MAKE_DIRECTORY ${WDIR})
  file(GLOB SS "${SDIR}/*.xml.in")
  foreach(S IN LISTS SS)
    get_filename_component(FN ${S} NAME_WE)
    # update PP_DIR and H5_DIR
    configure_file(${S} "${WDIR}/${FN}.xml")
  endforeach()

  if(NOT QMC_SYMLINK_TEST_FILES)
   configure_file("${QMC_DATA}/NiO/${H5_FILE}" "${WDIR}/../${H5_FILE}" COPYONLY)
  endif()

  separate_arguments(ADJUST_INPUT)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} ${qmcpack_SOURCE_DIR}/tests/performance/adjust_qmcpack_input.py ${ADJUST_INPUT}
            ${TEST_DIR}/${INPUT_FILE} WORKING_DIRECTORY "${qmcpack_BINARY_DIR}/tests/performance/NiO")

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
elseif(NOT EXISTS ${QMC_DATA}/NiO)
  message("NiO directory under QMC_DATA does not exist. NiO performance tests not added.")
else()

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
  set(NIO_SIZES
      4
      8
      16
      32
      64
      96
      128
      192
      256
      512
      1024)
  set(DELAY_RANKS
      8
      16
      16
      32
      32
      32
      32
      32
      64
      64
      64)
  # walker count preset for NIO_SIZES problems
  # These numbers are not optimal but should be indicative of performance
  # Walker count scan for each problem size is still required to maximize throughput
  # for a given hardware if its memory capacity allows.
  set(WALKER_COUNTS
      1024
      512
      256
      128
      64
      32
      16
      16
      16
      16
      16)

  if(QMC_PERFORMANCE_NIO_MAX_ATOMS)
    foreach(SIZE IN LISTS NIO_SIZES)
      if(SIZE GREATER QMC_PERFORMANCE_NIO_MAX_ATOMS)
        list(REMOVE_ITEM NIO_SIZES ${SIZE})
      endif()
    endforeach()
    message("NiO sizes to benchmark: ${NIO_SIZES}")
  endif()

  if(QMC_NIO_MAX_SIZE)
    message(FATAL_ERROR "Variable QMC_NIO_MAX_SIZE has been renamed to QMC_PERFORMANCE_NIO_MAX_ATOMS")
  endif()

  # set up PP_DIR which will be used to update test input files.
  set(PP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
  if(QMC_SYMLINK_TEST_FILES)
    set(PP_DIR ${PP_SOURCE_DIR})
  else()
    # copy to build directory
    foreach(F Ni.opt.xml O.xml)
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${F} ${CMAKE_CURRENT_BINARY_DIR}/${F} COPYONLY)
    endforeach()
    set(PP_DIR ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  # set up H5_DIR which will be used to update test input files.
  if(QMC_SYMLINK_TEST_FILES)
    set(H5_DIR ${QMC_DATA}/NiO)
  else()
    set(H5_DIR ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/process_perf.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  list(LENGTH NIO_SIZES LENGTH_MAX)
  math(EXPR LENGTH_MAX "${LENGTH_MAX} - 1")
  foreach(INDEX RANGE ${LENGTH_MAX})
    list(GET NIO_SIZES ${INDEX} ATOM_COUNT)
    math(EXPR SIZE "${ATOM_COUNT} / 4")
    math(EXPR ELECTRON_COUNT "${ATOM_COUNT} * 12")
    set(TEST_SOURCE_DIR dmc-a${ATOM_COUNT}-e${ELECTRON_COUNT}-cpu)
    set(INPUT_FILE NiO-fcc-S${SIZE}-dmc.xml)
    set(H5_FILE NiO-fcc-supertwist111-supershift000-S${SIZE}.h5)
    set(H5_FULL_PATH "${QMC_DATA}/NiO/${H5_FILE}")
    if(EXISTS ${H5_FULL_PATH})
      list(GET WALKER_COUNTS ${INDEX} WALKER_COUNT)
      set(ADJUST_INPUT "-i")
      # CPU driver
      set(DRIVER_TYPE cpu_driver)
      set(PERF_TEST performance-NiO-a${ATOM_COUNT}-e${ELECTRON_COUNT}-SM1-${DRIVER_TYPE}-w${WALKER_COUNT})
      set(TEST_DIR dmc-a${ATOM_COUNT}-e${ELECTRON_COUNT}-SM1-${DRIVER_TYPE})
      add_nio_test(
        ${PERF_TEST}
        1
        16
        ${TEST_DIR}
        ${TEST_SOURCE_DIR}
        ${INPUT_FILE}
        ${H5_FILE}
        "${ADJUST_INPUT} -w ${WALKER_COUNT} -d 1")
      ## delayed update
      list(GET DELAY_RANKS ${INDEX} NDELAY)
      set(PERF_TEST performance-NiO-a${ATOM_COUNT}-e${ELECTRON_COUNT}-DU${NDELAY}-${DRIVER_TYPE}-w${WALKER_COUNT})
      set(TEST_DIR dmc-a${ATOM_COUNT}-e${ELECTRON_COUNT}-DU${NDELAY}-${DRIVER_TYPE})
      add_nio_test(
        ${PERF_TEST}
        1
        16
        ${TEST_DIR}
        ${TEST_SOURCE_DIR}
        ${INPUT_FILE}
        ${H5_FILE}
        "${ADJUST_INPUT} -w ${WALKER_COUNT} -d ${NDELAY}")
      ## J3
      set(PERF_TEST performance-NiO-a${ATOM_COUNT}-e${ELECTRON_COUNT}-DU${NDELAY}-J3-${DRIVER_TYPE}-w${WALKER_COUNT})
      set(TEST_DIR dmc-a${ATOM_COUNT}-e${ELECTRON_COUNT}-DU${NDELAY}-J3-${DRIVER_TYPE})
      add_nio_test(
        ${PERF_TEST}
        1
        16
        ${TEST_DIR}
        ${TEST_SOURCE_DIR}
        ${INPUT_FILE}
        ${H5_FILE}
        "${ADJUST_INPUT} -w ${WALKER_COUNT} -d ${NDELAY} -j ${CMAKE_CURRENT_SOURCE_DIR}/sample/${TEST_SOURCE_DIR}/J123.xml")

      # Batched driver
      set(DRIVER_TYPE batched_driver)
      ## delayed update
      set(PERF_TEST performance-NiO-a${ATOM_COUNT}-e${ELECTRON_COUNT}-DU${NDELAY}-${DRIVER_TYPE}-w${WALKER_COUNT})
      set(TEST_DIR dmc-a${ATOM_COUNT}-e${ELECTRON_COUNT}-DU${NDELAY}-${DRIVER_TYPE})
      add_nio_test(
        ${PERF_TEST}
        1
        4
        ${TEST_DIR}
        ${TEST_SOURCE_DIR}
        ${INPUT_FILE}
        ${H5_FILE}
        "${ADJUST_INPUT} -w ${WALKER_COUNT} -d ${NDELAY} -u --detbatched")
      # To help analyze multithread offload implementations, add a single threaded batched test with same total workload as 4 thread case
      if(ENABLE_OFFLOAD)
        add_nio_test(
          ${PERF_TEST}
          1
          1
          ${TEST_DIR}
          ${TEST_SOURCE_DIR}
          ${INPUT_FILE}
          ${H5_FILE}
          "${ADJUST_INPUT} -w ${WALKER_COUNT} -d ${NDELAY} -u --detbatched")
      endif(ENABLE_OFFLOAD)

    else()
      message(
        "NiO-a${ATOM_COUNT}-e${ELECTRON_COUNT} performance test not added because the corresponding h5 file not found: ${H5_FULL_PATH}"
      )
    endif()
  endforeach()

  # Example with SIZE and DRIVER_TYPE expanded
  #SET(PERF_TEST NiO-cpu-a32-e384)
  #SET(TEST_DIR dmc-a32-e384-cpu)
  #SET(INPUT_FILE NiO-fcc-S8-dmc.xml)
  #SET(H5_FILE NiO-fcc-supertwist111-supershift000-S8.h5)
  #ADD_NIO_TEST(${PERF_TEST} ${TEST_DIR} ${INPUT_FILE} ${H5_FILE})
endif()
