# functions for RMG workflow test

if(QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)
  function(ADD_RMG_TEST)

  endfunction()
  function(RUN_RMG_TEST)

  endfunction()
else(QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)

  function(
    ADD_RMG_TEST
    TESTNAME
    NPROCS
    NTHREADS
    TEST_BINARY
    WORKDIR
    TEST_INPUT)
    #if(HAVE_MPI)
    #  add_test(NAME ${TESTNAME} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NPROCS} ${MPIEXEC_PREFLAGS}
    #                                    ${TEST_BINARY} ${TEST_INPUT})
    #else(HAVE_MPI)
    add_test(NAME ${TESTNAME} COMMAND ${TEST_BINARY} ${TEST_INPUT})
    #endif(HAVE_MPI)
    set_tests_properties(
      ${TESTNAME}
      PROPERTIES ENVIRONMENT
                 "OMP_NUM_THREADS=${NTHREADS};RMG_NUM_THREADS=${NTHREADS}"
                 PROCESSORS
                 ${NPROCS}
                 PROCESSOR_AFFINITY
                 TRUE
                 WORKING_DIRECTORY
                 ${WORKDIR})
    set_property(
      TEST ${TESTNAME}
      APPEND
      PROPERTY LABELS "converter;rmg")
  endfunction()

  function(ADD_RMG_CONVERT_TEST TESTNAME PREFIX WORKDIR TEST_INPUT)
    add_test(NAME ${TESTNAME} COMMAND $<TARGET_FILE:convert4qmc> -nojastrow -rmg ${TEST_INPUT} -prefix ${PREFIX})
    set_tests_properties(${TESTNAME} PROPERTIES WORKING_DIRECTORY ${WORKDIR})
    set_property(TEST ${TESTNAME} APPEND PROPERTY LABELS "converter;rmg")
  endfunction()

  function(RUN_RMG_TEST BASE_NAME SRC_DIR NPROCS NTHREADS TEST_NAME)
    set(FULL_NAME ${BASE_NAME}-np-${NPROCS})
    set(${TEST_NAME}
        ${FULL_NAME}
        PARENT_SCOPE)
    set(MY_WORKDIR ${CMAKE_CURRENT_BINARY_DIR}/${FULL_NAME})
    message(VERBOSE "Adding test ${FULL_NAME}")
    copy_directory("${SRC_DIR}" "${MY_WORKDIR}")
    message("workdir: ${MY_WORKDIR}")
    add_rmg_test(${FULL_NAME}-scf ${NPROCS} ${NTHREADS} ${RMG_CPU_EXE} ${MY_WORKDIR} input)
    softlink_h5_rmg_waves(${FULL_NAME} ${BASE_NAME})
    add_rmg_convert_test(${FULL_NAME}-rmg2qmc ${BASE_NAME} ${MY_WORKDIR} ${BASE_NAME}.h5)
    # rmg2qmc needs h5 soft linked in work directory
    set_tests_properties(${FULL_NAME}-rmg2qmc PROPERTIES DEPENDS LINK_${FULL_NAME}_h5_Waves)
  endfunction()

endif(QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)

function(SOFTLINK_H5_RMG_WAVES SOURCE PREFIX)
  # set(${TEST_NAME}
  #     "LINK_${SOURCE}_h5_Waves"
  #     PARENT_SCOPE)
  add_test(NAME LINK_${SOURCE}_h5_Waves COMMAND ${qmcpack_SOURCE_DIR}/tests/scripts/clean_and_link_h5.sh
                                                ${SOURCE}/Waves/wave.out.h5 ${SOURCE}/${PREFIX}.h5)
  # creating soft link needs the scf run to generate h5.
  set_tests_properties(LINK_${SOURCE}_h5_Waves PROPERTIES DEPENDS ${SOURCE}-scf)
  set_property(TEST LINK_${SOURCE}_h5_Waves APPEND PROPERTY LABELS "rmg")
endfunction()

function(SOFTLINK_RMG_INPUT SOURCE TARGET PREFIX TEST_NAME)
  add_test(NAME LINK_${SOURCE}_TO_${TARGET} COMMAND ${qmcpack_SOURCE_DIR}/tests/scripts/clean_and_link_h5.sh
                                                    ${SOURCE}/${PREFIX}.h5 ${SOURCE}-${TARGET}/${PREFIX}.h5)
  # creating soft link needs LINK_${SOURCE}_h5_Waves because the link source is a softlink created by LINK_${SOURCE}_h5_Waves
  set_tests_properties(LINK_${SOURCE}_TO_${TARGET} PROPERTIES DEPENDS LINK_${SOURCE}_h5_Waves)
  set_property(TEST LINK_${SOURCE}_TO_${TARGET} APPEND PROPERTY LABELS "rmg")
  add_test(
    NAME COPY_${SOURCE}_XML_TO_${TARGET}
    COMMAND
      bash -c "mkdir -p ${SOURCE}-${TARGET}; \
                cp ${SOURCE}/${PREFIX}.structure.xml ${SOURCE}-${TARGET}/${PREFIX}.structure.xml ; \
                cp ${SOURCE}/${PREFIX}.wfnoj.xml ${SOURCE}-${TARGET}/${PREFIX}.wfnoj.xml ; \
                cp ${SOURCE}/*.qmcpp.xml ${SOURCE}-${TARGET}/")
  # there is only one dependency output, thus LINK_${SOURCE}_TO_${TARGET} is added as a dependency.
  # xml files require ${SOURCE}-rmg2qmc
  set_tests_properties(COPY_${SOURCE}_XML_TO_${TARGET} PROPERTIES DEPENDS "${SOURCE}-rmg2qmc;LINK_${SOURCE}_TO_${TARGET}")
  set_property(TEST COPY_${SOURCE}_XML_TO_${TARGET} APPEND PROPERTY LABELS "rmg")
  # the final and only dependency out
  set(${TEST_NAME}
      "COPY_${SOURCE}_XML_TO_${TARGET}"
      PARENT_SCOPE)
endfunction()
