# functions for QE workflow test

if(QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)
  function(ADD_QE_TEST)

  endfunction()
  function(RUN_QE_TEST)

  endfunction()
else(QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)

  function(
    ADD_QE_TEST
    TESTNAME
    PROCS
    TEST_BINARY
    NPOOL
    WORKDIR
    TEST_INPUT)
    if(HAVE_MPI)
      add_test(NAME ${TESTNAME} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
                                        ${TEST_BINARY} -npool ${NPOOL} -inp ${TEST_INPUT})
    else()
      add_test(NAME ${TESTNAME} COMMAND ${TEST_BINARY} -npool 1 -inp ${TEST_INPUT})
    endif()
    set_tests_properties(
      ${TESTNAME}
      PROPERTIES ENVIRONMENT
                 OMP_NUM_THREADS=1
                 PROCESSORS
                 ${PROCS}
                 PROCESSOR_AFFINITY
                 TRUE
                 WORKING_DIRECTORY
                 ${WORKDIR})
    set_property(
      TEST ${TESTNAME}
      APPEND
      PROPERTY LABELS "converter")
  endfunction()

  function(
    RUN_QE_TEST
    BASE_NAME
    SRC_DIR
    PROCS1
    PROCS2
    PROCS3
    NPOOL1
    NPOOL2
    NPOOL3
    TEST_INPUT_PREFIX
    TEST_NAME)
    set(FULL_NAME ${BASE_NAME}-np-${PROCS1}-${PROCS2}-${PROCS3}-nk-${NPOOL1}-${NPOOL2}-${NPOOL3})
    set(${TEST_NAME}
        ${FULL_NAME}
        PARENT_SCOPE)
    set(MY_WORKDIR ${CMAKE_CURRENT_BINARY_DIR}/${FULL_NAME})
    message(VERBOSE "Adding test ${FULL_NAME}")
    copy_directory("${SRC_DIR}" "${MY_WORKDIR}")
    add_qe_test(${FULL_NAME}-scf ${PROCS1} ${QE_PW_DIR}/pw.x ${NPOOL1} ${MY_WORKDIR} ${TEST_INPUT_PREFIX}-scf.in)
    if(PROCS2 EQUAL 0)
      add_qe_test(${FULL_NAME}-pw2x ${PROCS3} ${QE_PW2Q_DIR}/pw2qmcpack.x ${NPOOL3} ${MY_WORKDIR}
                  ${TEST_INPUT_PREFIX}-pw2x.in)
      set_tests_properties(${FULL_NAME}-pw2x PROPERTIES DEPENDS ${FULL_NAME}-scf)
    else(PROCS2 EQUAL 0)
      add_qe_test(${FULL_NAME}-nscf ${PROCS2} ${QE_PW_DIR}/pw.x ${NPOOL2} ${MY_WORKDIR} ${TEST_INPUT_PREFIX}-nscf.in)
      set_tests_properties(${FULL_NAME}-nscf PROPERTIES DEPENDS ${FULL_NAME}-scf)
      add_qe_test(${FULL_NAME}-pw2x ${PROCS3} ${QE_PW2Q_DIR}/pw2qmcpack.x ${NPOOL3} ${MY_WORKDIR}
                  ${TEST_INPUT_PREFIX}-pw2x.in)
      set_tests_properties(${FULL_NAME}-pw2x PROPERTIES DEPENDS ${FULL_NAME}-nscf)
    endif(PROCS2 EQUAL 0)
  endfunction()

endif(QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)

function(SOFTLINK_H5 SOURCE TARGET PREFIX FILENAME TEST_NAME)
  set(${TEST_NAME}
      "LINK_${SOURCE}_TO_${TARGET}"
      PARENT_SCOPE)
  add_test(NAME LINK_${SOURCE}_TO_${TARGET} COMMAND ${qmcpack_SOURCE_DIR}/tests/scripts/clean_and_link_h5.sh
                                                    ${SOURCE}/out/${PREFIX}.pwscf.h5 ${SOURCE}-${TARGET}/${FILENAME})
  set_tests_properties(LINK_${SOURCE}_TO_${TARGET} PROPERTIES DEPENDS ${SOURCE}-pw2x)
  set_property(TEST LINK_${SOURCE}_TO_${TARGET} APPEND PROPERTY LABELS "converter")
endfunction()
