# Runs unit tests
FUNCTION( ADD_QE_TEST TESTNAME PROCS TEST_BINARY NPOOL WORKDIR TEST_INPUT)
    IF ( USE_MPI )
        ADD_TEST( NAME ${TESTNAME} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${TEST_BINARY} -npool ${NPOOL} -inp ${TEST_INPUT} )
    ELSE()
        ADD_TEST( NAME ${TESTNAME} COMMAND ${TEST_BINARY} -npool 1 ${TEST_INPUT} )
    ENDIF()
    SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1 WORKING_DIRECTORY ${WORKDIR} )
    SET_PROPERTY( TEST ${TESTNAME} APPEND PROPERTY LABELS "Quantum ESPRESSO" )
ENDFUNCTION()

FUNCTION( RUN_QE_TEST BASE_NAME SRC_DIR PROCS1 PROCS2 PROCS3 NPOOL1 NPOOL2 NPOOL3 TEST_INPUT_PREFIX TEST_NAME)
    SET( FULL_NAME ${BASE_NAME}-np-${PROCS1}-${PROCS2}-${PROCS3}-nk-${NPOOL1}-${NPOOL2}-${NPOOL3} )
    SET( ${TEST_NAME} ${FULL_NAME} PARENT_SCOPE)
    SET( MY_WORKDIR ${CMAKE_CURRENT_BINARY_DIR}/${FULL_NAME} )
    MESSAGE("Adding test ${FULL_NAME}")
    COPY_DIRECTORY( "${SRC_DIR}" "${MY_WORKDIR}" )
    ADD_QE_TEST(${FULL_NAME}-scf  ${PROCS1} ${QE_BIN}/pw.x         ${NPOOL1} ${MY_WORKDIR} ${TEST_INPUT_PREFIX}-scf.in )
    IF(PROCS2 EQUAL 0)
        ADD_QE_TEST(${FULL_NAME}-pw2x ${PROCS3} ${QE_BIN}/pw2qmcpack.x ${NPOOL3} ${MY_WORKDIR} ${TEST_INPUT_PREFIX}-pw2x.in )
        SET_TESTS_PROPERTIES(${FULL_NAME}-pw2x PROPERTIES DEPENDS ${FULL_NAME}-scf)
    ELSE(PROCS2 EQUAL 0)
        ADD_QE_TEST(${FULL_NAME}-nscf ${PROCS2} ${QE_BIN}/pw.x         ${NPOOL2} ${MY_WORKDIR} ${TEST_INPUT_PREFIX}-nscf.in )
        SET_TESTS_PROPERTIES(${FULL_NAME}-nscf PROPERTIES DEPENDS ${FULL_NAME}-scf)
        ADD_QE_TEST(${FULL_NAME}-pw2x ${PROCS3} ${QE_BIN}/pw2qmcpack.x ${NPOOL3} ${MY_WORKDIR} ${TEST_INPUT_PREFIX}-pw2x.in )
        SET_TESTS_PROPERTIES(${FULL_NAME}-pw2x PROPERTIES DEPENDS ${FULL_NAME}-nscf)
    ENDIF(PROCS2 EQUAL 0)
ENDFUNCTION()

FUNCTION( SOFTLINK_H5 SOURCE TARGET PREFIX FILENAME TEST_NAME)
    SET(${TEST_NAME} "LINK_${SOURCE}_TO_${TARGET}" PARENT_SCOPE)
    ADD_TEST( NAME LINK_${SOURCE}_TO_${TARGET} COMMAND ${qmcpack_SOURCE_DIR}/utils/clean_and_link_h5.sh ${SOURCE}/out/${PREFIX}.pwscf.h5 ${SOURCE}-${TARGET}/${FILENAME} )
    SET_TESTS_PROPERTIES(LINK_${SOURCE}_TO_${TARGET} PROPERTIES DEPENDS ${SOURCE}-pw2x)
    SET_PROPERTY( TEST LINK_${SOURCE}_TO_${TARGET} APPEND PROPERTY LABELS "Link h5" )
ENDFUNCTION()
