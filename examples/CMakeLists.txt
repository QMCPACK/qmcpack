# Ensure that examples run

# This function was copied and modified from src/QMCApp/test/CMakeLists.txt
# Changes: add 'example' label, and add SERIES parameter

# Add a test run and associated scalar checks
# BASE_NAME - name of test (number of MPI processes, number of threads, and value to check (if applicable)
#             will be appended to get the full test name)
# BASE_DIR - source location of test input files
# PREFIX - prefix for output files
# INPUT_FILE - XML input file to QMCPACK
# PROCS - number of MPI processes
# THREADS - number of OpenMP threads
# SCALAR_VALUES - list of output values to check with check_scalars.py
#                 The list entries alternate between the value name and the value (usually a string with the
#                 both the average and error).
# SERIES - series index to compute
# SHOULD_SUCCEED - whether the test is expected to pass or fail.  Expected failing tests will not have
#                  the scalar tests added.

FUNCTION(QMC_RUN_AND_CHECK BASE_NAME BASE_DIR PREFIX INPUT_FILE PROCS THREADS SCALAR_VALUES SERIES SHOULD_SUCCEED)
    # Map from name of check to appropriate flag for check_scalars.py
    LIST(APPEND SCALAR_CHECK_TYPE "kinetic" "totenergy" "eeenergy" "samples" "potential" "localpotential" "flux")
    LIST(APPEND CHECK_SCALAR_FLAG "--ke"    "--le"      "--ee"     "--ts"    "--lp"      "--lpp"          "--fl")

    SET( TEST_ADDED FALSE )
    SET( FULL_NAME "${BASE_NAME}-${PROCS}-${THREADS}" )
    MESSAGE("Adding test ${FULL_NAME}")
    RUN_QMC_APP(${FULL_NAME} ${BASE_DIR} ${PROCS} ${THREADS} TEST_ADDED ${INPUT_FILE})
    SET_PROPERTY(TEST ${FULL_NAME} APPEND PROPERTY LABELS "example")


    IF ( TEST_ADDED AND NOT SHOULD_SUCCEED)
        SET_PROPERTY(TEST ${FULL_NAME} APPEND PROPERTY WILL_FAIL TRUE)
        #MESSAGE("Test ${FULL_NAME} should fail")
    ENDIF()

    IF ( TEST_ADDED AND SHOULD_SUCCEED)
        FOREACH(SCALAR_CHECK IN LISTS SCALAR_CHECK_TYPE)
            LIST(FIND ${SCALAR_VALUES} ${SCALAR_CHECK} IDX1)
            IF (IDX1 GREATER -1)
                LIST(FIND SCALAR_CHECK_TYPE ${SCALAR_CHECK} IDX)
                LIST(GET CHECK_SCALAR_FLAG ${IDX} FLAG)

                MATH( EXPR IDX2 "${IDX1} + 1")
                LIST(GET ${SCALAR_VALUES} ${IDX2} VALUE)

                SET( TEST_NAME "${FULL_NAME}-${SCALAR_CHECK}" )
                #MESSAGE("Adding scalar check ${TEST_NAME}")
                SET(CHECK_CMD ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 --series ${SERIES} -p ${PREFIX} -e 2 ${FLAG} ${VALUE})
                #MESSAGE("check command = ${CHECK_CMD}")
                ADD_TEST( NAME ${TEST_NAME}
                    COMMAND ${CHECK_CMD}
                    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${FULL_NAME}"
                )
                SET_PROPERTY( TEST ${TEST_NAME} APPEND PROPERTY DEPENDS ${FULL_NAME} )
                SET_PROPERTY( TEST ${TEST_NAME} APPEND PROPERTY LABELS "example" )
            ENDIF()
        ENDFOREACH(SCALAR_CHECK)
    ENDIF()
ENDFUNCTION()

SUBDIRS(molecules/H2O)
SUBDIRS(molecules/He)
