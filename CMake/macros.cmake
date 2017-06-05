IF ( ${CMAKE_MAJOR_VERSION} EQUAL 3 )
    CMAKE_POLICY(SET CMP0026 OLD)
ENDIF()


# Function to copy a directory
FUNCTION( COPY_DIRECTORY SRC_DIR DST_DIR )
    EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_directory "${SRC_DIR}" "${DST_DIR}" )
ENDFUNCTION()

# Function to copy a directory using symlinks for the files. This saves storage
# space with large test files.
# SRC_DIR must be an absolute path
# The -s flag copies using symlinks
# The -T ${DST_DIR} ensures the destination is copied as the directory, and not
#  placed as a subdirectory if the destination already exists.
FUNCTION( COPY_DIRECTORY_USING_SYMLINK SRC_DIR DST_DIR )
    EXECUTE_PROCESS( COMMAND cp -as --remove-destination "${SRC_DIR}" -T "${DST_DIR}" )
ENDFUNCTION()

# Copy files, but symlink the *.h5 files (which are the large ones)
FUNCTION( COPY_DIRECTORY_SYMLINK_H5 SRC_DIR DST_DIR)
    # Copy everything but *.h5 files and pseudopotential files
    FILE(COPY "${SRC_DIR}/" DESTINATION "${DST_DIR}"
         PATTERN "*.h5" EXCLUDE
         PATTERN "*.BFD.xml" EXCLUDE)

    # Now find and symlink the *.h5 files and psuedopotential files
    FILE(GLOB_RECURSE H5 "${SRC_DIR}/*.h5" "${SRC_DIR}/*.BFD.xml")
    FOREACH(F IN LISTS H5)
      FILE(RELATIVE_PATH R "${SRC_DIR}" "${F}")
      #MESSAGE("Creating symlink from  ${SRC_DIR}/${R} to ${DST_DIR}/${R}")
      EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink "${SRC_DIR}/${R}" "${DST_DIR}/${R}")
    ENDFOREACH()
ENDFUNCTION()

# Control copy vs. symlink with top-level variable
FUNCTION( COPY_DIRECTORY_MAYBE_USING_SYMLINK SRC_DIR DST_DIR )
  IF (QMC_SYMLINK_TEST_FILES)
    #COPY_DIRECTORY_USING_SYMLINK("${SRC_DIR}" "${DST_DIR}")
    COPY_DIRECTORY_SYMLINK_H5("${SRC_DIR}" "${DST_DIR}" )
  ELSE()
    COPY_DIRECTORY("${SRC_DIR}" "${DST_DIR}")
  ENDIF()
ENDFUNCTION()

# Symlink or copy an individual file
FUNCTION(MAYBE_SYMLINK SRC_DIR DST_DIR)
  IF (QMC_SYMLINK_TEST_FILES)
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink "${SRC_DIR}" "${DST_DIR}")
  ELSE()
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy "${SRC_DIR}" "${DST_DIR}")
  ENDIF()
ENDFUNCTION()


# Macro to add the dependencies and libraries to an executable
MACRO( ADD_QMC_EXE_DEP EXE )
    # Add the package dependencies
    TARGET_LINK_LIBRARIES(${EXE} qmc qmcdriver qmcham qmcwfs qmcbase qmcutil adios_config)
    FOREACH(l ${QMC_UTIL_LIBS})
        TARGET_LINK_LIBRARIES(${EXE} ${l})
    ENDFOREACH(l ${QMC_UTIL_LIBS})
    IF(ENABLE_TAU_PROFILE)
        TARGET_LINK_LIBRARIES(${EXE} tau)
    ENDIF(ENABLE_TAU_PROFILE)
    IF(MPI_LIBRARY)
        TARGET_LINK_LIBRARIES(${EXE} ${MPI_LIBRARY})
    ENDIF(MPI_LIBRARY)
ENDMACRO()

#############################################################
# Add tests to ctest
#############################################################
# Useful macros to compile and run an executable:
#   ADD_QMC_PROVISIONAL_TEST( EXECFILE )
#       Add a provisional test that will be compiled 
#       but not executed
#   ADD_QMC_TEST( EXECFILE ARGS )
#       Add a serial test passing the given args to the test
#   ADD_QMC_TEST_PARALLEL( EXECFILE PROCS ARGS )
#       Add a parallel test with the given number of
#       processors and arguments
#   ADD_QMC_TEST_1_2_4( EXECFILE args )
#       Add a test that will run on 1, 2, and 4 processors
#   ADD_QMC_WEEKLY_TEST( EXECFILE PROCS ARGS )
#       Add a "WEEKLY" parallel test with the given number 
#       of processors and arguments
#   ADD_QMC_TEST_THREAD_MPI( EXECFILE PROCS THREADS ARGS )
#       Add a parallel test with the given number of threads
#       per processes
# Useful macros to run an existing executable:
#   RUN_QMC_APP( TESTNAME SRC_DIR PROCS THREADS ARGS )
#############################################################

# Add a provisional test
FUNCTION( ADD_QMC_PROVISIONAL_TEST EXEFILE )
    # Change the output directory so we do not install in bin
    SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
    # Check if we actually want to add the test
    # SET( EXCLUDE_TESTS_FROM_ALL 1 )
    # Check if test has already been added
    SET( tmp )
    IF ( TARGET ${EXEFILE} )
        GET_TARGET_PROPERTY(tmp ${EXEFILE} LOCATION)
        STRING(REGEX REPLACE "//" "/" tmp "${tmp}" )        
    ENDIF()
    IF ( NOT tmp )
        # The target has not been added
        SET( CXXFILE ${EXEFILE}.cpp )
        SET( TESTS_SO_FAR ${TESTS_SO_FAR} ${EXEFILE} )
        IF ( NOT EXCLUDE_TESTS_FROM_ALL )
            ADD_EXECUTABLE( ${EXEFILE} ${CXXFILE} )
        ELSE()
            ADD_EXECUTABLE( ${EXEFILE} EXCLUDE_FROM_ALL ${CXXFILE} )
        ENDIF()
        ADD_QMC_EXE_DEP( ${EXEFILE} )
    ELSEIF( ${tmp} STREQUAL "${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}" )
        # The correct target has already been added
    ELSEIF( ${tmp} STREQUAL "${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}.exe" )
        # The correct target has already been added
    ELSEIF( ${tmp} STREQUAL "${CMAKE_CURRENT_BINARY_DIR}/$(Configuration)/${EXEFILE}.exe" )
        # The correct target has already been added
    ELSEIF( ${tmp} STREQUAL "${CMAKE_CURRENT_BINARY_DIR}/$(OutDir)/${EXEFILE}.exe" )
        # The correct target has already been added
    ELSE()
        # We are trying to add 2 different tests with the same name
        MESSAGE ( "Existing test: ${tmp}" )
        MESSAGE ( "New test:      ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}" )
        MESSAGE ( FATAL_ERROR "Trying to add 2 different tests with the same name" )
    ENDIF()
ENDFUNCTION()


# Macro to create the test name
MACRO( CREATE_TEST_NAME TEST ${ARGN} )
    SET( TESTNAME "${TEST}" )
    FOREACH( tmp ${ARGN} )
        SET( TESTNAME "${TESTNAME}--${tmp}")
    endforeach()
    # STRING(REGEX REPLACE "--" "-" TESTNAME ${TESTNAME} )
ENDMACRO()


# Add a executable as a test
FUNCTION( ADD_QMC_TEST EXEFILE ${ARGN} )
    ADD_QMC_PROVISIONAL_TEST ( ${EXEFILE} )
    CREATE_TEST_NAME( ${EXEFILE} ${ARGN} )
    GET_TARGET_PROPERTY(EXE ${EXEFILE} LOCATION)
    STRING(REGEX REPLACE "\\$\\(Configuration\\)" "${CONFIGURATION}" EXE "${EXE}" )
    IF ( USE_MPI_FOR_SERIAL_TESTS )
        ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${EXE} ${ARGN} )
    ELSE()
        ADD_TEST( ${TESTNAME} ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE} ${ARGN} )
    ENDIF()
    SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS 1 )
ENDFUNCTION()


# Add a executable as a weekly test
FUNCTION( ADD_QMC_WEEKLY_TEST EXEFILE PROCS ${ARGN} )
    ADD_QMC_PROVISIONAL_TEST ( ${EXEFILE} )
    GET_TARGET_PROPERTY(EXE ${EXEFILE} LOCATION)
    STRING(REGEX REPLACE "\\$\\(Configuration\\)" "${CONFIGURATION}" EXE "${EXE}" )
    IF( ${PROCS} STREQUAL "1" )
        CREATE_TEST_NAME( "${EXEFILE}_WEEKLY" ${ARGN} )
    ELSEIF( USE_MPI AND NOT (${PROCS} GREATER ${TEST_MAX_PROCS}) )
        CREATE_TEST_NAME( "${EXEFILE}_${PROCS}procs_WEEKLY" ${ARGN} )
    ENDIF()
    IF ( ${PROCS} GREATER ${TEST_MAX_PROCS} )
        MESSAGE("Disabling test ${TESTNAME} (exceeds maximum number of processors ${TEST_MAX_PROCS})")
    ELSEIF( ${PROCS} STREQUAL "1" )
        CREATE_TEST_NAME( "${EXEFILE}_WEEKLY" ${ARGN} )
        IF ( USE_MPI_FOR_SERIAL_TESTS )
            ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${EXE} ${ARGN} )
        ELSE()
            ADD_TEST( ${TESTNAME} ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE} ${ARGN} )
        ENDIF()
        SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS 1 )
    ELSEIF( USE_MPI AND NOT (${PROCS} GREATER ${TEST_MAX_PROCS}) )
        CREATE_TEST_NAME( "${EXEFILE}_${PROCS}procs_WEEKLY" ${ARGN} )
        ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${EXE} ${ARGN} )
        SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS ${PROCS} )
    ENDIF()
ENDFUNCTION()


# Add a executable as a parallel test
FUNCTION( ADD_QMC_TEST_PARALLEL EXEFILE PROCS ${ARGN} )
    ADD_QMC_PROVISIONAL_TEST ( ${EXEFILE} )
    GET_TARGET_PROPERTY(EXE ${EXEFILE} LOCATION)
    STRING(REGEX REPLACE "\\$\\(Configuration\\)" "${CONFIGURATION}" EXE "${EXE}" )
    CREATE_TEST_NAME( "${EXEFILE}_${PROCS}procs" ${ARGN} )
    IF ( NOT USE_MPI )
        MESSAGE("Disabling test ${TESTNAME} (configured without MPI)")
    ELSEIF ( ${PROCS} GREATER ${TEST_MAX_PROCS} )
        MESSAGE("Disabling test ${TESTNAME} (exceeds maximum number of processors ${TEST_MAX_PROCS})")
    ELSE()
        ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${EXE} ${ARGN} )
        SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS ${PROCS} )
    ENDIF()
ENDFUNCTION()


# Add a executable as a parallel 1, 2, 4 processor test
MACRO( ADD_QMC_TEST_1_2_4 EXENAME ${ARGN} )
    ADD_QMC_TEST ( ${EXENAME} ${ARGN} )
    ADD_QMC_TEST_PARALLEL( ${EXENAME} 2 ${ARGN} )
    ADD_QMC_TEST_PARALLEL( ${EXENAME} 4 ${ARGN} )
ENDMACRO()


# Add a executable as a parallel 8, 12, 16 processor test
MACRO( ADD_QMC_TEST_8_12_16 EXENAME ${ARGN} )
    ADD_QMC_TEST_PARALLEL( ${EXENAME} 8 ${ARGN} )
    ADD_QMC_TEST_PARALLEL( ${EXENAME} 12 ${ARGN} )
    ADD_QMC_TEST_PARALLEL( ${EXENAME} 16 ${ARGN} )
ENDMACRO()


# Add a parallel test that may use both MPI and threads
# This allows us to correctly compute the number of processors used by the test
MACRO( ADD_QMC_TEST_THREAD_MPI EXEFILE PROCS THREADS ${ARGN} )
    ADD_QMC_PROVISIONAL_TEST( ${EXEFILE} )
    GET_TARGET_PROPERTY(EXE ${EXEFILE} LOCATION)
    STRING(REGEX REPLACE "\\$\\(Configuration\\)" "${CONFIGURATION}" EXE "${EXE}" )
    CREATE_TEST_NAME( "${EXEFILE}_${PROCS}procs_${THREADS}threads" ${ARGN} )
    MATH( EXPR TOT_PROCS "${PROCS} * ${THREADS}" )
    IF ( ${TOT_PROCS} GREATER ${TEST_MAX_PROCS} )
        MESSAGE("Disabling test ${TESTNAME} (exceeds maximum number of processors ${TEST_MAX_PROCS})")
    ELSEIF ( ( ${PROCS} STREQUAL "1" ) AND NOT USE_EXT_MPI_FOR_SERIAL_TESTS )
        ADD_TEST( ${TESTNAME} ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE} ${ARGN} )
        SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS ${TOT_PROCS} )
    ELSEIF ( USE_MPI )
        ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${EXE} ${ARGN} )
        SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS ${TOT_PROCS} )
    ENDIF()
ENDMACRO()


# Runs qmcpack
#  Note that TEST_ADDED is an output variable
FUNCTION( RUN_QMC_APP_NO_COPY TESTNAME WORKDIR PROCS THREADS TEST_ADDED ${ARGN} )
    MATH( EXPR TOT_PROCS "${PROCS} * ${THREADS}" )
    SET( QMC_APP "${qmcpack_BINARY_DIR}/bin/qmcpack" )
    SET( ${TEST_ADDED} FALSE PARENT_SCOPE )
    IF ( USE_MPI )
        IF ( ${TOT_PROCS} GREATER ${TEST_MAX_PROCS} )
            MESSAGE("Disabling test ${TESTNAME} (exceeds maximum number of processors ${TEST_MAX_PROCS})")
        ELSEIF ( USE_MPI )
            ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${QMC_APP} ${ARGN} )
            SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" 
                PROCESSORS ${TOT_PROCS} WORKING_DIRECTORY ${WORKDIR}
                ENVIRONMENT OMP_NUM_THREADS=${THREADS} )
            SET( ${TEST_ADDED} TRUE PARENT_SCOPE )
        ENDIF()
    ELSE()
        IF ( ( ${PROCS} STREQUAL "1" ) )
            ADD_TEST( ${TESTNAME} ${QMC_APP} ${ARGN} )
            SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" 
                PROCESSORS ${TOT_PROCS} WORKING_DIRECTORY ${WORKDIR}
                ENVIRONMENT OMP_NUM_THREADS=${THREADS} )
            SET( ${TEST_ADDED} TRUE PARENT_SCOPE )
        ENDIF()
    ENDIF()
ENDFUNCTION()

# Runs qmcpack
#  Note that TEST_ADDED is an output variable
FUNCTION( RUN_QMC_APP TESTNAME SRC_DIR PROCS THREADS TEST_ADDED ${ARGN} )
    COPY_DIRECTORY_MAYBE_USING_SYMLINK( "${SRC_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/${TESTNAME}" )
    SET( TEST_ADDED_TEMP FALSE )
    RUN_QMC_APP_NO_COPY( ${TESTNAME} ${CMAKE_CURRENT_BINARY_DIR}/${TESTNAME} ${PROCS} ${THREADS} TEST_ADDED_TEMP ${ARGN} )
    SET( ${TEST_ADDED} ${TEST_ADDED_TEMP} PARENT_SCOPE )
ENDFUNCTION()


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
    LIST(APPEND SCALAR_CHECK_TYPE "kinetic" "totenergy" "eeenergy" "samples" "potential" "ionion" "localecp" "nonlocalecp" "flux" "kinetic_mixed" "kinetic_pure" "eeenergy_mixed" "eeenergy_pure" "potential_pure")
    LIST(APPEND CHECK_SCALAR_FLAG "--ke"    "--le"      "--ee"     "--ts"    "--lp"      "--ii"       "--lpp"    "--nlpp" "--fl" "--ke_m" "--ke_p" "--ee_m" "--ee_p" "--lp_p")

    SET( TEST_ADDED FALSE )
    SET( FULL_NAME "${BASE_NAME}-${PROCS}-${THREADS}" )
    MESSAGE("Adding test ${FULL_NAME}")
    RUN_QMC_APP(${FULL_NAME} ${BASE_DIR} ${PROCS} ${THREADS} TEST_ADDED ${INPUT_FILE})
    IF ( TEST_ADDED )
        SET_PROPERTY(TEST ${FULL_NAME} APPEND PROPERTY LABELS "QMCPACK")
    ENDIF()


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
                SET_PROPERTY( TEST ${TEST_NAME} APPEND PROPERTY LABELS "QMCPACK-checking-results" )
            ENDIF()
        ENDFOREACH(SCALAR_CHECK)
    ENDIF()
ENDFUNCTION()

function(SIMPLE_RUN_AND_CHECK base_name base_dir input_file procs threads check_script)
  
  # "simple run and check" function does 2 things:
  #  1. run qmcpack executable on $input_file located in $base_dir
  #  2. run $check_script located in the same folder ($base_dir)
  # note: NAME, COMMAND, and WORKING_DIRECTORY must be upper case in add_test!

  # build test name
  set(full_name "${base_name}-${procs}-${threads}")
  message("Adding test ${full_name}")

  # add run (task 1)
  set (test_added false)
  RUN_QMC_APP(${full_name} ${base_dir} ${procs} ${threads} test_added ${input_file})
  if ( NOT test_added)
    message(FATAL_ERROR "test ${full_name} cannot be added")
  endif()

  # set up command to run check, assume check_script is in the same folder as input
  set(check_cmd ${CMAKE_CURRENT_BINARY_DIR}/${full_name}/${check_script})
  #message(${check_cmd})

  # add test (task 2)
  set(test_name "${full_name}-check") # hard-code for single test
  set(work_dir "${CMAKE_CURRENT_BINARY_DIR}/${full_name}")
  #message(${work_dir})
  add_test(NAME "${test_name}"
    COMMAND "${check_cmd}"
    WORKING_DIRECTORY "${work_dir}"
  )

  # make test depend on the run
  set_property(TEST ${test_name} APPEND PROPERTY DEPENDS ${full_name})

endfunction()
