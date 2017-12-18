
#############################################################
# Functions for adding tests to ctest
#############################################################
# Useful macros to run an existing executable:
#   RUN_QMC_APP
#     Run QMCPACK with a given number of threads and MPI processes
#
#   QMC_RUN_AND_CHECK
#     Run QMCPACK and check scalar output values.  This is the
#     primary function used for system tests.
#
#   SIMPLE_RUN_AND_CHECK
#     Run QMCPACK on the given input file and check output
#     using a specified script
#############################################################

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



# Macro to create the test name
MACRO( CREATE_TEST_NAME TEST ${ARGN} )
    SET( TESTNAME "${TEST}" )
    FOREACH( tmp ${ARGN} )
        SET( TESTNAME "${TESTNAME}--${tmp}")
    endforeach()
    # STRING(REGEX REPLACE "--" "-" TESTNAME ${TESTNAME} )
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
# ---required inputs---
# BASE_NAME - name of test (number of MPI processes, number of threads, and value to check (if applicable)
#             will be appended to get the full test name)
# BASE_DIR - source location of test input files
# PREFIX - prefix for output files
# INPUT_FILE - XML input file to QMCPACK
# PROCS - number of MPI processes
# THREADS - number of OpenMP threads
# SHOULD_SUCCEED - whether the test is expected to pass or fail.  Expected failing tests will not have
#                  the scalar tests added.
# ---optional inputs---
# ---any number of SERIES/SCALAR_VALUES list pairs can be provided
# ---input pairs beyond the first result in the series number being added to the test name
# ---support for this functionality is provided through the ARGN catch-all input list
# SERIES - series index to compute
# SCALAR_VALUES - list of output values to check with check_scalars.py
#                 The list entries alternate between the value name and the value (usually a string with the
#                 both the average and error).

FUNCTION(QMC_RUN_AND_CHECK BASE_NAME BASE_DIR PREFIX INPUT_FILE PROCS THREADS SHOULD_SUCCEED)
    # Map from name of check to appropriate flag for check_scalars.py
    LIST(APPEND SCALAR_CHECK_TYPE "kinetic" "totenergy" "variance" "eeenergy" "samples" "potential" "ionion" "localecp" "nonlocalecp" "flux" "kinetic_mixed" "kinetic_pure" "eeenergy_mixed" "eeenergy_pure" "potential_pure" "totenergy_A" "totenergy_B" "dtotenergy_AB" "ionion_A" "ionion_B" "dionion_AB" "eeenergy_A" "eeenergy_B" "deeenergy_AB")
    LIST(APPEND CHECK_SCALAR_FLAG "--ke"    "--le"    "--va"    "--ee"     "--ts"    "--lp"      "--ii"       "--lpp"    "--nlpp" "--fl" "--ke_m" "--ke_p" "--ee_m" "--ee_p" "--lp_p" "--le_A" "--le_B" "--dle_AB" "--ii_A" "--ii_B" "--dii_AB" "--ee_A" "--ee_B" "--dee_AB")

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
        SET(IDX0 0)
        FOREACH(V ${ARGN})
            MATH(EXPR MOD_IDX0 "${IDX0}%2")
            IF(MOD_IDX0 EQUAL 0)
                #MESSAGE("   SERIES   : ${V}")
                SET(SERIES ${V})
            ENDIF()
            IF(MOD_IDX0 EQUAL 1)
                #MESSAGE("   CHECKLIST: ${V}")
                SET(SCALAR_VALUES ${V})
                FOREACH(SCALAR_CHECK IN LISTS SCALAR_CHECK_TYPE)
                    LIST(FIND ${SCALAR_VALUES} ${SCALAR_CHECK} IDX1)
                    IF (IDX1 GREATER -1)
                        LIST(FIND SCALAR_CHECK_TYPE ${SCALAR_CHECK} IDX)
                        LIST(GET CHECK_SCALAR_FLAG ${IDX} FLAG)
                
                        MATH( EXPR IDX2 "${IDX1} + 1")
                        LIST(GET ${SCALAR_VALUES} ${IDX2} VALUE)

                        IF(IDX0 LESS 2)
                            SET( TEST_NAME "${FULL_NAME}-${SCALAR_CHECK}" )
                        ELSE()
                            SET( TEST_NAME "${FULL_NAME}-${SERIES}-${SCALAR_CHECK}" )
                        ENDIF()
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
            MATH(EXPR IDX0 "${IDX0}+1")
        ENDFOREACH(V)
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
