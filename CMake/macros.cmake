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

include(test_labels)

# Function to copy a directory
function(COPY_DIRECTORY SRC_DIR DST_DIR)
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${SRC_DIR}" "${DST_DIR}")
endfunction()

# Create symlinks for a list of files.
function(SYMLINK_LIST_OF_FILES FILENAMES DST_DIR)
  foreach(F IN LISTS FILENAMES)
    get_filename_component(NAME_ONLY ${F} NAME)
    file(CREATE_LINK ${F} "${DST_DIR}/${NAME_ONLY}" SYMBOLIC)
  endforeach()
endfunction()

# Function to copy a directory using symlinks for the files to save storage space.
# Subdirectories are ignored.
# SRC_DIR must be an absolute path
# The -s flag copies using symlinks
# The -t ${DST_DIR} ensures the destination must be a directory
function(COPY_DIRECTORY_USING_SYMLINK SRC_DIR DST_DIR)
  file(MAKE_DIRECTORY "${DST_DIR}")
  # Find all the files but not subdirectories
  file(
    GLOB FILE_ONLY_NAMES
    LIST_DIRECTORIES TRUE
    "${SRC_DIR}/*")
  symlink_list_of_files("${FILE_ONLY_NAMES}" "${DST_DIR}")
endfunction()

# Copy selected files only. h5, pseudopotentials, wavefunction, structure and the used one input file are copied.
function(COPY_DIRECTORY_USING_SYMLINK_LIMITED SRC_DIR DST_DIR ${ARGN})
  file(MAKE_DIRECTORY "${DST_DIR}")
  # Find all the files but not subdirectories
  file(
    GLOB FILE_FOLDER_NAMES
    LIST_DIRECTORIES TRUE
    "${SRC_DIR}/qmc_ref"
    "${SRC_DIR}/qmc-ref"
    "${SRC_DIR}/*.h5"
    "${SRC_DIR}/*.opt.xml"
    "${SRC_DIR}/*.ncpp.xml"
    "${SRC_DIR}/*.BFD.xml"
    "${SRC_DIR}/*.ccECP.xml"
    "${SRC_DIR}/*.py"
    "${SRC_DIR}/*.sh"
    "${SRC_DIR}/*.restart.xml"
    "${SRC_DIR}/Li.xml"
    "${SRC_DIR}/H.xml"
    "${SRC_DIR}/*.L2_test.xml"
    "${SRC_DIR}/*.opt_L2.xml"
    "${SRC_DIR}/*.wfnoj.xml"
    "${SRC_DIR}/*.wfj*.xml"
    "${SRC_DIR}/*.wfs*.xml"
    "${SRC_DIR}/*.wfn*.xml"
    "${SRC_DIR}/*.cuspInfo.xml"
    "${SRC_DIR}/*.H*.xml"
    "${SRC_DIR}/*.structure.xml"
    "${SRC_DIR}/*ptcl.xml")
  symlink_list_of_files("${FILE_FOLDER_NAMES}" "${DST_DIR}")
  list(TRANSFORM ARGN PREPEND "${SRC_DIR}/")
  symlink_list_of_files("${ARGN}" "${DST_DIR}")
endfunction()

# Control copy vs. symlink with top-level variable
function(COPY_DIRECTORY_MAYBE_USING_SYMLINK SRC_DIR DST_DIR ${ARGN})
  if(QMC_SYMLINK_TEST_FILES)
    copy_directory_using_symlink_limited("${SRC_DIR}" "${DST_DIR}" ${ARGN})
  else()
    copy_directory("${SRC_DIR}" "${DST_DIR}")
  endif()
endfunction()

# Symlink or copy an individual file
function(MAYBE_SYMLINK SRC_DIR DST_DIR)
  if(QMC_SYMLINK_TEST_FILES)
    file(CREATE_LINK ${SRC_DIR} ${DST_DIR} SYMBOLIC)
  else()
    file(COPY ${SRC_DIR} DESTINATION ${DST_DIR})
  endif()
endfunction()

# Macro to create the test name
macro(CREATE_TEST_NAME TEST ${ARGN})
  set(TESTNAME "${TEST}")
  foreach(tmp ${ARGN})
    set(TESTNAME "${TESTNAME}--${tmp}")
  endforeach()
  # STRING(REGEX REPLACE "--" "-" TESTNAME ${TESTNAME} )
endmacro()

# Runs qmcpack
#  Note that TEST_ADDED is an output variable
function(
  RUN_QMC_APP_NO_COPY
  TESTNAME
  WORKDIR
  PROCS
  THREADS
  TEST_ADDED
  TEST_LABELS
  ${ARGN})
  math(EXPR TOT_PROCS "${PROCS} * ${THREADS}")
  set(QMC_APP $<TARGET_FILE:qmcpack>)
  set(TEST_ADDED_TEMP FALSE)

  if(NOT QMC_OMP)
    if(${THREADS} GREATER 1)
      message(VERBOSE
              "Disabling test ${TESTNAME} (exceeds maximum number of threads=1 if OpenMP is disabled -DQMC_OMP=0)")
      return()
    endif()
  endif()

  if(HAVE_MPI)
    if(${TOT_PROCS} GREATER ${TEST_MAX_PROCS})
      message(VERBOSE "Disabling test ${TESTNAME} (exceeds maximum number of processors ${TEST_MAX_PROCS})")
    else()
      add_test(NAME ${TESTNAME} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
                                        ${QMC_APP} ${ARGN})
      set_tests_properties(
        ${TESTNAME}
        PROPERTIES FAIL_REGULAR_EXPRESSION
                   "ERROR"
                   PASS_REGULAR_EXPRESSION
                   "QMCPACK execution completed successfully"
                   PROCESSORS
                   ${TOT_PROCS}
                   PROCESSOR_AFFINITY
                   TRUE
                   WORKING_DIRECTORY
                   ${WORKDIR}
                   ENVIRONMENT
                   OMP_NUM_THREADS=${THREADS})
      set(TEST_ADDED_TEMP TRUE)
    endif()
  else()
    if((${PROCS} STREQUAL "1"))
      add_test(NAME ${TESTNAME} COMMAND ${QMC_APP} ${ARGN})
      set_tests_properties(
        ${TESTNAME}
        PROPERTIES FAIL_REGULAR_EXPRESSION
                   "ERROR"
                   PASS_REGULAR_EXPRESSION
                   "QMCPACK execution completed successfully"
                   PROCESSORS
                   ${TOT_PROCS}
                   PROCESSOR_AFFINITY
                   TRUE
                   WORKING_DIRECTORY
                   ${WORKDIR}
                   ENVIRONMENT
                   OMP_NUM_THREADS=${THREADS})
      set(TEST_ADDED_TEMP TRUE)
    else()
      message(VERBOSE "Disabling test ${TESTNAME} (building without MPI)")
    endif()
  endif()

  if(TEST_ADDED_TEMP
     AND (QMC_CUDA
          OR ENABLE_CUDA
          OR ENABLE_ROCM
          OR ENABLE_OFFLOAD
         ))
    set_tests_properties(${TESTNAME} PROPERTIES RESOURCE_LOCK exclusively_owned_gpus)
  endif()

  set(TEST_LABELS_TEMP "")
  if(TEST_ADDED_TEMP)
    add_test_labels(${TESTNAME} TEST_LABELS_TEMP)
  endif()
  set(${TEST_ADDED}
      ${TEST_ADDED_TEMP}
      PARENT_SCOPE)
  set(${TEST_LABELS}
      ${TEST_LABELS_TEMP}
      PARENT_SCOPE)
endfunction()

# Runs qmcpack
#  Note that TEST_ADDED is an output variable
function(
  RUN_QMC_APP
  TESTNAME
  SRC_DIR
  PROCS
  THREADS
  TEST_ADDED
  TEST_LABELS
  ${ARGN})
  # restrict ARGN to only one file or empty
  list(LENGTH ARGN INPUT_FILE_LENGTH)
  if(INPUT_FILE_LENGTH GREATER 1)
    message(FATAL_ERROR "Incorrect invocation of RUN_QMC_APP by ${TESTNAME}. ARGN value is \"${ARGN}\"")
  endif()

  copy_directory_maybe_using_symlink("${SRC_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/${TESTNAME}" "${ARGN}")
  set(TEST_ADDED_TEMP FALSE)
  set(TEST_LABELS_TEMP "")
  run_qmc_app_no_copy(
    ${TESTNAME}
    ${CMAKE_CURRENT_BINARY_DIR}/${TESTNAME}
    ${PROCS}
    ${THREADS}
    TEST_ADDED_TEMP
    TEST_LABELS_TEMP
    ${ARGN})
  set(${TEST_ADDED}
      ${TEST_ADDED_TEMP}
      PARENT_SCOPE)
  set(${TEST_LABELS}
      ${TEST_LABELS_TEMP}
      PARENT_SCOPE)
endfunction()

if(QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)
  function(QMC_RUN_AND_CHECK)

  endfunction()
  function(QMC_RUN_AND_CHECK_CUSTOM_SCALAR)

  endfunction()
  function(SIMPLE_RUN_AND_CHECK)

  endfunction()
else(QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)

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

  function(
    QMC_RUN_AND_CHECK
    BASE_NAME
    BASE_DIR
    PREFIX
    INPUT_FILE
    PROCS
    THREADS
    SHOULD_SUCCEED)
    # Map from name of check to appropriate flag for check_scalars.py
    list(
      APPEND
      SCALAR_CHECK_TYPE
      "kinetic"
      "totenergy"
      "variance"
      "eeenergy"
      "samples"
      "potential"
      "ionion"
      "localecp"
      "nonlocalecp"
      "flux"
      "kinetic_mixed"
      "kinetic_pure"
      "eeenergy_mixed"
      "eeenergy_pure"
      "potential_pure"
      "totenergy_A"
      "totenergy_B"
      "dtotenergy_AB"
      "ionion_A"
      "ionion_B"
      "dionion_AB"
      "eeenergy_A"
      "eeenergy_B"
      "deeenergy_AB"
      "Eloc"
      "ElocEstim"
      "latdev"
      "EnergyEstim__nume_real"
      "kecorr"
      "mpc")
    list(
      APPEND
      CHECK_SCALAR_FLAG
      "--ke"
      "--le"
      "--va"
      "--ee"
      "--ts"
      "--lp"
      "--ii"
      "--lpp"
      "--nlpp"
      "--fl"
      "--ke_m"
      "--ke_p"
      "--ee_m"
      "--ee_p"
      "--lp_p"
      "--le_A"
      "--le_B"
      "--dle_AB"
      "--ii_A"
      "--ii_B"
      "--dii_AB"
      "--ee_A"
      "--ee_B"
      "--dee_AB"
      "--eloc"
      "--elocest"
      "--latdev"
      "--el"
      "--kec"
      "--mpc")

    set(TEST_ADDED FALSE)
    set(TEST_LABELS "")
    set(FULL_NAME "${BASE_NAME}-${PROCS}-${THREADS}")
    message(VERBOSE "Adding test ${FULL_NAME}")
    run_qmc_app(
      ${FULL_NAME}
      ${BASE_DIR}
      ${PROCS}
      ${THREADS}
      TEST_ADDED
      TEST_LABELS
      ${INPUT_FILE})
    if(TEST_ADDED)
      set_property(
        TEST ${FULL_NAME}
        APPEND
        PROPERTY LABELS "QMCPACK")
    endif()

    if(TEST_ADDED AND NOT SHOULD_SUCCEED)
      set_property(TEST ${FULL_NAME} APPEND PROPERTY WILL_FAIL TRUE)
      #MESSAGE("Test ${FULL_NAME} should fail")
    endif()

    if(TEST_ADDED AND SHOULD_SUCCEED)
      set(IDX0 0)
      foreach(V ${ARGN})
        math(EXPR MOD_IDX0 "${IDX0}%2")
        if(MOD_IDX0 EQUAL 0)
          #MESSAGE("   SERIES   : ${V}")
          set(SERIES ${V})
        endif()
        if(MOD_IDX0 EQUAL 1)
          #MESSAGE("   CHECKLIST: ${V}")
          set(SCALAR_VALUES ${V})
          set(SCALAR_VALUE_FOUND FALSE)
          if(NOT ${SCALAR_VALUES})
            message(FATAL_ERROR "Scalar values not found in variable ${SCALAR_VALUES}")
          endif()
          foreach(SCALAR_CHECK IN LISTS SCALAR_CHECK_TYPE)
            list(FIND ${SCALAR_VALUES} ${SCALAR_CHECK} IDX1)
            if(IDX1 GREATER -1)
              set(SCALAR_VALUE_FOUND TRUE)
              list(FIND SCALAR_CHECK_TYPE ${SCALAR_CHECK} IDX)
              list(GET CHECK_SCALAR_FLAG ${IDX} FLAG)

              math(EXPR IDX2 "${IDX1} + 1")
              list(GET ${SCALAR_VALUES} ${IDX2} VALUE)

              if(IDX0 LESS 2)
                set(TEST_NAME "${FULL_NAME}-${SCALAR_CHECK}")
              else()
                set(TEST_NAME "${FULL_NAME}-${SERIES}-${SCALAR_CHECK}")
              endif()
              #MESSAGE("Adding scalar check ${TEST_NAME}")
              set(CHECK_CMD
                  ${qmcpack_SOURCE_DIR}/tests/scripts/check_scalars.py
                  --ns
                  3
                  --series
                  ${SERIES}
                  -p
                  ${PREFIX}
                  -e
                  2
                  ${FLAG}
                  ${VALUE})
              #MESSAGE("check command = ${CHECK_CMD}")
              add_test(
                NAME ${TEST_NAME}
                COMMAND ${Python3_EXECUTABLE} ${CHECK_CMD}
                WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${FULL_NAME}")
              set_property(TEST ${TEST_NAME} APPEND PROPERTY DEPENDS ${FULL_NAME})
              set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS "QMCPACK-checking-results")
              set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS ${TEST_LABELS})
            endif()
          endforeach(SCALAR_CHECK)
          if(NOT SCALAR_VALUE_FOUND)
            message(FATAL_ERROR "Unknown scalar value to check in ${${SCALAR_VALUES}}")
          endif()
        endif()
        math(EXPR IDX0 "${IDX0}+1")
      endforeach(V)
    endif()
  endfunction()

  # Add a test run and associated scalar checks for a custom named scalar
  # Arguments
  # BASE_NAME - name of test (number of MPI processes, number of threads, and value to check (if applicable)
  #             will be appended to get the full test name)
  # BASE_DIR - source location of test input files
  # PREFIX - prefix for output files
  # INPUT_FILE - XML input file to QMCPACK
  # PROCS - number of MPI processes (default: 1)
  # THREADS - number of OpenMP threads (default: 1)
  # SERIES - series index to compute (default: 0)
  # SCALAR_VALUES - name of list of output values to check with check_scalars.py
  #                 The list entries contain consecutively the name, the value, and the error.
  #                 The name of the variable is passed (instead of the value) in case future support
  #                 for multiple SERIES/SCALAR_VALUES pairs is added

  function(QMC_RUN_AND_CHECK_CUSTOM_SCALAR)
    set(OPTIONS SHOULD_FAIL)
    set(ONE_VALUE_ARGS
        BASE_NAME
        BASE_DIR
        PREFIX
        INPUT_FILE
        PROCS
        THREADS
        SERIES
        SCALAR_VALUES
        EQUILIBRATION)
    # Eventually many want to support multiple SERIES/SCALAR_VALUES pairs
    #SET(MULTI_VALUE_ARGS SERIES SCALAR_VALUES)

    cmake_parse_arguments(QRC "${options}" "${ONE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}" ${ARGN})

    set(PROCS 1)
    if(QRC_PROCS)
      set(PROCS ${QRC_PROCS})
    endif()

    set(THREADS 1)
    if(QRC_THREADS)
      set(THREADS ${QRC_THREADS})
    endif()

    set(BASE_NAME ${QRC_BASE_NAME})
    set(BASE_DIR ${QRC_BASE_DIR})
    set(PREFIX ${QRC_PREFIX})
    set(INPUT_FILE ${QRC_INPUT_FILE})

    set(EQUIL 2)
    if(DEFINED QRC_EQUILIBRATION)
      set(EQUIL ${QRC_EQUILIBRATION})
    endif()

    set(TEST_ADDED FALSE)
    set(TEST_LABELS "")
    set(FULL_NAME "${BASE_NAME}-${PROCS}-${THREADS}")
    message(VERBOSE "Adding test ${FULL_NAME}")
    run_qmc_app(
      ${FULL_NAME}
      ${BASE_DIR}
      ${PROCS}
      ${THREADS}
      TEST_ADDED
      TEST_LABELS
      ${INPUT_FILE})
    if(TEST_ADDED)
      set_property(TEST ${FULL_NAME} APPEND PROPERTY LABELS "QMCPACK")
    endif()

    if(TEST_ADDED AND SHOULD_FAIL)
      set_property(TEST ${FULL_NAME} APPEND PROPERTY WILL_FAIL TRUE)
    endif()

    if(TEST_ADDED AND NOT SHOULD_FAIL)
      # Derefence the list of scalar values by variable name
      set(SCALAR_VALUES "${${QRC_SCALAR_VALUES}}")

      list(LENGTH SCALAR_VALUES listlen)
      math(EXPR listlen2 "${listlen}-1")
      foreach(sv_idx RANGE 0 ${listlen2} 3)

        math(EXPR sv_idx_p1 "${sv_idx}+1")
        math(EXPR sv_idx_p2 "${sv_idx}+2")

        list(GET SCALAR_VALUES ${sv_idx} SCALAR_NAME)
        list(GET SCALAR_VALUES ${sv_idx_p1} SCALAR_VALUE)
        list(GET SCALAR_VALUES ${sv_idx_p2} SCALAR_ERROR)

        set(SERIES 0)
        if(QRC_SERIES)
          set(SERIES ${QRC_SERIES})
          set(TEST_NAME "${FULL_NAME}-${SERIES}-${SCALAR_NAME}")
        else()
          set(TEST_NAME "${FULL_NAME}-${SCALAR_NAME}")
        endif()
        set(CHECK_CMD
            ${qmcpack_SOURCE_DIR}/tests/scripts/check_scalars.py
            --ns
            3
            --series
            ${SERIES}
            -p
            ${PREFIX}
            -e
            ${EQUIL}
            --name
            ${SCALAR_NAME}
            --ref-value
            ${SCALAR_VALUE}
            --ref-error
            ${SCALAR_ERROR})
        add_test(
          NAME ${TEST_NAME}
          COMMAND ${Python3_EXECUTABLE} ${CHECK_CMD}
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${FULL_NAME}")
        set_property(TEST ${TEST_NAME} APPEND PROPERTY DEPENDS ${FULL_NAME})
        set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS "QMCPACK-checking-results")
        set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS ${TEST_LABELS})
      endforeach()
    endif()
  endfunction()

  function(
    SIMPLE_RUN_AND_CHECK
    base_name
    base_dir
    input_file
    procs
    threads
    check_script)

    # "simple run and check" function does 2 things:
    #  1. run qmcpack executable on $input_file located in $base_dir
    #  2. run $check_script located in the same folder ($base_dir)
    # note: NAME, COMMAND, and WORKING_DIRECTORY must be upper case in add_test!

    # build test name
    set(full_name "${base_name}-${procs}-${threads}")
    message(VERBOSE "Adding test ${full_name}")

    # add run (task 1)
    set(test_added false)
    set(test_labels "")
    run_qmc_app(
      ${full_name}
      ${base_dir}
      ${procs}
      ${threads}
      test_added
      test_labels
      ${input_file})
    if(NOT test_added)
      return()
    endif()

    # set up command to run check, assume check_script is in the same folder as input
    if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/${full_name}/${check_script}")
      set(check_cmd "${CMAKE_CURRENT_BINARY_DIR}/${full_name}/${check_script}")
    elseif(EXISTS "${qmcpack_SOURCE_DIR}/tests/scripts/${check_script}")
      set(check_cmd "${qmcpack_SOURCE_DIR}/tests/scripts/${check_script}")
    else()
      message(FATAL_ERROR "Check script not found: ${check_script}")
    endif()
    #message(${check_cmd})

    # add test (task 2)
    set(test_name "${full_name}-check") # hard-code for single test
    set(work_dir "${CMAKE_CURRENT_BINARY_DIR}/${full_name}")
    #message(${work_dir})

    add_test(
      NAME "${test_name}"
      COMMAND ${Python3_EXECUTABLE} ${check_cmd} ${ARGN}
      WORKING_DIRECTORY "${work_dir}")

    # make test depend on the run
    set_property(TEST ${test_name} APPEND PROPERTY DEPENDS ${full_name})
    set_property(TEST ${test_name} APPEND PROPERTY LABELS ${test_labels})

  endfunction()

endif(QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)

function(COVERAGE_RUN TESTNAME SRC_DIR PROCS THREADS ${ARGN})
  set(FULLNAME "coverage-${TESTNAME}")
  set(TEST_ADDED FALSE)
  set(TEST_LABELS "")
  run_qmc_app(
    ${FULLNAME}
    ${SRC_DIR}
    ${PROCS}
    ${THREADS}
    TEST_ADDED
    TEST_LABELS
    ${ARGN})
  if(TEST_ADDED)
    set_property(TEST ${FULLNAME} APPEND PROPERTY LABELS "coverage")
  endif()
endfunction()

function(
  CPU_LIMIT_RUN
  TESTNAME
  SRC_DIR
  PROCS
  THREADS
  TIME
  ${ARGN})
  set(FULLNAME "cpu_limit-${TESTNAME}")
  set(TEST_ADDED FALSE)
  set(TEST_LABELS "")
  run_qmc_app(
    ${FULLNAME}
    ${SRC_DIR}
    ${PROCS}
    ${THREADS}
    TEST_ADDED
    TEST_LABELS
    ${ARGN})
  if(TEST_ADDED)
    set_property(TEST ${FULLNAME} APPEND PROPERTY TIMEOUT ${TIME})
    set_property(TEST ${FULLNAME} APPEND PROPERTY PASS_REGULAR_EXPRESSION "Time limit reached for")
  endif()
endfunction()
