include(test_labels)

#[=======================================================================[.rst:
ADD_UNIT_TEST
-------------

Registers a unit test with CTest, handling MPI execution, thread counts,
working directory isolation, and symlinking of required input files.

.. command:: ADD_UNIT_TEST

  .. code-block:: cmake

    ADD_UNIT_TEST(
      BASENAME <name>
      TEST_BINARY <executable_target>
      [PROCS <num_mpi_ranks>]
      [THREADS <num_omp_threads>]
      [INPUT_FILES <file1> [<file2> ...]]
      [ARGS <arg1> [<arg2> ...]]
    )

  ``BASENAME``
    Base name for the test. The registered CTest name will be appended
    with ``-r<PROCS>-t<THREADS>``.

  ``TEST_BINARY``
    The executable to run. Usually a generator expression like
    ``$<TARGET_FILE:my_test_exe>``.

  ``PROCS``
    Number of MPI ranks. Defaults to 1. If greater than 1 and MPI is
    disabled in the build, the test is not registered.

  ``THREADS``
    Number of OpenMP threads (sets ``OMP_NUM_THREADS``). Defaults to 1.

  ``INPUT_FILES``
    List of files to symlink or copy into the test's isolated working
    directory before execution. Relative paths are resolved against
    ``CMAKE_CURRENT_SOURCE_DIR``.

  ``ARGS``
    Extra arguments to pass to the test executable.

  The function sets the variable ``ADDED_UNIT_TEST_NAME`` in the parent
  scope containing the fully registered CTest name (or empty if the test
  was not added due to configuration, e.g., MPI missing).
#]=======================================================================]
function(ADD_UNIT_TEST)
  set(options "")
  set(oneValueArgs BASENAME PROCS THREADS TEST_BINARY)
  set(multiValueArgs ARGS INPUT_FILES)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT ARG_BASENAME)
    message(FATAL_ERROR "ADD_UNIT_TEST requires BASENAME")
  endif()
  if(NOT ARG_PROCS)
    set(ARG_PROCS 1)
  endif()
  if(NOT ARG_THREADS)
    set(ARG_THREADS 1)
  endif()
  if(NOT ARG_TEST_BINARY)
    message(FATAL_ERROR "ADD_UNIT_TEST requires TEST_BINARY")
  endif()

  set(TESTNAME "${ARG_BASENAME}-r${ARG_PROCS}-t${ARG_THREADS}")
  message(VERBOSE "Adding test ${TESTNAME}")
  math(EXPR TOT_PROCS "${ARG_PROCS} * ${ARG_THREADS}")
  if(HAVE_MPI)
    add_test(NAME ${TESTNAME} COMMAND ${QMC_GPU_TEST_LAUNCHER} ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ARG_PROCS}
                                      ${MPIEXEC_PREFLAGS} ${ARG_TEST_BINARY} ${ARG_ARGS})
    set(TEST_ADDED TRUE)
  else()
    if((${ARG_PROCS} STREQUAL "1"))
      add_test(NAME ${TESTNAME} COMMAND ${QMC_GPU_TEST_LAUNCHER} ${ARG_TEST_BINARY} ${ARG_ARGS})
      set(TEST_ADDED TRUE)
    else()
      message(VERBOSE "Disabling test ${TESTNAME} (building without MPI)")
    endif()
  endif()

  if(TEST_ADDED)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${TESTNAME}")

    foreach(file IN LISTS ARG_INPUT_FILES)
      cmake_path(GET file FILENAME fname)
      if(IS_ABSOLUTE "${file}")
        maybe_symlink("${file}" "${CMAKE_CURRENT_BINARY_DIR}/${TESTNAME}/${fname}")
      else()
        maybe_symlink("${CMAKE_CURRENT_SOURCE_DIR}/${file}" "${CMAKE_CURRENT_BINARY_DIR}/${TESTNAME}/${fname}")
      endif()
    endforeach()

    set_tests_properties(${TESTNAME} PROPERTIES PROCESSORS ${TOT_PROCS} ENVIRONMENT OMP_NUM_THREADS=${ARG_THREADS}
                                                PROCESSOR_AFFINITY TRUE WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${TESTNAME}")
    if("asan" IN_LIST ENABLE_SANITIZER)
      set_property(
        TEST ${TESTNAME}
        APPEND
        PROPERTY ENVIRONMENT LSAN_OPTIONS=${LSAN_OPTIONS})
    endif()

    set_test_gpu_resources(${TESTNAME})

    if(ENABLE_OFFLOAD)
      set_property(
        TEST ${TESTNAME}
        APPEND
        PROPERTY ENVIRONMENT "OMP_TARGET_OFFLOAD=mandatory")
    endif()

    set_property(
      TEST ${TESTNAME}
      APPEND
      PROPERTY LABELS "unit")

    set(ADDED_UNIT_TEST_NAME ${TESTNAME} PARENT_SCOPE)
  else()
    set(ADDED_UNIT_TEST_NAME "" PARENT_SCOPE)
  endif()
endfunction()

#[=======================================================================[.rst:
make_file_alias
---------------

Creates a symlink (or copies, depending on configuration) for a given file
into the current binary directory under a new alias name. Also appends the
resulting aliased file path to the ``ALIASED_FILES`` list variable.

.. command:: make_file_alias

  .. code-block:: cmake

    make_file_alias(<src_file> <dest_file_name>)
#]=======================================================================]
macro(make_file_alias file dst_fname)
  maybe_symlink("${file}" "${CMAKE_CURRENT_BINARY_DIR}/${dst_fname}")
  list(APPEND ALIASED_FILES "${CMAKE_CURRENT_BINARY_DIR}/${dst_fname}")
endmacro()

#[=======================================================================[.rst:
add_test_target_in_output_location
----------------------------------

Registers a simple test that checks if the output executable for a CMake
target exists in the expected output location (`qmcpack_BINARY_DIR/bin/`).

.. command:: add_test_target_in_output_location

  .. code-block:: cmake

    add_test_target_in_output_location(<target_name> <relative_exe_dir>)
#]=======================================================================]
function(add_test_target_in_output_location TARGET_NAME_TO_TEST EXE_DIR_RELATIVE_TO_BUILD)

  # obtain BASE_NAME
  get_target_property(BASE_NAME ${TARGET_NAME_TO_TEST} OUTPUT_NAME)
  if(NOT BASE_NAME)
    set(BASE_NAME ${TARGET_NAME_TO_TEST})
  endif()

  set(TESTNAME build_output_${TARGET_NAME_TO_TEST}_exists)
  add_test(NAME ${TESTNAME} COMMAND ls ${qmcpack_BINARY_DIR}/bin/${BASE_NAME})

  set_property(
    TEST ${TESTNAME}
    APPEND
    PROPERTY LABELS "unit;deterministic")
endfunction()
