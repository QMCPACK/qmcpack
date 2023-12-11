# ctest script for building, running, and submitting the test results
# Usage:  ctest -s script,build
#   build = debug / optimized / valgrind / coverage
# Note: this test will use use the number of processors defined in the variable N_PROCS,
#   the environment variables
#   N_PROCS, or the number of processors available (if not specified)
#   N_PROCS_BUILD, or N_PROCS (if not specified)
#   N_CONCURRENT_TESTS, or N_PROCS (if not specified)
#   TEST_SITE_NAME, or HOSTNAME (if not specified)

# Get the source directory based on the current directory
if(NOT DEFINED QMC_SOURCE_DIR)
  set(QMC_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/..")
endif()

# Check that we specified the build type to run
if(NOT DEFINED CTEST_BUILD_NAME)
  set(CTEST_BUILD_NAME "$ENV{QMCPACK_TEST_SUBMIT_NAME}")
endif()
if(NOT DEFINED CTEST_BUILD_NAME)
  set(CTEST_BUILD_NAME "QMCPACK_TEST_SUBMIT_NAME-unset")
endif()
if(NOT DEFINED CTEST_SCRIPT_ARG)
  message(FATAL_ERROR "No build specified: ctest -S /path/to/script,build (debug/optimized/valgrind")
elseif(${CTEST_SCRIPT_ARG} STREQUAL "debug")
  set(CMAKE_BUILD_TYPE "Debug")
  set(USE_VALGRIND FALSE)
elseif(
  (${CTEST_SCRIPT_ARG} STREQUAL "optimized")
  OR (${CTEST_SCRIPT_ARG} STREQUAL "opt")
  OR (${CTEST_SCRIPT_ARG} STREQUAL "release"))
  set(CMAKE_BUILD_TYPE "Release")
  set(CTEST_COVERAGE_COMMAND)
  set(ENABLE_GCOV "false")
  set(USE_VALGRIND FALSE)
elseif(${CTEST_SCRIPT_ARG} STREQUAL "valgrind")
  set(CMAKE_BUILD_TYPE "Debug")
  set(CTEST_COVERAGE_COMMAND)
  set(ENABLE_GCOV "false")
  set(USE_VALGRIND TRUE)
elseif(${CTEST_SCRIPT_ARG} STREQUAL "coverage")
  set(CMAKE_BUILD_TYPE "Debug")
  set(CTEST_COVERAGE_COMMAND "gcov")
  set(ENABLE_GCOV "true")
  set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-coverage")
else()
  message(FATAL_ERROR "Invalid build (${CTEST_SCRIPT_ARG}): ctest -S /path/to/script,build (debug/opt/valgrind")
endif()

# Set the number of processors
if(NOT DEFINED N_PROCS)
  set(N_PROCS $ENV{N_PROCS})
endif()
if(NOT DEFINED N_PROCS)
  set(N_PROCS 1)
  # Linux:
  set(cpuinfo_file "/proc/cpuinfo")
  if(EXISTS "${cpuinfo_file}")
    file(STRINGS "${cpuinfo_file}" procs REGEX "^processor.: [0-9]+$")
    list(LENGTH procs N_PROCS)
  endif()
  # Mac:
  if(APPLE)
    find_program(cmd_sys_pro "system_profiler")
    if(cmd_sys_pro)
      execute_process(COMMAND ${cmd_sys_pro} OUTPUT_VARIABLE info)
      string(REGEX REPLACE "^.*Total Number of Cores: ([0-9]+).*$" "\\1" N_PROCS "${info}")
    endif()
  endif()
  # Windows:
  if(WIN32)
    set(N_PROCS "$ENV{NUMBER_OF_PROCESSORS}")
  endif()
endif()

# Set the number of processors for compilation and running tests
if(NOT DEFINED N_PROCS_BUILD)
  if(DEFINED ENV{N_PROCS_BUILD})
    set(N_PROCS_BUILD $ENV{N_PROCS_BUILD})
  else()
    set(N_PROCS_BUILD ${N_PROCS})
  endif()
endif()
if(NOT DEFINED N_CONCURRENT_TESTS)
  if(DEFINED ENV{N_CONCURRENT_TESTS})
    set(N_CONCURRENT_TESTS $ENV{N_CONCURRENT_TESTS})
  else()
    set(N_CONCURRENT_TESTS ${N_PROCS})
  endif()
endif()

message("Testing with ${N_CONCURRENT_TESTS} processors")

# Set basic variables
set(CTEST_PROJECT_NAME "QMCPACK")
set(CTEST_SOURCE_DIRECTORY "${QMC_SOURCE_DIR}")
set(CTEST_BINARY_DIRECTORY ".")
set(CTEST_DASHBOARD "Nightly")
set(CTEST_TEST_TIMEOUT 900)
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS 500)
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 500)
set(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE 100000)
set(CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE 100000)
set(NIGHTLY_START_TIME "18:00:00 EST")
set(CTEST_NIGHTLY_START_TIME "22:00:00 EST")
set(CTEST_COMMAND "\"${CTEST_EXECUTABLE_NAME}\" -D ${CTEST_DASHBOARD}")

set(CTEST_USE_LAUNCHERS TRUE)
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
if(NOT DEFINED MAKE_CMD)
  set(MAKE_CMD make)
endif()
if(BUILD_SERIAL)
  set(CTEST_BUILD_COMMAND "${MAKE_CMD} -i")
else(BUILD_SERIAL)
  set(CTEST_BUILD_COMMAND "${MAKE_CMD} -i -j ${N_PROCS_BUILD}")
  message("Building with ${N_PROCS_BUILD} processors")
endif(BUILD_SERIAL)

# Set valgrind options
if(USE_VALGRIND)
  set(VALGRIND_COMMAND_OPTIONS
      "--tool=memcheck --leak-check=yes --track-fds=yes --num-callers=50 --show-reachable=yes --suppressions=${QMC_SOURCE_DIR}/src/ValgrindSuppresionFile"
  )
  set(MEMORYCHECK_COMMAND ${VALGRIND_COMMAND})
  set(MEMORYCHECKCOMMAND ${VALGRIND_COMMAND})
  set(CTEST_MEMORYCHECK_COMMAND ${VALGRIND_COMMAND})
  set(CTEST_MEMORYCHECKCOMMAND ${VALGRIND_COMMAND})
  set(CTEST_MEMORYCHECK_COMMAND_OPTIONS ${VALGRIND_COMMAND_OPTIONS})
  set(CTEST_MEMORYCHECKCOMMAND_OPTIONS ${VALGRIND_COMMAND_OPTIONS})
endif()

# Clear the binary directory and create an initial cache
ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
file(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" "CTEST_TEST_CTEST:BOOL=1")

# Set the configure options
set(CTEST_OPTIONS)

if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  if(DEFINED CMAKE_C_COMPILER)
    set(CTEST_OPTIONS "${CTEST_OPTIONS};-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}")
  endif()
  if(DEFINED CMAKE_CXX_COMPILER)
    set(CTEST_OPTIONS "${CTEST_OPTIONS};-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
  endif()
  if(DEFINED CMAKE_C_FLAGS)
    set(CTEST_OPTIONS "${CTEST_OPTIONS};-DCMAKE_C_FLAGS='${CMAKE_C_FLAGS}'")
  endif()
  if(DEFINED CMAKE_CXX_FLAGS)
    set(CTEST_OPTIONS "${CTEST_OPTIONS};-DCMAKE_CXX_FLAGS='${CMAKE_CXX_FLAGS}'")
  endif()
else()
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DCMAKE_TOOLCHAIN_FILE='${CMAKE_TOOLCHAIN_FILE}'")
endif()

if(DEFINED ENABLE_GCOV)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DENABLE_GCOV:BOOL=${ENABLE_GCOV}")
endif()

if(DEFINED CMAKE_BUILD_TYPE)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}")
endif()

if(DEFINED QMC_DATA)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DQMC_DATA='${QMC_DATA}'")
endif()

if(DEFINED QE_BIN)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DQE_BIN='${QE_BIN}'")
endif()

if(DEFINED RMG_BIN)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DRMG_BIN='${RMG_BIN}'")
endif()

if(DEFINED ENABLE_CUDA)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DENABLE_CUDA=${ENABLE_CUDA}")
endif()

if(DEFINED QMC_CUDA2HIP)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DQMC_CUDA2HIP=${QMC_CUDA2HIP}")
endif()

if(DEFINED QMC_COMPLEX)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DQMC_COMPLEX=${QMC_COMPLEX}")
endif()

if(DEFINED QMC_MPI)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DQMC_MPI=${QMC_MPI}")
endif()

if(DEFINED MPIEXEC_EXECUTABLE)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}")
endif()

if(DEFINED MPIEXEC_NUMPROC_FLAG)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}")
endif()

if(DEFINED MPIEXEC_PREFLAGS)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}")
endif()

if(DEFINED QMC_MIXED_PRECISION)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DQMC_MIXED_PRECISION=${QMC_MIXED_PRECISION}")
endif()

if(DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DCMAKE_CUDA_ARCHITECTURES='${CMAKE_CUDA_ARCHITECTURES}'")
endif()

if(DEFINED BUILD_AFQMC)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DBUILD_AFQMC=${BUILD_AFQMC}")
endif()

if(DEFINED BLA_VENDOR)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DBLA_VENDOR='${BLA_VENDOR}'")
endif()

if(DEFINED ENABLE_TIMERS)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DENABLE_TIMERS=${ENABLE_TIMERS}")
endif()

if(DEFINED CMAKE_PREFIX_PATH)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DCMAKE_PREFIX_PATH='${CMAKE_PREFIX_PATH}'")
endif()

if(DEFINED QMC_EXTRA_LIBS)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DQMC_EXTRA_LIBS:STRING='${QMC_EXTRA_LIBS}'")
endif()

if(DEFINED QMC_OPTIONS)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};${QMC_OPTIONS}")
endif()

if(DEFINED HDF5_ROOT)
  set(CTEST_OPTIONS "${CTEST_OPTIONS};-DHDF5_ROOT='${HDF5_ROOT}'")
endif()

message("Configure options:")
message("   ${CTEST_OPTIONS}")

if(NOT DEFINED CTEST_SITE)
  set(CTEST_SITE $ENV{TEST_SITE_NAME})
endif()
if(NOT DEFINED CTEST_SITE)
  site_name(HOSTNAME)
  set(CTEST_SITE ${HOSTNAME})
endif()

# Configure and run the tests
ctest_start("${CTEST_DASHBOARD}")
ctest_update()
ctest_configure(
  BUILD ${CTEST_BINARY_DIRECTORY}
  SOURCE ${CTEST_SOURCE_DIRECTORY}
  OPTIONS "${CTEST_OPTIONS}")

# Run the configure, build and tests
ctest_build()

# Submit the results to oblivion
set(CTEST_DROP_METHOD "https")
set(CTEST_DROP_SITE "cdash.qmcpack.org")
set(CTEST_DROP_LOCATION "/CDash/submit.php?project=QMCPACK")
set(CTEST_DROP_SITE_CDASH TRUE)
set(DROP_SITE_CDASH TRUE)
ctest_submit(PARTS Configure Build)

if(USE_VALGRIND)
  ctest_memcheck(EXCLUDE procs PARALLEL_LEVEL ${N_PROCS})
  ctest_submit(PARTS MemCheck)
elseif(CTEST_COVERAGE_COMMAND)
  # Skip the normal tests when doing coverage
else()
  #    CTEST_TEST( INCLUDE short PARALLEL_LEVEL ${N_PROCS} )
  # run and submit the classified tests to their corresponding track
  ctest_start("${CTEST_DASHBOARD}" TRACK "Deterministic" APPEND)
  ctest_test(INCLUDE_LABEL "deterministic" PARALLEL_LEVEL ${N_CONCURRENT_TESTS})
  ctest_submit(PARTS Test)
  ctest_start("${CTEST_DASHBOARD}" TRACK "Converter" APPEND)
  ctest_test(INCLUDE_LABEL "converter" PARALLEL_LEVEL ${N_CONCURRENT_TESTS})
  ctest_submit(PARTS Test)
  ctest_start("${CTEST_DASHBOARD}" TRACK "Performance" APPEND)
  ctest_test(INCLUDE_LABEL "performance" PARALLEL_LEVEL 16)
  ctest_submit(PARTS Test)
  # run and submit unclassified tests to the default track
  ctest_start("${CTEST_DASHBOARD}" TRACK "${CTEST_DASHBOARD}" APPEND)
  ctest_test(EXCLUDE_LABEL "deterministic|performance|converter|unstable" PARALLEL_LEVEL ${N_CONCURRENT_TESTS})
  ctest_submit(PARTS Test)
  # Only the result checking is placed in the unstable category and parent process must run before it
  # To fullfil this implicit dependency, the unstable category is placed last.
  # A better solution is needed to make the dependency explicit.
  ctest_start("${CTEST_DASHBOARD}" TRACK "Unstable" APPEND)
  ctest_test(INCLUDE_LABEL "unstable" PARALLEL_LEVEL ${N_CONCURRENT_TESTS})
  ctest_submit(PARTS Test)
endif()

if(CTEST_COVERAGE_COMMAND)

  # Path prefix to remove to shorten some file names.  The final SRC_ROOT Should not contain '..'
  get_filename_component(SRC_ROOT2 ${QMC_SOURCE_DIR} DIRECTORY)
  get_filename_component(SRC_ROOT ${SRC_ROOT2} DIRECTORY)

  execute_process(COMMAND "pwd" OUTPUT_VARIABLE CURRENT_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
  #MESSAGE("Using new code coverage path in ${CTEST_SOURCE_DIRECTORY}")
  #MESSAGE("Using new code coverage path bin: ${CTEST_BINARY_DIRECTORY}")
  include("${CTEST_SOURCE_DIRECTORY}/CMake/compareGCOV.cmake")
  # Base test
  clear_gcda(${CTEST_BINARY_DIRECTORY})
  ctest_test(INCLUDE_LABEL coverage)
  file(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY}/tgcov_base_raw)
  generate_gcov(${CTEST_BINARY_DIRECTORY} ${CTEST_BINARY_DIRECTORY}/tgcov_base_raw "USE_LONG_FILE_NAMES" ${SRC_ROOT})
  filter_gcov(${CTEST_BINARY_DIRECTORY}/tgcov_base_raw)

  file(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY}/tgcov_base)
  merge_gcov(${CTEST_BINARY_DIRECTORY}/tgcov_base_raw ${CTEST_BINARY_DIRECTORY}/tgcov_base ${SRC_ROOT})

  # Generate gcov
  clear_gcda(${CTEST_BINARY_DIRECTORY})
  # Remove gcda files
  ctest_test(INCLUDE_LABEL unit)
  # Generate gcov
  file(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY}/tgcov_unit_raw)
  generate_gcov(${CTEST_BINARY_DIRECTORY} ${CTEST_BINARY_DIRECTORY}/tgcov_unit_raw "USE_LONG_FILE_NAMES" ${SRC_ROOT})
  filter_gcov(${CTEST_BINARY_DIRECTORY}/tgcov_unit_raw)

  file(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY}/tgcov_unit)
  merge_gcov(${CTEST_BINARY_DIRECTORY}/tgcov_unit_raw ${CTEST_BINARY_DIRECTORY}/tgcov_unit ${SRC_ROOT})

  # Generate diff
  file(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY}/tgcov_diff)
  compare_gcov(${CTEST_BINARY_DIRECTORY}/tgcov_base ${CTEST_BINARY_DIRECTORY}/tgcov_unit
               ${CTEST_BINARY_DIRECTORY}/tgcov_diff tgcov_diff)

  # create tar file
  create_gcov_tar(${CTEST_BINARY_DIRECTORY} tgcov_unit)
  create_gcov_tar(${CTEST_BINARY_DIRECTORY} tgcov_base)

  file(GLOB DIFF_GCOV_FILES ${CTEST_BINARY_DIRECTORY}/tgcov_diff/*.gcov)

  set(CTEST_BUILD_NAME_ORIGINAL "${CTEST_BUILD_NAME}")
  set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME_ORIGINAL}-diff")
  ctest_start(${CTEST_DASHBOARD})
  if(EXISTS "${CTEST_BINARY_DIRECTORY}/gcov.tar")
    message("submitting ${CTEST_BINARY_DIRECTORY}/gcov.tar")
    ctest_submit(CDASH_UPLOAD "${CTEST_BINARY_DIRECTORY}/gcov.tar" CDASH_UPLOAD_TYPE GcovTar)
  endif()

  set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME_ORIGINAL}-base")
  ctest_start(${CTEST_DASHBOARD})

  if(EXISTS "${CTEST_BINARY_DIRECTORY}/gcov_tgcov_base.tar")
    message("submitting ${CTEST_BINARY_DIRECTORY}/gcov_tgcov_base.tar")
    ctest_submit(CDASH_UPLOAD "${CTEST_BINARY_DIRECTORY}/gcov_tgcov_base.tar" CDASH_UPLOAD_TYPE GcovTar)
  endif()

  set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME_ORIGINAL}-unit")
  ctest_start(${CTEST_DASHBOARD})

  if(EXISTS "${CTEST_BINARY_DIRECTORY}/gcov_tgcov_unit.tar")
    message("submitting ${CTEST_BINARY_DIRECTORY}/gcov_tgcov_unit.tar")
    ctest_submit(CDASH_UPLOAD "${CTEST_BINARY_DIRECTORY}/gcov_tgcov_unit.tar" CDASH_UPLOAD_TYPE GcovTar)
  endif()

endif()

# Clean up
# exec_program("make distclean")
