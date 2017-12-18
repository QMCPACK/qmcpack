# ctest script for building, running, and submitting the test results 
# Usage:  ctest -s script,build
#   build = debug / optimized / valgrind / coverage
# Note: this test will use use the number of processors defined in the variable N_PROCS,
#   the enviornmental variable N_PROCS, or the number of processors availible (if not specified)

# Set platform specific variables
SITE_NAME( HOSTNAME )
IF( ${HOSTNAME} STREQUAL "mbt01" )
    SET( CC mpicc )
    SET( CXX mpicxx )
    SET( CTEST_CMAKE_GENERATOR "Unix Makefiles")
    SET( QMC_DATA /projects/QMCPACK/qmc-data )
    SET( ENV{HDF5_ROOT} /packages/TPLs/install/debug/hdf5 )
    SET( ENV{BOOST_ROOT} /packages/TPLs/install/debug/boost )
    SET( QMC_OPTIONS )
    SET( QMC_OPTIONS "${QMC_OPTIONS};-DHAVE_ACML=1;-DACML_HOME=/packages/acml-5.3.1/gfortran64" )
    SET( QMC_OPTIONS "${QMC_OPTIONS};-DLibxml2_INCLUDE_DIRS=/usr/include/libxml2" )
    SET( QMC_OPTIONS "${QMC_OPTIONS};-DLibxml2_LIBRARY_DIRS=/usr/lib/x86_64-linux-gnu" )
    SET( QMC_OPTIONS "${QMC_OPTIONS};-DFFTW_INCLUDE_DIRS=/usr/include" )
    SET( QMC_OPTIONS "${QMC_OPTIONS};-DFFTW_LIBRARY_DIRS=/usr/lib/x86_64-linux-gnu" )
    SET( QMC_EXTRA_LIBS "-ldl /packages/acml-5.3.1/gfortran64/lib/libacml.a -lgfortran" )
ELSEIF( ${HOSTNAME} MATCHES "oxygen" )
    SET( CC ${CMAKE_C_COMPILER} )
    SET( CXX ${CMAKE_CXX_COMPILER} )
    SET( CTEST_CMAKE_GENERATOR "Unix Makefiles")
ELSEIF( ${HOSTNAME} MATCHES "hyperion" )
# Setup for hyperion.alcf.anl.gov
    SET( CC ${CMAKE_C_COMPILER} )
    SET( CXX ${CMAKE_CXX_COMPILER} )
    SET( CTEST_CMAKE_GENERATOR "Unix Makefiles")
    SET( CTEST_SITE "hyperion.alcf.anl.gov" )
    SET( N_PROCS 64)
ELSEIF( ${HOSTNAME} MATCHES "cori" )
# Setup for cori.nersc.gov phase 1
    SET( CC ${CMAKE_C_COMPILER} )
    SET( CXX ${CMAKE_CXX_COMPILER} )
    SET( CTEST_CMAKE_GENERATOR "Unix Makefiles")
    SET( CTEST_SITE "cori.nersc.gov" )
    SET( N_PROCS 32)
ELSEIF( ${HOSTNAME} MATCHES "eos" )
# Setup for eos.ccs.ornl.gov Cray XC30 Intel E5-2670 Aries interconnect.
# N_NPROCS and N_PROCS_BUILD should be set to respect current OLCF limits
    SET( CC ${CMAKE_C_COMPILER} )
    SET( CXX ${CMAKE_CXX_COMPILER} )
    SET( CTEST_CMAKE_GENERATOR "Unix Makefiles")
    SET( CTEST_SITE "eos.ccs.ornl.gov" )
    SET( N_PROCS_BUILD 8 )
ELSEIF( ${HOSTNAME} MATCHES "titan" )
# Setup for titan.ccs.ornl.gov Cray XK7
# N_NPROCS and N_PROCS_BUILD should be set to respect current OLCF limits
    SET( CC ${CMAKE_C_COMPILER} )
    SET( CXX ${CMAKE_CXX_COMPILER} )
    SET( CTEST_CMAKE_GENERATOR "Unix Makefiles")
    SET( CTEST_SITE "titan.ccs.ornl.gov" )
    SET( N_PROCS_BUILD 8 )
ELSEIF( ${HOSTNAME} MATCHES "cetus" )
# Setup for cetus.alcf.anl.gov BlueGene/Q
    SET( CTEST_CMAKE_GENERATOR "Unix Makefiles")
    SET( CTEST_SITE "cetus.alcf.anl.gov" )
    SET( N_PROCS 24)
    SET( TEST_PARALLEL_LEVEL 1 ) 
ELSEIF( ${HOSTNAME} MATCHES "billmp1" )
# Setup for billmp1.ornl.gov Mac
    SET( CC /opt/mpich/install/mpich-3.2/bin/mpicc )
    SET( CXX /opt/mpich/install/mpich-3.2/bin/mpic++ )
    SET( MPIEXEC /opt/mpich/install/mpich-3.2/bin/mpirun )
    SET( ENV{HDF5_ROOT} /opt/TPLs/install/opt/hdf5 )
    SET( ENV{BOOST_ROOT} /opt/TPLs/install/opt/boost )
    SET( QMC_OPTIONS "${QMC_OPTIONS};-DMPIEXEC=${MPIEXEC}" )
    SET( CTEST_CMAKE_GENERATOR "Unix Makefiles")
    SET( CTEST_SITE "billmp1.ornl.gov" )
    SET( N_PROCS 8)
    IF( ${CTEST_SCRIPT_ARG} STREQUAL "debug" )
        SET( CTEST_BUILD_NAME Clang-Debug-Mac )
    ELSE()
        SET( CTEST_BUILD_NAME Clang-Release-Mac )
    ENDIF()
ELSEIF( ${HOSTNAME} MATCHES "ryzen-box" )
    SET( CC mpicc )
    SET( CXX mpicxx )
    SET( CTEST_CMAKE_GENERATOR "Unix Makefiles")
    SET( QMC_OPTIONS "${QMC_OPTIONS};-DCMAKE_PREFIX_PATH=/opt/OpenBLAS" )
    SET( N_PROCS 16)
ELSE()
    MESSAGE( FATAL_ERROR "Unknown host: ${HOSTNAME}" )
ENDIF()

# Get the source directory based on the current directory
IF ( NOT QMC_SOURCE_DIR )
    SET( QMC_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/.." )
ENDIF()
IF ( NOT MAKE_CMD )
    SET( MAKE_CMD make )
ENDIF()


# Check that we specified the build type to run
SET( USE_MATLAB 1 )
IF ( NOT CTEST_BUILD_NAME )
    SET( CTEST_BUILD_NAME "$ENV{QMCPACK_TEST_SUBMIT_NAME}" )
ENDIF()
IF ( NOT CTEST_BUILD_NAME )
    SET( CTEST_BUILD_NAME "QMCPACK_TEST_SUBMIT_NAME-unset" )
ENDIF()
IF( NOT CTEST_SCRIPT_ARG )
    MESSAGE(FATAL_ERROR "No build specified: ctest -S /path/to/script,build (debug/optimized/valgrind")
ELSEIF( ${CTEST_SCRIPT_ARG} STREQUAL "debug" )
    SET( CMAKE_BUILD_TYPE "Debug" )
    SET( USE_VALGRIND FALSE )
    SET( USE_VALGRIND_MATLAB FALSE )
ELSEIF( (${CTEST_SCRIPT_ARG} STREQUAL "optimized") OR (${CTEST_SCRIPT_ARG} STREQUAL "opt") OR (${CTEST_SCRIPT_ARG} STREQUAL "release") )
    SET( CMAKE_BUILD_TYPE "Release" )
    SET( CTEST_COVERAGE_COMMAND )
    SET( ENABLE_GCOV "false" )
    SET( USE_VALGRIND FALSE )
    SET( USE_VALGRIND_MATLAB FALSE )
ELSEIF( ${CTEST_SCRIPT_ARG} STREQUAL "valgrind" )
    SET( CMAKE_BUILD_TYPE "Debug" )
    SET( CTEST_COVERAGE_COMMAND )
    SET( ENABLE_GCOV "false" )
    SET( USE_VALGRIND TRUE )
    SET( USE_VALGRIND_MATLAB FALSE )
    SET( USE_MATLAB 0 )
ELSEIF( ${CTEST_SCRIPT_ARG} STREQUAL "coverage" )
    SET( CMAKE_BUILD_TYPE "Debug" )
    SET( CTEST_COVERAGE_COMMAND "gcov" )
    SET( ENABLE_GCOV "true" )
    SET( CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-coverage" )
ELSE()
    MESSAGE(FATAL_ERROR "Invalid build (${CTEST_SCRIPT_ARG}): ctest -S /path/to/script,build (debug/opt/valgrind")
ENDIF()


# Set the number of processors
IF( NOT DEFINED N_PROCS )
    SET( N_PROCS $ENV{N_PROCS} )
ENDIF()
IF( NOT DEFINED N_PROCS )
    SET(N_PROCS 1)
    # Linux:
    SET(cpuinfo_file "/proc/cpuinfo")
    IF(EXISTS "${cpuinfo_file}")
        FILE(STRINGS "${cpuinfo_file}" procs REGEX "^processor.: [0-9]+$")
        list(LENGTH procs N_PROCS)
    ENDIF()
    # Mac:
    IF(APPLE)
        find_program(cmd_sys_pro "system_profiler")
        if(cmd_sys_pro)
            execute_process(COMMAND ${cmd_sys_pro} OUTPUT_VARIABLE info)
            STRING(REGEX REPLACE "^.*Total Number of Cores: ([0-9]+).*$" "\\1" N_PROCS "${info}")
        ENDIF()
    ENDIF()
    # Windows:
    IF(WIN32)
        SET(N_PROCS "$ENV{NUMBER_OF_PROCESSORS}")
    ENDIF()
ENDIF()

# Set the number of processors
IF( NOT DEFINED N_PROCS_BUILD )
     IF ( DEFINED $ENV{N_PROCS_BUILD} )
        SET( N_PROCS_BUILD $ENV{N_PROCS_BUILD} )
     ENDIF()
ENDIF()	

# Set basic variables
SET( CTEST_PROJECT_NAME "QMCPACK" )
SET( CTEST_SOURCE_DIRECTORY "${QMC_SOURCE_DIR}" )
SET( CTEST_BINARY_DIRECTORY "." )
SET( CTEST_DASHBOARD "Nightly" )
SET( CTEST_TEST_TIMEOUT 900 )
SET( CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS 500 )
SET( CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 500 )
SET( CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE 100000 )
SET( CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE 100000 )
SET( NIGHTLY_START_TIME "18:00:00 EST" )
SET( CTEST_NIGHTLY_START_TIME "22:00:00 EST" )
SET( CTEST_COMMAND "\"${CTEST_EXECUTABLE_NAME}\" -D ${CTEST_DASHBOARD}" )
IF ( BUILD_SERIAL )
    SET( CTEST_BUILD_COMMAND "${MAKE_CMD} -i" )
ELSEIF ( DEFINED N_PROCS_BUILD )
    SET( CTEST_BUILD_COMMAND "${MAKE_CMD} -i -j ${N_PROCS_BUILD}" )
ELSE()
    SET( CTEST_BUILD_COMMAND "${MAKE_CMD} -i -j ${N_PROCS}" )
ENDIF()

MESSAGE("Building with ${N_PROCS_BUILD} processors")
MESSAGE("Testing with ${N_PROCS} processors")

# Set valgrind options
SET( VALGRIND_COMMAND_OPTIONS  "--tool=memcheck --leak-check=yes --track-fds=yes --num-callers=50 --show-reachable=yes --suppressions=${QMC_SOURCE_DIR}/src/ValgrindSuppresionFile" )
IF ( USE_VALGRIND )
    SET( MEMORYCHECK_COMMAND ${VALGRIND_COMMAND} )
    SET( MEMORYCHECKCOMMAND ${VALGRIND_COMMAND} )
    SET( CTEST_MEMORYCHECK_COMMAND ${VALGRIND_COMMAND} )
    SET( CTEST_MEMORYCHECKCOMMAND ${VALGRIND_COMMAND} )
    SET( CTEST_MEMORYCHECK_COMMAND_OPTIONS ${VALGRIND_COMMAND_OPTIONS} )
    SET( CTEST_MEMORYCHECKCOMMAND_OPTIONS  ${VALGRIND_COMMAND_OPTIONS} )
ENDIF()


# Clear the binary directory and create an initial cache
CTEST_EMPTY_BINARY_DIRECTORY( ${CTEST_BINARY_DIRECTORY} )
FILE(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" "CTEST_TEST_CTEST:BOOL=1")


# Set the configure options
IF ( NOT DEFINED CMAKE_TOOLCHAIN_FILE )
  SET( CTEST_OPTIONS )
  SET( CTEST_OPTIONS "-DCMAKE_C_COMPILER=${CC};-DCMAKE_CXX_COMPILER=${CXX}" )
  SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DCMAKE_C_FLAGS='${C_FLAGS}';-DCMAKE_CXX_FLAGS='${CXX_FLAGS}'" )
  SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DENABLE_GCOV:BOOL=${ENABLE_GCOV}" )
  SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}" )
  SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DQMC_DATA='${QMC_DATA}'" )
  SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DQMC_CUDA='${QMC_CUDA}'" )
  SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DQE_BIN='${QE_BIN}'" )
ELSE()
  SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DCMAKE_TOOLCHAIN_FILE='${CMAKE_TOOLCHAIN_FILE}'" )
ENDIF()

IF ( QMC_COMPLEX )
  SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DQMC_COMPLEX='${QMC_COMPLEX}'" )
ENDIF()

IF ( DEFINED QMC_MPI )
  SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DQMC_MPI='${QMC_MPI}'" )
ENDIF()

IF ( QMC_MIXED_PRECISION )
  SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DQMC_MIXED_PRECISION='${QMC_MIXED_PRECISION}'" )
ENDIF()

IF ( ENABLE_SOA )
  SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DENABLE_SOA='${ENABLE_SOA}'" )
ENDIF()

IF ( BUILD_AFQMC )
  SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DBUILD_AFQMC='${BUILD_AFQMC}'" )
ENDIF()

IF ( BLA_VENDOR )
  SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DBLA_VENDOR=${BLA_VENDOR}" )
ENDIF()

IF ( ENABLE_TIMERS )
  SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DENABLE_TIMERS=${ENABLE_TIMERS}" )
ENDIF()

IF ( CMAKE_PREFIX_PATH )
  SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}" )
ENDIF()

IF ( QMC_EXTRA_LIBS )
    SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DQMC_EXTRA_LIBS:STRING='${QMC_EXTRA_LIBS}'" )
ENDIF()

IF ( QMC_OPTIONS )
    SET( CTEST_OPTIONS "${CTEST_OPTIONS};${QMC_OPTIONS}" )
ENDIF()

IF ( HDF5_ROOT )
    SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DHDF5_ROOT=${HDF5_ROOT}" )
ENDIF()

MESSAGE("Configure options:")
MESSAGE("   ${CTEST_OPTIONS}")


# Configure and run the tests
IF ( NOT DEFINED CTEST_SITE )
     SET( CTEST_SITE ${HOSTNAME} )
ENDIF()
CTEST_START("${CTEST_DASHBOARD}")
CTEST_UPDATE()
CTEST_CONFIGURE(
    BUILD   ${CTEST_BINARY_DIRECTORY}
    SOURCE  ${CTEST_SOURCE_DIRECTORY}
    OPTIONS "${CTEST_OPTIONS}"
)


IF ( DEFINED CMAKE_TOOLCHAIN_FILE )
#need to trigger the cmake configuration twice
CTEST_CONFIGURE(
    BUILD   ${CTEST_BINARY_DIRECTORY}
    SOURCE  ${CTEST_SOURCE_DIRECTORY}
    OPTIONS "${CTEST_OPTIONS}"
)
ENDIF()

# Run the configure, build and tests
CTEST_BUILD()
IF ( USE_VALGRIND )
    CTEST_MEMCHECK( EXCLUDE procs   PARALLEL_LEVEL ${N_PROCS} )
ELSEIF (CTEST_COVERAGE_COMMAND)
  # Skip the normal tests when doing coverage
ELSE()
#    CTEST_TEST( INCLUDE short PARALLEL_LEVEL ${N_PROCS} )
    IF( DEFINED TEST_PARALLEL_LEVEL )
         CTEST_TEST( PARALLEL_LEVEL ${TEST_PARALLEL_LEVEL} )
    ELSE()
         CTEST_TEST( PARALLEL_LEVEL ${N_PROCS} )
    ENDIF()
ENDIF()


# Submit the results to oblivion
SET( CTEST_DROP_METHOD "https" )
SET( CTEST_DROP_SITE "cdash.qmcpack.org" )
SET( CTEST_DROP_LOCATION "/CDash/submit.php?project=QMCPACK" )
SET( CTEST_DROP_SITE_CDASH TRUE )
SET( DROP_SITE_CDASH TRUE )
CTEST_SUBMIT()

IF( CTEST_COVERAGE_COMMAND )

  # Path prefix to remove to shorten some file names.  The final SRC_ROOT Should not contain '..'
  GET_FILENAME_COMPONENT(SRC_ROOT2 ${QMC_SOURCE_DIR} DIRECTORY)
  GET_FILENAME_COMPONENT(SRC_ROOT ${SRC_ROOT2} DIRECTORY)

  EXECUTE_PROCESS(COMMAND "pwd" OUTPUT_VARIABLE CURRENT_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
  #MESSAGE("Using new code coverage path in ${CTEST_SOURCE_DIRECTORY}")
  #MESSAGE("Using new code coverage path bin: ${CTEST_BINARY_DIRECTORY}")
  INCLUDE("${CTEST_SOURCE_DIRECTORY}/CMake/compareGCOV.cmake")
  # Base test
  CLEAR_GCDA(${CTEST_BINARY_DIRECTORY})
  CTEST_TEST(INCLUDE_LABEL coverage)
  FILE(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY}/tgcov_base_raw)
  GENERATE_GCOV(${CTEST_BINARY_DIRECTORY} ${CTEST_BINARY_DIRECTORY}/tgcov_base_raw "USE_LONG_FILE_NAMES" ${SRC_ROOT})
  FILTER_GCOV(${CTEST_BINARY_DIRECTORY}/tgcov_base_raw)

  FILE(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY}/tgcov_base)
  MERGE_GCOV(${CTEST_BINARY_DIRECTORY}/tgcov_base_raw ${CTEST_BINARY_DIRECTORY}/tgcov_base ${SRC_ROOT})

  # Generate gcov
  CLEAR_GCDA(${CTEST_BINARY_DIRECTORY})
  # Remove gcda files
  CTEST_TEST(INCLUDE_LABEL unit)
  # Generate gcov
  FILE(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY}/tgcov_unit_raw)
  GENERATE_GCOV(${CTEST_BINARY_DIRECTORY} ${CTEST_BINARY_DIRECTORY}/tgcov_unit_raw "USE_LONG_FILE_NAMES" ${SRC_ROOT})
  FILTER_GCOV(${CTEST_BINARY_DIRECTORY}/tgcov_unit_raw)

  FILE(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY}/tgcov_unit)
  MERGE_GCOV(${CTEST_BINARY_DIRECTORY}/tgcov_unit_raw ${CTEST_BINARY_DIRECTORY}/tgcov_unit ${SRC_ROOT})

  # Generate diff
  FILE(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY}/tgcov_diff)
  COMPARE_GCOV(${CTEST_BINARY_DIRECTORY}/tgcov_base
               ${CTEST_BINARY_DIRECTORY}/tgcov_unit
               ${CTEST_BINARY_DIRECTORY}/tgcov_diff
               tgcov_diff)

  # create tar file
  CREATE_GCOV_TAR(${CTEST_BINARY_DIRECTORY} tgcov_unit)
  CREATE_GCOV_TAR(${CTEST_BINARY_DIRECTORY} tgcov_base)

  FILE(GLOB DIFF_GCOV_FILES ${CTEST_BINARY_DIRECTORY}/tgcov_diff/*.gcov)


  SET( CTEST_BUILD_NAME_ORIGINAL "${CTEST_BUILD_NAME}")
  SET( CTEST_BUILD_NAME "${CTEST_BUILD_NAME_ORIGINAL}-diff" )
  CTEST_START(${CTEST_DASHBOARD})
  IF(EXISTS "${CTEST_BINARY_DIRECTORY}/gcov.tar")
    MESSAGE("submitting ${CTEST_BINARY_DIRECTORY}/gcov.tar")
    CTEST_SUBMIT(CDASH_UPLOAD "${CTEST_BINARY_DIRECTORY}/gcov.tar"
      CDASH_UPLOAD_TYPE GcovTar)
  ENDIF()

  SET( CTEST_BUILD_NAME "${CTEST_BUILD_NAME_ORIGINAL}-base" )
  CTEST_START(${CTEST_DASHBOARD})

  IF(EXISTS "${CTEST_BINARY_DIRECTORY}/gcov_tgcov_base.tar")
    MESSAGE("submitting ${CTEST_BINARY_DIRECTORY}/gcov_tgcov_base.tar")
    CTEST_SUBMIT(CDASH_UPLOAD "${CTEST_BINARY_DIRECTORY}/gcov_tgcov_base.tar"
      CDASH_UPLOAD_TYPE GcovTar)
  ENDIF()

  SET( CTEST_BUILD_NAME "${CTEST_BUILD_NAME_ORIGINAL}-unit" )
  CTEST_START(${CTEST_DASHBOARD})


  IF(EXISTS "${CTEST_BINARY_DIRECTORY}/gcov_tgcov_unit.tar")
    MESSAGE("submitting ${CTEST_BINARY_DIRECTORY}/gcov_tgcov_unit.tar")
    CTEST_SUBMIT(CDASH_UPLOAD "${CTEST_BINARY_DIRECTORY}/gcov_tgcov_unit.tar"
      CDASH_UPLOAD_TYPE GcovTar)
  ENDIF()

ENDIF()

# Clean up
# exec_program("make distclean")


