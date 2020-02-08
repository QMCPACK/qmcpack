######################################################################
# CMake version and policies
######################################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.10.0)

cmake_policy(SET CMP0003 NEW)
IF(NOT $ENV{CRAYPE_VERSION} MATCHES ".")
  cmake_policy(SET CMP0056 NEW)
  # This policy insures that the CMAKE_EXE_LINKER_FLAGS of the calling project
  # are used in try_compile test cmake projects.
  # CHECK_CXX_SOURCE_COMPILES and others depend on try_compile
endif(NOT $ENV{CRAYPE_VERSION} MATCHES ".")
IF(CMAKE_GENERATOR MATCHES "Ninja")
  cmake_policy(SET CMP0058 NEW)
ENDIF(CMAKE_GENERATOR MATCHES "Ninja")
# CMP0074: CMake find_package will use <PackageName>_ROOT CMake variable
# and environment variable in search path.
# various Find<Package>.cmake modules may not follow this policy
# CMP0075: check_include_file_cxx and other follow CMAKE_REQUIRED_LIBRARIES
# when doing link check
IF(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")
  cmake_policy(SET CMP0074 NEW)
  cmake_policy(SET CMP0075 NEW)
ENDIF(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")

######################################################################
# QMCPACK
######################################################################
PROJECT(qmcpack)

SET(QMCPACK_VERSION_MAJOR 3)
SET(QMCPACK_VERSION_MINOR 8)
SET(QMCPACK_VERSION_PATCH 9)
SET(QMCPACK_VERSION "${QMCPACK_VERSION_MAJOR}.${QMCPACK_VERSION_MINOR}.${QMCPACK_VERSION_PATCH}")

######################################################
# Directory where customize cmake files reside
######################################################
SET (PROJECT_CMAKE ${qmcpack_SOURCE_DIR}/CMake)

SET(QMCPACK_UNIT_TEST_DIR ${qmcpack_BINARY_DIR}/tests/bin)

######################################################################
# Verify Python3 available
######################################################################
INCLUDE(CMake/python.cmake)
IF(${CMAKE_VERSION} VERSION_LESS "3.12.0")
  TEST_PYTHON_MODULE( "os"  Python3_FOUND )
ELSE()
  find_package(Python3)
ENDIF()
IF ( NOT Python3_FOUND )
  MESSAGE( FATAL_ERROR "Could not find required python3" )
ENDIF ( NOT Python3_FOUND )
######################################################################
# Verify QE executables present if QE_BIN specified
######################################################################
IF ( DEFINED QE_BIN )
  INCLUDE(CMake/run_qe.cmake)
  VERIFY_QE_PRESENT( ${QE_BIN} )
ENDIF()

######################################################################
# Verify PYSCF package is present
######################################################################
IF ( DEFINED PYSCF_BIN )
  TEST_PYTHON_MODULE(pyscf HAVE_PYSCF)
  IF (NOT HAVE_PYSCF)
    MESSAGE( FATAL_ERROR "PYSCF_BIN was specified but could not import pyscf." )
  ENDIF()
ENDIF()

######################################################################
# CTest
######################################################################
INCLUDE( "${qmcpack_SOURCE_DIR}/CMake/macros.cmake" )
SET( DROP_METHOD "http" )
SET( DROP_SITE "cdash.qmcpack.org" )
SET( DROP_LOCATION "/CDash/submit.php?project=QMCPACK" )
SET( TRIGGER_SITE "" )
SET( DROP_SITE_CDASH TRUE )
# Increase timeout per test over the default of 1500 seconds (25 minutes)
SET( DART_TESTING_TIMEOUT 3600 CACHE STRING "Maximum time for one test")
ENABLE_TESTING()
INCLUDE( CTest )

SET( QMC_SYMLINK_TEST_FILES TRUE CACHE BOOL "Use symbolic links for test files to save space.  Set to false to copy files instead.")


IF (QMC_SYMLINK_TEST_FILES)
  SET(SYMLINK_MSG "Using symbolic links for large test files may cause test failures if the build is installed on a separate filesystem from the source.  For example, Titan at OLCF.")
ELSE()
  SET(SYMLINK_MSG "Copying large test files uses more disk space than using symbolic links.")
ENDIF()
MESSAGE(STATUS "QMC_SYMLINK_TEST_FILES = ${QMC_SYMLINK_TEST_FILES}.  ${SYMLINK_MSG}")

# The following option disables the extremely slow setup of full system run tests
# This is a developer option
# It offers a more reasonable compile debug loop if CMakeLists.txt files are being changed
# These constom commands execution time far exceeds even the feature and compiler checks for time consumed
# when rerunning cmake.
#
# On leconte: rerun of cmake with
# slow custom commands 27 secs
# without 8 secs
SET( QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS FALSE CACHE BOOL "Disable the slow cmake custom commands for integration tests.")
MARK_AS_ADVANCED(QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)


######################################################################
# Build level
######################################################################

SET(QMC_BUILD_SANDBOX_ONLY OFF CACHE BOOL "Build only applications in Sandbox directory")
IF ( NOT CMAKE_BUILD_TYPE AND NOT CMAKE_TOOLCHAIN_FILE)
  SET( CMAKE_BUILD_TYPE Release )
ENDIF()
STRING( TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE )
MESSAGE( STATUS "CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")

######################################################################
#define the paths for library and executable for out-of-soruce compilation
######################################################################
SET (LIBRARY_OUTPUT_PATH ${qmcpack_BINARY_DIR}/lib CACHE PATH "Single output directory for building all libraries.")
SET (EXECUTABLE_OUTPUT_PATH ${qmcpack_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")

######################################################################
# build include and library paths required for the projects
# For XYZ library, set an environment as
# export XYZ_HOME='root-of-the-library'
# XYZ_HOME/include
# XYZ_HOME/lib
# IF XYZ library is found, XYZ_FOUND is set to true
#
# find_package is used with cmake 2.8 and higher
######################################################################
#if(APPLE)
#  set(QMC_INCLUDE_PATHS ${QMC_INCLUDE_PATHS} /sw/include)
#  set(QMC_LIBRARY_PATHS ${QMC_LIBRARY_PATHS} /sw/lib)
#endif(APPLE)

######################################################################
# Set the compiler-time parameters
# WALKER_MAX_PROPERTIES max number of observables + 12 or so standard
#   props. Things like forces are per particle so you may need this to
#   rather large.  Each property increases walker size by
#   sizeof(FULLPRECREALTYPE)
# OHMMS_DIM =  dimension of the problem
# OHMMS_INDEXTYPE = type of index
# OHMMS_PRECISION  = base precision, float, double etc
# OHMMS_PRECISION_FULL  = full precision, double etc
# QMC_COMPLEX = true if using complex wavefunctions
# QMC_MPI =  enable MPI
# QMC_OMP = enable OMP
######################################################################
SET(WALKER_MAX_PROPERTIES 256 CACHE STRING "Maximum number of properties tracked by walkers")
MARK_AS_ADVANCED(WALKER_MAX_PROPERTIES)
SET(OHMMS_DIM 3 CACHE STRING "Select physical dimension")
SET(OHMMS_INDEXTYPE int)
MESSAGE(STATUS "defining the float point precision")
SET(OHMMS_PRECISION_FULL double)
IF(QMC_CUDA)
  SET(QMC_MIXED_PRECISION 1 CACHE BOOL "Enable/disable mixed precision")
  SET(OHMMS_PRECISION double)
  SET(CUDA_PRECISION_FULL double)
  IF(QMC_MIXED_PRECISION)
    SET(CUDA_PRECISION float)
  ELSE(QMC_MIXED_PRECISION)
    SET(CUDA_PRECISION double)
  ENDIF(QMC_MIXED_PRECISION)
  MESSAGE("   Base precision = ${OHMMS_PRECISION}")
  MESSAGE("   Full precision = ${OHMMS_PRECISION_FULL}")
  MESSAGE("   CUDA base precision = ${CUDA_PRECISION}")
  MESSAGE("   CUDA full precision = ${CUDA_PRECISION_FULL}")
ELSE(QMC_CUDA)
  SET(QMC_MIXED_PRECISION 0 CACHE BOOL "Enable/disable mixed precision")
  IF(QMC_MIXED_PRECISION)
    SET(OHMMS_PRECISION float)
    SET(MIXED_PRECISION 1)
  ELSE(QMC_MIXED_PRECISION)
    SET(OHMMS_PRECISION double)
  ENDIF(QMC_MIXED_PRECISION)
  MESSAGE("   Base precision = ${OHMMS_PRECISION}")
  MESSAGE("   Full precision = ${OHMMS_PRECISION_FULL}")
ENDIF(QMC_CUDA)

# Code coverage
SET(GCOV_SUPPORTED FALSE)
SET(ENABLE_GCOV FALSE CACHE BOOL "Enable code coverage")

######################################################################
# enable MPI and OPNEMP  if possible
######################################################################
SET(QMC_MPI 1 CACHE BOOL "Enable/disable MPI")
SET(QMC_OMP 1 CACHE BOOL "Enable/disable OpenMP")
SET(QMC_COMPLEX 0 CACHE BOOL "Build for complex binary")
SET(PRINT_DEBUG 0 CACHE BOOL "Enable/disable debug printing")
SET(QMC_CUDA 0 CACHE BOOL "Build with GPU support through CUDA")
SET(ENABLE_CUDA 0 CACHE BOOL "Build with the second generation of GPU support through CUDA (production quality for AFQMC, experimental for real space)")
SET(ENABLE_SOA 1 CACHE BOOL "Enable/disable SoA optimization")
SET(ENABLE_OFFLOAD 0 CACHE BOOL "Enable OpenMP offload")

IF(ENABLE_CUDA AND QMC_CUDA)
  MESSAGE(FATAL_ERROR "ENABLE_CUDA=1 and QMC_CUDA=1 can not be set together!")
ENDIF(ENABLE_CUDA AND QMC_CUDA)

SET(QMC_VERBOSE_CONFIGURATION 0 CACHE BOOL "Give details why specific tests are/are not added")

######################################################################
# set debug printout
# DEBUG_PSIBUFFER_ON keep track buffer's locator
######################################################################
SET(DEBUG_PSIBUFFER_ON 0 CACHE BOOL "Build with GPU support through CUDA")

######################################################################
# Standard test
######################################################################
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})


######################################################################
# add macros to compile einspline if necessary
######################################################################
INCLUDE(${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckCXXSymbolExists.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckIncludeFile.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckTypeSize.cmake)

CHECK_FUNCTION_EXISTS(posix_memalign HAVE_POSIX_MEMALIGN)

CHECK_INCLUDE_FILE(dlfcn.h HAVE_DLFCN_H)
CHECK_INCLUDE_FILE(inttypes.h HAVE_INTTYPES_H)
CHECK_INCLUDE_FILE(memory.h HAVE_MEMORY_H)
CHECK_INCLUDE_FILE(sys/stat.h HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILE(sys/time.h HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILE(sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILE(unistd.h HAVE_UNISTD_H)

######################################################################
# Build options
# QMC_BUILD_STATIC build static/dynamic  libraries
# BUILD_QMCTOOLS   build utility programs
# MPIP_PROFILE     profile mpi performance
######################################################################
SET(BUILD_UNIT_TESTS 1 CACHE BOOL "Build unit tests")
SET(BUILD_LMYENGINE_INTERFACE 1 CACHE BOOL "Build LMY engine")
IF (QMC_CUDA AND BUILD_LMYENGINE_INTERFACE)
  MESSAGE(STATUS "LMY engine is not compatible with CUDA build! Disabling LMY engine")
  SET(BUILD_LMYENGINE_INTERFACE 0)
ENDIF()
IF(MIXED_PRECISION AND BUILD_LMYENGINE_INTERFACE)
  MESSAGE(STATUS "LMY engine is not compatible with CPU mixed precision build! Disabling LMY engine")
  SET(BUILD_LMYENGINE_INTERFACE 0)
ENDIF()
SET(BUILD_AFQMC 0 CACHE BOOL "Build with AFQMC")
SET(BUILD_AFQMC_WITH_NCCL 0 CACHE BOOL "Build AFQMC with NCCL library.")
# AFQMC requires MPI.
If (BUILD_AFQMC AND NOT QMC_MPI)
  MESSAGE(FATAL_ERROR "AFQMC requires building with MPI (QMC_MPI=1). Set BUILD_AFQMC=0 or configure MPI.")
ENDIF()
SET(BUILD_FCIQMC 0 CACHE BOOL "Build with FCIQMC")
#SET(BUILD_QMCTOOLS 1 CACHE BOOL "Build tools for QMCPACK")
#SET(MPIP_PROFILE 0 CACHE BOOL "Build with mpip for mpi profile")
#SET(BUILD_QMCTOOLS 1)
#SET(MPIP_PROFILE 0)
SET(QMC_BUILD_STATIC 0 CACHE BOOL "Link to static libraries")
SET(ENABLE_TIMERS 0 CACHE BOOL "Enable internal timers")

######################################################################
# Install options
# INSTALL_NEXUS install Nexus alongside QMCPACK
######################################################################
SET(INSTALL_NEXUS 1 CACHE BOOL "Install Nexus alongside QMCPACK")

######################################################################
#list of options to determine libraries. First, try to use libraries available
#FindXYZ.cmake will set the values if successful
######################################################################
SET(HAVE_EINSPLINE 1) # to be removed
SET(HAVE_LIBXML2 1)
SET(HAVE_LIBFFTW 0)
SET(HAVE_LIBXMLPP 0)
SET(HAVE_LIBSPRNG 0)
SET(HAVE_LIBBLITZ 0)
SET(HAVE_LIBGSL 0)
SET(MAC_VECLIB 0)

######################################################
# Git information (see src/CMakeLists.txt for the
# command to get repository information)
######################################################
FIND_PACKAGE(Git)

IF (GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  SET(IS_GIT_PROJECT 1)
ENDIF()


######################################################
#COMPILER choose one of the cmake files to customize the compiler options
#If nothing is chosen, default settings by cmake will  be used.
#If the automatic detection does not work, comment out everything
#upto COMPILER MANUAL SELECTION and use one of the customized cmake file.
######################################################
INCLUDE(${CMAKE_ROOT}/Modules/CheckCCompilerFlag.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckCXXCompilerFlag.cmake)

######################################################
# DO NOT MODIFY BELOW
# If you like to customize your build, create a toolchain as instructed
######################################################
######################################################
# now set the compiler flags and the macros in src/config.h
######################################################
IF(CMAKE_TOOLCHAIN_FILE)
  MESSAGE(STATUS "Using ${CMAKE_TOOLCHAIN_FILE} toolchain ")
ELSE(CMAKE_TOOLCHAIN_FILE)
  MESSAGE(STATUS "Trying to figure out .....")
  SET(ENABLE_OPENMP 0)
  SET(HAVE_MPI 0)
  SET(INTEL_COMPILER 0)
  SET(HAVE_MKL 0)
  SET(HAVE_MKL_VML 0)
  SET(HAVE_ACML 0)

  ######################################################
  # Try to identify CPU identity
  ######################################################
  SET(CPU_IDENTITY "generic")
  #INCLUDE(${PROJECT_CMAKE}/CheckProcessorID.cmake)

  #------------------------
  # On Cray's machine
  #------------------------
  IF($ENV{CRAYPE_VERSION} MATCHES ".")
    MESSAGE(STATUS "Running on a Cray machine.")
    SET(CMAKE_SKIP_RPATH TRUE)

    IF(NOT $ENV{CRAYPE_LINK_TYPE} STREQUAL "dynamic")
      MESSAGE(WARNING "Cray Programming Environment uses static linking by default. "
                      "We prefer dynamic linking which makes library searching easier. "
                      "Set environment variable CRAYPE_LINK_TYPE to \"dynamic\" enables dynamic linking.\n")
    ENDIF()

    IF ( NOT DEFINED MPIEXEC_EXECUTABLE )
      FIND_PROGRAM(MPIEXEC_EXECUTABLE aprun)
      IF(NOT MPIEXEC_EXECUTABLE)
        # For slurm use srun
        FIND_PROGRAM(MPIEXEC_EXECUTABLE srun)
      ENDIF(NOT MPIEXEC_EXECUTABLE)
    ENDIF()

    IF(QMC_CUDA)
      IF(QMC_BUILD_STATIC)
        MESSAGE(FATAL_ERROR "Static linking is not supported with CUDA build. Set QMC_BUILD_STATIC=0 (default)!")
      ENDIF(QMC_BUILD_STATIC)
      MESSAGE(STATUS "Building GPU code. Linking dynamically.")
    ELSE(QMC_CUDA)
      IF(QMC_BUILD_STATIC)
        MESSAGE(STATUS "Building CPU code. Linking statically.")
        SET(BUILD_SHARED_LIBS off)
        SET(CMAKE_FIND_LIBRARY_SUFFIXES .a)
        FOREACH(type SHARED_LIBRARY SHARED_MODULE EXE)
          SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-Wl,-Bstatic")
          SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-static")
          SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-Wl,-Bstatic")
          SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-static")
        ENDFOREACH(type)
      ELSE(QMC_BUILD_STATIC)
        MESSAGE(STATUS "Building CPU code. Linking dynamically.")
      ENDIF(QMC_BUILD_STATIC)
    ENDIF(QMC_CUDA)
  ENDIF($ENV{CRAYPE_VERSION} MATCHES ".")

  #--------------------------------
  #  C++ 14 is the minimum standard
  #--------------------------------
  set(CMAKE_CXX_STANDARD 14)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)

  #------------------------------------
  # Identify the compiler -- This serves only to deal with case where only C or CXX is set, refactor? 
  #------------------------------------
  IF ( CMAKE_C_COMPILER_WORKS OR CMAKE_C_COMPILER_WORKS )
    IF( (CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR ( CMAKE_C_COMPILER_ID MATCHES "GNU") )
      SET( COMPILER GNU )
    ELSEIF( (CMAKE_CXX_COMPILER_ID MATCHES "XL") OR ( CMAKE_C_COMPILER_ID MATCHES "XL") )
      SET( COMPILER IBM )
    ELSEIF( (CMAKE_CXX_COMPILER_ID MATCHES "Intel") OR (CMAKE_C_COMPILER_ID MATCHES "Intel") )
      SET( COMPILER Intel )
    ELSEIF( (CMAKE_CXX_COMPILER_ID MATCHES "PGI") OR (CMAKE_C_COMPILER_ID MATCHES "PGI") )
      SET( COMPILER PGI )
    ELSEIF( (CMAKE_CXX_COMPILER_ID MATCHES "Cray") OR (CMAKE_C_COMPILER_ID MATCHES "Cray") )
      SET( COMPILER Cray )
    ELSEIF( (CMAKE_CXX_COMPILER_ID MATCHES "Clang") OR (CMAKE_C_COMPILER_ID MATCHES "Clang") )
      SET( COMPILER Clang )
    ELSE()
        SET(USING_DEFAULT TRUE)
        MESSAGE("${CMAKE_C_COMPILER_ID}")
        MESSAGE(WARNING "Unknown C/C++ compiler, default flags will be used")
    ENDIF()
  ELSE()
    MESSAGE(WARNING "No compiler identified")
  ENDIF()
  MESSAGE(STATUS "Compiler: ${COMPILER}")

  #------------------------------------
  # Include compiler-specific cmake file
  #------------------------------------
  IF( ${COMPILER} MATCHES "IBM" )
    INCLUDE(${PROJECT_CMAKE}/IBMCompilers.cmake)
  ELSEIF( ${COMPILER} MATCHES "Intel" )
    INCLUDE(${PROJECT_CMAKE}/IntelCompilers.cmake)
  ELSEIF( ${COMPILER} MATCHES "GNU" )
    INCLUDE(${PROJECT_CMAKE}/GNUCompilers.cmake)
  ELSEIF( ${COMPILER} MATCHES "Clang" )
    INCLUDE(${PROJECT_CMAKE}/ClangCompilers.cmake)
  ELSEIF( ${COMPILER} MATCHES "PGI" )
    INCLUDE(${PROJECT_CMAKE}/PGICompilers.cmake)
  ELSE()
    MESSAGE(WARNING "No default file for compiler (${COMPILER})")
  ENDIF()


  #--------------------------------------------------------------------------
  # Check that a C++ compiler is compatiable with the underlying libstdc++
  #--------------------------------------------------------------------------
  include(CMake/Testlibstdc++.cmake)

  #---------------------------------------------------
  # Check that a C++ 14 standard library is configured
  #---------------------------------------------------
  include(CMake/TestCxx14Library.cmake)

  #-------------------------------------------------------------------
  # check lapack/blas
  #-------------------------------------------------------------------
  # Quick search for MKL (this will need to be improved).
  # This reproduces the old intel -mkl behavior or
  # if ENABLE_MKL is set to 1 supports other compilers
  IF (INTEL_COMPILER)
    SET(ENABLE_MKL TRUE CACHE BOOL "ENABLE MKL for Intel compilers")
  ELSE (INTEL_COMPILER)
    SET(ENABLE_MKL FALSE CACHE BOOL "ENABLE MKL for non-intel compilers (Experimental)")
  ENDIF (INTEL_COMPILER)
  IF (ENABLE_MKL)
   INCLUDE(CMake/FindMKL.cmake)
  ENDIF(ENABLE_MKL)
  IF ( HAVE_MKL )
    # Skip adding Math::BLAS_LAPACK target which has been imported by FindMKL.
  ELSEIF ( $ENV{CRAYPE_VERSION} MATCHES "." )
    # Cray specific environment variable indicates we are on a cray. Assume libsci will be linked
    # preferable to attempt a compile/link
    MESSAGE(STATUS "Assuming Cray libsci or other BLAS/LAPACK will be linked by compiler wrappers")
    ADD_LIBRARY(Math::BLAS_LAPACK INTERFACE IMPORTED)
  ELSE()
    # The last search for Blas/Lapack
    find_package(LAPACK REQUIRED)
    IF("${LAPACK_LIBRARIES}" MATCHES "mkl")
      MESSAGE(FATAL_ERROR "MKL found via LAPACK/BLAS. Please pass -DENABLE_MKL=1 to cmake for full Intel MKL utilization. -DMKL_ROOT may also be needed for the include files to be found.")
    ENDIF()
    MESSAGE(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")
    MESSAGE(STATUS "LAPACK linker flags: ${LAPACK_LINKER_FLAGS}")
    SET(LAPACK_LIBRARY ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
    # Add Math::BLAS_LAPACK target
    ADD_LIBRARY(Math::BLAS_LAPACK INTERFACE IMPORTED)
    SET_TARGET_PROPERTIES(Math::BLAS_LAPACK PROPERTIES INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARY}")
  ENDIF()
ENDIF(CMAKE_TOOLCHAIN_FILE)

IF (ENABLE_GCOV)
  IF (NOT GCOV_SUPPORTED)
    MESSAGE(FATAL_ERROR "ENABLE_GCOV was specified but compiler does not support GCC coverage flag")
  ENDIF()
ENDIF(ENABLE_GCOV)

# AFQMC requires MKL sparse for good performance (roughly a factor of 2x)
IF (BUILD_AFQMC AND NOT MKL_FOUND)
  MESSAGE(WARNING "AFQMC - MKL not found, using simple sparse matrix routines.  Link with MKL sparse libraries for better performance.")
ENDIF()

IF (BUILD_AFQMC AND NOT APPLE)
  LINK_LIBRARIES("rt")
ENDIF()

#-------------------------------------------------------------------
#  check MPI installation
#-------------------------------------------------------------------
IF(QMC_MPI)
  SET(MPI_CXX_SKIP_MPICXX TRUE)
  SET(MPI_DETERMINE_LIBRARY_VERSION TRUE)
  FIND_PACKAGE(MPI COMPONENTS CXX)

  IF(NOT MPI_FOUND)
    MESSAGE(FATAL_ERROR "MPI support not found! Provide MPI compiler wrappers or build without MPI by passing '-DQMC_MPI=0' to cmake.")
  ENDIF(NOT MPI_FOUND)

  IF(${MPI_CXX_LIBRARY_VERSION_STRING} MATCHES "MVAPICH2" AND NOT ${MPI_CXX_LIBRARY_VERSION_STRING} MATCHES "disable-registration-cache")
    MESSAGE(FATAL_ERROR "MVAPICH2 with registration cache enabled breaks QMCPACK. "
                        "Use a different MPI library or build MVAPICH2 with --disable-registration-cache configure option.")
  ENDIF()

  IF(NOT ${CMAKE_CXX_COMPILER} STREQUAL ${MPI_CXX_COMPILER})
    SET(MPI_WARNING_LIST
        "Building MPI version without using MPI compiler wrappers.\n"
        "This may not build qmcpack correctly. To ensure the correct version, specify the compiler wrappers to cmake.\n"
        "For example: cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpic++\n"
        "To build without MPI, pass '-DQMC_MPI=0' to cmake")
    MESSAGE(WARNING ${MPI_WARNING_LIST})
  ENDIF()

  IF ( NOT TEST_MAX_PROCS )
    SET( TEST_MAX_PROCS 100 )
  ENDIF()

  MESSAGE(STATUS "MPI runner MPIEXEC_EXECUTABLE : ${MPIEXEC_EXECUTABLE}")
  MESSAGE(STATUS "MPIEXEC_NUMPROC_FLAG : ${MPIEXEC_NUMPROC_FLAG}")
  MESSAGE(STATUS "MPIEXEC_PREFLAGS : ${MPIEXEC_PREFLAGS}")
  MESSAGE(STATUS "Tests run as : ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} NUM_PROCS ${MPIEXEC_PREFLAGS} EXECUTABLE")

  SET(HAVE_MPI 1)
  MESSAGE(STATUS "MPI is enabled")
ELSE(QMC_MPI)
  SET(HAVE_MPI 0)
  MESSAGE(STATUS "MPI is disabled")
ENDIF(QMC_MPI)

####################################################################
#First check the required libraries. Abort if these are not found.
#Changed the fall-back versions to use similar conventions to 2.8.0 and higher
####################################################################
# use customized cmake file for boost so that we can deal with multiple
# installataions
# list of libraries to be linked with the main applications using I/O

#-------------------------------------------------------------------
#  set up scalar/vector math libraries
#-------------------------------------------------------------------
SET(ENABLE_MASS FALSE CACHE BOOL "ENABLE MASS scalar and vector math libraries for Power architecture")

# This needs to go before HAVE_SINCOS
IF(ENABLE_MASS)
  INCLUDE(CMake/FindIBMMASS.cmake)
ENDIF(ENABLE_MASS)

IF(NOT HAVE_MASS AND NOT HAVE_MKL_VML)
  MESSAGE(STATUS "No usable vector math library detected.")
  ADD_LIBRARY(Math::scalar_vector_functions INTERFACE IMPORTED)
ENDIF(NOT HAVE_MASS AND NOT HAVE_MKL_VML)

# CheckSincos relies on SINCOS_INCLUDE which may be modified based on vendor libraries
IF(NOT SINCOS_INCLUDE)
  SET(SINCOS_INCLUDE cmath)
ENDIF(NOT SINCOS_INCLUDE)
MESSAGE(STATUS "SINCOS_INCLUDE : ${SINCOS_INCLUDE}")
INCLUDE(${PROJECT_CMAKE}/CheckSincos.cmake)

#-------------------------------------------------------------------
#  set up FFTW3 library
#-------------------------------------------------------------------
SET( FFTW_FOUND 0 )
IF ( HAVE_MKL )
  IF ( HAVE_MKL_FFTW3 )
    SET( FFTW_FOUND 1 )
  ENDIF()
ENDIF()
IF ( NOT FFTW_FOUND )
  find_package(FFTW)
ENDIF()
IF(FFTW_FOUND)
  SET(HAVE_LIBFFTW 1)
  MESSAGE(STATUS "Found FFTW")
ELSE(FFTW_FOUND)
#  MESSAGE(STATUS "Disabling FFTW")
  MESSAGE(FATAL_ERROR "Could not find required library FFTW. Verify FFTW_HOME")
ENDIF(FFTW_FOUND)

SET(QMC_UTIL_LIBS ${QMC_UTIL_LIBS} ${FORTRAN_LIBRARIES})

find_package(ZLIB)
find_package(LibXml2)

IF(LIBXML2_FOUND)
  INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR})
  SET(HAVE_LIBXML2 1)
  SET(QMC_UTIL_LIBS ${QMC_UTIL_LIBS} ${LIBXML2_LIBRARIES})
ELSE()
  MESSAGE(FATAL_ERROR "Could not find required library Libxml2")
ENDIF(LIBXML2_FOUND)

#-------------------------------------------------------------------
#  set up HDF5 library
#-------------------------------------------------------------------
IF(HAVE_MPI)
  SET(HDF5_PREFER_PARALLEL 1 CACHE BOOL "Request parallel/serial HDF5 library")
ELSE(HAVE_MPI)
  SET(HDF5_PREFER_PARALLEL 0 CACHE BOOL "Request parallel/serial HDF5 library")
  IF(HDF5_PREFER_PARALLEL)
    MESSAGE(FATAL_ERROR "Parallel HDF5 library cannot be selected with QMCPACK non-MPI build. "
                        "Please set HDF5_PREFER_PARALLEL=0.")
  ENDIF(HDF5_PREFER_PARALLEL)
ENDIF(HAVE_MPI)

if(QMC_BUILD_STATIC)
  MESSAGE(STATUS "Linking static HDF5 library")
  set(HDF5_USE_STATIC_LIBRARIES on)
else()
  MESSAGE(STATUS "Linking dynamic HDF5 library")
  set(HDF5_USE_STATIC_LIBRARIES off)
endif()

find_package(HDF5 COMPONENTS C)

IF(HDF5_FOUND)
  SET(HAVE_LIBHDF5 1)

  IF(HDF5_IS_PARALLEL)
    MESSAGE(STATUS "Parallel HDF5 library found")
    SET(ENABLE_PHDF5 1 CACHE BOOL "Enable code paths using parallel HDF5")
  ELSE(HDF5_IS_PARALLEL)
    MESSAGE(STATUS "Serial HDF5 library found")
    SET(ENABLE_PHDF5 0 CACHE BOOL "Enable code paths using parallel HDF5")
    IF(ENABLE_PHDF5)
      IF(HAVE_MPI)
        MESSAGE(FATAL_ERROR "Parallel HDF5 code paths requested but serial HDF5 library found! "
                            "Please either provide parallel HDF5 library or set ENABLE_PHDF5=0.")
      ELSE(HAVE_MPI)
        MESSAGE(FATAL_ERROR "Parallel HDF5 code paths cannot be enabled on non-MPI builds! Please set ENABLE_PHDF5=0.")
      ENDIF(HAVE_MPI)
    ENDIF(ENABLE_PHDF5)
  ENDIF(HDF5_IS_PARALLEL)

  IF(ENABLE_PHDF5)
    MESSAGE(STATUS "Using HDF5 parallel collective I/O code paths")
  ELSE(ENABLE_PHDF5)
    MESSAGE(STATUS "Using HDF5 non-scalable serial I/O code paths")
  ENDIF(ENABLE_PHDF5)

  IF(HAVE_MPI AND NOT ENABLE_PHDF5)
    MESSAGE(WARNING "MPI builds may have potential performance loss not using parallel HDF5!")
  ENDIF()

  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
  IF ( CMAKE_BUILD_TYPE AND HDF5_LIBRARIES_DEBUG )
    IF ( CMAKE_BUILD_TYPE MATCHES DEBUG )
      SET( HDF5_LIBRARIES ${HDF5_LIBRARIES_DEBUG} )
    ELSE()
      SET( HDF5_LIBRARIES ${HDF5_LIBRARIES_RELEASE} )
    ENDIF()
  ENDIF()
  SET(QMC_UTIL_LIBS ${QMC_UTIL_LIBS} ${HDF5_LIBRARIES})
ELSE(HDF5_FOUND)
  MESSAGE(FATAL_ERROR "HDF5 not found. Set HDF5_ROOT")
ENDIF(HDF5_FOUND)

#make sure we can find boost if it's not in /usr
set(Boost_NO_BOOST_CMAKE on)
find_package(Boost 1.61.0 REQUIRED)
IF(Boost_FOUND)
  SET(HAVE_LIBBOOST 1)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  MESSAGE(STATUS "Setting Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIRS}")
ELSE()
  MESSAGE(FATAL_ERROR "Need boost 1.61.0 or higher.  Set BOOST_ROOT")
ENDIF()

SET(HAVE_CUDA 0)
IF(QMC_CUDA OR ENABLE_CUDA)
  IF(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CUDA_PROPAGATE_HOST_FLAGS ON)
  ELSE(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  ENDIF(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  FIND_PACKAGE(CUDA REQUIRED)
  SET(CUDA_LINK_LIBRARIES_KEYWORD PRIVATE)
  #set(CUDA_NVCC_FLAGS
  #  "-arch=sm_20;-Drestrict=__restrict__;-DNO_CUDA_MAIN;-O3;-use_fast_math")
  if(CUDA_NVCC_FLAGS MATCHES "arch")
    # User defined NVCC flags
    MESSAGE(STATUS "Setting CUDA FLAGS=${CUDA_NVCC_FLAGS}")
  else(CUDA_NVCC_FLAGS MATCHES "arch")
    # Automatically set the default NVCC flags
    SET(CUDA_NVCC_FLAGS "-Drestrict=__restrict__;-DNO_CUDA_MAIN")
    if (QMC_COMPLEX)
      SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-DQMC_COMPLEX")
    endif()
    if ( CMAKE_BUILD_TYPE STREQUAL "DEBUG" )
      SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-g;-G")
    else()
      # Temporarily disable fast_math because it causes multiple test failures
      # SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-O3;-use_fast_math")
      SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-O3")
    ENDIF()
    SET( CUDA_ARCH sm_35 CACHE STRING "CUDA architecture sm_XX")
    SET(CUDA_NVCC_FLAGS "-arch=${CUDA_ARCH};${CUDA_NVCC_FLAGS}")
  endif(CUDA_NVCC_FLAGS MATCHES "arch")
  INCLUDE_DIRECTORIES(${CUDA_INCLUDE_DIRS})
  set(HAVE_CUDA 1)
  MESSAGE("   CUDA_NVCC_FLAGS=${CUDA_NVCC_FLAGS}")
ELSE(QMC_CUDA OR ENABLE_CUDA)
  MESSAGE(STATUS "Disabling CUDA")
ENDIF(QMC_CUDA OR ENABLE_CUDA)

SET(USE_NVTX_API 0 CACHE BOOL "Enable/disable NVTX regions in CUDA code.")
IF(USE_NVTX_API)
  IF(HAVE_CUDA)
    FIND_LIBRARY(NVTX_API_LIB
      NAME nvToolsExt
      HINTS ${CUDA_TOOLKIT_ROOT_DIR}
      PATH_SUFFIXES lib lib64)
    IF(NOT NVTX_API_LIB)
      MESSAGE(FATAL_ERROR "USE_NVTX_API set but NVTX_API_LIB not found")
    ENDIF(NOT NVTX_API_LIB)
    MESSAGE("CUDA nvToolsExt library: ${NVTX_API_LIB}")
    LINK_LIBRARIES(${NVTX_API_LIB})
  ENDIF(HAVE_CUDA)
ENDIF(USE_NVTX_API)

######################################################
# Add optional macros
######################################################

#include qmcpack/src build/src
INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/src ${PROJECT_BINARY_DIR}/src)

IF (USE_VTUNE_TASKS)
  IF (NOT ENABLE_TIMERS)
    MESSAGE(FATAL_ERROR "USE_VTUNE_TASKS is set, but timers are not enabled.  Set ENABLE_TIMERS=1.")
  ENDIF()
  SET(USE_VTUNE_API 1)
ENDIF()

IF (USE_VTUNE_API)
  include(CheckIncludeFileCXX)
  CHECK_INCLUDE_FILE_CXX(ittnotify.h HAVE_ITTNOTIFY_H)
  IF (NOT HAVE_ITTNOTIFY_H)
    MESSAGE(FATAL_ERROR "USE_VTUNE_API is defined, but the ittnotify.h include file is not found.  Check that the correct include directory is present in CMAKE_CXX_FLAGS.")
  ENDIF()

  FIND_LIBRARY(VTUNE_ITTNOTIFY_LIBRARY ittnotify)
  IF (NOT VTUNE_ITTNOTIFY_LIBRARY)
    MESSAGE(FATAL_ERROR "USE_VTUNE_API is defined, but the ittnotify library is not found.  Check that correct library path is present in CMAKE_LIBRARY_PATH.")
  ENDIF()
  LINK_LIBRARIES("${VTUNE_ITTNOTIFY_LIBRARY}")
ENDIF()

SET(QMC_EXP_THREADING FALSE CACHE BOOL "Experimental non openmp threading models")
MARK_AS_ADVANCED(QMC_EXP_THREADING)
IF(QMC_EXP_THREADING)
  include(CheckAtomic)
  IF(HAVE_LIBATOMIC)
    LINK_LIBRARIES("${HAVE_LIBATOMIC}")
  ENDIF(HAVE_LIBATOMIC)
  ADD_DEFINITIONS(-DQMC_EXP_THREADING)
ENDIF(QMC_EXP_THREADING)

#include(ExternalProject)
#  set(einspline_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/einspline")
#  set(einspline_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/einspline")
#  set(einspline_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${einspline_INSTALL_DIR})
#ExternalProject_Add(
#   einspline
#   SVN_REPOSITORY http://qmctools.googlecode.com/svn/ext-packages/einspline
#   )
#    INSTALL_DIR ${einspline_INSTALL_DIR}
#    CMAKE_ARGS ${einspline_CMAKE_ARGS}
#    )
#  find_package(einspline REQUIRED)
#  set(HAVE_EINSPLINE 1)


###############################################################
# for specific linker in the case of profiling like hpctoolkit.
###############################################################
if(CMAKE_CXX_LINKER)
  set(CMAKE_CXX_LINK_EXECUTABLE
  "${CMAKE_CXX_LINKER} <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS>  <OBJECTS> -o <TARGET> <LINK_LIBRARIES>")
  MESSAGE(STATUS "Linker is specified as '${CMAKE_CXX_LINKER}' instead of '${CMAKE_CXX_COMPILER}'")
endif()

###################################################
# Add user extra flags
###################################################
IF(QMC_INCLUDE)
  FOREACH ( tmp ${QMC_INCLUDE} )
    MESSAGE(STATUS "Adding '${tmp}' to include directories")
    INCLUDE_DIRECTORIES( "${tmp}" )
  ENDFOREACH()
ENDIF()
IF(QMC_EXTRA_LIBS)
  string(REPLACE " " ";" QMC_EXTRA_LIBS ${QMC_EXTRA_LIBS})
  SET(QMC_UTIL_LIBS ${QMC_UTIL_LIBS} ${QMC_EXTRA_LIBS})
  MESSAGE("QMC_UTIL_LIBS=${QMC_UTIL_LIBS}")
ENDIF()

# Avoid unused variables warnings (for now) by using for information purposes
IF(QMC_DATA)
MESSAGE(STATUS "QMC_DATA is ${QMC_DATA}")
ENDIF()
IF(ENABLE_GCOV)
MESSAGE(STATUS "GCOV is enabled")
ENDIF()

############################################################
# FLAGS at the project level
############################################################
#this might be redundant but maybe not in all CMake versions.
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" THIS_CONFIG)

FOREACH(lang IN ITEMS C CXX)
  SET(PROJECT_CMAKE_${lang}_FLAGS ${CMAKE_${lang}_FLAGS})
  #pre 3.0 cmake does not have string CONCAT
  SET(TMP_PROJECT_CMAKE_${lang}_FLAGS "${PROJECT_CMAKE_${lang}_FLAGS} ${CMAKE_${lang}_FLAGS_${THIS_CONFIG}}")
  MESSAGE("Project ${lang}_FLAGS: ${TMP_PROJECT_CMAKE_${lang}_FLAGS}")
ENDFOREACH(lang IN ITEMS C CXX)

get_directory_property(TMP_PROJECT_INCLUDE_DIRECTORIES INCLUDE_DIRECTORIES)
MESSAGE("Project INCLUDE_DIRECTORIES: ${TMP_PROJECT_INCLUDE_DIRECTORIES}")

FOREACH(target_type IN ITEMS EXE SHARED)
  SET(PROJECT_CMAKE_${target_type}_LINKER_FLAGS ${CMAKE_${target_type}_LINKER_FLAGS})
  SET(TMP_PROJECT_CMAKE_${target_type}_LINKER_FLAGS "${PROJECT_CMAKE_${target_type}_LINKER_FLAGS} ${CMAKE_${target_type}_LINKER_FLAGS_${THIS_CONFIG}}")
  MESSAGE("Project ${target_type}_LINKER_FLAGS: ${TMP_PROJECT_CMAKE_${target_type}_LINKER_FLAGS}")
ENDFOREACH(target_type IN ITEMS EXE SHARED)


###################################################
# Done searching libraries. Start building qmcpack
###################################################

MESSAGE(STATUS "Ready to build qmcpack")

####################################################################
# add definitions
#  -DH5_USE_16_API for HDF5
####################################################################

ADD_DEFINITIONS(-DHAVE_CONFIG_H -DH5_USE_16_API -DUSE_REAL_STRUCT_FACTOR -DINLINE_ALL=inline -DADD_)

EXECUTE_PROCESS(COMMAND hostname OUTPUT_VARIABLE QMC_HOSTNAME OUTPUT_STRIP_TRAILING_WHITESPACE)
CONFIGURE_FILE(${qmcpack_SOURCE_DIR}/src/qmcpack.settings ${qmcpack_BINARY_DIR}/bin/qmcpack.settings)
INSTALL(FILES ${qmcpack_BINARY_DIR}/bin/qmcpack.settings DESTINATION bin)
IF(INSTALL_NEXUS)
  INSTALL(CODE "EXECUTE_PROCESS(COMMAND ${qmcpack_SOURCE_DIR}/nexus/install --leave_paths ${CMAKE_INSTALL_PREFIX}/bin)")
ENDIF()

SUBDIRS(src)

IF(NOT QMC_BUILD_SANDBOX_ONLY)
  IF( NOT QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS )
    SUBDIRS(tests)
    SUBDIRS(examples)
  ENDIF()
  IF(IS_DIRECTORY "${PROJECT_SOURCE_DIR}/nexus" AND INSTALL_NEXUS)
    SUBDIRS(nexus/tests)
  ENDIF()
ENDIF()

#ADD_CUSTOM_TARGET(print_settings ALL
#      ${CMAKE_COMMAND} -E echo "Print build settings"
#          VERBATIM)
#
#add_custom_command(TARGET print_settings
#  POST_BUILD
#  COMMAND svn info ${qmcpack_SOURCE_DIR})
#
#file(APPEND ${qmcpack_BINARY_DIR}/bin/qmcpack.settings ${print_settings})
