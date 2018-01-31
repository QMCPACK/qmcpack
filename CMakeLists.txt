PROJECT(qmcpack)


#####################################################
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.10)
IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)


######################################################################
# Version 
######################################################################
SET(QMCPACK_VERSION_MAJOR 3)
SET(QMCPACK_VERSION_MINOR 4)
SET(QMCPACK_VERSION_PATCH 0)
SET(QMCPACK_VERSION "${QMCPACK_VERSION_MAJOR}.${QMCPACK_VERSION_MINOR}.${QMCPACK_VERSION_PATCH}")

######################################################################
# CTest 
######################################################################
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

######################################################################
# Build level
######################################################################
SET(QMC_BUILD_LEVEL 3 CACHE INTEGER 
  "QMC Build Level: 1=bare, 2=developer, 3=experimental, 4=minimal")
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
# OHMMS_DIM =  dimension of the problem
# OHMMS_INDEXTYPE = type of index
# OHMMS_PRECISION  = base precision, float, double etc
# OHMMS_PRECISION_FULL  = full precision, double etc
# QMC_COMPLEX = true if using complex wavefunctions
# QMC_MPI =  enable MPI 
# QMC_OMP = enable OMP 
######################################################################
SET(OHMMS_DIM 3 CACHE INTEGER "Select physical dimension")
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
SET(QMC_COMPLEX 0 CACHE INTEGER "Build for complex binary")
SET(PRINT_DEBUG 0 CACHE BOOL "Enable/disable debug printing")
SET(QMC_CUDA 0 CACHE BOOL "Build with GPU support through CUDA")

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
INCLUDE(${CMAKE_ROOT}/Modules/CheckIncludeFile.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckTypeSize.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckCXXSourceCompiles.cmake)

CHECK_FUNCTION_EXISTS(sincos HAVE_SINCOS)
CHECK_FUNCTION_EXISTS(posix_memalign HAVE_POSIX_MEMALIGN)
SET(CHECK_ISZERO_SOURCE "#include <math.h>
                         int main() {iszero(1.0); return 0;}")
CHECK_CXX_SOURCE_COMPILES("${CHECK_ISZERO_SOURCE}" HAVE_ISZERO)

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
# BUILD_SANDBOX    build test programs
# BUILD_MINIAPPS   build miniapps (Incomplete)
# MPIP_PROFILE     profile mpi performance
######################################################################
SET(QMC_ADIOS 0 CACHE BOOL "Build with ADIOS")
SET(BUILD_UNIT_TESTS 1 CACHE BOOL "Build unit tests")
SET(BUILD_MINIAPPS 1 CACHE BOOL "Build miniapps")
SET(BUILD_LMYENGINE_INTERFACE 1 CACHE BOOL "Build LMY engine")
IF (QMC_CUDA AND BUILD_LMYENGINE_INTERFACE)
  MESSAGE(STATUS "LMY engine is not compatiable with CUDA build! Disabling LMY engine")
  SET(BUILD_LMYENGINE_INTERFACE 0)
ENDIF()
IF(MIXED_PRECISION AND BUILD_LMYENGINE_INTERFACE)
  MESSAGE(STATUS "LMY engine is not compatiable with CPU mixed precision build! Disabling LMY engine")
  SET(BUILD_LMYENGINE_INTERFACE 0)
ENDIF()
SET(BUILD_AFQMC 0 CACHE BOOL "Build with AFQMC")
SET(BUILD_FCIQMC 0 CACHE BOOL "Build with FCIQMC")
#SET(BUILD_QMCTOOLS 1 CACHE BOOL "Build tools for QMCPACK")
#SET(BUILD_SANDBOX 0 CACHE BOOL "Build snadbox for testing")
#SET(MPIP_PROFILE 0 CACHE BOOL "Build with mpip for mpi profile")
#SET(BUILD_QMCTOOLS 1)
#SET(BUILD_SANDBOX 0)
#SET(MPIP_PROFILE 0)
SET(QMC_BUILD_STATIC 0 CACHE BOOL "Link to static libraries")
SET(ENABLE_TIMERS 0 CACHE BOOL "Enable internal timers")

######################################################################
# Performance-related macros
# QMC_SK_USE_RECURSIVE enable/disable recursive evalaution of SK
######################################################################
SET(QMC_SK_USE_RECURSIVE 0 CACHE BOOL "Use recursive evalaution of SK")
IF(QMC_SK_USE_RECURSIVE)
  MESSAGE(STATUS "SK structure factor uses a recursive algorithm.")
ENDIF(QMC_SK_USE_RECURSIVE)

######################################################################
# FIXED PARAMETERS for test and legacy reasons
# DO NOT TOUCH THESE
# QMC_ORTHO  = 0 for general unit cell, 1 for orthogonal cells
# QMC_FASTWALKER  = true, when operation is faster for each walker
######################################################################
SET(QMC_FASTWALKER 1)
SET(OHMMS_ORTHO 0)

######################################################################
#list of options to determine libraries. First, try to use libraries available
#FindXYZ.cmake will set the values if successful
######################################################################
SET(ENABLE_PHDF5 1 CACHE BOOL "Enable phdf5 for output")
SET(HAVE_LIBXML2 1)
SET(HAVE_LIBFFTW 0)
SET(HAVE_LIBXMLPP 0)
SET(HAVE_LIBHDF5 0)
SET(HAVE_LIBSPRNG 0)
SET(HAVE_LIBBLITZ 0)
SET(HAVE_LIBGSL 0)
SET(MAC_VECLIB 0)
SET(ENABLE_GTK_GUI 0)
set(ENABLE_TAU_PROFILE 0 CACHE BOOL "Enable tau for profiling")
set(TAU_HOME $ENV{TAU_HOME})

######################################################
# Directory where customize cmake files reside
######################################################
SET (PROJECT_CMAKE ${qmcpack_SOURCE_DIR}/CMake)

######################################################
# set ADIOS flags 
######################################################
IF(QMC_ADIOS)
  IF(NOT HAVE_ADIOS)
    SET (ADIOS_FOUND 0)
    INCLUDE(${PROJECT_CMAKE}/FindADIOS.cmake)
    IF (${ADIOS_FOUND} EQUAL 1)
      INCLUDE_DIRECTORIES(${ADIOS_INCLUDES})
      MESSAGE(STATUS "Found ADIOS")
      LINK_LIBRARIES(${ADIOS_LIBRARIES})
      ADD_DEFINITIONS(-DHAVE_ADIOS)
    ELSE()
      MESSAGE(STATUS "Couldn't find required library ADIOS")
      MESSAGE(FATAL_ERROR "ADIOS_HOME: export ADIOS_HOME=/path/to/adios")
    ENDIF()
  ELSE()
    ADD_DEFINITIONS(-DHAVE_ADIOS)
  ENDIF()
ENDIF(QMC_ADIOS)

SET(QMCPACK_UNIT_TEST_DIR ${qmcpack_BINARY_DIR}/tests/bin)

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
  SET(USE_PREFETCH 0)
  SET(PREFETCH_AHEAD 12)
  SET(HAVE_C_VARARRAYS 1)

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

    # Flags for ctest
    SET(USE_EXT_MPI_FOR_SERIAL_TESTS 1)
    FIND_PROGRAM(MPIEXEC aprun)
    IF(NOT MPIEXEC)
    # For slurm use srun
      FIND_PROGRAM(MPIEXEC srun)
    ENDIF(NOT MPIEXEC)
    SET(MPIEXEC_NUMPROC_FLAG "-n")

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

  #------------------------------------
  # Identify the compiler
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
  ELSE()
    MESSAGE(WARNING "No default file for compiler (${COMPILER})")
  ENDIF()

  #-------------------------------------------------------------------
  #  check MPI installation
  #-------------------------------------------------------------------
  IF(QMC_MPI)
    ## mpi compilers
    GET_FILENAME_COMPONENT(BASE_CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)
    if(${BASE_CXX_COMPILER_NAME} MATCHES "^mp")
      SET(HAVE_MPI 1)
      SET(HAVE_OOMPI 1)
      SET(MPI_FOUND TRUE)
    elseif($ENV{CRAYPE_VERSION} MATCHES ".")
      SET(HAVE_MPI 1)
      SET(HAVE_OOMPI 1)
      SET(MPI_FOUND TRUE)
    else()
      include(FindMPI)
      IF(MPI_FOUND)
        MESSAGE(STATUS "Found mpi.h and libraries. Check for working mpi compiler ")
        INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
        FILE(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testmpi.cxx 
          "#include <mpi.h>"
          "int main(int argc, char** argv){MPI_Init(&argc,&argv); MPI_Finalize();return 0;}\n")
        #TRY_COMPILE(MPI_WORKS ${CMAKE_BINARY_DIR} 
        #${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testmpi.cxx
        #OUTPUT_VARIABLE OUTPUT)
        set(MPI_WORKS 1)
        SET(MPI_WARNING_LIST
            "Building MPI version without using MPI compiler wrappers.\n"
            "This may not build qmcpack correctly. To ensure the correct version, specify the compiler wrappers to cmake.\n"
            "For example: cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpic++\n"
            "To build without MPI, pass '-DQMC_MPI=0' to cmake")
        MESSAGE(WARNING ${MPI_WARNING_LIST})

        IF(MPI_WORKS)
          MESSAGE(STATUS "MPI is enabled")
          SET(HAVE_MPI 1)
          SET(HAVE_OOMPI 1)
          LINK_LIBRARIES(${MPI_CXX_LIBRARIES})
        ELSE(MPI_WORKS)
          MESSAGE(STATUS "MPI is disabled")
          SET(HAVE_MPI 0)
          SET(HAVE_OOMPI 0)
          SET(MPI_LIBRARY "")
        ENDIF(MPI_WORKS)
      ELSE(MPI_FOUND)
        MESSAGE(STATUS "MPI not found")
        SET(MPI_LIBRARY "")
      ENDIF(MPI_FOUND)
    ENDIF()
  ENDIF(QMC_MPI)

  #-------------------------------------------------------------------
  # check lapack/blas
  #-------------------------------------------------------------------
  # Quick search for MKL (this will need to be improved). Currently only for Intel compilers
  IF (INTEL_COMPILER)
   INCLUDE(CMake/FindMKL.cmake)
  ENDIF(INTEL_COMPILER)
  IF ( HAVE_MKL )
    # We found MKL, add the link library
    SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${MKL_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MKL_FLAGS}")
    SET(LAPACK_FOUND 1)
    SET(LAPACK_LIBRARIES "${MKL_LIBRARIES}")
    SET(LAPACK_LINKER_FLAGS "${MKL_LINKER_FLAGS}")
  ELSEIF ( $ENV{CRAYPE_VERSION} MATCHES "." )
    # Cray specific environment variable indicates we are on a cray. Assume libsci will be linked
    # preferable to attempt a compile/link 
    SET(LAPACK_FOUND 1)
    MESSAGE(STATUS "Assuming Cray libsci or other BLAS/LAPACK will be linked by compiler wrappers")
  ELSE()	 
    # Search for Blas/Lapack
    find_package(LAPACK)
    # When cmake finds and sets up the MKL libraries with a non-Intel compiler, set the flag
    # This should be tidied up when FindMKL.cmake supports non-Intel compilers.
    IF("${LAPACK_LIBRARIES}" MATCHES "mkl")
      MESSAGE(STATUS "MKL found via LAPACK/BLAS")
      SET( MKL_FOUND 1 )
      SET( HAVE_MKL 1 )
    ENDIF()
  ENDIF()
  IF(LAPACK_FOUND)
    MESSAGE(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")
    MESSAGE(STATUS "LAPACK linker flags: ${LAPACK_LINKER_FLAGS}")
    SET(LAPACK_LIBRARY ${LAPACK_LIBRARIES} ${LAPACK_LINKER_FLAGS})
  ELSE()
    MESSAGE(FATAL_ERROR "Could not find required libraries LAPACK &/or BLAS") 
  ENDIF()

  INCLUDE(CMake/FindVectorMath.cmake)

ENDIF(CMAKE_TOOLCHAIN_FILE)

IF (ENABLE_GCOV)
  IF (NOT GCOV_SUPPORTED)
    MESSAGE(FATAL_ERROR "ENABLE_GCOV was specified but compiler does not support GCC coverage flag")
  ENDIF()
ENDIF(ENABLE_GCOV)

# Require C++11 while leaving possibility to fall back by modifying this file
SET(CXX11_NEEDED TRUE)
IF(CXX11_NEEDED)
  SET(CXX11_FLAG "-std=c++11")
  #check if the CXX compiler supports -std=c++11 option
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG(${CXX11_FLAG} CXX_COMPILER_SUPPORT_CXX11)
  # Turn on C++11 for this code

  # Force the flag on Cray with Intel compiler, because the Cray wrapper
  #  prints an warning that interferes with the flag detection code
  #  with older versions of CMake.
  IF($ENV{CRAYPE_VERSION} MATCHES ".")
    IF( ${COMPILER} MATCHES "Intel" AND NOT CXX_COMPILER_SUPPORT_CXX11)
      SET(CXX_COMPILER_SUPPORT_CXX11 TRUE)
      MESSAGE(STATUS "Forcing C++11 support on Cray with Intel")
    ENDIF()
  ENDIF()

  IF (CXX_COMPILER_SUPPORT_CXX11)
    # avoid repeated -std=c++11 flag
    STRING(REPLACE "++" "\\+\\+" CXX11_FLAG_MATCH ${CXX11_FLAG})
    IF(NOT CMAKE_CXX_FLAGS MATCHES ${CXX11_FLAG_MATCH})
      SET (CMAKE_CXX_FLAGS "${CXX11_FLAG} ${CMAKE_CXX_FLAGS}")
    ENDIF()
  ELSE()
    MESSAGE(FATAL_ERROR "A compiler supporting C++11 is required. Use a newer C++ compiler.")
  ENDIF()
ENDIF(CXX11_NEEDED)

# AFQMC requires MKL sparse for good performance (roughly a factor of 2x)
IF (BUILD_AFQMC AND NOT MKL_FOUND)
  MESSAGE(WARNING "AFQMC - MKL not found, using simple sparse matrix routines.  Link with MKL sparse libraries for better performance.")
ENDIF()

# AFQMC requires MPI
If (BUILD_AFQMC AND NOT QMC_MPI)
  MESSAGE(FATAL_ERROR "AFQMC requires building with MPI (QMC_MPI=1). Set BUILD_AFQMC=0 or configure MPI.")
ENDIF()

IF (BUILD_AFQMC AND NOT APPLE)
  LINK_LIBRARIES("rt")
ENDIF()
  
# setup ctest variables
IF ( HAVE_MPI )
  SET( USE_MPI 1 )
  IF ( NOT MPIEXEC )
    SET( MPIEXEC "mpirun" )
  ENDIF()

  IF ( NOT DEFINED MPIEXEC_NUMPROC_FLAG )
    SET( MPIEXEC_NUMPROC_FLAG "-np" )
  ENDIF()

  IF ( NOT TEST_MAX_PROCS )
    SET( TEST_MAX_PROCS 100 )
  ENDIF()

  MESSAGE(STATUS "MPI runner is ${MPIEXEC}")
ENDIF(HAVE_MPI)

MESSAGE(STATUS "CXX compiler flags: ${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "CC  compiler flags: ${CMAKE_C_FLAGS}")

SET(HAVE_OOMPI ${HAVE_MPI})
INCLUDE(${CMAKE_ROOT}/Modules/FindThreads.cmake)

####################################################################
#First check the required libraries. Abort if these are not found.
#Changed the fall-back versions to use similar conventions to 2.8.0 and higher
####################################################################
# use customized cmake file for boost so that we can deal with multiple
# installataions
# list of libraries to be linked with the main applications using I/O

SET(QMC_UTIL_LIBS ${LAPACK_LIBRARY} ${BLAS_LIBRARY})

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
  INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIR})
  SET(QMC_UTIL_LIBS ${QMC_UTIL_LIBS} ${FFTW_LIBRARIES})
  MESSAGE(STATUS "Found FFTW")
ELSE(FFTW_FOUND)
#  MESSAGE(STATUS "Disabling FFTW")
  MESSAGE(FATAL_ERROR "Could not find required library FFTW. Verify FFTW_HOME")
ENDIF(FFTW_FOUND)
 
SET(QMC_UTIL_LIBS ${QMC_UTIL_LIBS} ${FORTRAN_LIBRARIES})

find_package(ZLIB)
#find_package(LibXml2) 
INCLUDE(CMake/FindLibxml2QMC.cmake)
   
IF(LIBXML2_FOUND)
  INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR})
  SET(HAVE_LIBXML2 1)
  SET(QMC_UTIL_LIBS ${QMC_UTIL_LIBS} ${LIBXML2_LIBRARIES})
ELSE()
  MESSAGE(FATAL_ERROR "Could not find required library Libxml2")
ENDIF(LIBXML2_FOUND)

set(HDF5_USE_STATIC_LIBRARIES off)
if(QMC_BUILD_STATIC)
  set(HDF5_USE_STATIC_LIBRARIES on)
endif()

SET(HDF5_PREFER_PARALLEL 0)
IF(HAVE_MPI)
 SET(HDF5_PREFER_PARALLEL 1)
ENDIF()
find_package(HDF5 COMPONENTS C)
IF(HDF5_FOUND)
  SET(HAVE_LIBHDF5 1)
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
  IF ( CMAKE_BUILD_TYPE AND HDF5_LIBRARIES_DEBUG )
    IF ( CMAKE_BUILD_TYPE MATCHES DEBUG )
      SET( HDF5_LIBRARIES ${HDF5_LIBRARIES_DEBUG} )
    ELSE()
      SET( HDF5_LIBRARIES ${HDF5_LIBRARIES_RELEASE} )
    ENDIF()
  ENDIF()
  SET(QMC_UTIL_LIBS ${QMC_UTIL_LIBS} ${HDF5_LIBRARIES})
  if(ZLIB_LIBRARY)
    SET(QMC_UTIL_LIBS ${QMC_UTIL_LIBS} ${ZLIB_LIBRARY})
  endif(ZLIB_LIBRARY)
  #IF(SZLIB_FOUND)
  #  SET(QMC_UTIL_LIBS ${QMC_UTIL_LIBS} ${SZLIB_LIBRARIES})
  #ENDIF(SZLIB_FOUND)
  IF(HDF5_IS_PARALLEL)
    IF(ENABLE_PHDF5)
      MESSAGE(STATUS "Using HDF5 parallel collective I/O")
    ELSE(ENABLE_PHDF5)
      MESSAGE(STATUS "Using HDF5 non-scalable serial I/O")
    ENDIF(ENABLE_PHDF5)
  ELSE(HDF5_IS_PARALLEL)
    IF(ENABLE_PHDF5)
      MESSAGE(STATUS "Using HDF5 non-scalable serial I/O due to the lack of library support for parallel")
      SET(ENABLE_PHDF5 0)
    ELSE(ENABLE_PHDF5)
      MESSAGE(STATUS "Using HDF5 non-scalable serial I/O")
    ENDIF(ENABLE_PHDF5)
  ENDIF(HDF5_IS_PARALLEL)
ELSE(HDF5_FOUND)
  if(NOT QMC_PHI)
    MESSAGE(FATAL_ERROR "Require hdf5 1.6.4 or higher. Set HDF5_ROOT")
  endif()
ENDIF(HDF5_FOUND)

#make sure we can find boost if it's not in /usr
set(Boost_NO_BOOST_CMAKE on)
find_package(Boost REQUIRED)
IF(Boost_FOUND)
  SET(HAVE_LIBBOOST 1)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  MESSAGE(STATUS "Setting Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIRS}")
ELSE()
  MESSAGE(STATUS "Using MT engine in Utilities/MersenneTwister.h")
  #MESSAGE(FATAL_ERROR "Require boost 1.33.x or higher. Set BOOST_ROOT")
ENDIF()

SET(HAVE_CUDA 0)
IF(QMC_CUDA)
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  FIND_PACKAGE(CUDA REQUIRED)
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
  set(HAVE_CUDA 1)
  MESSAGE("   CUDA_NVCC_FLAGS=${CUDA_NVCC_FLAGS}")
ELSE(QMC_CUDA)
  MESSAGE(STATUS "Disabling CUDA")
ENDIF(QMC_CUDA)

#INCLUDE(${PROJECT_CMAKE}/FindPkgConfig.cmake)
##################################################################
# TODO:use profile tools
##################################################################
#IF(MPIP_PROFILE)
#  LINK_LIBRARIES("/usr/apps/tools/mpip/lib/cmpi/libmpiP.a")
#  LINK_LIBRARIES("/usr/lib/libbfd.a")
#  LINK_LIBRARIES("/usr/lib/libiberty.a")
#ENDIF(MPIP_PROFILE)
#if(ENABLE_TAU_PROFILE)
#  INCLUDE(${PROJECT_CMAKE}/FindTau.cmake)
#endif(ENABLE_TAU_PROFILE)
#    #add tau stuff
#    IF(FOUND_TAU)
#      ADD_DEFINITIONS(
#        -DTAU_MPI
#        -DTAU_MPI_THREADED
#        -DTAU_MPIGREQUEST
#        -DTAU_MPIDATAREP
#        -DTAU_MPIERRHANDLER
#        -DTAU_MPIATTRFUNCTION
#        -DTAU_MPITYPEEX
#        -DTAU_MPIADDERROR
#      )
#      #build libtau
#      subdirs(src/tau)
#      LINK_LIBRARIES(tau)
#    ENDIF(FOUND_TAU)
#IF(FOUND_TAU)
#  INCLUDE_DIRECTORIES(${TAU_INCLUDE_DIR})
#  ADD_DEFINITIONS (-DPROFILING_ON
#      -DTAU_STDCXXLIB
#      -DTAU_DOT_H_LESS_HEADERS
#      -DTAU_CALLPATH
#      -DTAU_LARGEFILE64_SOURCE
#      -DTAU_BFD
#      )
#  if(ENABLE_OPENMP)
#    ADD_DEFINITIONS(-DTAU_OPENMP)
#  endif(ENABLE_OPENMP)
#ENDIF(FOUND_TAU)

######################################################
# Add optional macros
######################################################

#include qmcpack/src build/src
INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/src ${PROJECT_BINARY_DIR}/src)

######################################################################
# Manage Einspline library
# - toolchain file can set 
# - check if libeinspline can be found in $EINSPLINE_HOME/include and
# $EINSPLINE_HOME/lib
# - if not found, build it from the source directory in $EINSPLINE_HOME
# - if everything fails, do not use it
######################################################################
SET(HAVE_EINSPLINE 1)
set(build_einspline 1)
#if(HAVE_EINSPLINE_EXT)
#  set(HAVE_EINSPLINE 0)
#  INCLUDE(${PROJECT_CMAKE}/FindEinspline.cmake)
#  IF(EINSPLINE_FOUND)
#    INCLUDE_DIRECTORIES(${EINSPLINE_INCLUDE_DIR})
#    SET(QMC_UTIL_LIBS ${EINSPLINE_LIBRARIES} ${QMC_UTIL_LIBS})
#    set(build_einspline 0)
#  endif()
#endif()
if(build_einspline)
  SUBDIRS(src/einspline)
  INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/include)
  SET(QMC_UTIL_LIBS ${QMC_UTIL_LIBS} einspline)
endif()

IF(BUILD_LMYENGINE_INTERFACE)
  ADD_DEFINITIONS(-DHAVE_LMY_ENGINE)
  SET(FORMIC_BINARY_DIR ${qmcpack_BINARY_DIR})
  SUBDIRS(src/formic/utils)
  SET(QMC_UTIL_LIBS formic_utils ${QMC_UTIL_LIBS})
ENDIF(BUILD_LMYENGINE_INTERFACE)

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
#######################################################################
## Manage GTK+ and GTKMM packages for GUI
## search headers and libraries
##   1) -lXYZ-${XYZ_VERSION} 
##   2) -lXYZ
## example: -lgtk-${GTK_VERSION} | -lgtk
#######################################################################
#SET(GTK_VERSION 2.0)
#SET(PANGO_VERSION 1.0)
#
#SET(GTKMM_VERSION 2.4)
#SET(PANGOMM_VERSION 1.4)
#SET(SIGCPP_VERSION 2.0)
#IF(ENABLE_GTK_GUI)
#  INCLUDE(${PROJECT_CMAKE}/FindGTK.cmake)
#  #THIS SHOULD GO TO the directory where GUI is built
#  IF(GTKMM_FOUND)
#   INCLUDE_DIRECTORIES(${GTKMM_INCLUDE_DIR})
#   LINK_LIBRARIES(${GTKMM_LIBRARIES})
#  ENDIF(GTKMM_FOUND)
#ENDIF(ENABLE_GTK_GUI)

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

###################################################
# Done searching libraries. Start building qmcpack
###################################################

MESSAGE(STATUS "Ready to build qmcpack")

####################################################################
# add definitions
#  -DH5_USE_16_API for HDF5
#  -DMPICH_SKIP_MPICXX : disable C++ binding of MPICH
####################################################################

ADD_DEFINITIONS(-DHAVE_CONFIG_H -DH5_USE_16_API -DUSE_REAL_STRUCT_FACTOR -DINLINE_ALL=inline -DADD_)

if(HAVE_MPI)
  ADD_DEFINITIONS(-DMPICH_SKIP_MPICXX)
endif()
#ADD_DEFINITIONS(-DDISABLE_WS_CELL)

EXECUTE_PROCESS(COMMAND hostname OUTPUT_VARIABLE QMC_HOSTNAME OUTPUT_STRIP_TRAILING_WHITESPACE)
CONFIGURE_FILE(${qmcpack_SOURCE_DIR}/src/qmcpack.settings ${qmcpack_BINARY_DIR}/bin/qmcpack.settings)

#########################
# add definitions for adios profile
# add definitions for adios verification
########################
IF(QMC_ADIOS)
  #ADD_DEFINITIONS(-DIO_PROFILE)
  #ADD_DEFINITIONS(-DADIOS_VERIFY)
  ADD_DEFINITIONS(-DWRITE_GLOBAL)
ENDIF(QMC_ADIOS)

SUBDIRS(src)

IF(NOT(QMC_BUILD_LEVEL GREATER 4))
  SUBDIRS(tests)
  SUBDIRS(examples)
  IF(IS_DIRECTORY nexus)
    SUBDIRS(nexus/examples)
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
