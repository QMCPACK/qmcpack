######################################################################
# CMake version and policies
######################################################################
cmake_minimum_required(VERSION 3.14.0)

# CMP0074: CMake find_package will use <PackageName>_ROOT CMake variable
# and environment variable in search path.
# various Find<Package>.cmake modules may not follow this policy
# CMP0075: check_include_file_cxx and other follow CMAKE_REQUIRED_LIBRARIES
# when doing link check
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0075 NEW)
include(CMakeDependentOption)

######################################################################
# QMCPACK
######################################################################
project(
  qmcpack
  VERSION 3.11.9
  LANGUAGES C CXX)

######################################################
# Directory where customize cmake files reside
######################################################
set(PROJECT_CMAKE ${qmcpack_SOURCE_DIR}/CMake)

set(QMCPACK_UNIT_TEST_DIR ${qmcpack_BINARY_DIR}/tests/bin)

include(CMake/macros.cmake)

######################################################################
# CTest
######################################################################
# The following option disables the extremely slow setup of full system run tests
# This is a developer option
# It offers a more reasonable compile debug loop if CMakeLists.txt files are being changed
# These constom commands execution time far exceeds even the feature and compiler checks for time consumed
# when rerunning cmake.
#
# On leconte: rerun of cmake with
# slow custom commands 27 secs
# without 8 secs
option(QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS "Disable the slow cmake custom commands for integration tests." OFF)
mark_as_advanced(QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)

set(DROP_METHOD "http")
set(DROP_SITE "cdash.qmcpack.org")
set(DROP_LOCATION "/CDash/submit.php?project=QMCPACK")
set(TRIGGER_SITE "")
set(DROP_SITE_CDASH TRUE)
# Increase timeout per test over the default of 1500 seconds (25 minutes)
set(DART_TESTING_TIMEOUT
    3600
    CACHE STRING "Maximum time for one test")
enable_testing()
include(CTest)

option(QMC_SYMLINK_TEST_FILES "Use symbolic links for test files to save space.  Set to false to copy files instead."
       ON)

if(QMC_SYMLINK_TEST_FILES)
  set(SYMLINK_MSG
      "Using symbolic links for large test files may cause test failures if the build is installed on a separate filesystem from the source."
  )
else()
  set(SYMLINK_MSG "Copying large test files uses more disk space than using symbolic links.")
endif()
message(STATUS "QMC_SYMLINK_TEST_FILES = ${QMC_SYMLINK_TEST_FILES}.  ${SYMLINK_MSG}")

######################################################################
# Verify Python3 available
######################################################################
include(CMake/python.cmake)
find_package(Python3)
if(NOT Python3_FOUND)
  message(FATAL_ERROR "Could not find required python3")
endif(NOT Python3_FOUND)

######################################################################
# Check if PySCF is available
######################################################################
if(NOT QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)
  test_python_module(pyscf HAVE_PYSCF)
  if(NOT HAVE_PYSCF)
    message(STATUS "Unable to import PySCF python module. PySCF tests will not be run.")
  else()
    message(STATUS "Successfully imported PySCF python module.")
  endif()
endif()

######################################################################
# Check if QE executables available at QE_BIN or on the PATH
######################################################################
if(NOT QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)
  include("CMake/FindQE.cmake")
  if(QE_FOUND)
    message(STATUS "Found and using patched Quantum ESPRESSO (QE) with pw2qmcpack.x at ${QE_PW_DIR}")
  else(QE_FOUND)
    message(
      STATUS "Did not find a patched Quantum ESPRESSO (QE) distribution with pw2qmcpack.x. QE tests will not be run.")
  endif(QE_FOUND)
endif()

######################################################################
# Build level
######################################################################

option(QMC_BUILD_SANDBOX_ONLY "Build only applications in Sandbox directory" OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_BUILD_TYPE Release)
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)
message(STATUS "CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")

######################################################################
#define the paths for library and executable for out-of-soruce compilation
######################################################################
set(LIBRARY_OUTPUT_PATH
    ${qmcpack_BINARY_DIR}/lib
    CACHE PATH "Single output directory for building all libraries.")
set(EXECUTABLE_OUTPUT_PATH
    ${qmcpack_BINARY_DIR}/bin
    CACHE PATH "Single output directory for building all executables.")

##################################
# Set C++ standard used by QMCPACK
##################################
set(QMC_CXX_STANDARD
    14
    CACHE STRING "QMCPACK C++ language standard")
if(NOT QMC_CXX_STANDARD EQUAL 14)
  message(WARNING "C++14 is the only language standard officially supported by this QMCPACK version. "
                  "Using other versions of the C++ standard is unsupported and done entirely at user's own risk.")
endif()

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
set(WALKER_MAX_PROPERTIES
    2048
    CACHE STRING "Maximum number of properties tracked by walkers")
mark_as_advanced(WALKER_MAX_PROPERTIES)
set(OHMMS_DIM
    3
    CACHE STRING "Select physical dimension")
set(OHMMS_INDEXTYPE int)
message(STATUS "defining the float point precision")
set(OHMMS_PRECISION_FULL double)

if(QMC_CUDA)
  set(QMC_MIXED_PRECISION
      1
      CACHE BOOL "Enable/disable mixed precision")
  set(OHMMS_PRECISION double)
  set(CUDA_PRECISION_FULL double)
  if(QMC_MIXED_PRECISION)
    set(CUDA_PRECISION float)
  else(QMC_MIXED_PRECISION)
    set(CUDA_PRECISION double)
  endif(QMC_MIXED_PRECISION)
  message("   Base precision = ${OHMMS_PRECISION}")
  message("   Full precision = ${OHMMS_PRECISION_FULL}")
  message("   CUDA base precision = ${CUDA_PRECISION}")
  message("   CUDA full precision = ${CUDA_PRECISION_FULL}")
else(QMC_CUDA)
  set(QMC_MIXED_PRECISION
      0
      CACHE BOOL "Enable/disable mixed precision")
  if(QMC_MIXED_PRECISION)
    set(OHMMS_PRECISION float)
    set(MIXED_PRECISION 1)
  else(QMC_MIXED_PRECISION)
    set(OHMMS_PRECISION double)
  endif(QMC_MIXED_PRECISION)
  message("   Base precision = ${OHMMS_PRECISION}")
  message("   Full precision = ${OHMMS_PRECISION_FULL}")
endif(QMC_CUDA)

# Code coverage
set(GCOV_SUPPORTED FALSE)
option(ENABLE_GCOV "Enable code coverage" OFF)

######################################################################
# enable MPI and OPENMP  if possible
######################################################################
option(QMC_MPI "Enable/disable MPI" ON)
option(QMC_OMP "Enable/disable OpenMP" ON)
option(QMC_COMPLEX "Build for complex binary" OFF)
set(QMC_CUDA
    0
    CACHE BOOL "Build with GPU support through CUDA")
option(
  ENABLE_CUDA
  "Build with the second generation of GPU support through CUDA (production quality for AFQMC, experimental for real space)"
  OFF)
option(ENABLE_HIP "Build with with GPU support through HIP" OFF)
option(ENABLE_ROCM "Build with with GPU support through ROCM" OFF)
option(ENABLE_OFFLOAD "Enable OpenMP offload" OFF)
# Use CMake object library targets to workaround clang linker not being able to handle fat
# binary archives which contain both host and device codes, for example OpenMP offload regions.
# CMake does not propagate indirect object files by design.
# So the dependency must be explicitly specified with USE_OBJECT_TARGET.
cmake_dependent_option(USE_OBJECT_TARGET "Use CMake object library target" OFF "ENABLE_OFFLOAD" OFF)

if(ENABLE_CUDA AND QMC_CUDA)
  message(FATAL_ERROR "ENABLE_CUDA=1 and QMC_CUDA=1 can not be set together!")
endif(ENABLE_CUDA AND QMC_CUDA)

option(QMC_VERBOSE_CONFIGURATION "Give details why specific tests are/are not added" OFF)

######################################################################
# set debug printout
# DEBUG_PSIBUFFER_ON keep track buffer's locator
######################################################################
set(PRINT_DEBUG
    0
    CACHE BOOL "Enable/disable debug printing")
set(DEBUG_PSIBUFFER_ON
    0
    CACHE BOOL "Debug messages for PSIBUFFER locator")
mark_as_advanced(DEBUG_PSIBUFFER_ON)
set(DEBUG_PER_STEP_ACCEPT_REJECT
    0
    CACHE BOOL "Print accepts and rejects at each step")
mark_as_advanced(DEBUG_PER_STEP_ACCEPT_REJECT)

######################################################################
# Standard test
######################################################################
set(CMAKE_MODULE_PATH "${qmcpack_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})

######################################################################
# add macros to compile einspline if necessary
######################################################################
include(${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
include(${CMAKE_ROOT}/Modules/CheckCXXSymbolExists.cmake)
include(${CMAKE_ROOT}/Modules/CheckIncludeFile.cmake)
include(${CMAKE_ROOT}/Modules/CheckTypeSize.cmake)

check_function_exists(posix_memalign HAVE_POSIX_MEMALIGN)

######################################################################
# Build options
# QMC_BUILD_STATIC build static/dynamic  libraries
# BUILD_QMCTOOLS   build utility programs
# MPIP_PROFILE     profile mpi performance
######################################################################
option(BUILD_UNIT_TESTS "Build unit tests" ON)
option(BUILD_MICRO_BENCHMARKS "Build micro benchmarks" ON)
option(BUILD_LMYENGINE_INTERFACE "Build LMY engine" ON)
if(QMC_CUDA AND BUILD_LMYENGINE_INTERFACE)
  message(STATUS "LMY engine is not compatible with CUDA build! Disabling LMY engine")
  set(BUILD_LMYENGINE_INTERFACE OFF)
endif()
if(MIXED_PRECISION AND BUILD_LMYENGINE_INTERFACE)
  message(STATUS "LMY engine is not compatible with CPU mixed precision build! Disabling LMY engine")
  set(BUILD_LMYENGINE_INTERFACE OFF)
endif()
option(BUILD_AFQMC "Build with AFQMC" OFF)
option(BUILD_AFQMC_WITH_NCCL "Build AFQMC with NCCL library." OFF)
# AFQMC requires MPI.
if(BUILD_AFQMC AND NOT QMC_MPI)
  message(FATAL_ERROR "AFQMC requires building with MPI (QMC_MPI=1). Set BUILD_AFQMC=0 or configure MPI.")
endif()
option(BUILD_FCIQMC "Build with FCIQMC" OFF)
option(QMC_BUILD_STATIC "Link to static libraries" OFF)
option(ENABLE_TIMERS "Enable internal timers" ON)
option(ENABLE_STACKTRACE "Enable use of boost::stacktrace" OFF)
option(USE_VTUNE_API "Enable use of VTune ittnotify APIs" OFF)
cmake_dependent_option(USE_VTUNE_TASKS "USE VTune ittnotify task annotation" OFF "ENABLE_TIMERS AND USE_VTUNE_API" OFF)
# CMake note - complex conditionals in cmake_dependent_option must have spaces around parentheses
cmake_dependent_option(USE_NVTX_API "Enable/disable NVTX regions in CUDA code." OFF
                       "ENABLE_TIMERS AND ( QMC_CUDA OR ENABLE_CUDA )" OFF)

######################################################################
# Sanitizer options
######################################################################

# Add optional sanitizers ASAN, UBSAN, MSAN
set(VALID_SANITIZERS "none" "asan" "ubsan" "tsan" "msan")
set(ENABLE_SANITIZER
    "none"
    CACHE STRING "none,asan,ubsan,tsan,msan")
set_property(CACHE ENABLE_SANITIZER PROPERTY STRINGS ${VALID_SANITIZERS})

# Perform sanitizer option check, only works in debug mode
if(NOT ENABLE_SANITIZER IN_LIST VALID_SANITIZERS)
  message(FATAL_ERROR "Invalid -DENABLE_SANITIZER=${ENABLE_SANITIZER}, value must be one of ${VALID_SANITIZERS}")
else()
  message(STATUS "Enable sanitizer ENABLE_SANITIZER=${ENABLE_SANITIZER}")
endif()

######################################################################
# Install options
# INSTALL_NEXUS install Nexus alongside QMCPACK
######################################################################
set(INSTALL_NEXUS
    1
    CACHE BOOL "Install Nexus alongside QMCPACK")

######################################################################
#list of options to determine libraries. First, try to use libraries available
#FindXYZ.cmake will set the values if successful
######################################################################
set(HAVE_EINSPLINE 1) # to be removed

######################################################
# Git information (see src/CMakeLists.txt for the
# command to get repository information)
######################################################
find_package(Git)

if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  set(IS_GIT_PROJECT 1)
endif()

######################################################
#COMPILER choose one of the cmake files to customize the compiler options
#If nothing is chosen, default settings by cmake will  be used.
#If the automatic detection does not work, comment out everything
#upto COMPILER MANUAL SELECTION and use one of the customized cmake file.
######################################################
include(${CMAKE_ROOT}/Modules/CheckCCompilerFlag.cmake)
include(${CMAKE_ROOT}/Modules/CheckCXXCompilerFlag.cmake)

######################################################
# DO NOT MODIFY BELOW
# If you like to customize your build, create a toolchain as instructed
######################################################
######################################################
# now set the compiler flags and the macros in src/config.h
######################################################
if(CMAKE_TOOLCHAIN_FILE)
  message(STATUS "Using ${CMAKE_TOOLCHAIN_FILE} toolchain ")
else(CMAKE_TOOLCHAIN_FILE)
  message(STATUS "Trying to figure out .....")
  set(ENABLE_OPENMP 0)
  set(HAVE_ACML 0)

  #------------------------
  # On Cray's machine
  #------------------------
  get_filename_component(BASE_CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)
  if($ENV{CRAYPE_VERSION} MATCHES "." AND NOT BASE_CXX_COMPILER_NAME STREQUAL "CC")
    message(WARNING "Cray Programming Environment has been detected but C++ compiler wrapper CC is not used! "
                    "Please use cc/CC as the C/C++ compiler to ensure correct cross-compiling for the compute nodes "
                    "unless you understand cross-compiling and intend not to use Cray compiler wrappers.")
  endif()
  if(CMAKE_SYSTEM_NAME STREQUAL "CrayLinuxEnvironment")
    message(STATUS "Running on a Cray machine.")

    if(NOT $ENV{CRAYPE_LINK_TYPE} STREQUAL "dynamic")
      message(WARNING "Cray Programming Environment uses static linking by default. "
                      "We prefer dynamic linking which makes library searching easier. "
                      "Set environment variable CRAYPE_LINK_TYPE to \"dynamic\" enables dynamic linking.\n")
    endif()

    if(NOT DEFINED MPIEXEC_EXECUTABLE)
      find_program(MPIEXEC_EXECUTABLE aprun)
      if(NOT MPIEXEC_EXECUTABLE)
        # For slurm use srun
        find_program(MPIEXEC_EXECUTABLE srun)
      endif(NOT MPIEXEC_EXECUTABLE)
    endif()
  else()
    if($ENV{CRAYPE_VERSION} MATCHES "." AND BASE_CXX_COMPILER_NAME STREQUAL "CC")
      message(
        FATAL_ERROR
          "Cray compiler wrapper detected. Empty the build folder and rerun cmake with -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment added."
      )
    endif()
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(COMPILER GNU)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "XL")
    set(COMPILER IBM)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
    set(COMPILER Clang)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(COMPILER Intel)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "PGI" OR CMAKE_CXX_COMPILER_ID MATCHES "NVHPC")
    set(COMPILER PGI)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Cray")
    set(COMPILER Cray)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(COMPILER Clang)
  else()
    message("${CMAKE_CXX_COMPILER_ID}")
    message(WARNING "Unknown C/C++ compiler ${CMAKE_CXX_COMPILER_ID}, default flags will be used")
  endif()
  message(STATUS "C++ Compiler is identified by QMCPACK as : ${COMPILER}")

  #------------------------------------
  # Include compiler-specific cmake file
  #------------------------------------
  if(COMPILER MATCHES "IBM")
    include(${PROJECT_CMAKE}/IBMCompilers.cmake)
  elseif(COMPILER MATCHES "Intel")
    include(${PROJECT_CMAKE}/IntelCompilers.cmake)
  elseif(COMPILER MATCHES "GNU")
    include(${PROJECT_CMAKE}/GNUCompilers.cmake)
  elseif(COMPILER MATCHES "Clang")
    include(${PROJECT_CMAKE}/ClangCompilers.cmake)
  elseif(COMPILER MATCHES "PGI")
    include(${PROJECT_CMAKE}/PGICompilers.cmake)
  else()
    message(WARNING "No default file for compiler (${COMPILER})")
  endif()

endif(CMAKE_TOOLCHAIN_FILE)

if(NOT CMAKE_CXX_COMPILER_ID STREQUAL CMAKE_C_COMPILER_ID)
  message(FATAL_ERROR "Mixing C and C++ compilers from different vendors is not permitted. Current "
                      "C compiler is ${CMAKE_C_COMPILER_ID} but C++ compiler is ${CMAKE_CXX_COMPILER_ID}.")
endif()

###############################################
# Set C++ minimum standard and run basic checks
###############################################
set(CMAKE_CXX_STANDARD ${QMC_CXX_STANDARD})
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#-----------------------------------------------------------------------
# Check that a C++ compiler is compatiable with the underlying libstdc++
#-----------------------------------------------------------------------
include(CMake/Testlibstdc++.cmake)

#---------------------------------------------------------
# Check that a C++ standard library is configured properly
#---------------------------------------------------------
if(QMC_CXX_STANDARD GREATER_EQUAL 17)
  include(CMake/TestCxx17Library.cmake)
else()
  include(CMake/TestCxx14Library.cmake)
endif()

if(ENABLE_GCOV)
  if(NOT GCOV_SUPPORTED)
    message(FATAL_ERROR "ENABLE_GCOV was specified but compiler does not support GCC coverage flag")
  endif()
endif(ENABLE_GCOV)

if(BUILD_AFQMC AND NOT APPLE)
  link_libraries("rt")
endif()

#-------------------------------------------------------------------------------
# Check SIMD alignment for CPU only
#-------------------------------------------------------------------------------
# This is intentionally placed before adding OpenMP offload compile options
# to avoid contamination from device compilation pass.
# When '-march=skylake-avx512 -fopenmp-targets=nvptx64 -march=sm_70' is added
# for Clang, the source code is parsed twice for both host and offload targets.
# A trap for macro __AVX512F__ always fails because the offload pass doesn't
# carry '-march=skylake-avx512' but only takes '-march=sm_70'.
#-------------------------------------------------------------------------------
include(CMake/CheckSIMDAlignment.cmake)
set(QMC_SIMD_ALIGNMENT
    ${default_alignment}
    CACHE STRING "Cache/SIMD alignment in bytes")
math(EXPR alignment_remainder "${QMC_SIMD_ALIGNMENT} % ${default_alignment}")
if(alignment_remainder)
  message(
    FATAL_ERROR
      "QMC_SIMD_ALIGNMENT must be a multiple of ${default_alignment}. Bad cached value is ${QMC_SIMD_ALIGNMENT}")
endif()
message(STATUS "QMC_SIMD_ALIGNMENT is set to ${QMC_SIMD_ALIGNMENT}")

#---------------------------------------------------------
# Determine if OpenMP taskloop works with the CXX compiler
#---------------------------------------------------------
include(CMake/TestOpenMPtaskloop.cmake)
option(ENABLE_OMP_TASKLOOP "Enable OpenMP taskloop" ${OMP_TASKLOOP_OKAY})
message(STATUS "ENABLE_OMP_TASKLOOP is set to ${ENABLE_OMP_TASKLOOP}")

#---------------------------------------------------------
# Set up OpenMP offload compile options
#---------------------------------------------------------
if(ENABLE_OFFLOAD AND DEFINED OPENMP_OFFLOAD_COMPILE_OPTIONS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMP_OFFLOAD_COMPILE_OPTIONS}")
endif()

#-------------------------------------------------------------------
#  check MPI installation
#-------------------------------------------------------------------
if(QMC_MPI)
  # for backward compatibility with MPIEXEC
  if(DEFINED MPIEXEC AND NOT DEFINED MPIEXEC_EXECUTABLE)
    message(
      WARNING
        "Setting MPIEXEC_EXECUTABLE based on MPIEXEC. MPIEXEC has been deprecated! Use MPIEXEC_EXECUTABLE instead!")
    set(MPIEXEC_EXECUTABLE ${MPIEXEC})
  endif(DEFINED MPIEXEC AND NOT DEFINED MPIEXEC_EXECUTABLE)
  set(MPI_CXX_SKIP_MPICXX TRUE)
  if(NOT CMAKE_CROSSCOMPILING)
    set(MPI_DETERMINE_LIBRARY_VERSION TRUE)
  endif()
  find_package(MPI COMPONENTS CXX)

  if(NOT MPI_FOUND)
    message(
      FATAL_ERROR
        "MPI support not found! Provide MPI compiler wrappers or build without MPI by passing '-DQMC_MPI=0' to cmake.")
  endif(NOT MPI_FOUND)

  if(${MPI_CXX_LIBRARY_VERSION_STRING} MATCHES "MVAPICH2")
    string(REGEX REPLACE "\n" ";" ONE_LINE "${MPI_CXX_LIBRARY_VERSION_STRING}")
    string(REGEX REPLACE " +|\t" ";" ONE_LINE "${ONE_LINE}")
    list(GET ONE_LINE 3 MVAPICH2_VERSION)
    message(STATUS "MVAPICH2 version ${MVAPICH2_VERSION}")
    if(${MVAPICH2_VERSION} VERSION_LESS "2.4" AND NOT ${MPI_CXX_LIBRARY_VERSION_STRING} MATCHES
                                                  "disable-registration-cache")
      message(
        FATAL_ERROR
          "MVAPICH2 (version < 2.4) with registration cache enabled breaks QMCPACK. "
          "Use a different MPI library or a MVAPICH2 >=2.4 "
          "or build MVAPICH2 with --disable-registration-cache configure option.")
    endif()
  endif()

  if(${MPI_CXX_LIBRARY_VERSION_STRING} MATCHES "Open MPI")
    if(NOT MPIEXEC_PREFLAGS)
      set(MPIEXEC_PREFLAGS "--bind-to;none")
      message(STATUS "Detected Open MPI. Setting bind-to options for thread affinity in MPIEXEC_PREFLAGS.")
    else(NOT MPIEXEC_PREFLAGS)
      if(NOT "${MPIEXEC_PREFLAGS}" MATCHES "--bind-to")
        message(
          WARNING "Default Open MPI thread affinity policy may bind all the threads to a single core and "
                  "significantly slow down testing. Add proper --bind-to options via MPIEXEC_PREFLAGS. "
                  "If you don't know exactly what affinity to add, "
                  "add '--bind-to;none' to your current MPIEXEC_PREFLAGS to stop this warning.")
      endif()
    endif(NOT MPIEXEC_PREFLAGS)
  endif()

  if(NOT ${CMAKE_CXX_COMPILER} STREQUAL ${MPI_CXX_COMPILER})
    set(MPI_WARNING_LIST
        "Building MPI version without using MPI compiler wrappers.\n"
        "This may not build qmcpack correctly. To ensure the correct version, specify the compiler wrappers to cmake.\n"
        "For example: cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpic++\n"
        "To build without MPI, pass '-DQMC_MPI=0' to cmake")
    message(WARNING ${MPI_WARNING_LIST})
  endif()

  if(NOT TEST_MAX_PROCS)
    set(TEST_MAX_PROCS 100)
  endif()

  message(STATUS "MPI runner MPIEXEC_EXECUTABLE : ${MPIEXEC_EXECUTABLE}")
  message(STATUS "MPIEXEC_NUMPROC_FLAG : ${MPIEXEC_NUMPROC_FLAG}")
  message(STATUS "MPIEXEC_PREFLAGS : ${MPIEXEC_PREFLAGS}")
  string(REPLACE ";" " " MPIEXEC_PREFLAGS_PRINT "${MPIEXEC_PREFLAGS}")
  message(
    STATUS "Tests run as : ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} NUM_PROCS ${MPIEXEC_PREFLAGS_PRINT} EXECUTABLE"
  )

  set(HAVE_MPI 1)
  message(STATUS "MPI is enabled")
else(QMC_MPI)
  set(HAVE_MPI 0)
  message(STATUS "MPI is disabled")
endif(QMC_MPI)

####################################################################
#First check the required libraries. Abort if these are not found.
#Changed the fall-back versions to use similar conventions to 2.8.0 and higher
####################################################################
# use customized cmake file for boost so that we can deal with multiple
# installataions
# list of libraries to be linked with the main applications using I/O

#-------------------------------------------------------------------
# Check LAPACK/BLAS
#-------------------------------------------------------------------
# Add Math::BLAS_LAPACK target
add_library(Math::BLAS_LAPACK INTERFACE IMPORTED)

foreach(MKL_LIB_PATH IN ITEMS ${MKL_ROOT} $ENV{MKL_ROOT} $ENV{MKLROOT} $ENV{MKL_HOME})
  set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${MKL_LIB_PATH}/lib/intel64)
endforeach()
if(NOT BLA_VENDOR)
  message(STATUS "Trying to find LAPACK from Intel MKL")
  if(QMC_OMP AND CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(BLA_VENDOR Intel10_64lp)
  else()
    set(BLA_VENDOR Intel10_64lp_seq)
  endif()
  find_package(LAPACK)
  if(NOT LAPACK_FOUND)
    message(STATUS "Intel MKL library files not found via FindLAPACK.")
    message(STATUS "Trying to find alternative LAPACK libraries")
    set(BLA_VENDOR All)
    find_package(LAPACK REQUIRED)
  endif()
else()
  message(STATUS "Trying to find LAPACK library as requested BLA_VENDOR=${BLA_VENDOR}")
  find_package(LAPACK REQUIRED)
endif()

message(STATUS "LAPACK linker flags: ${LAPACK_LINKER_FLAGS}")
message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")
set_target_properties(Math::BLAS_LAPACK PROPERTIES INTERFACE_LINK_LIBRARIES
                                                   "${LAPACK_LINKER_FLAGS};${LAPACK_LIBRARIES}")

# Detects MKL header files and other components.
if(LAPACK_LIBRARIES MATCHES "mkl_core")
  include(CMake/FindMKL.cmake)
endif()

# AFQMC requires MKL sparse for good performance (roughly a factor of 2x)
if(BUILD_AFQMC AND NOT MKL_FOUND)
  message(
    WARNING
      "AFQMC - MKL not found, using simple sparse matrix routines.  Link with MKL sparse libraries for better performance."
  )
endif()

#-------------------------------------------------------------------
#  set up scalar/vector math libraries
#-------------------------------------------------------------------
set(ENABLE_MASS
    FALSE
    CACHE BOOL "ENABLE MASS scalar and vector math libraries for Power architecture")

# This needs to go before HAVE_SINCOS
if(ENABLE_MASS)
  include(CMake/FindIBMMASS.cmake)
endif(ENABLE_MASS)

if(NOT MASS_FOUND AND NOT HAVE_MKL_VML)
  message(STATUS "No usable vector math library detected.")
  add_library(Math::scalar_vector_functions INTERFACE IMPORTED)
endif(NOT MASS_FOUND AND NOT HAVE_MKL_VML)

# CheckSincos relies on SINCOS_INCLUDE which may be modified based on vendor libraries
if(NOT SINCOS_INCLUDE)
  set(SINCOS_INCLUDE cmath)
endif(NOT SINCOS_INCLUDE)
message(STATUS "SINCOS_INCLUDE : ${SINCOS_INCLUDE}")
include(${PROJECT_CMAKE}/CheckSincos.cmake)

#-------------------------------------------------------------------
#  set up FFTW3 library
#-------------------------------------------------------------------
set(FFTW_FOUND 0)
if(HAVE_MKL)
  if(HAVE_MKL_FFTW3)
    set(FFTW_FOUND 1)
  endif()
endif()
if(NOT FFTW_FOUND)
  find_package(FFTW)
endif()
if(FFTW_FOUND)
  set(HAVE_LIBFFTW 1)
  message(STATUS "Found FFTW")
else(FFTW_FOUND)
  #  MESSAGE(STATUS "Disabling FFTW")
  message(FATAL_ERROR "Could not find required library FFTW. Verify FFTW_HOME")
endif(FFTW_FOUND)

set(QMC_UTIL_LIBS ${QMC_UTIL_LIBS} ${FORTRAN_LIBRARIES})

#-------------------------------------------------------------------
#  set up Libxml2 library
#-------------------------------------------------------------------
find_package(ZLIB)
find_package(LibXml2 REQUIRED)

#-------------------------------------------------------------------
#  set up HDF5 library
#-------------------------------------------------------------------
if(HAVE_MPI)
  set(HDF5_PREFER_PARALLEL
      TRUE
      CACHE BOOL "Request parallel/serial HDF5 library")
else(HAVE_MPI)
  set(HDF5_PREFER_PARALLEL
      FALSE
      CACHE BOOL "Request parallel/serial HDF5 library")
  if(HDF5_PREFER_PARALLEL)
    message(FATAL_ERROR "Parallel HDF5 library cannot be selected with QMCPACK non-MPI build. "
                        "Please set HDF5_PREFER_PARALLEL=0.")
  endif(HDF5_PREFER_PARALLEL)
endif(HAVE_MPI)

if(QMC_BUILD_STATIC)
  message(STATUS "Linking static HDF5 library")
  set(HDF5_USE_STATIC_LIBRARIES on)
else()
  message(STATUS "Linking dynamic HDF5 library")
  set(HDF5_USE_STATIC_LIBRARIES off)
endif()

find_package(HDF5 COMPONENTS C)

if(HDF5_FOUND)
  set(HAVE_LIBHDF5 1)

  if(HDF5_IS_PARALLEL)
    message(STATUS "Parallel HDF5 library found")
    option(ENABLE_PHDF5 "Enable code paths using parallel HDF5" ON)
  else(HDF5_IS_PARALLEL)
    message(STATUS "Serial HDF5 library found")
    option(ENABLE_PHDF5 "Enable code paths using parallel HDF5" OFF)
    if(ENABLE_PHDF5)
      if(HAVE_MPI)
        message(FATAL_ERROR "Parallel HDF5 code paths requested but serial HDF5 library found! "
                            "Please either provide parallel HDF5 library or set ENABLE_PHDF5=0.")
      else(HAVE_MPI)
        message(FATAL_ERROR "Parallel HDF5 code paths cannot be enabled on non-MPI builds! Please set ENABLE_PHDF5=0.")
      endif(HAVE_MPI)
    endif(ENABLE_PHDF5)
  endif(HDF5_IS_PARALLEL)

  if(ENABLE_PHDF5)
    message(STATUS "Using HDF5 parallel collective I/O code paths")
  else(ENABLE_PHDF5)
    message(STATUS "Using HDF5 non-scalable serial I/O code paths")
  endif(ENABLE_PHDF5)

  if(HAVE_MPI AND NOT ENABLE_PHDF5)
    message(
      WARNING
        "MPI builds may have performance loss by not using parallel HDF5! (Safe to ignore for workstation builds).")
  endif()

  if(CMAKE_BUILD_TYPE AND HDF5_LIBRARIES_DEBUG)
    if(CMAKE_BUILD_TYPE MATCHES DEBUG)
      set(HDF5_LIBRARIES ${HDF5_LIBRARIES_DEBUG})
    else()
      set(HDF5_LIBRARIES ${HDF5_LIBRARIES_RELEASE})
    endif()
  endif()

  add_library(IO::HDF5 INTERFACE IMPORTED)
  set_target_properties(
    IO::HDF5
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIR}" INTERFACE_COMPILE_DEFINITIONS
                                                                   "HAVE_LIBHDF5;H5_USE_16_API" INTERFACE_LINK_LIBRARIES
                                                                                                "${HDF5_LIBRARIES}")
  if(ENABLE_PHDF5)
    set_property(TARGET IO::HDF5 APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS "ENABLE_PHDF5")
  endif(ENABLE_PHDF5)
else(HDF5_FOUND)
  message(FATAL_ERROR "HDF5 not found. Set HDF5_ROOT")
endif(HDF5_FOUND)

#-------------------------------------------------------------------
#  set up libboost library, header only
#-------------------------------------------------------------------
#make sure we can find boost if it's not in /usr
set(Boost_NO_BOOST_CMAKE on)
find_package(Boost 1.61.0 REQUIRED)
if(Boost_FOUND)
  set_property(TARGET Boost::boost APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS "HAVE_LIBBOOST")
  message(STATUS "Setting Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIRS}")
else()
  message(FATAL_ERROR "Need boost 1.61.0 or higher.  Set BOOST_ROOT")
endif()

#-------------------------------------------------------------------
#  set up CUDA compiler options and libraries
#-------------------------------------------------------------------
set(HAVE_CUDA 0)
if(QMC_CUDA OR ENABLE_CUDA)
  if(NOT QMC_CUDA2HIP)
    # FindCUDA default CUDA_PROPAGATE_HOST_FLAGS to ON but we prefer OFF
    # It happened -ffast-math from host caused numerical issue in CUDA kernels.
    option(CUDA_PROPAGATE_HOST_FLAGS "Propagate C/CXX_FLAGS and friends to the host compiler via -Xcompile" OFF)
    find_package(CUDA REQUIRED)
    set(CUDA_LINK_LIBRARIES_KEYWORD PRIVATE)
    #set(CUDA_NVCC_FLAGS
    #  "-arch=sm_20;-Drestrict=__restrict__;-DNO_CUDA_MAIN;-O3;-use_fast_math")
    if(CUDA_NVCC_FLAGS MATCHES "arch")
      # User defined NVCC flags
      message(STATUS "Setting CUDA FLAGS=${CUDA_NVCC_FLAGS}")
    else(CUDA_NVCC_FLAGS MATCHES "arch")
      # Automatically set the default NVCC flags
      set(CUDA_NVCC_FLAGS "-Drestrict=__restrict__;-DNO_CUDA_MAIN;-std=c++14")
      if(QMC_COMPLEX)
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-DQMC_COMPLEX")
      endif()
      if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-g;-G")
      else()
        # Temporarily disable fast_math because it causes multiple test failures
        # SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-O3;-use_fast_math")
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-O3")
      endif()
      set(CUDA_ARCH
          sm_70
          CACHE STRING "CUDA architecture sm_XX")
      set(CUDA_NVCC_FLAGS "-arch=${CUDA_ARCH};${CUDA_NVCC_FLAGS}")
    endif(CUDA_NVCC_FLAGS MATCHES "arch")
    include_directories(${CUDA_INCLUDE_DIRS})
    set(HAVE_CUDA 1)
    message("   CUDA_NVCC_FLAGS=${CUDA_NVCC_FLAGS}")
  else()
    set(ROCM_ROOT "/opt/rocm" CACHE PATH "Root directory of ROCM")
    message(STATUS "ROCM_ROOT: ${ROCM_ROOT}")
    list(APPEND CMAKE_MODULE_PATH ${ROCM_ROOT}/hip/cmake)
    list(APPEND CMAKE_PREFIX_PATH ${ROCM_ROOT})
    find_package(HIP REQUIRED)
    find_package(hipblas REQUIRED)
    # architecture flags
    set(HIP_ARCH
        "gfx906,gfx908"
        CACHE STRING "HIP architecture gfxXXX")
    list(APPEND HIP_HIPCC_FLAGS "-fPIC -ffast-math -O3")
    list(APPEND HIP_HIPCC_FLAGS "--amdgpu-target=${HIP_ARCH}")
    list(APPEND HIP_HIPCC_FLAGS "--gpu-max-threads-per-block=256")
    # warning suppression
    list(APPEND HIP_HIPCC_FLAGS "-Wno-vla")
    list(APPEND HIP_HIPCC_FLAGS "-Wno-deprecated-declarations")
    list(APPEND HIP_HIPCC_FLAGS "-Wno-unused-command-line-argument")
  endif(NOT QMC_CUDA2HIP)
else(QMC_CUDA OR ENABLE_CUDA)
  message(STATUS "CUDA disabled")
endif(QMC_CUDA OR ENABLE_CUDA)

if(USE_NVTX_API)
  message(STATUS "Enabling use of CUDA NVTX APIs")
  find_library(
    NVTX_API_LIB NAME nvToolsExt
    HINTS ${CUDA_TOOLKIT_ROOT_DIR}
    PATH_SUFFIXES lib lib64)
  if(NOT NVTX_API_LIB)
    message(FATAL_ERROR "USE_NVTX_API set but NVTX_API_LIB not found")
  endif(NOT NVTX_API_LIB)
  message("CUDA nvToolsExt library: ${NVTX_API_LIB}")
  link_libraries(${NVTX_API_LIB})
else()
  message(STATUS "CUDA NVTX APIs disabled")
endif(USE_NVTX_API)

#-------------------------------------------------------------------
#  set up ROCM compiler options and libraries
#-------------------------------------------------------------------
if(ENABLE_ROCM)
  message(STATUS "ROCM_ROOT: ${ROCM_ROOT}")
  add_library(ROCM::libraries INTERFACE IMPORTED)
  # temporarily put rocsolver rocrand here for convenience, should be moved to Platforms.
  set_target_properties(ROCM::libraries PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${ROCM_ROOT}/include"
                                                   INTERFACE_LINK_LIBRARIES "-L${ROCM_ROOT}/lib;-lrocsolver;-lrocrand")
endif(ENABLE_ROCM)

#-------------------------------------------------------------------
#  set up HIP compiler options
#-------------------------------------------------------------------
if(ENABLE_HIP)
  if(NOT ENABLE_ROCM)
    message(FATAL_ERROR "ROCM is required to use HIP. Please set ENABLE_ROCM=ON.")
  endif()
  set(CMAKE_MODULE_PATH "${ROCM_ROOT}/hip/cmake" ${CMAKE_MODULE_PATH})
  find_package(HIP REQUIRED)

  add_library(HIP::HIP INTERFACE IMPORTED)
  # temporarily put hipsparse hipblas here for convenience, should be moved to Platforms.
  set_target_properties(
    HIP::HIP PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${ROCM_ROOT}/include" INTERFACE_COMPILE_DEFINITIONS "ENABLE_HIP"
                        INTERFACE_LINK_LIBRARIES "-L${ROCM_ROOT}/lib;-lhipsparse;-lhipblas;-lamdhip64")
endif(ENABLE_HIP)

######################################################
# Add optional macros
######################################################

#include qmcpack/src build/src
include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_BINARY_DIR}/src)

if(USE_VTUNE_API)
  message(STATUS "Enabling use of VTune ittnotify APIs")
  find_path(
    VTUNE_ITTNOTIFY_INCLUDE_DIR ittnotify.h
    HINTS ${VTUNE_ROOT} $ENV{VTUNE_ROOT}
    PATH_SUFFIXES include REQUIRED)
  message(STATUS "Found VTUNE_ITTNOTIFY_INCLUDE_DIR ${VTUNE_ITTNOTIFY_INCLUDE_DIR}")
  find_library(
    VTUNE_ITTNOTIFY_LIBRARY ittnotify
    HINTS ${VTUNE_ROOT} $ENV{VTUNE_ROOT}
    PATH_SUFFIXES lib64 lib REQUIRED)
  message(STATUS "Found VTUNE_ITTNOTIFY_LIBRARY ${VTUNE_ITTNOTIFY_LIBRARY}")

  include_directories(${VTUNE_ITTNOTIFY_INCLUDE_DIR})
  link_libraries(${VTUNE_ITTNOTIFY_LIBRARY})
  if(USE_VTUNE_TASKS)
    message(STATUS "VTune ittnotify tasks enabled")
  endif()
else()
  message(STATUS "VTune ittnotify APIs disabled")
endif()

option(QMC_EXP_THREADING "Experimental non openmp threading models" OFF)
mark_as_advanced(QMC_EXP_THREADING)
if(QMC_EXP_THREADING)
  include(CheckAtomic)
  if(HAVE_LIBATOMIC)
    link_libraries("${HAVE_LIBATOMIC}")
  endif(HAVE_LIBATOMIC)
  add_definitions(-DQMC_EXP_THREADING)
endif(QMC_EXP_THREADING)

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
  message(STATUS "Linker is specified as '${CMAKE_CXX_LINKER}' instead of '${CMAKE_CXX_COMPILER}'")
endif()

###################################################
# Add user extra flags
###################################################
if(QMC_INCLUDE)
  foreach(tmp ${QMC_INCLUDE})
    message(STATUS "Adding '${tmp}' to include directories")
    include_directories("${tmp}")
  endforeach()
endif()
if(QMC_EXTRA_LIBS)
  string(REPLACE " " ";" QMC_EXTRA_LIBS ${QMC_EXTRA_LIBS})
  set(QMC_UTIL_LIBS ${QMC_UTIL_LIBS} ${QMC_EXTRA_LIBS})
  message("QMC_UTIL_LIBS=${QMC_UTIL_LIBS}")
endif()

# Avoid unused variables warnings (for now) by using for information purposes
if(QMC_DATA)
  message(STATUS "QMC_DATA is ${QMC_DATA}")
endif()
if(ENABLE_GCOV)
  message(STATUS "GCOV is enabled")
endif()

# SETUP SANITIZERS FLAGS
if(NOT "${ENABLE_SANITIZER}" STREQUAL "none")
  if(NOT ${COMPILER} MATCHES "GNU" AND NOT ${COMPILER} MATCHES "Clang")
    message(FATAL_ERROR "-DENABLE_SANITIZER=${ENABLE_SANITIZER} only works with GNU or Clang compilers")
  endif()

  if("${ENABLE_SANITIZER}" STREQUAL "asan")
    set(CMAKE_CXX_FLAGS_SAN
        "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer"
        CACHE STRING "AddressSanitizer C++ compiler builds." FORCE)
  elseif("${ENABLE_SANITIZER}" STREQUAL "ubsan")
    set(CMAKE_CXX_FLAGS_SAN
        "-fsanitize=undefined"
        CACHE STRING "UndefinedBehaviorSanitizer C++ compiler builds." FORCE)
  elseif("${ENABLE_SANITIZER}" STREQUAL "msan")
    set(CMAKE_CXX_FLAGS_SAN
        "-fsanitize=memory"
        CACHE STRING "MemorySanitizer C++ compiler builds." FORCE)
  elseif("${ENABLE_SANITIZER}" STREQUAL "tsan")
    set(CMAKE_CXX_FLAGS_SAN
        "-fsanitize=thread"
        CACHE STRING "ThreadSanitizer C++ compiler builds." FORCE)
  endif()

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_SAN}")

endif()

############################################################
# FLAGS at the project level
############################################################
#this might be redundant but maybe not in all CMake versions.
string(TOUPPER "${CMAKE_BUILD_TYPE}" THIS_CONFIG)

foreach(lang IN ITEMS C CXX)
  set(PROJECT_CMAKE_${lang}_FLAGS ${CMAKE_${lang}_FLAGS})
  #pre 3.0 cmake does not have string CONCAT
  set(TMP_PROJECT_CMAKE_${lang}_FLAGS "${PROJECT_CMAKE_${lang}_FLAGS} ${CMAKE_${lang}_FLAGS_${THIS_CONFIG}}")
  message("Project ${lang}_FLAGS: ${TMP_PROJECT_CMAKE_${lang}_FLAGS}")
endforeach(lang IN ITEMS C CXX)

get_directory_property(TMP_PROJECT_INCLUDE_DIRECTORIES INCLUDE_DIRECTORIES)
message("Project INCLUDE_DIRECTORIES: ${TMP_PROJECT_INCLUDE_DIRECTORIES}")

foreach(target_type IN ITEMS EXE SHARED)
  set(PROJECT_CMAKE_${target_type}_LINKER_FLAGS ${CMAKE_${target_type}_LINKER_FLAGS})
  set(TMP_PROJECT_CMAKE_${target_type}_LINKER_FLAGS
      "${PROJECT_CMAKE_${target_type}_LINKER_FLAGS} ${CMAKE_${target_type}_LINKER_FLAGS_${THIS_CONFIG}}")
  message("Project ${target_type}_LINKER_FLAGS: ${TMP_PROJECT_CMAKE_${target_type}_LINKER_FLAGS}")
endforeach(target_type IN ITEMS EXE SHARED)

###################################################
# Done searching libraries. Start building qmcpack
###################################################

message(STATUS "Ready to build qmcpack")

####################################################################
# add definitions
####################################################################

add_definitions(-DHAVE_CONFIG_H -DUSE_REAL_STRUCT_FACTOR -DADD_)

execute_process(
  COMMAND hostname
  OUTPUT_VARIABLE QMC_HOSTNAME
  OUTPUT_STRIP_TRAILING_WHITESPACE)
configure_file(${qmcpack_SOURCE_DIR}/src/qmcpack.settings ${qmcpack_BINARY_DIR}/bin/qmcpack.settings)
install(FILES ${qmcpack_BINARY_DIR}/bin/qmcpack.settings DESTINATION bin)
if(INSTALL_NEXUS)
  install(CODE "EXECUTE_PROCESS(COMMAND ${qmcpack_SOURCE_DIR}/nexus/install --leave_paths ${CMAKE_INSTALL_PREFIX}/bin)")
endif()

add_subdirectory(src)

if(NOT QMC_BUILD_SANDBOX_ONLY)
  if(NOT QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)
    add_subdirectory(tests)
    add_subdirectory(examples)
  endif()
  if(IS_DIRECTORY "${PROJECT_SOURCE_DIR}/nexus" AND INSTALL_NEXUS)
    add_subdirectory(nexus/tests)
  endif()
endif()
