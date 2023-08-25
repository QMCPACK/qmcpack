######################################################################
# CMake version and policies
######################################################################
cmake_minimum_required(VERSION 3.17.0)
# Note that cmake_minimum_required affects policy defaults.
# All policies known to the running version of CMake and introduced in
# cmake_minimum_required version or earlier will be set to use NEW behavior

if(POLICY CMP0127)
  cmake_policy(SET CMP0127 NEW) # Condition evaluation v3.22+ https://cmake.org/cmake/help/latest/policy/CMP0127.html
endif()

######################################################################
# QMCPACK project
######################################################################
project(
  qmcpack
  VERSION 3.17.1
  LANGUAGES C CXX)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#--------------------------------------------------------------------
# Directory where customize cmake files reside
#--------------------------------------------------------------------
set(PROJECT_CMAKE ${qmcpack_SOURCE_DIR}/CMake)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_CMAKE})
include(macros)

######################################################################
# Build and install options
######################################################################

#--------------------------------------------------------------------
# Set C++ standard
#--------------------------------------------------------------------
set(QMC_CXX_STANDARD
    17
    CACHE STRING "QMCPACK C++ language standard")
if(NOT QMC_CXX_STANDARD EQUAL 17)
  message(WARNING "C++17 is the only language standard officially supported by this QMCPACK version. "
                  "Using other versions of the C++ standard is unsupported and done entirely at user's own risk.")
endif()

#--------------------------------------------------------------------
# Programmind model related build options
# MPI, OpenMP, GPU acceleration
#--------------------------------------------------------------------
include(CMakeDependentOption)
option(QMC_MPI "Enable/disable MPI" ON)
option(QMC_OMP "Enable/disable OpenMP" ON)
option(QMC_COMPLEX "Build for complex binary" OFF)
option(ENABLE_CUDA "Build with GPU support through CUDA" OFF)
option(QMC_CUDA2HIP "Map all CUDA kernels and library calls to HIP" OFF)
if(QMC_CUDA2HIP OR ENABLE_HIP)
  set(ENABLE_ROCM ON) # option(ENABLE_ROCM) will be no-op
endif()
# explicit HIP source codes only exist in AFQMC.
option(ENABLE_HIP "Build with with GPU support through explicit HIP source code" OFF)
option(ENABLE_ROCM "Build with with GPU support through ROCM libraries" OFF)
option(ENABLE_OFFLOAD "Enable OpenMP offload" OFF)
option(ENABLE_SYCL "Enable SYCL offload" OFF)
# Use CMake object library targets to workaround clang linker not being able to handle fat
# binary archives which contain both host and device codes, for example OpenMP offload regions.
# CMake does not propagate indirect object files by design.
# So the dependency must be explicitly specified with USE_OBJECT_TARGET.
cmake_dependent_option(USE_OBJECT_TARGET "Use CMake object library target" OFF "ENABLE_OFFLOAD" OFF)

if(QMC_CUDA)
  message(
    FATAL_ERROR
      "Legacy GPU implementation using CUDA has been removed! QMC_CUDA cannot be set ON. "
      "For NVIDIA GPU support use ENABLE_OFFLOAD with ENABLE_CUDA. Use v3.16.0 or earlier if you need QMC_CUDA.")
endif(QMC_CUDA)

if(ENABLE_CUDA AND ENABLE_SYCL)
  message(FATAL_ERROR "ENABLE_CUDA=ON and ENABLE_SYCL=ON can not be set together!")
endif(ENABLE_CUDA AND ENABLE_SYCL)

if(ENABLE_CUDA)
  if(QMC_CUDA2HIP)
    message(STATUS "ENABLE_CUDA enabled, QMC_CUDA2HIP enabled")
    if(CMAKE_VERSION VERSION_LESS 3.21.0)
      message(FATAL_ERROR "ENABLE_ROCM or QMC_CUDA2HIP require CMake 3.21.0 or later")
    endif()
  else(QMC_CUDA2HIP)
    message(STATUS "ENABLE_CUDA enabled, QMC_CUDA2HIP disabled")
    if(CMAKE_VERSION VERSION_LESS 3.18.0)
      message(FATAL_ERROR "ENABLE_CUDA require CMake 3.18.0 or later")
    endif()
  endif()
else(ENABLE_CUDA)
  if(QMC_CUDA2HIP)
    message(FATAL_ERROR "QMC_CUDA2HIP requires ENABLE_CUDA=ON.")
  endif()
  message(STATUS "ENABLE_CUDA disabled")
endif(ENABLE_CUDA)

# stop CUDA_ARCH
if(ENABLE_CUDA AND NOT QMC_CUDA2HIP)
  if(DEFINED CUDA_ARCH)
    unset(CUDA_ARCH CACHE)
    message(
      FATAL_ERROR "CUDA_ARCH option has been removed. Use -DCMAKE_CUDA_ARCHITECTURES=80 if -DCUDA_ARCH=sm_80 was used.")
  endif()
endif()

include(DetermineDeviceArchitectures)
message(STATUS "GPU device architectures: ${QMC_GPU_ARCHS}")

#--------------------------------------------------------------------
# Set compiler-time parameters
# WALKER_MAX_PROPERTIES max number of observables + 12 or so standard
#   properties. Properties such forces are per particle so you may need this
#   to be rather large.  Each property increases walker size by
#   sizeof(FULLPRECREALTYPE)
# OHMMS_DIM =  dimension of the problem
# OHMMS_INDEXTYPE = type of index
# OHMMS_PRECISION  = base precision, float, double etc
# OHMMS_PRECISION_FULL  = full precision, double etc
#--------------------------------------------------------------------
set(WALKER_MAX_PROPERTIES
    2048
    CACHE STRING "Maximum number of properties tracked by walkers")
mark_as_advanced(WALKER_MAX_PROPERTIES)
set(OHMMS_DIM
    3
    CACHE STRING "Select physical dimension")
set(OHMMS_INDEXTYPE int)
set(OHMMS_PRECISION_FULL double)

if(QMC_VERBOSE_CONFIGURATION)
  message(WARNING "QMC_VERBOSE_CONFIGURATION option has been removed. Use --log-level=VERBOSE cmake option instead.")
endif()

#--------------------------------------------------------------------
message(STATUS "Defining the float point precision")
#--------------------------------------------------------------------
option(QMC_MIXED_PRECISION "Enable/disable mixed precision" OFF)
if(QMC_MIXED_PRECISION)
  set(OHMMS_PRECISION float)
  set(MIXED_PRECISION ON)
else(QMC_MIXED_PRECISION)
  set(OHMMS_PRECISION double)
endif(QMC_MIXED_PRECISION)
message("   Base precision = ${OHMMS_PRECISION}")
message("   Full precision = ${OHMMS_PRECISION_FULL}")

#--------------------------------------------------------------------
# Set debug printout
# DEBUG_PSIBUFFER_ON keep track buffer's locator
#--------------------------------------------------------------------
option(PRINT_DEBUG "Enable/disable debug printing" OFF)
option(DEBUG_PSIBUFFER_ON "Debug messages for PSIBUFFER locator" OFF)
mark_as_advanced(DEBUG_PSIBUFFER_ON)
option(DEBUG_PER_STEP_ACCEPT_REJECT "Print accepts and rejects at each step" OFF)
mark_as_advanced(DEBUG_PER_STEP_ACCEPT_REJECT)

#--------------------------------------------------------------------
# Set build feature options
#--------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_BUILD_TYPE Release)
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)
message(STATUS "CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")

option(QMC_BUILD_SANDBOX_ONLY "Build only applications in Sandbox directory" OFF)
option(ENABLE_GCOV "Enable code coverage" OFF)
option(BUILD_UNIT_TESTS "Build unit tests" ON)
option(BUILD_MICRO_BENCHMARKS "Build micro benchmarks" ON)
option(BUILD_LMYENGINE_INTERFACE "Build LMY engine" ON)
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
cmake_dependent_option(USE_NVTX_API "Enable/disable NVTX regions in CUDA code." OFF "ENABLE_TIMERS AND ENABLE_CUDA" OFF)
set(HAVE_EINSPLINE 1) # to be removed
option(QMC_EXP_THREADING "Experimental non openmp threading models" OFF)
mark_as_advanced(QMC_EXP_THREADING)

#--------------------------------------------------------------------
# Sanitizer options
#--------------------------------------------------------------------

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

#--------------------------------------------------------------------
# INSTALL_NEXUS install Nexus alongside QMCPACK
#--------------------------------------------------------------------
option(INSTALL_NEXUS "Install Nexus alongside QMCPACK" ON)
if(INSTALL_NEXUS)
  install(
    CODE "EXECUTE_PROCESS(COMMAND ${qmcpack_SOURCE_DIR}/nexus/install --leave_paths \$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/bin)"
  )
endif()

######################################################################
# Set compiler specific options/flags
######################################################################
#-------------------------------------------------------------------------------
# Check for OpenMP compiler support
#
# The CMake OpenMP package is *only* used to to verify compiler support.
# No feature of the package is used and we customize OpenMP options directly
# on a per vendor/compiler basis.
#-------------------------------------------------------------------------------
if(QMC_OMP)
  find_package(OpenMP)
  if(NOT OpenMP_FOUND)
    message(
      FATAL_ERROR "No compiler support for OpenMP found. Switching to a compiler with OpenMP support is recommended."
                  "Alternatively, you will need to run CMake configure with -DQMC_OMP=OFF")
  endif()
endif(QMC_OMP)

#-------------------------------------------------------------------------------
# Set vendor specific compiler options
#-------------------------------------------------------------------------------
if(CMAKE_TOOLCHAIN_FILE)
  message(STATUS "Using ${CMAKE_TOOLCHAIN_FILE} toolchain ")
else(CMAKE_TOOLCHAIN_FILE)
  message(STATUS "Trying to figure out compiler options ....")

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
                      "Setting environment variable CRAYPE_LINK_TYPE to \"dynamic\" enables dynamic linking.\n")
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

  include(inspectCompiler)

endif(CMAKE_TOOLCHAIN_FILE)

if(NOT CMAKE_CXX_COMPILER_ID STREQUAL CMAKE_C_COMPILER_ID)
  message(FATAL_ERROR "Mixing C and C++ compilers from different vendors is not permitted. Current "
                      "C compiler is ${CMAKE_C_COMPILER_ID} but C++ compiler is ${CMAKE_CXX_COMPILER_ID}.")
endif()

#-----------------------------------------------------------------------
# Set C++ minimum standard and run basic checks
#-----------------------------------------------------------------------
set(CMAKE_CXX_STANDARD ${QMC_CXX_STANDARD})
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Check that a C++ compiler is compatible with the underlying libstdc++
include(Testlibstdc++)

# Check that a C++ standard library is configured properly
include(TestCxx17Library)

#-----------------------------------------------------------------------
# SETUP SANITIZERS FLAGS
#-----------------------------------------------------------------------
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

#-----------------------------------------------------------------------
# verify code coverage compiler support
#-----------------------------------------------------------------------
if(ENABLE_GCOV)
  if(GCOV_SUPPORTED)
    message(STATUS "GCOV is enabled")
  else()
    message(FATAL_ERROR "ENABLE_GCOV was specified but compiler does not support GCC coverage flag")
  endif()
endif(ENABLE_GCOV)

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
include(CheckSIMDAlignment)
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
include(TestOpenMPtaskloop)
option(ENABLE_OMP_TASKLOOP "Enable OpenMP taskloop" ${OMP_TASKLOOP_OKAY})
message(STATUS "ENABLE_OMP_TASKLOOP is set to ${ENABLE_OMP_TASKLOOP}")

#---------------------------------------------------------
# Set up OpenMP offload compile options
#---------------------------------------------------------
set(QMC_OFFLOAD_MEM_ASSOCIATED_DEFAULT OFF)
if(ENABLE_OFFLOAD AND DEFINED OPENMP_OFFLOAD_COMPILE_OPTIONS)
  message(STATUS "OpenMP offload CXX flags: ${OPENMP_OFFLOAD_COMPILE_OPTIONS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMP_OFFLOAD_COMPILE_OPTIONS}")
  if(${COMPILER} MATCHES "Clang"
     AND OPENMP_OFFLOAD_COMPILE_OPTIONS MATCHES "gfx"
     AND QMC_CUDA2HIP)
    # As of 11/2021, QMC_OFFLOAD_MEM_ASSOCIATED=ON is needed for AMD and mainline LLVM compilers
    # when using OpenMP offload to AMD GPU together with HIP.
    set(QMC_OFFLOAD_MEM_ASSOCIATED_DEFAULT ON)
  endif()
endif()
# Some OpenMP offload runtime libraries have composibility issue with a vendor native runtime.
# A workaround is making the vendor native runtime responsible for memory allocations and OpenMP associate/disassocate them.
cmake_dependent_option(QMC_OFFLOAD_MEM_ASSOCIATED "Manage OpenMP memory allocations via the vendor runtime"
                       ${QMC_OFFLOAD_MEM_ASSOCIATED_DEFAULT} "ENABLE_OFFLOAD;ENABLE_CUDA" OFF)

#-------------------------------------------------------------------------------------
# consider making this always on if OpenMP is no longer UB with Thread Support Library
# We should test out abstractions consistently and c++17 stdlib's all have support for
# std::thread and std::atomic
#-------------------------------------------------------------------------------------
if(QMC_EXP_THREADING)
  include(CheckAtomic)
  if(HAVE_LIBATOMIC)
    link_libraries("${HAVE_LIBATOMIC}")
  endif(HAVE_LIBATOMIC)
  add_definitions(-DQMC_EXP_THREADING)
endif(QMC_EXP_THREADING)

#-------------------------------------------------------------------
# Check MPI installation. MPI is treated as a part of compiler
#-------------------------------------------------------------------
if(QMC_MPI)
  # for backward compatibility with MPIEXEC
  if(DEFINED MPIEXEC AND NOT DEFINED MPIEXEC_EXECUTABLE)
    message(
      WARNING
        "Setting MPIEXEC_EXECUTABLE based on MPIEXEC. MPIEXEC has been deprecated! Use MPIEXEC_EXECUTABLE instead!")
    set(MPIEXEC_EXECUTABLE ${MPIEXEC})
  endif(DEFINED MPIEXEC AND NOT DEFINED MPIEXEC_EXECUTABLE)
  set(MPI_CXX_SKIP_MPICXX
      TRUE
      CACHE BOOL "Must be TRUE for QMCPACK. Cached for working around CMake FindMPI issue.")
  if(NOT CMAKE_CROSSCOMPILING)
    set(MPI_DETERMINE_LIBRARY_VERSION TRUE)
  endif()
  find_package(MPI COMPONENTS CXX)

  if(NOT MPI_FOUND)
    message(
      FATAL_ERROR
        "MPI support not found! Provide MPI compiler wrappers or build without MPI by passing '-DQMC_MPI=OFF' to cmake."
    )
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
        "To build without MPI, pass '-DQMC_MPI=OFF' to cmake")
    message(WARNING ${MPI_WARNING_LIST})
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

#--------------------------------------------------------------
# final test for compile options before searching for libraries
#--------------------------------------------------------------
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("int main(){}" TEST_CXX_COMPILE_MAIN)
if(NOT TEST_CXX_COMPILE_MAIN)
  unset(TEST_CXX_COMPILE_MAIN CACHE)
  message(FATAL_ERROR "Failed in compiling a main() function likely due to incorrect compile options. "
                      "Check error in \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeError.log\".")
endif()

#--------------------------------------------------------------
# Workaround for breakage of OMP Target kernels at -O0 by Clang
#--------------------------------------------------------------
# When the build type is Debug, default -O0 compilation is broken by
# Clang 15 as of 5d2ce7663b10c107328a4ae0c678165209e64619.
# Previous compilers are not suggested for clang offload builds.
#
# You can set this option to false on the command line to check
# if this problem has been fixed.

if((${COMPILER} MATCHES "Clang") AND ENABLE_OFFLOAD)
  option(ENABLE_OFFLOAD_CLANG_DEBUG_O3 "build OMP target kernels with -O3 in the build type Debug" ON)
endif()
mark_as_advanced(ENABLE_OFFLOAD_CLANG_DEBUG_O3)

######################################################################
# Check external libraries.
######################################################################

#-------------------------------------------------------------------
# check OS related libraries
#-------------------------------------------------------------------
include(CheckFunctionExists)
check_function_exists(posix_memalign HAVE_POSIX_MEMALIGN)

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
  if(QMC_OMP AND COMPILER STREQUAL "Intel")
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
target_link_libraries(Math::BLAS_LAPACK INTERFACE "${LAPACK_LINKER_FLAGS};${LAPACK_LIBRARIES}")

# Detects MKL header files and other components.
if(LAPACK_LIBRARIES MATCHES "mkl_core")
  include(FindMKL)
endif()

if(MKL_FOUND)
  target_link_libraries(Math::BLAS_LAPACK INTERFACE MKL::BLAS_LAPACK)
  set(MATH_VENDOR_GUESS "INTEL_VML")
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
add_library(Math::scalar_vector_functions INTERFACE IMPORTED)

set(VALID_MATH_VENDORS "GENERIC" "INTEL_VML" "IBM_MASS" "AMD_LIBM")
if(NOT DEFINED MATH_VENDOR_GUESS)
  set(MATH_VENDOR_GUESS "GENERIC")
endif()
set(QMC_MATH_VENDOR
    ${MATH_VENDOR_GUESS}
    CACHE STRING "Vendor optimized libraries for scalar and vector math functions")

# Perform math vendor option check
if(NOT QMC_MATH_VENDOR IN_LIST VALID_MATH_VENDORS)
  message(FATAL_ERROR "Invalid vendor math library ${QMC_MATH_VENDOR}, value must be one of ${VALID_MATH_VENDORS}")
else()
  message(STATUS "Selected vendor math library ${QMC_MATH_VENDOR}")
endif()

# This needs to go before HAVE_SINCOS
if(QMC_MATH_VENDOR STREQUAL "INTEL_VML")
  include(FindMKLVML)
elseif(QMC_MATH_VENDOR STREQUAL "IBM_MASS")
  include(FindIBMMASS)
elseif(QMC_MATH_VENDOR STREQUAL "AMD_LIBM")
  include(FindAMDLibM)
endif()

# CheckSincos relies on SINCOS_INCLUDE which may be modified based on vendor libraries
if(NOT SINCOS_INCLUDE)
  set(SINCOS_INCLUDE cmath)
endif(NOT SINCOS_INCLUDE)
message(STATUS "SINCOS_INCLUDE : ${SINCOS_INCLUDE}")
include(CheckSincos)

#-------------------------------------------------------------------
#  set up FFTW3 library
#-------------------------------------------------------------------
add_library(Math::FFTW3 INTERFACE IMPORTED)
set(FFTW_FOUND 0)
if(HAVE_MKL)
  if(HAVE_MKL_FFTW3)
    set(FFTW_FOUND 1)
    target_link_libraries(Math::FFTW3 INTERFACE MKL::FFTW3)
  endif()
endif()
if(NOT FFTW_FOUND)
  find_package(FFTW REQUIRED)
  target_link_libraries(Math::FFTW3 INTERFACE FFTW::FFTW3)
endif()

set(HAVE_LIBFFTW 1)
target_compile_definitions(Math::FFTW3 INTERFACE "HAVE_LIBFFTW")

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
  option(HDF5_PREFER_PARALLEL "Request parallel/serial HDF5 library" ON)
else(HAVE_MPI)
  option(HDF5_PREFER_PARALLEL "Request parallel/serial HDF5 library" OFF)
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

find_package(HDF5 COMPONENTS C) # Note: minimum version check is done below to bypass find_package
                                # and HDF5 version compatibility subtleties

if(HDF5_FOUND)
  if(HDF5_VERSION)
    if (HDF5_VERSION VERSION_LESS 1.10.0)
      message(FATAL_ERROR "QMCPACK requires HDF5 version >= 1.10.0")
    endif()
  endif(HDF5_VERSION)
  if(HDF5_IS_PARALLEL)
    if(HAVE_MPI)
      message(STATUS "Parallel HDF5 library found")
      option(ENABLE_PHDF5 "Enable code paths using parallel HDF5" ON)
    else(HAVE_MPI)
      message(FATAL_ERROR "Parallel HDF5 library found but cannot be used with QMCPACK non-MPI build. "
                          "Please provide a serial HDF5 library or switch to building QMCPACK with MPI.")
    endif(HAVE_MPI)
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
  target_include_directories(IO::HDF5 INTERFACE "${HDF5_INCLUDE_DIR}")
  target_compile_definitions(IO::HDF5 INTERFACE "H5_USE_110_API")
  target_link_libraries(IO::HDF5 INTERFACE "${HDF5_LIBRARIES}")
  if(ENABLE_PHDF5)
    target_compile_definitions(IO::HDF5 INTERFACE "ENABLE_PHDF5")
  endif(ENABLE_PHDF5)
else(HDF5_FOUND)
  message(FATAL_ERROR "HDF5 not found. Set HDF5_ROOT")
endif(HDF5_FOUND)

#-------------------------------------------------------------------
#  set up libboost library, header only
#-------------------------------------------------------------------
set(Boost_NO_BOOST_CMAKE on)
find_package(Boost 1.61.0 REQUIRED)
if(Boost_FOUND)
  target_compile_definitions(Boost::boost INTERFACE "HAVE_LIBBOOST")
  message(STATUS "Setting Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIRS}")
else()
  message(FATAL_ERROR "Need boost 1.61.0 or higher.  Set BOOST_ROOT")
endif()

#-------------------------------------------------------------------
#  set up CUDA compiler options and libraries
#-------------------------------------------------------------------
set(HAVE_CUDA 0)
if(ENABLE_CUDA AND NOT QMC_CUDA2HIP)
  if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD 17)
  endif()
  set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)
  set(CMAKE_CUDA_EXTENSIONS OFF)
  # CMAKE_CUDA_STANDARD doesn't affect nvcc detection and using unsupported compilers may choke the detection.
  # This happens when Clang default to C++17. Setting CMAKE_CUDA_FLAGS directly works around the issue.
  # but CUDA_FLAGS in flags.make may contain -std=c++17 twice with no harm.
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++${CMAKE_CUDA_STANDARD}")

  set(CMAKE_CUDA_HOST_COMPILER
      ${CMAKE_CXX_COMPILER}
      CACHE STRING "nvcc host compiler passed via -ccbin")
  if(NOT CMAKE_CUDA_FLAGS MATCHES "allow-unsupported-compiler")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --allow-unsupported-compiler")
  endif()

  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    string(REPLACE "sm_" "" CUDA_ARCH_NUMBERS "${QMC_GPU_ARCHS}")
    set(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_NUMBERS})
  endif()

  enable_language(CUDA)
  if(ENABLE_CUDA)
    include(TestCUDAHostCompatibility)
  endif()
  find_package(CUDAToolkit 11.0 REQUIRED)
  if(NOT TARGET CUDA::cublas)
    message(
      FATAL_ERROR
        "Found an incomplete CUDA toolkit installation, target CUDA::cublas not found. "
        "This often happens when CMake failed in recognizing the NVHPC internal CUDA toolkit with libcublas. "
        "Set CMAKE_CUDA_COMPILER to the full path of nvcc from a complete CUDA toolkit installation.")
  endif()
  # Automatically set the default NVCC flags
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Drestrict=__restrict__ -DNO_CUDA_MAIN")
  if(QMC_COMPLEX)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -DQMC_COMPLEX=${QMC_COMPLEX}")
  endif()
  set(HAVE_CUDA 1)
  message("Project CUDA_FLAGS: ${CMAKE_CUDA_FLAGS}")
endif()

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
  enable_language(HIP)
  set(CMAKE_HIP_LINKER_PREFERENCE_PROPAGATES 0)
  if(DEFINED ROCM_ROOT)
    message(STATUS "Using ROCM_ROOT: ${ROCM_ROOT}")
    # explicit ROCM_ROOT takes precedence
    list(PREPEND CMAKE_PREFIX_PATH ${ROCM_ROOT})
  elseif(CMAKE_HIP_COMPILER_ROCM_ROOT)
    list(PREPEND CMAKE_PREFIX_PATH ${CMAKE_HIP_COMPILER_ROCM_ROOT})
  endif()

  # use hip package to supply the hip runtime for CXX source code.
  find_package(hip CONFIG)
  if(NOT hip_FOUND)
    message(
      FATAL_ERROR
        "The HIP runtime library cmake configure file hip-config.cmake cannot be found! "
        "This happens when the ROCm installation cannot be found by CMake. "
        "Rerun cmake with -DROCM_ROOT=<your ROCm installation root directory> added.")
  endif()
  find_package(hipblas CONFIG REQUIRED)
  find_package(rocsolver CONFIG REQUIRED)
  find_package(rocthrust CONFIG REQUIRED)

  # architecture flags
  if(DEFINED HIP_ARCH)
    unset(HIP_ARCH CACHE)
    message(
      FATAL_ERROR
        "HIP_ARCH option has been removed. Use -DCMAKE_HIP_ARCHITECTURES=gfx90a if -DHIP_ARCH=gfx90a was used.")
  endif()
  if(NOT DEFINED CMAKE_HIP_ARCHITECTURES)
    set(CMAKE_HIP_ARCHITECTURES ${QMC_GPU_ARCHS})
  endif()

  set(CMAKE_HIP_FLAGS "${CMAKE_HIP_FLAGS} -fPIC -ffast-math --gpu-max-threads-per-block=256")
  # warning suppression
  set(CMAKE_HIP_FLAGS "${CMAKE_HIP_FLAGS} -Wno-vla -Wno-deprecated-declarations -Wno-unused-command-line-argument")

  if(NOT DEFINED ROCM_ROOT)
    list(PREPEND CMAKE_PREFIX_PATH ${CMAKE_HIP_COMPILER_ROCM_ROOT})
  endif()

  if(BUILD_AFQMC)
    find_package(rocrand CONFIG REQUIRED)
    # connect ROCm shipped libraries
    add_library(ROCM::libraries INTERFACE IMPORTED)
    # temporarily put rocsolver rocrand here for convenience, should be moved to Platforms.
    target_link_libraries(ROCM::libraries INTERFACE roc::rocsolver roc::rocrand)
  endif()
  message("Project HIP_FLAGS: ${CMAKE_HIP_FLAGS}")
  option(QMC_DISABLE_HIP_HOST_REGISTER "Disable hipHostRegister for pinning host memory" ON)
endif(ENABLE_ROCM)

#-------------------------------------------------------------------
#  set up HIP compiler options
#-------------------------------------------------------------------
if(ENABLE_HIP)
  if(NOT ENABLE_ROCM)
    message(FATAL_ERROR "ROCM is required to use HIP. Please set ENABLE_ROCM=ON.")
  endif()

  find_package(hipsparse CONFIG REQUIRED)
  add_library(HIP::HIP INTERFACE IMPORTED)
  # temporarily put hipsparse hipblas here for convenience, should be moved to Platforms.
  target_link_libraries(HIP::HIP INTERFACE roc::hipsparse roc::hipblas)
  target_compile_definitions(HIP::HIP INTERFACE "ENABLE_HIP")

  # FindHIP sets HIP_PLATFORM to amd or nvcc
  message(STATUS "FindHIP determined HIP_PLATFORM is ${HIP_PLATFORM}.")
  if(HIP_PLATFORM STREQUAL "amd")
    target_compile_definitions(HIP::HIP INTERFACE "__HIP_PLATFORM_AMD__")
  elseif(HIP_PLATFORM STREQUAL "nvcc")
    target_compile_definitions(HIP::HIP INTERFACE "__HIP_PLATFORM_NVIDIA__")
  else()
    message(FATAL_ERROR "Unknown HIP platform ${HIP_PLATFORM} only support 'amd' or 'nvcc'")
  endif()
endif(ENABLE_HIP)

#-------------------------------------------------------------------
#  set up SYCL compiler options and libraries
#-------------------------------------------------------------------
if(ENABLE_SYCL)
  if(NOT "x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xIntelLLVM")
    message(FATAL_ERROR "QMCPACK only supports SYCL with LLVM-based Intel oneAPI compiler (icpx).")
  endif()
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 2023.0)
    message(FATAL_ERROR "Requires Intel oneAPI 2023.0 or higher!")
  endif()
  add_library(SYCL::host INTERFACE IMPORTED)
  add_library(SYCL::device INTERFACE IMPORTED)
  find_package(
    IntelDPCPP
    REQUIRED
    CONFIGS
    IntelDPCPPConfig-modified.cmake
    PATHS
    ${PROJECT_CMAKE}
    NO_DEFAULT_PATH)
  target_link_libraries(SYCL::host INTERFACE OneAPI::DPCPP-host)
  target_link_libraries(SYCL::device INTERFACE OneAPI::DPCPP-device)
  if(TARGET MKL::sycl)
    target_link_libraries(MKL::sycl INTERFACE OneAPI::DPCPP-host)
  endif()
endif(ENABLE_SYCL)

#-------------------------------------------------------------------
# set up VTune ittnotify library
#-------------------------------------------------------------------

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

#-------------------------------------------------------------------
# Add user extra flags
#-------------------------------------------------------------------
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

if(BUILD_AFQMC AND NOT APPLE)
  link_libraries("rt")
endif()

######################################################################
# All the settings related to building QMCPACK are complete
######################################################################

#include build/src
include_directories(${qmcpack_BINARY_DIR}/src)
#include qmcpack/src
include_directories(${PROJECT_SOURCE_DIR}/src)

#-------------------------------------------------------------------
# add definitions DO NOT EXPAND
#-------------------------------------------------------------------
add_definitions(-DHAVE_CONFIG_H -DADD_)

#-------------------------------------------------------------------
# qmcpack.settings
#-------------------------------------------------------------------
execute_process(
  COMMAND hostname
  OUTPUT_VARIABLE QMC_HOSTNAME
  OUTPUT_STRIP_TRAILING_WHITESPACE)
configure_file(${qmcpack_SOURCE_DIR}/src/qmcpack.settings ${qmcpack_BINARY_DIR}/bin/qmcpack.settings)
install(FILES ${qmcpack_BINARY_DIR}/bin/qmcpack.settings DESTINATION bin)

#-------------------------------------------------------------------
# FLAGS at the project level
#-------------------------------------------------------------------
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

#-------------------------------------------------------------------
# Git information (see src/CMakeLists.txt for the
# command to get repository information)
#-------------------------------------------------------------------
find_package(Git)

if(GIT_FOUND AND EXISTS "${qmcpack_SOURCE_DIR}/.git")
  set(IS_GIT_PROJECT 1)
endif()

# set up config.h
configure_file(${qmcpack_SOURCE_DIR}/src/config.h.cmake.in ${qmcpack_BINARY_DIR}/src/config.h)

######################################################################
# Testing related settings
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

if(NOT TEST_MAX_PROCS)
  set(TEST_MAX_PROCS 100)
endif()

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

if(QMC_DATA)
  message(STATUS "QMC_DATA is ${QMC_DATA}")
endif()

#-------------------------------------------------------------------
# CTest
#-------------------------------------------------------------------
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

#-------------------------------------------------------------------
# Verify Python3 available
#-------------------------------------------------------------------
find_package(Python3)
if(NOT Python3_FOUND)
  message(FATAL_ERROR "Could not find required python3")
endif(NOT Python3_FOUND)
include(python)

#-------------------------------------------------------------------
# Check if PySCF is available
#-------------------------------------------------------------------
if(NOT QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)
  check_python_reqs(pyscf "" HAVE_PYSCF)
  if(NOT ${HAVE_PYSCF})
    message(STATUS "Unable to import PySCF python module. PySCF tests will not be run.")
  else()
    message(STATUS "Successfully imported PySCF python module.")
  endif()
endif()

#-------------------------------------------------------------------
# Check if GPAW is available
#-------------------------------------------------------------------
if(NOT QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)
  check_python_reqs(gpaw "" HAVE_GPAW)
  if(NOT ${HAVE_GPAW})
    message(STATUS "Unable to import GPAW python module. GPAW converter tests will not be run.")
  else()
    message(STATUS "Successfully imported GPAW python modules.")
  endif()
endif()

#-------------------------------------------------------------------
# Check if QE executables available at QE_BIN or on the PATH
#-------------------------------------------------------------------
if(NOT QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)
  include("CMake/FindQE.cmake")
  if(QE_FOUND)
    message(STATUS "Found and using patched Quantum ESPRESSO (QE) with pw2qmcpack.x at ${QE_PW_DIR}")
  else(QE_FOUND)
    message(
      STATUS "Did not find a patched Quantum ESPRESSO (QE) distribution with pw2qmcpack.x. QE tests will not be run.")
  endif(QE_FOUND)
endif()

#-------------------------------------------------------------------
# Check if RMG executable available at RMG_BIN or on the PATH
#-------------------------------------------------------------------
if(NOT QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)
  include("CMake/FindRMG.cmake")
  if(RMG_FOUND)
    message(STATUS "Found and using RMG (rmg-cpu) at ${RMG_CPU_EXE}")
  else(RMG_FOUND)
    message(STATUS "Did not find RMG (rmg-cpu). RMG tests will not be run.")
  endif(RMG_FOUND)
endif()

#-------------------------------------------------------------------
# If the timestamp of test_labels.py changes, unset any cached values
#-------------------------------------------------------------------

file(TIMESTAMP ${qmcpack_SOURCE_DIR}/tests/scripts/test_labels.py TIMESTAMP_TEST_LABELS)
if(NOT DEFINED CACHE_TIMESTAMP_TEST_LABELS)
  set(CACHE_TIMESTAMP_TEST_LABELS
      ${TIMESTAMP_TEST_LABELS}
      CACHE INTERNAL "Timestamp used to validate TEST_LABELS_ cache variables")
endif()
if(NOT ${TIMESTAMP_TEST_LABELS} STREQUAL ${CACHE_TIMESTAMP_TEST_LABELS})
  get_cmake_property(_variableNames VARIABLES)
  list(FILTER _variableNames INCLUDE REGEX "^TEST_LABELS_*")
  foreach(_variableName ${_variableNames})
    unset(${_variableName} CACHE)
  endforeach()
  set(CACHE_TIMESTAMP_TEST_LABELS
      ${TIMESTAMP_TEST_LABELS}
      CACHE INTERNAL "Timestamp used to validate TEST_LABELS_ cache variables" FORCE)
endif()

######################################################################
message(STATUS "Ready to parse QMCPACK source tree")
######################################################################

add_subdirectory(external_codes)
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
