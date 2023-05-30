# Check compiler version
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.0)
  message(STATUS "Compiler Version ${CMAKE_CXX_COMPILER_VERSION}")
  message(FATAL_ERROR "Requires clang 7.0 or higher ")
endif()

if(CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 11.0.0
   AND QMC_CXX_STANDARD EQUAL 17
   AND BUILD_AFQMC)
  message(FATAL_ERROR "Avoid Clang 11.0.0 which cannot compile AFQMC properly with C++17!")
endif()

# Enable OpenMP
if(QMC_OMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

  if(ENABLE_OFFLOAD)
    if(DEFINED OFFLOAD_TARGET)
      set(OPENMP_OFFLOAD_COMPILE_OPTIONS "-fopenmp-targets=${OFFLOAD_TARGET}")
      if(DEFINED OFFLOAD_ARCH)
        set(OPENMP_OFFLOAD_COMPILE_OPTIONS
            "${OPENMP_OFFLOAD_COMPILE_OPTIONS} -Xopenmp-target=${OFFLOAD_TARGET} -march=${OFFLOAD_ARCH}")
      endif()
    elseif(QMC_GPU_ARCHS)
      string(REGEX REPLACE ";" "," QMC_GPU_ARCHS_COMMA_SEPARATED "${QMC_GPU_ARCHS}")
      set(OPENMP_OFFLOAD_COMPILE_OPTIONS "--offload-arch=${QMC_GPU_ARCHS_COMMA_SEPARATED}")
    else()
      message(FATAL_ERROR "Require QMC_GPU_ARCHS or OFFLOAD_TARGET set for OpenMP offload using Clang.")
    endif()

    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-Wno-linker-warnings" LINKER_WARNING_SUPPORTED)
    if(LINKER_WARNING_SUPPORTED)
      set(OPENMP_OFFLOAD_COMPILE_OPTIONS "${OPENMP_OFFLOAD_COMPILE_OPTIONS} -Wno-linker-warnings")
    endif()

    # additional customization for NVIDIA toolchain
    if(OPENMP_OFFLOAD_COMPILE_OPTIONS MATCHES "sm_")
      set(OPENMP_OFFLOAD_COMPILE_OPTIONS "${OPENMP_OFFLOAD_COMPILE_OPTIONS} -Wno-unknown-cuda-version")
      # unfortunately this removes standalone-debug altogether for offload builds
      # but until we discover how to use the ${OPENMP_OFFLOAD_COMPILE_OPTIONS} more selectively
      # this is the only way to avoid a warning per compilation unit that contains an omp symbol.
      message(STATUS "QMCPACK adds -fstandalone-debug for Debug builds")
      set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fstandalone-debug")
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fstandalone-debug")
    endif()

    check_cxx_compiler_flag(-fopenmp-assume-no-thread-state OPENMP_ASSUME_NO_THREAD_STATE_WORKS)
    option(OMPTARGET_ASSUME_NO_THREAD_STATE "Use -fopenmp-assume-no-thread-state compile option" ${OPENMP_ASSUME_NO_THREAD_STATE_WORKS})
    if(OMPTARGET_ASSUME_NO_THREAD_STATE)
      set(OPENMP_OFFLOAD_COMPILE_OPTIONS "${OPENMP_OFFLOAD_COMPILE_OPTIONS} -fopenmp-assume-no-thread-state")
    endif()

    check_cxx_compiler_flag(-fopenmp-assume-no-nested-parallelism OPENMP_ASSUME_NO_NESTED_PAR_WORKS)
    option(OMPTARGET_ASSUME_NO_NESTED_PAR "Use -fopenmp-assume-no-nested-parallelism compile option" ${OPENMP_ASSUME_NO_NESTED_PAR_WORKS})
    if(OMPTARGET_ASSUME_NO_NESTED_PAR)
      set(OPENMP_OFFLOAD_COMPILE_OPTIONS "${OPENMP_OFFLOAD_COMPILE_OPTIONS} -fopenmp-assume-no-nested-parallelism")
    endif()

    # AOMP/ROCM special option -fdisable-host-devmem to disable unecessary data transfers
    # The down side of using this option is that it disables printf from offload regions.
    # see https://github.com/ROCm-Developer-Tools/aomp/issues/526
    check_cxx_compiler_flag(-fdisable-host-devmem DISABLE_HOST_DEVMEM_WORKS)
    option(AMDGPU_DISABLE_HOST_DEVMEM "Use -fdisable-host-devmem link option" ${DISABLE_HOST_DEVMEM_WORKS})
    if(AMDGPU_DISABLE_HOST_DEVMEM)
      string(APPEND CMAKE_EXE_LINKER_FLAGS " -fdisable-host-devmem")
    endif()
  endif()
endif(QMC_OMP)

# Set clang specific flags (which we always want)
add_definitions(-Drestrict=__restrict__)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fstrict-aliasing")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fstrict-aliasing")

# treat VLA as error
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror=vla")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wvla")

# set compiler warnings
string(APPEND CMAKE_CXX_FLAGS
       " -Wall -Wno-unused-variable -Wno-overloaded-virtual -Wno-unused-private-field -Wno-unused-local-typedef")

if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 11.0)
  string(APPEND CMAKE_CXX_FLAGS " -Wsuggest-override")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 10.0)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wmisleading-indentation")
endif()

# Set extra optimization specific flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffast-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")

# Set extra debug flags
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer")

#--------------------------------------
# Special architectural flags
#--------------------------------------
# case arch
#     x86_64: -march
#     powerpc: -mpcu
#     arm: -mpcu
#     default or cray: none
#--------------------------------------
if(CMAKE_SYSTEM_NAME STREQUAL "CrayLinuxEnvironment")
  # It's a cray machine. Don't do anything
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
  # the case for x86_64
  # check if the user has already specified -march=XXXX option for cross-compiling.
  if(CMAKE_CXX_FLAGS MATCHES "-march=" OR CMAKE_C_FLAGS MATCHES "-march=")
    # make sure that the user specifies -march= for both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS.
    if(CMAKE_CXX_FLAGS MATCHES "-march=" AND CMAKE_C_FLAGS MATCHES "-march=")

    else() #(CMAKE_CXX_FLAGS MATCHES "-march=" AND CMAKE_C_FLAGS MATCHES "-march=")
      message(
        FATAL_ERROR
          "if -march=ARCH is specified by the user, it should be added in both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS!")
    endif() #(CMAKE_CXX_FLAGS MATCHES "-march=" AND CMAKE_C_FLAGS MATCHES "-march=")
  else() #(CMAKE_CXX_FLAGS MATCHES "-march=" OR CMAKE_C_FLAGS MATCHES "-march=")
    # use -march=native
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  endif() #(CMAKE_CXX_FLAGS MATCHES "-march=" OR CMAKE_C_FLAGS MATCHES "-march=")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "ppc64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64")
  # the case for PowerPC and ARM
  # check if the user has already specified -mcpu=XXXX option for cross-compiling.
  if(CMAKE_CXX_FLAGS MATCHES "-mcpu=" OR CMAKE_C_FLAGS MATCHES "-mcpu=")
    # make sure that the user specifies -mcpu= for both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS.
    if(CMAKE_CXX_FLAGS MATCHES "-mcpu=" AND CMAKE_C_FLAGS MATCHES "-mcpu=")

    else() #(CMAKE_CXX_FLAGS MATCHES "-mcpu=" AND CMAKE_C_FLAGS MATCHES "-mcpu=")
      message(
        FATAL_ERROR
          "if -mcpu=ARCH is specified by the user, it should be added in both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS!")
    endif() #(CMAKE_CXX_FLAGS MATCHES "-mcpu=" AND CMAKE_C_FLAGS MATCHES "-mcpu=")
  else() #(CMAKE_CXX_FLAGS MATCHES "-mcpu=" OR CMAKE_C_FLAGS MATCHES "-mcpu=")
    # use -mcpu=native
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mcpu=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=native")
  endif() #(CMAKE_CXX_FLAGS MATCHES "-mcpu=" OR CMAKE_C_FLAGS MATCHES "-mcpu=")
endif()

# Add static flags if necessary
if(QMC_BUILD_STATIC)
  set(CMAKE_CXX_LINK_FLAGS " -static")
endif(QMC_BUILD_STATIC)

# Coverage
if(ENABLE_GCOV)
  set(GCOV_SUPPORTED TRUE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
endif(ENABLE_GCOV)

set(XRAY_PROFILE
    FALSE
    CACHE BOOL "Use llvm xray profiling")
set(XRAY_INSTRUCTION_THRESHOLD
    200
    CACHE STRING "Instruction threshold for xray instrumentation")

if(XRAY_PROFILE)
  set(XRAY_FLAGS "-fxray-instrument -fxray-instruction-threshold=${XRAY_INSTRUCTION_THRESHOLD}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${XRAY_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${XRAY_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${XRAY_FLAGS}")
endif(XRAY_PROFILE)
