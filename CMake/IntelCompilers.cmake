include(CheckCXXCompilerFlag)

# Check compiler version
message(DEBUG "CMAKE_CXX_COMPILER_VERSION ${CMAKE_CXX_COMPILER_VERSION}")
if(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 2021.3)
    message(FATAL_ERROR "Requires Intel oneAPI 2021.3 or higher!")
  endif()
else()
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 2021.1)
    message(FATAL_ERROR "Requires Intel classic compiler 2021.1 or higher!")
  endif()
endif()

# Enable OpenMP
if(QMC_OMP)
  if(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
    if(ENABLE_OFFLOAD)
      if(OFFLOAD_ARCH OR QMC_GPU_ARCHS)
        # for ahead-of-time compilation and linking
        # if OFFLOAD_ARCH not defined, overwrite it with QMC_GPU_ARCHS
        if(NOT OFFLOAD_ARCH AND QMC_GPU_ARCHS)
          set(OFFLOAD_ARCH ${QMC_GPU_ARCHS})
        endif()
        if(OFFLOAD_ARCH MATCHES "^intel_gpu_")
          set(OPENMP_OFFLOAD_COMPILE_OPTIONS "-fopenmp-targets=spir64_gen")
	  string(REGEX REPLACE "^intel_gpu_" "" INTEL_GPU_ARCH "${OFFLOAD_ARCH}")
          set(OpenMP_OFFLOAD_LINKER_FLAGS "-Xs \"-device ${INTEL_GPU_ARCH}\"")
        else()
          message(FATAL_ERROR "Invalid Intel GPU architecture \"${OFFLOAD_ARCH}\"! Did you miss \"intel_gpu_\" prefix?")
        endif()
      else()
        set(OFFLOAD_TARGET
            "spir64"
            CACHE STRING "Offload target architecture")
        set(OPENMP_OFFLOAD_COMPILE_OPTIONS "-fopenmp-targets=${OFFLOAD_TARGET}")
      endif()
      # Select the intra-team reduction implementation using shared local memory.
      set(OPENMP_OFFLOAD_COMPILE_OPTIONS "${OPENMP_OFFLOAD_COMPILE_OPTIONS} -mllvm -vpo-paropt-atomic-free-reduction-slm=true" )
    endif(ENABLE_OFFLOAD)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fiopenmp")
  else()
    if(ENABLE_OFFLOAD)
      message(FATAL_ERROR "OpenMP offload requires using Intel oneAPI compilers.")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
  endif()
endif(QMC_OMP)

if(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
  # oneAPI compiler options

  # Set extra optimization specific flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")

  # Set clang specific flags (which we always want)
  add_compile_definitions(restrict=__restrict__)

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fstrict-aliasing")

  # Force frame-pointer kept in DEBUG build.
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer")

  # The internal device-linking step of oneapi clutters link.d with temporary files on /tmp and results in re-linking every time.
  if(QMC_GPU AND NOT DEFINED CMAKE_LINK_DEPENDS_USE_LINKER)
    set(CMAKE_LINK_DEPENDS_USE_LINKER FALSE)
  endif()
else()
  # classic compiler options

  # Suppress compile warnings
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")

  # Set extra optimization specific flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -restrict -unroll -ip")

  # Set prefetch flag
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-prefetch")

  if(MIXED_PRECISION)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -prec-sqrt")
  endif()

  #check if -ftz is accepted
  check_cxx_compiler_flag("${CMAKE_CXX_FLAGS} -ftz" INTEL_FTZ)
  if(INTEL_FTZ)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftz")
  endif(INTEL_FTZ)
endif()

#------------------------
# Not on Cray's machine
#------------------------
if(NOT CMAKE_SYSTEM_NAME STREQUAL "CrayLinuxEnvironment")

  # use -x for classic compiler only. this option is not robust with oneAPI compiler as 2021.3 release
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
    set(X_OPTION "^-x| -x")
    set(AX_OPTION "^-ax| -ax")
    #check if the user has already specified -x option for cross-compiling.
    if(NOT (CMAKE_CXX_FLAGS MATCHES ${X_OPTION} OR CMAKE_CXX_FLAGS MATCHES ${AX_OPTION}))
      #check if -xHost is accepted
      check_cxx_compiler_flag("-xHost" INTEL_CXX_FLAGS)
      if(INTEL_CXX_FLAGS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost")
      endif(INTEL_CXX_FLAGS)
    endif()
  else()
    # check if the user has already specified -march=XXXX option for cross-compiling.
    if(NOT CMAKE_CXX_FLAGS MATCHES "-march=")
      # use -march=native
      # skipped in OneAPI 2022.0 when using SYCL which caused linking failure.
      if (NOT (CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 2022.0 AND ENABLE_SYCL))
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
      endif()
    endif()
  endif()

endif()
