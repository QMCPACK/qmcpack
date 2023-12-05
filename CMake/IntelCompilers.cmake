include(CheckCCompilerFlag)
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
      set(OFFLOAD_TARGET
          "spir64"
          CACHE STRING "Offload target architecture")
      set(OPENMP_OFFLOAD_COMPILE_OPTIONS "-fopenmp-targets=${OFFLOAD_TARGET}")
    endif(ENABLE_OFFLOAD)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fiopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fiopenmp")
  else()
    if(ENABLE_OFFLOAD)
      message(FATAL_ERROR "OpenMP offload requires using Intel oneAPI compilers.")
    endif()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
  endif()
endif(QMC_OMP)

if(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
  # oneAPI compiler options

  # Set extra optimization specific flags
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffast-math")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")

  # Set clang specific flags (which we always want)
  add_compile_definitions(restrict=__restrict__)

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fstrict-aliasing")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fstrict-aliasing")

  # Force frame-pointer kept in DEBUG build.
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fno-omit-frame-pointer")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer")

  if(MIXED_PRECISION)
    # without this test_longrange went off in the mixed precsion build
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffp-model=precise")
  endif()
else()
  # classic compiler options

  # Suppress compile warnings
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-deprecated")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")

  # Set extra optimization specific flags
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -restrict -unroll -ip")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -restrict -unroll -ip")

  # Set prefetch flag
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-prefetch")

  if(MIXED_PRECISION)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -prec-sqrt")
  endif()

  #check if -ftz is accepted
  check_cxx_compiler_flag("${CMAKE_CXX_FLAGS} -ftz" INTEL_FTZ)
  if(INTEL_FTZ)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftz")
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
    if(CMAKE_CXX_FLAGS MATCHES ${X_OPTION}
       OR CMAKE_C_FLAGS MATCHES ${X_OPTION}
       OR CMAKE_CXX_FLAGS MATCHES ${AX_OPTION}
       OR CMAKE_C_FLAGS MATCHES ${AX_OPTION})
      # make sure that the user specifies -x for both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS.
      if(CMAKE_CXX_FLAGS MATCHES ${X_OPTION} AND CMAKE_C_FLAGS MATCHES ${X_OPTION})

      else() #(CMAKE_CXX_FLAGS MATCHES "-x" AND CMAKE_C_FLAGS MATCHES "-x")
        if(CMAKE_CXX_FLAGS MATCHES ${AX_OPTION} AND CMAKE_C_FLAGS MATCHES ${AX_OPTION})

        else()
          message(
            FATAL_ERROR
              "if -xcode or -axcode is specified by the user, it should be added in both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS!"
          )
        endif()
      endif() #(CMAKE_CXX_FLAGS MATCHES "-x" AND CMAKE_C_FLAGS MATCHES "-x")
    else() #(CMAKE_CXX_FLAGS MATCHES "-x" OR CMAKE_C_FLAGS MATCHES "-x")
      #check if -xHost is accepted
      check_c_compiler_flag("-xHost" INTEL_CC_FLAGS)
      if(INTEL_CC_FLAGS)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xHost")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost")
      endif(INTEL_CC_FLAGS)
    endif() #(CMAKE_CXX_FLAGS MATCHES "-x" OR CMAKE_C_FLAGS MATCHES "-x")
  else()
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
      # skipped in OneAPI 2022.0 when using SYCL which caused linking failure.
      if (NOT (CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 2022.0 AND ENABLE_SYCL))
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
      endif()
    endif() #(CMAKE_CXX_FLAGS MATCHES "-march=" OR CMAKE_C_FLAGS MATCHES "-march=")
  endif()

endif()
