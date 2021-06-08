# Check compiler version
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 19.0.0.20190206)
  message(FATAL_ERROR "Requires Intel 19 update 3 (19.0.0.20190206) or higher!")
endif()

# Enable OpenMP
if(QMC_OMP)
  set(ENABLE_OPENMP 1)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qopenmp")
  if(ENABLE_OFFLOAD)
    set(OFFLOAD_TARGET
        "host"
        CACHE STRING "Offload target architecture")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp -qopenmp-offload=${OFFLOAD_TARGET}")
  else(ENABLE_OFFLOAD)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
  endif(ENABLE_OFFLOAD)
endif(QMC_OMP)

# Suppress compile warnings
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-deprecated")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")

# Set extra optimization specific flags
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -restrict -unroll -ip")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -restrict -unroll -ip")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -restrict -unroll -ip")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -restrict -unroll -ip")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -restrict -unroll -ip")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -restrict -unroll -ip")

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

#------------------------
# Not on Cray's machine
#------------------------
if(NOT CMAKE_SYSTEM_NAME STREQUAL "CrayLinuxEnvironment")

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

endif()
