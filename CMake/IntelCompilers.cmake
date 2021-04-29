# Check compiler version
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 19.0.0.20190206)
message(FATAL_ERROR "Requires Intel 19 update 3 (19.0.0.20190206) or higher!")
endif()

# Set the std
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -std=c99")

# Enable OpenMP
IF(QMC_OMP)
  SET(ENABLE_OPENMP 1)
  SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -qopenmp")
  IF(ENABLE_OFFLOAD)
    SET(OFFLOAD_TARGET "host" CACHE STRING "Offload target architecture")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp -qopenmp-offload=${OFFLOAD_TARGET}")
  ELSE(ENABLE_OFFLOAD)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
  ENDIF(ENABLE_OFFLOAD)
ENDIF(QMC_OMP)

# Suppress compile warnings
SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -Wno-deprecated")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")

# Set extra optimization specific flags
SET( CMAKE_C_FLAGS_DEBUG     "${CMAKE_C_FLAGS_DEBUG} -restrict -unroll -ip" )
SET( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -restrict -unroll -ip" )
SET( CMAKE_C_FLAGS_RELEASE     "${CMAKE_C_FLAGS_RELEASE} -restrict -unroll -ip" )
SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -restrict -unroll -ip" )
SET( CMAKE_C_FLAGS_RELWITHDEBINFO     "${CMAKE_C_FLAGS_RELWITHDEBINFO} -restrict -unroll -ip" )
SET( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -restrict -unroll -ip" )

# Set prefetch flag
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-prefetch" )

IF(MIXED_PRECISION)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -prec-sqrt" )
ENDIF()
#check if -ftz is accepted
CHECK_CXX_COMPILER_FLAG( "${CMAKE_CXX_FLAGS} -ftz" INTEL_FTZ )
IF( INTEL_FTZ)
SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -ftz" )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftz" )
ENDIF( INTEL_FTZ)

#------------------------
# Not on Cray's machine
#------------------------
IF(NOT CMAKE_SYSTEM_NAME STREQUAL "CrayLinuxEnvironment")

SET(X_OPTION "^-x| -x")
SET(AX_OPTION "^-ax| -ax")
#check if the user has already specified -x option for cross-compiling.
if(CMAKE_CXX_FLAGS MATCHES ${X_OPTION} OR CMAKE_C_FLAGS MATCHES ${X_OPTION} OR
    CMAKE_CXX_FLAGS MATCHES ${AX_OPTION} OR CMAKE_C_FLAGS MATCHES ${AX_OPTION})
  # make sure that the user specifies -x for both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS.
  if(CMAKE_CXX_FLAGS MATCHES ${X_OPTION} AND CMAKE_C_FLAGS MATCHES ${X_OPTION})
  else() #(CMAKE_CXX_FLAGS MATCHES "-x" AND CMAKE_C_FLAGS MATCHES "-x")
    if(CMAKE_CXX_FLAGS MATCHES ${AX_OPTION} AND CMAKE_C_FLAGS MATCHES ${AX_OPTION})
    else()
      MESSAGE(FATAL_ERROR "if -xcode or -axcode is specified by the user, it should be added in both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS!")
    endif()
  endif() #(CMAKE_CXX_FLAGS MATCHES "-x" AND CMAKE_C_FLAGS MATCHES "-x")
else() #(CMAKE_CXX_FLAGS MATCHES "-x" OR CMAKE_C_FLAGS MATCHES "-x")
  #check if -xHost is accepted
  CHECK_C_COMPILER_FLAG( "-xHost" INTEL_CC_FLAGS )
  IF(INTEL_CC_FLAGS)
    SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -xHost")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost")
  ENDIF(INTEL_CC_FLAGS)
endif() #(CMAKE_CXX_FLAGS MATCHES "-x" OR CMAKE_C_FLAGS MATCHES "-x")

ENDIF()
