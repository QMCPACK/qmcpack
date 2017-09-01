# Check compiler version
SET(INTEL_COMPILER 1)
IF ( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0 )
MESSAGE(FATAL_ERROR "Requires Intel 15.0 or higher ")
ENDIF()

# Set the std
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -std=c99")

# Enable OpenMP
SET(ENABLE_OPENMP 1)
IF ( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 16 )
SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -openmp")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")
ELSE()
SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -qopenmp")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
ENDIF()

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

# Use deprecated options prior to 11.1
SET(ICC_DEPRECATED_OPTS FALSE)
IF ( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11.1 )
SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -prefetch ")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -prefetch" )
ELSEIF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 16 )  
SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -opt-prefetch" )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -opt-prefetch" )
ELSE()
SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -qopt-prefetch" )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-prefetch" )
ENDIF()

#check if -ftz is accepted
CHECK_C_COMPILER_FLAG( "${CMAKE_CXX_FLAGS} -ftz" INTEL_FTZ )
IF( INTEL_FTZ)
SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -ftz" )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftz" )
ENDIF( INTEL_FTZ)

#------------------------
# Not on Cray's machine
#------------------------
IF(NOT $ENV{CRAYPE_VERSION} MATCHES ".")

#check if the user has already specified -x option for cross-compiling.
if(CMAKE_CXX_FLAGS MATCHES "-x" OR CMAKE_C_FLAGS MATCHES "-x" OR
    CMAKE_CXX_FLAGS MATCHES "-ax" OR CMAKE_C_FLAGS MATCHES "-ax")
  # make sure that the user specifies -x for both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS.
  if(CMAKE_CXX_FLAGS MATCHES "-x" AND CMAKE_C_FLAGS MATCHES "-x")
  else() #(CMAKE_CXX_FLAGS MATCHES "-x" AND CMAKE_C_FLAGS MATCHES "-x")
    if(CMAKE_CXX_FLAGS MATCHES "-ax" AND CMAKE_C_FLAGS MATCHES "-ax")
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

ENDIF(NOT $ENV{CRAYPE_VERSION} MATCHES ".")
