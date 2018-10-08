# Check compiler version
IF ( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.3 )
  MESSAGE(STATUS "Compiler Version ${CMAKE_CXX_COMPILER_VERSION}")
  MESSAGE(FATAL_ERROR "Requires clang 3.3 or higher ")
ENDIF()

# Set the std
SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -std=c99")

# Enable OpenMP
IF(QMC_OMP)
  SET(ENABLE_OPENMP 1)
  SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -fopenmp")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
ENDIF(QMC_OMP)

# Set clang specfic flags (which we always want)
ADD_DEFINITIONS( -Drestrict=__restrict__ )

SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -fomit-frame-pointer -fstrict-aliasing")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fomit-frame-pointer -fstrict-aliasing -D__forceinline=inline")

# Suppress compile warnings
SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -Wno-deprecated -Wno-unused-value")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wno-unused-value")
IF ( CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "3.8.0" )
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-undefined-var-template")
ENDIF()

# Set extra optimization specific flags
SET( CMAKE_C_FLAGS_RELEASE     "${CMAKE_C_FLAGS_RELEASE} -ffast-math" )
SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math" )
SET( CMAKE_C_FLAGS_RELWITHDEBINFO     "${CMAKE_C_FLAGS_RELWITHDEBINFO} -ffast-math" )
SET( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -ffast-math" )

#--------------------------------------
# Neither on Cray's machine nor PowerPC
#--------------------------------------
IF((NOT $ENV{CRAYPE_VERSION} MATCHES ".") AND (NOT CMAKE_SYSTEM_PROCESSOR MATCHES "ppc64"))

#check if the user has already specified -march=XXXX option for cross-compiling.
if(CMAKE_CXX_FLAGS MATCHES "-march=" OR CMAKE_C_FLAGS MATCHES "-march=")
  # make sure that the user specifies -march= for both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS.
  if(CMAKE_CXX_FLAGS MATCHES "-march=" AND CMAKE_C_FLAGS MATCHES "-march=")
  else() #(CMAKE_CXX_FLAGS MATCHES "-march=" AND CMAKE_C_FLAGS MATCHES "-march=")
    MESSAGE(FATAL_ERROR "if -march=ARCH is specified by the user, it should be added in both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS!")
  endif() #(CMAKE_CXX_FLAGS MATCHES "-march=" AND CMAKE_C_FLAGS MATCHES "-march=")
else() #(CMAKE_CXX_FLAGS MATCHES "-march=" OR CMAKE_C_FLAGS MATCHES "-march=")
  # use -march=native
  SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -march=native")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif() #(CMAKE_CXX_FLAGS MATCHES "-march=" OR CMAKE_C_FLAGS MATCHES "-march=")

ENDIF((NOT $ENV{CRAYPE_VERSION} MATCHES ".") AND (NOT CMAKE_SYSTEM_PROCESSOR MATCHES "ppc64"))

# Add static flags if necessary
IF(QMC_BUILD_STATIC)
    SET(CMAKE_CXX_LINK_FLAGS " -static")
ENDIF(QMC_BUILD_STATIC)

# Coverage
IF (ENABLE_GCOV)
  SET(GCOV_COVERAGE TRUE)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
ENDIF(ENABLE_GCOV)

SET(XRAY_PROFILE FALSE CACHE BOOL "Use llvm xray profiling")
SET(XRAY_INSTRUCTION_THRESHOLD 200 CACHE INT "Instruction threshold for xray instrumentation")
SET(XRAY_GPU_MOST FALSE CACHE BOOL "Manually instrument almost all of the CUDA implementations calls on the CPU")

IF(XRAY_PROFILE)
  set(XRAY_FLAGS "-fxray-instrument -fxray-instruction-threshold=${XRAY_INSTRUCTION_THRESHOLD}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${XRAY_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${XRAY_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${XRAY_FLAGS}")
  if(XRAY_GPU_MOST)
    set_property(DIRECTORY ${CMAKE_SOURCE_DIR} APPEND PROPERTY COMPILE_DEFINITIONS GPU_XRAY_TRACE_ON)
  endif(XRAY_GPU_MOST)
ENDIF(XRAY_PROFILE)
