#compiler flags for intel 8.x, 9.x, 10.x
#check the version and warn users if older compilers are being used
exec_program(icpc  
  ARGS -v
  OUTPUT_VARIABLE ICPC_OUTPUT
  RETURN_VALUE ICPC_INFO
  )
STRING(REGEX MATCH "[0-9]+" ICPC_VERSION "${ICPC_OUTPUT}")
IF(ICPC_VERSION LESS 10)
  MESSAGE("-- Found Intel Compilers Version ${ICPC_OUTPUT}. Please upgrade to 10.0.3 or higher")
ENDIF(ICPC_VERSION LESS 10)

SET(INTEL_COMPILER 1)
#ADD_DEFINITIONS(-DADD_ -DINLINE_ALL=inline)
ADD_DEFINITIONS(-DADD_ -DINLINE_ALL=inline -DMPICH_SKIP_MPICXX)
#enable Interprocedural (IP) Optimizations
#-ipo_obj force generation of real object files (requires -ipo)
#SET(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -Ob=1 -ansi -ipo -ipo_obj -cxxlib-icc")
#SET(CMAKE_CC_FLAGS "-restrict -unroll -fno-alias -O3 -Ob=1 -ansi -ipo -ipo_obj")
#SET(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -g -ansi")
#SET(CMAKE_CC_FLAGS "-restrict -unroll -fno-alias -g -ansi")
#SET(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -cxxlib-icc")
#SET(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -ansi -fno-fnalias -ivdep_parallel -Ob=2")
#SET(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -ivdep_parallel -Ob=2")
#SET(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -Ob=2 -qp")
#SET(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -Ob=2 -cxxlib-icc")
#SET(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3")

# common options for intel compilers
SET(INTEL_OPTS "-g -restrict -unroll -fno-alias -O3 -ip")

IF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "i386")
  SET(CMAKE_CXX_FLAGS "${INTEL_OPTS}")
  SET(CMAKE_C_FLAGS "${INTEL_OPTS}")

  IF(CPU_IDENTITY MATCHES "E4")
    SET(CMAKE_CXX_FLAGS "${INTEL_OPTS} -xP ")
    SET(CMAKE_C_FLAGS "${INTEL_OPTS} -xP")
  ENDIF(CPU_IDENTITY MATCHES "E4")

ENDIF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "i386")

IF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "ia64")
  ADD_DEFINITIONS(-DSGI_IA64)
  SET(CMAKE_CXX_FLAGS "${INTEL_OPTS} -ftz -prefetch")
  SET(CMAKE_C_FLAGS "${INTEL_OPTS} -ftz -prefetch")
ENDIF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "ia64")

IF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
  IF(CPU_IDENTITY MATCHES "E5")
    SET(CMAKE_CXX_FLAGS "${INTEL_OPTS} -xT ")
    SET(CMAKE_C_FLAGS "${INTEL_OPTS} -xT")
  ENDIF(CPU_IDENTITY MATCHES "E5")
  SET(HAVE_SSE 1)
  SET(HAVE_SSE2 1)
  SET(HAVE_SSE3 1)
  SET(HAVE_SSSE3 1)
ENDIF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")

IF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "unknown")
  IF(CPU_IDENTITY MATCHES "E5")
    SET(CMAKE_CXX_FLAGS "${INTEL_OPTS} -xT ")
    SET(CMAKE_C_FLAGS "${INTEL_OPTS} -xT")
  ELSE(CPU_IDENTITY MATCHES "E5")
  SET(CMAKE_CXX_FLAGS "${INTEL_OPTS} -xP ")
  SET(CMAKE_C_FLAGS "${INTEL_OPTS} -xP")
  ENDIF(CPU_IDENTITY MATCHES "E5")
  SET(HAVE_SSE 1)
  SET(HAVE_SSE2 1)
  SET(HAVE_SSE3 1)
  SET(HAVE_SSSE3 1)
  #SET(USE_PREFETCH 1)
  #SET(PREFETCH_AHEAD 12)
ENDIF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "unknown")

SET(ENABLE_OPENMP 1)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -openmp")
SET(CMAKE_CXX_FLAGS "$ENV{CXX_FLAGS} ${CMAKE_CXX_FLAGS} -Wno-deprecated")
SET(CMAKE_C_FLAGS "$ENV{CC_FLAGS} ${CMAKE_C_FLAGS} -Wno-deprecated")

#IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
#ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#  IF(QMC_BITS MATCHES 64)
#    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftz -cxxlib-icc")
#    IF(HAVE_MPI)
#      LINK_LIBRARIES(-lmpi)
#    ENDIF(HAVE_MPI)
#  ELSE(QMC_BITS MATCHES 64)
#    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xW -cxxlib-icc")
#  ENDIF(QMC_BITS MATCHES 64)
#ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

IF($ENV{CXX} MATCHES cmpi)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ccl icpc")
ENDIF($ENV{CXX} MATCHES cmpi)

#ifc -> ifort
#SET(FORTRAN_LIBS " -lifcore -lifport")
SET(F77 ifort)
SET(F77OPTFLAGS  -fpp2 -O3)
SET(F77FLAGS ${F77OPTFLAGS})

######################################################
#KCC needs to be used to build static libraries
######################################################
#SET(CMAKE_AR xild) 
#SET(CMAKE_CXX_CREATE_STATIC_LIBRARY "<CMAKE_AR> -lib -o <TARGET> <OBJECTS>")
