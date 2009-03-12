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

set(INTEL_COMPILER 1)
add_definitions(-DADD_ -DINLINE_ALL=inline)
# common options for intel compilers
set(ENABLE_OPENMP 1)
set(INTEL_OPTS "-g -restrict -unroll -fno-alias -O3 -ip -openmp")

#enable Interprocedural (IP) Optimizations
#-ipo_obj force generation of real object files (requires -ipo)
#set(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -Ob=1 -ansi -ipo -ipo_obj -cxxlib-icc")
#set(CMAKE_CC_FLAGS "-restrict -unroll -fno-alias -O3 -Ob=1 -ansi -ipo -ipo_obj")
#set(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -g -ansi")
#set(CMAKE_CC_FLAGS "-restrict -unroll -fno-alias -g -ansi")
#set(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -cxxlib-icc")
#set(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -ansi -fno-fnalias -ivdep_parallel -Ob=2")
#set(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -ivdep_parallel -Ob=2")
#set(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -Ob=2 -qp")
#set(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -Ob=2 -cxxlib-icc")
#set(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3")
MESSAGE("-- CPU_IDENTITY ${CPU_IDENTITY}")

if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "i386")
  if(CPU_IDENTITY MATCHES "E4")
    set(INTEL_OPTS "${INTEL_OPTS} -xP ")
  endif(CPU_IDENTITY MATCHES "E4")
elseif(${CMAKE_SYSTEM_PROCESSOR} MATCHES "ia64")
  add_definitions(-DSGI_IA64)
  set(INTEL_OPTS "${INTEL_OPTS} -ftz -prefetch")
else(${CMAKE_SYSTEM_PROCESSOR} MATCHES "i386")
  #if(CPU_IDENTITY MATCHES "E5")
  #  set(INTEL_OPTS "${INTEL_OPTS} -xT ")
  #else(CPU_IDENTITY MATCHES "E5")
  #  set(INTEL_OPTS "${INTEL_OPTS} -xP ")
  #endif(CPU_IDENTITY MATCHES "E5")
  set(INTEL_OPTS "${INTEL_OPTS} -xT")
  set(HAVE_SSE 1)
  set(HAVE_SSE2 1)
  set(HAVE_SSE3 1)
  set(HAVE_SSSE3 1)
endif(${CMAKE_SYSTEM_PROCESSOR} MATCHES "i386")

#if($ENV{CXX} matches cmpi)
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ccl icpc")
#endif($ENV{CXX} matches cmpi)

#IF(${CMAKE_SYSTEM_PROCESSOR} matches "unknown")
#  IF(CPU_IDENTITY matches "E5")
#    set(CMAKE_CXX_FLAGS "${INTEL_OPTS} -xT ")
#    set(CMAKE_C_FLAGS "${INTEL_OPTS} -xT")
#  ELSE(CPU_IDENTITY matches "E5")
#  set(CMAKE_CXX_FLAGS "${INTEL_OPTS} -xP ")
#  set(CMAKE_C_FLAGS "${INTEL_OPTS} -xP")
#  endif(CPU_IDENTITY matches "E5")
#  set(HAVE_SSE 1)
#  set(HAVE_SSE2 1)
#  set(HAVE_SSE3 1)
#  set(HAVE_SSSE3 1)
#  #set(USE_PREFETCH 1)
#  #set(PREFETCH_AHEAD 12)
#endif(${CMAKE_SYSTEM_PROCESSOR} matches "unknown")

set(CMAKE_CXX_FLAGS "$ENV{CXX_FLAGS} ${INTEL_OPTS} -Wno-deprecated")
set(CMAKE_C_FLAGS "$ENV{CC_FLAGS} ${INTEL_OPTS} -std=c99 -Wno-deprecated")

#ifc -> ifort
#set(FORTRAN_LIBS " -lifcore -lifport")
set(F77 ifort)
set(F77OPTFLAGS  -fpp2 -O3)
set(F77FLAGS ${F77OPTFLAGS})

######################################################
#KCC needs to be used to build static libraries
######################################################
#set(CMAKE_AR xild) 
#set(CMAKE_CXX_CREATE_STATIC_LIBRARY "<CMAKE_AR> -lib -o <TARGET> <OBJECTS>")
