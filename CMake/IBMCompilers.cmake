#IBM Visual Age C/C++ compilers
#
#check guidec++ and overwrites the BITS
IF(CMAKE_SYSTEM_NAME MATCHES "AIX")
MESSAGE(STATUS "AIX system using xlC/xlc/xlf")

######################################################################
#AIX_BIT_FLAGS is the option to switch between 32/64
#If BITS=64, compiler and ar options have to have -q64 or -X64
######################################################################
SET(AIX_BIT_FLAGS "")
IF(QMC_BITS MATCHES 64)
  SET(AIX_BIT_FLAGS " -q64")
  SET(AR_OPTIONS "-X64")
ENDIF(QMC_BITS MATCHES 64)

######################################################################
#AIX_CXX_COMMON_FLAGS for the flags required for any level of optimizations
#AIX_CXX_OPT_FLAGS    for aggressive optimization
#AIX_CXX_COMMON_FLAGS for basic optimization level
######################################################################
#SET(AIX_CXX_COMMON_FLAGS "-qarch=pwr4 -qtune=pwr4 -qkeyword=restrict -qstrict -qansialias -qalias=typeptr -qnoeh -qrtti=dyna -qtemplateregistry=${LIBRARY_OUTPUT_PATH}/tempegistry")
#SET(AIX_CXX_OPT_FLAGS "-O3 -qipa -qinline -qspill=2000")
SET(AIX_CXX_COMMON_FLAGS "-qarch=pwr4 -qtune=pwr4 -qkeyword=restrict -qstrict -qansialias -qrtti=dyna ")
#SET(AIX_CXX_OPT_FLAGS "-O3 -DINLINE_ALL=  -Q -qipa=inline -qinline ")
#SET(AIX_CXX_OPT_FLAGS "-O3 -DINLINE_ALL=  -Q ")
SET(AIX_CXX_OPT_FLAGS "-DINLINE_ALL=  -Q ") # Optimization settings will be picked by build type
SET(AIX_CXX_FLAGS "-O")

######################################################################
#set the CXX flags: common + opt + bit
######################################################################
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${AIX_CXX_COMMON_FLAGS} ${AIX_CXX_OPT_FLAGS} ${AIX_BIT_FLAGS}")

######################################################################
#add openmp-related flags: common + opt + bit + thread
######################################################################
IF(QMC_OMP)
  # check if the compiler is correct xlC_r xlc_r mpCC_r mpcc_r
  IF($ENV{CXX} MATCHES "_r")
    SET(ENABLE_OPENMP 1)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qsmp=omp")
    SET(ENABLE_OMP 1 CACHE BOOL "OpenMP is enabled")
  ELSE($ENV{CXX} MATCHES "_r")
    MESSAGE(STATUS "OpenMP is not enabled. Change CXX to xlC_r/mpCC_r")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qnothreaded")
  ENDIF($ENV{CXX} MATCHES "_r")
ELSE(QMC_OMP)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qnothreaded")
ENDIF(QMC_OMP)

IF(QMC_BITS MATCHES 64)
  #SET(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_FLAGS})
  SET(CMAKE_CXX_CREATE_STATIC_LIBRARY
        "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS> " 
        "<CMAKE_RANLIB> <TARGET> ")
ENDIF(QMC_BITS MATCHES 64)
ENDIF(CMAKE_SYSTEM_NAME MATCHES "AIX")

IF(CMAKE_SYSTEM_PROCESSOR MATCHES "ppc64le")
MESSAGE(STATUS "Power8+ system using xlC/xlc/xlf")

ADD_DEFINITIONS( -Drestrict=__restrict__ )
ADD_DEFINITIONS( -DADD_ )
ADD_DEFINITIONS( -DINLINE_ALL=inline )

# Clean up flags

IF(CMAKE_C_FLAGS MATCHES "-qhalt=e")
  SET(CMAKE_C_FLAGS "")
ENDIF()
IF(CMAKE_CXX_FLAGS MATCHES "-qhalt=e")
  SET(CMAKE_CXX_FLAGS "")
ENDIF()

# Suppress compile warnings
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -Wno-deprecated -Wno-unused-value")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wno-unused-value")

# Set standard and other optimization settings
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -std=gnu11 -O3")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -qnoxlcompatmacros -O3")

IF(QMC_OMP)
  SET(ENABLE_OPENMP 1)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qsmp=omp")
  SET(ENABLE_OMP 1 CACHE BOOL "OpenMP is enabled")
ELSE(QMC_OMP)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qnothreaded")
ENDIF(QMC_OMP)

# Add static flags if necessary
IF(QMC_BUILD_STATIC)
    SET(CMAKE_CXX_LINK_FLAGS " -static")
ENDIF(QMC_BUILD_STATIC)

ENDIF(CMAKE_SYSTEM_PROCESSOR MATCHES "ppc64le")
