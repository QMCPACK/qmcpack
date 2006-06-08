#
#adding compiler flags
#
#Since F77FLAGS is used with ADD_CUSTOM_COMMAND variable F77FLAGS needs to be comma separated.
# incorrect:   SET(F77OPTFLAGS  "-fpp2 -O3 -unrol -fno-alias -ftz")
# correct:     SET(F77OPTFLAGS  -fpp2 -O3 -unrol -fno-alias -ftz)

MESSAGE(STATUS "Configuring compiler flags")

IF($ENV{CXX} MATCHES icc)
  SET(INTEL_COMPILER 1)
  ADD_DEFINITIONS(-DADD_)
  SET(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -xW")
#  SET(FORTRAN_LIBS "-lPEPCF90 -lCEPCF90 -lF90 -lintrins")
  SET(FORTRAN_LIBS "-lifcore")
  SET(F77 ifc)
  SET(F77OPTFLAGS  -fpp2 -O3)
  SET(F77FLAGS ${F77OPTFLAGS})
ENDIF($ENV{CXX} MATCHES icc)

IF($ENV{CXX} MATCHES ecc)
  SET(INTEL_COMPILER 1)
  ADD_DEFINITIONS(-DADD_)
  SET(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -ftz")
  SET(FORTRAN_LIBS "-lPEPCF90 -lCEPCF90 -lF90 -lintrins")
  SET(F77 efc)
  SET(F77OPTFLAGS  -fpp2 -O3 -unrol -fno-alias -ftz)
  SET(F77FLAGS ${F77OPTFLAGS})
ENDIF($ENV{CXX} MATCHES ecc)

IF(INTEL_COMPILER)
  MESSAGE(STATUS "Using intel compiler")
  IF(OHMMS_OMP)
    ADD_DEFINITIONS(-DUSE_OPENMP)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")
    SET(ENABLE_OMP 1 CACHE BOOL "OpenMP is enabled")
  ENDIF(OHMMS_OMP)
ENDIF(INTEL_COMPILER)

IF(${CMAKE_SYSTEM_NAME} MATCHES "AIX")

  # set the customized compiler flags
  SET(CMAKE_CXX_FLAGS "-O3  -qstrict -qarch=auto -qtune=auto -qansialias -qkeepinlines -qinline=1000 -qspill=2000  -qalias=typeptr -qnoeh -qnothreaded -qrtti=dyna")

  # set the fortran libraries for linker
  SET(FORTRAN_LIBS " -lxlf90 -lxlf")

  # check if the compiler is correct xlC_r xlc_r mpCC_r mpcc_r
  IF(OHMMS_OMP)
    IF($ENV{CXX} MATCHES "_r")
      ADD_DEFINITIONS(-DUSE_OPENMP)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qsmp=omp")
      SET(F77 xlf_r)
      SET(ENABLE_OMP 1 CACHE BOOL "OpenMP is enabled")
    ELSE($ENV{CXX} MATCHES "_r")
      MESSAGE(STATUS "OpenMP is not enabled. Change CXX to xlC_r/mpCC_r")
      SET(F77 xlf)
    ENDIF($ENV{CXX} MATCHES "_r")
  ELSE(OHMMS_OMP)
    SET(F77 xlf)
  ENDIF(OHMMS_OMP)

  SET(F77OPTFLAGS  -O3 -qstrict -qarch=auto -qtune=auto)
  SET(F77FLAGS ${F77OPTFLAGS})
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "AIX")

IF(NOT INTEL_COMPILER)
  IF(CMAKE_COMPILER_IS_GNUCXX) 
    ADD_DEFINITIONS(-Drestrict=__restrict__ -DADD_)
    SET(CMAKE_CXX_FLAGS "-O3 -ftemplate-depth-60 -Drestrict=__restrict__ -fstrict-aliasing -funroll-all-loops   -finline-limit=1000 -ffast-math -Wno-deprecated")
#    SET(CMAKE_CXX_FLAGS "-g")
#    SET(FORTRAN_LIBS "-lg2c")
    SET(F77 g77)
    SET(F77FLAGS  -funroll-loops -O3)
  ENDIF(CMAKE_COMPILER_IS_GNUCXX) 
ENDIF(NOT INTEL_COMPILER)
  
IF(APPLE)
  INCLUDE_DIRECTORIES(/sw/include)
ENDIF(APPLE)
