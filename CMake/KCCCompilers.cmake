#
#KCC or guidec++ with IBM Visual Age C/C++ backends
#JK cannot access other platforms with KCC or guidec++
#
MESSAGE(STATUS "AIX system with KCC or KAP")

SET(F77 xlf)
#####################################################
#guidec++ only works with -q64
#####################################################
IF(OHMMS_OMP)
  MESSAGE(STATUS "OpenMp is enabled. Make sure CXX=guidec++.")
  SET(QMC_BITS 64)
  SET(F77 guidef77)
ENDIF(OHMMS_OMP)

SET(KCC_BASE_FLAGS "+K3 -O --inline_keyword_space_time=10000 --backend -qassert=typ --backend -qspill=2000 -qmaxmem=-1 --strict --restrict --no_exceptions --one_per --no_implicit_include")


IF(QMC_BITS MATCHES 64)
  ADD_DEFINITIONS(-D_LONG_LONG)
  SET(F77OPTFLAGS  -O3 -qstrict -qarch=pwr4 -qtune=pwr4 -q64)
  SET(CMAKE_CXX_FLAGS "-q64 -qlonglong ${KCC_BASE_FLAGS} --diag_suppress 450")
  SET(AR_OPTIONS "-X64")
ELSE(QMC_BITS MATCHES 64)
  SET(F77OPTFLAGS  -O3 -qstrict -qarch=pwr4 -qtune=pwr4)
  SET(CMAKE_CXX_FLAGS ${KCC_BASE_FLAGS})
ENDIF(QMC_BITS MATCHES 64)

SET(F77FLAGS ${F77OPTFLAGS})
SET(FORTRAN_LIBS " -lxlf90 -lxlf")

######################################################
#KCC needs to be used to build static libraries
######################################################
SET(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_FLAGS})
SET(CMAKE_CXX_CREATE_STATIC_LIBRARY
   "<CMAKE_CXX_COMPILER> <CMAKE_CXX_LINK_FLAGS> -o <TARGET> <OBJECTS>")
