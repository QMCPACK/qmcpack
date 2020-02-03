# FFTW_INCLUDE_DIR = fftw3.h
# FFTW_LIBRARIES = libfftw3.a
# FFTW_FOUND = true if FFTW3 is found

SET(Libfftw fftw3)
IF(QMC_BUILD_STATIC)
  SET(Libfftw libfftw3.a)
ENDIF(QMC_BUILD_STATIC)

IF(FFTW_INCLUDE_DIRS)
  FIND_PATH(FFTW_INCLUDE_DIR fftw3.h  ${FFTW_INCLUDE_DIRS})
  FIND_LIBRARY(FFTW_LIBRARIES ${Libfftw} ${FFTW_LIBRARY_DIRS})
ELSE(FFTW_INCLUDE_DIRS)
  FIND_PATH(FFTW_INCLUDE_DIR fftw3.h ${FFTW_HOME}/include $ENV{FFTW_HOME}/include)
  FIND_LIBRARY(FFTW_LIBRARIES ${Libfftw} ${FFTW_HOME}/lib $ENV{FFTW_HOME}/lib) 
ENDIF(FFTW_INCLUDE_DIRS)

SET(FFTW_FOUND FALSE)
IF(FFTW_INCLUDE_DIR AND FFTW_LIBRARIES)
  MESSAGE(STATUS "FFTW_INCLUDE_DIR=${FFTW_INCLUDE_DIR}")
  MESSAGE(STATUS "FFTW_LIBRARIES=${FFTW_LIBRARIES}")
  SET(FFTW_FOUND TRUE)
  #create FFTW3 target
  ADD_LIBRARY(Math::FFTW3 INTERFACE IMPORTED)
  SET_TARGET_PROPERTIES(Math::FFTW3 PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIR}"
                                              INTERFACE_LINK_LIBRARIES "${FFTW_LIBRARIES}")
ENDIF()

MARK_AS_ADVANCED(
   FFTW_INCLUDE_DIR
   FFTW_LIBRARIES
   FFTW_FOUND
)
