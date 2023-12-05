# FFTW_INCLUDE_DIR = fftw3.h
# FFTW_LIBRARIES = libfftw3.a
# FFTW_FOUND = true if FFTW3 is found

set(Libfftw fftw3)
if(QMC_BUILD_STATIC)
  set(Libfftw libfftw3.a)
endif(QMC_BUILD_STATIC)

if(FFTW_INCLUDE_DIRS)
  find_path(FFTW_INCLUDE_DIR fftw3.h ${FFTW_INCLUDE_DIRS})
  find_library(FFTW_LIBRARIES ${Libfftw} ${FFTW_LIBRARY_DIRS})
else(FFTW_INCLUDE_DIRS)
  find_path(FFTW_INCLUDE_DIR fftw3.h ${FFTW_HOME}/include $ENV{FFTW_HOME}/include)
  find_library(FFTW_LIBRARIES ${Libfftw} ${FFTW_HOME}/lib $ENV{FFTW_HOME}/lib)
endif(FFTW_INCLUDE_DIRS)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW REQUIRED_VARS FFTW_LIBRARIES FFTW_INCLUDE_DIR)

if(FFTW_FOUND)
  message(STATUS "FFTW_INCLUDE_DIR=${FFTW_INCLUDE_DIR}")
  message(STATUS "FFTW_LIBRARIES=${FFTW_LIBRARIES}")
  set(FFTW_FOUND TRUE)
  add_library(FFTW::FFTW3 INTERFACE IMPORTED)
  target_include_directories(FFTW::FFTW3 INTERFACE "${FFTW_INCLUDE_DIR}")
  target_link_libraries(FFTW::FFTW3 INTERFACE "${FFTW_LIBRARIES}")
endif()

mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIBRARIES FFTW_FOUND)
