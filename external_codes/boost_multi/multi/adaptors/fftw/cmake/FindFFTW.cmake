
# ==================================================================================================
# This file is part of the CodeVault project. The project is licensed under Apache Version 2.0.
# CodeVault is part of the EU-project PRACE-4IP (WP7.3.C).
#
# Author(s):
#   Cedric Nugteren <cedric.nugteren@surfsara.nl>
#
# ==================================================================================================
#
# Defines the following variables:
#   FFTW_FOUND          Boolean holding whether or not the FFTW3 library was found
#   FFTW_INCLUDE_DIRS   The FFTW3 include directory
#   FFTW_LIBRARIES      The FFTW3 library
#
# In case FFTW3 is not installed in the default directory, set the FFTW_ROOT variable to point to
# the root of FFTW3, such that 'fftw3.h' can be found in $FFTW_ROOT/include. This can either be done
# using an environmental variable (e.g. export FFTW_ROOT=/path/to/fftw3) or using a CMake variable
# (e.g. cmake -DFFTW_ROOT=/path/to/fftw3 ..).
#
# ==================================================================================================

# Sets the possible install locations
set(FFTW_HINTS
  ${FFTW_ROOT}
  $ENV{FFTW_ROOT}
)
set(FFTW_PATHS
  /usr
  /usr/local
)

# Finds the include directories
find_path(FFTW_INCLUDE_DIRS
  NAMES fftw3.h
  HINTS ${FFTW_HINTS}
  PATH_SUFFIXES include api inc include/x86_64 include/x64
  PATHS ${FFTW_PATHS}
  DOC "FFTW3 include header fftw3.h"
)
mark_as_advanced(FFTW_INCLUDE_DIRS)

# Finds the library
find_library(FFTW_LIBRARIES
  NAMES fftw3
  HINTS ${FFTW_HINTS}
  PATH_SUFFIXES lib lib64 lib/x86_64 lib/x64 lib/x86 lib/Win32
  PATHS ${FFTW_PATHS}
  DOC "FFTW3 library"
)
mark_as_advanced(FFTW_LIBRARIES)

# ==================================================================================================

# Notification messages
if(NOT FFTW_INCLUDE_DIRS)
    message(STATUS "Could NOT find 'fftw3.h', install FFTW3 or set FFTW_ROOT")
endif()
if(NOT FFTW_LIBRARIES)
    message(STATUS "Could NOT find the FFTW3 library, install it or set FFTW_ROOT")
endif()

# Determines whether or not FFTW3 was found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_INCLUDE_DIRS FFTW_LIBRARIES)

# ==================================================================================================

