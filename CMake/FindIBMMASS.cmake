# Simple file to find IBM MASS (if available)
include(CheckCXXSourceCompiles)

message(STATUS "Looking for IBM MASS libraries")
# Finding and setting the MASS_INCLUDE_DIRECTORIES
set(SUFFIXES include)
find_path(
  MASS_INCLUDE_DIRECTORIES name "mass.h"
  HINTS ${MASS_ROOT}
  PATH_SUFFIXES ${SUFFIXES})
if(NOT MASS_INCLUDE_DIRECTORIES)
  message(FATAL_ERROR "MASS_INCLUDE_DIRECTORIES not set. \"mass.h\" not found in MASS_ROOT/(${SUFFIXES})")
endif(NOT MASS_INCLUDE_DIRECTORIES)
message(STATUS "MASS_INCLUDE_DIRECTORIES: ${MASS_INCLUDE_DIRECTORIES}")

# Finding and setting the MASS_LINK_DIRECTORIES
# the directory organization varies with platform and targets
# these suffixes are not exhaustive
set(MASS_FIND_LIB "libmass${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(SUFFIXES lib lib64)
find_path(
  MASS_LINK_DIRECTORIES name "${MASS_FIND_LIB}"
  HINTS ${MASS_ROOT}
  PATH_SUFFIXES ${SUFFIXES})
if(NOT MASS_LINK_DIRECTORIES)
  message(FATAL_ERROR "MASS_LINK_DIRECTORIES not set. ${MASS_FIND_LIB} " "not found in MASS_ROOT/(${SUFFIXES})")
endif(NOT MASS_LINK_DIRECTORIES)
message(STATUS "MASS_LINK_DIRECTORIES: ${MASS_LINK_DIRECTORIES}")

set(MASS_LINKER_FLAGS -L${MASS_LINK_DIRECTORIES} -Wl,-rpath,${MASS_LINK_DIRECTORIES})

set(MASS_LIBRARY "-lmass")
file(
  WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mass.cxx"
  "#include <cmath>
#include <mass.h>
#include <iostream>

int main(void) {
double input = 1.1;
double outsin, outcos;;
sincos(input, &outsin, &outcos);
}
")

try_compile(
  HAVE_MASS ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mass.cxx
  CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${MASS_INCLUDE_DIRECTORIES} " "-DLINK_DIRECTORIES=${MASS_LINK_DIRECTORIES}"
  LINK_LIBRARIES "${MASS_LIBRARY}"
  OUTPUT_VARIABLE MASS_OUT)
if(NOT HAVE_MASS)
  message("${MASS_OUT}")
endif(NOT HAVE_MASS)

set(MASSV_LIBRARY "-lmassv")
file(
  WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_massv.cxx"
  "#include <massv.h>
#include <iostream>

int main(void) {
int in_size = 10;
double inputv[in_size];
double resultv[in_size];
for( int i = 0; i < in_size; ++i)
  inputv[i] = i;
vlog10(resultv, inputv, &in_size);
}
")

try_compile(
  HAVE_MASSV ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_massv.cxx
  CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${MASS_INCLUDE_DIRECTORIES} " "-DLINK_DIRECTORIES=${MASS_LINK_DIRECTORIES}"
  LINK_LIBRARIES "${MASSV_LIBRARY}"
  OUTPUT_VARIABLE MASSV_OUT)
if(NOT HAVE_MASSV)
  message("${MASSV_OUT}")
endif(NOT HAVE_MASSV)

if(HAVE_MASS OR HAVE_MASSV)
  set(MASS_FOUND TRUE)
  message(STATUS "MASS found: HAVE_MASS=${HAVE_MASS}, HAVE_MASSV=${HAVE_MASSV}")

  target_include_directories(Math::scalar_vector_functions INTERFACE "${MASS_INCLUDE_DIRECTORIES}")
  target_link_options(Math::scalar_vector_functions INTERFACE "${MASS_LINKER_FLAGS}")

  if(HAVE_MASS)
    target_compile_definitions(Math::scalar_vector_functions INTERFACE "HAVE_MASS")
    target_link_libraries(Math::scalar_vector_functions INTERFACE "${MASS_LIBRARY}")
    set(SINCOS_INCLUDE mass.h)
  endif()

  if(HAVE_MASSV)
    target_compile_definitions(Math::scalar_vector_functions INTERFACE "HAVE_MASSV")
    target_link_libraries(Math::scalar_vector_functions INTERFACE "${MASSV_LIBRARY}")
  endif()

else()
  set(MASS_FOUND FALSE)
  message(STATUS "MASS not found")
endif()
