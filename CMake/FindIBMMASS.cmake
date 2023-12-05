# Simple file to find IBM MASS (if available)
include(CheckCXXSourceCompiles)

message(STATUS "Looking for IBM MASS libraries")

# Finding MASS_INCLUDE
find_path(
  MASS_INCLUDE mass.h
  HINTS ${MASS_ROOT}
  PATH_SUFFIXES include)
message(STATUS "MASS_INCLUDE: ${MASS_INCLUDE}")

# Finding and setting the MASS_LIBRARY
set(SUFFIXES lib lib64)
find_library(
  MASS_LIBRARY mass
  HINTS ${MASS_ROOT}
  PATH_SUFFIXES ${SUFFIXES})
message(STATUS "MASS_LIBRARY: ${MASS_LIBRARY}")

if(MASS_INCLUDE AND MASS_LIBRARY)
  add_library(IBM::MASS INTERFACE IMPORTED)
  target_compile_definitions(IBM::MASS INTERFACE "HAVE_MASS")
  target_include_directories(IBM::MASS INTERFACE ${MASS_INCLUDE})
  target_link_libraries(IBM::MASS INTERFACE ${MASS_LIBRARY})
  # Check if MASS works with the compiler
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
    LINK_LIBRARIES IBM::MASS
    OUTPUT_VARIABLE MASS_OUT)
  if(NOT HAVE_MASS)
    message("${MASS_OUT}")
  endif(NOT HAVE_MASS)
endif()

# Finding MASSV_INCLUDE
find_path(
  MASSV_INCLUDE massv.h
  HINTS ${MASS_ROOT}
  PATH_SUFFIXES include)
message(STATUS "MASSV_INCLUDE: ${MASSV_INCLUDE}")

# Finding and setting the MASSV_LIBRARY
set(SUFFIXES lib lib64)
find_library(
  MASSV_LIBRARY massv
  HINTS ${MASS_ROOT}
  PATH_SUFFIXES ${SUFFIXES})
message(STATUS "MASSV_LIBRARY: ${MASSV_LIBRARY}")

if(MASSV_INCLUDE AND MASSV_LIBRARY)
  add_library(IBM::MASSV INTERFACE IMPORTED)
  target_compile_definitions(IBM::MASSV INTERFACE "HAVE_MASSV")
  target_include_directories(IBM::MASSV INTERFACE ${MASSV_INCLUDE})
  target_link_libraries(IBM::MASSV INTERFACE ${MASSV_LIBRARY})
  # Check if MASSV works with the compiler
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
    LINK_LIBRARIES IBM::MASSV
    OUTPUT_VARIABLE MASSV_OUT)
  if(NOT HAVE_MASSV)
    message("${MASSV_OUT}")
  endif(NOT HAVE_MASSV)
endif()

if(HAVE_MASS OR HAVE_MASSV)
  set(MASS_FOUND TRUE)
  message(STATUS "MASS found: HAVE_MASS=${HAVE_MASS}, HAVE_MASSV=${HAVE_MASSV}")
  if(HAVE_MASS)
    target_link_libraries(Math::scalar_vector_functions INTERFACE IBM::MASS)
    set(SINCOS_INCLUDE mass.h)
  endif()
  if(HAVE_MASSV)
    target_link_libraries(Math::scalar_vector_functions INTERFACE IBM::MASSV)
  endif()
else()
  set(MASS_FOUND FALSE)
  message(STATUS "MASS not found")
endif()
