# Simple file to find IBM MASS (if available)
INCLUDE( CheckCXXSourceCompiles )

MESSAGE(STATUS "Looking for IBM MASS libraries")
# Finding and setting the MASS_INCLUDE_DIRECTORIES
set(SUFFIXES include)
find_path(MASS_INCLUDE_DIRECTORIES name "mass.h" HINTS ${MASS_ROOT}
  PATH_SUFFIXES ${SUFFIXES})
if (NOT MASS_INCLUDE_DIRECTORIES)
  message(FATAL_ERROR "MASS_INCLUDE_DIRECTORIES not set. \"mass.h\" not found in MASS_ROOT/(${SUFFIXES})")
endif (NOT MASS_INCLUDE_DIRECTORIES)
message("MASS_INCLUDE_DIRECTORIES: ${MASS_INCLUDE_DIRECTORIES}")

  # Finding and setting the MASS_LINK_DIRECTORIES
  # the directory organization varies with platform and targets
  # these suffixes are not exhaustive
set(MASS_FIND_LIB "libmass${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(SUFFIXES lib lib64)
find_path(MASS_LINK_DIRECTORIES name "${MASS_FIND_LIB}" HINTS ${MASS_ROOT}
  PATH_SUFFIXES ${SUFFIXES})
if (NOT MASS_LINK_DIRECTORIES)
  message(FATAL_ERROR "MASS_LINK_DIRECTORIES not set. ${MASS_FIND_LIB} "
    "not found in MASS_ROOT/(${SUFFIXES})")
endif (NOT MASS_LINK_DIRECTORIES)
message("MASS_LINK_DIRECTORIES: ${MASS_LINK_DIRECTORIES}")

set(MASS_LINKER_FLAGS -L${MASS_LINK_DIRECTORIES} -Wl,-rpath,${MASS_LINK_DIRECTORIES})
set(MASS_LIBRARIES "-lmass -lmassv")

FILE( WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mass.cxx"
  "#include <cmath>
#include <mass.h>
#include <massv.h>
#include <iostream>

int main(void) {
double input = 1.1;
double outsin, outcos;;
sincos(input,&outsin, &outcos);
int in_size = 10;
double inputv[in_size];
double resultv[in_size];
for( int i = 0; i < in_size; ++i)
  inputv[i] = i;
vlog10(resultv, inputv, &in_size);
}
")
try_compile(HAVE_MASS ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mass.cxx
  CMAKE_FLAGS
  "-DINCLUDE_DIRECTORIES=${MASS_INCLUDE_DIRECTORIES} "
  "-DLINK_DIRECTORIES=${MASS_LINK_DIRECTORIES}"
  LINK_LIBRARIES "${MASS_LIBRARIES}"
  COMPILE_DEFINITIONS "${MASS_COMPILE_DEFINITIONS}"
  OUTPUT_VARIABLE MASS_OUT)
if ( NOT HAVE_MASS )
  MESSAGE( "${MASS_OUT}" )
endif ( NOT HAVE_MASS )

IF ( HAVE_MASS )
  SET( MASS_FOUND 1 )
  SET( MASS_FLAGS ${MASS_COMPILE_DEFINITIONS} )
  include_directories( ${MASS_INCLUDE_DIRECTORIES} )
  set( HAVE_VECTOR_MATH 1 )
  set( HAVE_MASSV 1 )
  set( SINCOS_INCLUDE mass.h )
  MESSAGE(STATUS "MASS found: HAVE_MASS=${HAVE_MASS}, HAVE_MASSV=${HAVE_MASSV}")
ELSE( HAVE_MASS )
  SET( MASS_FOUND 0 )
  SET( MASS_FLAGS )
  SET( MASS_LIBRARIES )
  SET( MASS_LINKER_FLAGS )
  MESSAGE("MASS not found")
ENDIF( HAVE_MASS )
