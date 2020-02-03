# Simple file to find MKL (if available)
# This needs a lot of work to make it robust
INCLUDE( CheckCXXSourceCompiles )

MESSAGE(STATUS "Looking for Intel MKL libraries")

# Extremely Basic Support of common mkl module environment variables
# or -DMKLROOT/-DMKL_HOME instead of preferred -DMKL_ROOT
# Some Linux distributions (Ubuntu >= 18.10) provide MKL
# headers in /usr/include/mkl, define MKL_ROOT but not used for library searching.
if (NOT MKL_ROOT)
  # Finding and setting the MKL_INCLUDE_DIRECTORIES
  find_path(MKL_INCLUDE_DIRECTORIES "mkl.h"
    HINTS ${MKLROOT} ${MKL_HOME} $ENV{MKLROOT} $ENV{MKL_ROOT} $ENV{MKL_HOME}
    PATH_SUFFIXES include mkl)
  string(REPLACE "/include" "" MKL_ROOT ${MKL_INCLUDE_DIRECTORIES})
else (NOT MKL_ROOT)
  # Finding and setting the MKL_INCLUDE_DIRECTORIES
  set(SUFFIXES include)
  find_path(MKL_INCLUDE_DIRECTORIES "mkl.h" HINTS ${MKL_ROOT}
    PATH_SUFFIXES ${SUFFIXES})
  if (NOT MKL_INCLUDE_DIRECTORIES)
    MESSAGE(FATAL_ERROR "MKL_INCLUDE_DIRECTORIES not set. \"mkl.h\" not found in MKL_ROOT/(${SUFFIXES})")
  endif ()
endif (NOT MKL_ROOT)

if (NOT MKL_ROOT)
  if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    MESSAGE(FATAL_ERROR "Intel's standard compilervar.sh sets the env variable MKLROOT.\n"
            "If you are invoking icc without the customary environment\n"
            "you must set the the environment variable or pass cmake MKL_ROOT.")
  else(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    MESSAGE(FATAL_ERROR "ENABLE_MKL is TRUE and mkl directory not found. Set MKL_ROOT." )
  endif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
endif (NOT MKL_ROOT)

MESSAGE(STATUS "MKL_INCLUDE_DIRECTORIES: ${MKL_INCLUDE_DIRECTORIES}")

if ( NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel" )
  # Finding and setting the MKL_LINK_DIRECTORIES
  # the directory organization varies with platform and targets
  # these suffixes are not exhaustive
  set(MKL_FIND_LIB "libmkl_intel_lp64${CMAKE_SHARED_LIBRARY_SUFFIX}")
  set(SUFFIXES lib lib/intel64)
  find_library(MKL_LINK_DIRECTORIES "${MKL_FIND_LIB}" HINTS ${MKL_ROOT}
    PATH_SUFFIXES ${SUFFIXES})
  if (NOT MKL_LINK_DIRECTORIES)
    MESSAGE(FATAL_ERROR "MKL_LINK_DIRECTORIES not set. ${MKL_FIND_LIB} "
      "not found in MKL_ROOT/(${SUFFIXES})")
  endif ()
  string(REPLACE "/${MKL_FIND_LIB}" "" MKL_LINK_DIRECTORIES ${MKL_LINK_DIRECTORIES})
  MESSAGE(STATUS "MKL_LINK_DIRECTORIES: ${MKL_LINK_DIRECTORIES}")
  set(MKL_LINKER_FLAGS "-L${MKL_LINK_DIRECTORIES} -Wl,-rpath,${MKL_LINK_DIRECTORIES}")
  set(MKL_LIBRARIES "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")
else ( NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel" )
  # this takes away link control for intel but is convenient
  # perhaps we should consider dropping it since it will more or less
  # unify the MKL setup.
  # Note -mkl implicitly includes that icc's mkl/include
  set(MKL_COMPILE_DEFINITIONS "-mkl")
endif (NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel" )

# Check for mkl.h
FILE( WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl.cxx"
  "#include <iostream>\n #include <mkl.h>\n int main() { return 0; }\n" )
try_compile(HAVE_MKL ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl.cxx
  CMAKE_FLAGS
  "-DINCLUDE_DIRECTORIES=${MKL_INCLUDE_DIRECTORIES} "
  "-DLINK_DIRECTORIES=${MKL_LINK_DIRECTORIES}"
  LINK_LIBRARIES "${MKL_LIBRARIES}"
  COMPILE_DEFINITIONS "${MKL_COMPILE_DEFINITIONS}"
  OUTPUT_VARIABLE MKL_OUT)
if ( NOT HAVE_MKL )
  MESSAGE( "${MKL_OUT}" )
endif ( NOT HAVE_MKL )

# Check for mkl_vml_functions.h
FILE( WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_vml.cxx"
  "#include <iostream>\n #include <mkl_vml_functions.h>\n int main() { return 0; }\n" )
try_compile(HAVE_MKL_VML ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_vml.cxx
  CMAKE_FLAGS
  "-DINCLUDE_DIRECTORIES=${MKL_INCLUDE_DIRECTORIES} "
  "-DLINK_DIRECTORIES=${MKL_LINK_DIRECTORIES}"
  COMPILE_DEFINITIONS "${MKL_COMPILE_DEFINITIONS}"
  OUTPUT_VARIABLE MKL_OUT)

# Check for fftw3
find_path(MKL_FFTW3 "fftw3.h"
  HINTS ${MKL_INCLUDE_DIRECTORIES}
  PATH_SUFFIXES fftw)
if(MKL_FFTW3)
  list(APPEND MKL_INCLUDE_DIRECTORIES ${MKL_FFTW3})
  FILE( WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_fftw3.cxx"
    "#include <iostream>\n #include <fftw3.h>\n int main() { return 0; }\n" )
  try_compile(HAVE_MKL_FFTW3 ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_fftw3.cxx
    CMAKE_FLAGS
    "-DINCLUDE_DIRECTORIES=${MKL_INCLUDE_DIRECTORIES} "
    "-DLINK_DIRECTORIES=${MKL_LINK_DIRECTORIES}"
    LINK_LIBRARIES "${MKL_LIBRARIES}"
    COMPILE_DEFINITIONS "${MKL_COMPILE_DEFINITIONS}"
    OUTPUT_VARIABLE MKL_OUT)
else(MKL_FFTW3)
  UNSET(HAVE_MKL_FFTW3)
endif(MKL_FFTW3)

IF ( HAVE_MKL )
  SET( MKL_FOUND 1 )
  SET( MKL_FLAGS ${MKL_COMPILE_DEFINITIONS} )

  #create BLAS_LAPACK target
  ADD_LIBRARY(Math::BLAS_LAPACK INTERFACE IMPORTED)
  IF (NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel" )
    SET_TARGET_PROPERTIES(Math::BLAS_LAPACK PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIRECTORIES}"
                                                                INTERFACE_LINK_LIBRARIES "${MKL_LINKER_FLAGS} ${MKL_LIBRARIES}")
  ELSE (NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel" )
    SET_TARGET_PROPERTIES(Math::BLAS_LAPACK PROPERTIES INTERFACE_COMPILE_OPTIONS  "${MKL_COMPILE_DEFINITIONS}")
  ENDIF (NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel" )

  IF( HAVE_MKL_FFTW3 )
    #create FFTW3 target
    ADD_LIBRARY(Math::FFTW3 INTERFACE IMPORTED)
    SET_TARGET_PROPERTIES(Math::FFTW3 PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MKL_FFTW3}")
    IF (NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel" )
      SET_TARGET_PROPERTIES(Math::FFTW3 PROPERTIES INTERFACE_LINK_LIBRARIES "${MKL_LINKER_FLAGS} ${MKL_LIBRARIES}")
    ELSE (NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel" )
      SET_TARGET_PROPERTIES(Math::FFTW3 PROPERTIES INTERFACE_COMPILE_OPTIONS  "${MKL_COMPILE_DEFINITIONS}")
    ENDIF (NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel" )
  ENDIF( HAVE_MKL_FFTW3 )
  MESSAGE(STATUS "MKL found: HAVE_MKL=${HAVE_MKL}, HAVE_MKL_VML=${HAVE_MKL_VML}, HAVE_MKL_FFTW3=${HAVE_MKL_FFTW3}")
ELSE( HAVE_MKL )
  SET( MKL_FOUND 0 )
  SET( MKL_FLAGS )
  SET( MKL_LIBRARIES )
  SET( MKL_LINKER_FLAGS )
  MESSAGE("MKL not found")
ENDIF( HAVE_MKL )
