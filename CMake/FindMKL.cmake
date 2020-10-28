# Simple file to refine MKL search.
# It relies on FindLAPACK to locate MKL library files and set up linker options first.
INCLUDE( CheckCXXSourceCompiles )

SET(MKL_LIBRARIES ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})

MESSAGE(STATUS "Looking for Intel MKL library header files")

# Finding and setting the MKL_INCLUDE_DIRECTORIES based on MKL_ROOT, $ENV{MKLROOT}, $ENV{MKL_ROOT}, $ENV{MKL_HOME}
# Extremely Basic Support of common mkl module environment variables
FIND_PATH(MKL_INCLUDE_DIRECTORIES "mkl.h"
  HINTS ${MKL_ROOT} $ENV{MKLROOT} $ENV{MKL_ROOT} $ENV{MKL_HOME}
  PATH_SUFFIXES include)
IF(NOT MKL_INCLUDE_DIRECTORIES)
  # Finding MKL headers in the system
  FIND_PATH(MKL_INCLUDE_DIRECTORIES "mkl.h" PATH_SUFFIXES mkl)
ENDIF()

IF(MKL_INCLUDE_DIRECTORIES)
  MESSAGE(STATUS "MKL_INCLUDE_DIRECTORIES: ${MKL_INCLUDE_DIRECTORIES}")
ELSE(MKL_INCLUDE_DIRECTORIES)
  MESSAGE(STATUS "mkl.h cannot be found")
  IF(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    MESSAGE(FATAL_ERROR "Intel's standard compilervar.sh sets the env variable MKLROOT.\n"
            "If you are invoking icc without the customary environment\n"
            "you must set the the environment variable or pass cmake MKL_ROOT.")
  ELSE(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    MESSAGE(FATAL_ERROR "Pass mkl root directory to cmake via MKL_ROOT." )
  ENDIF(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
ENDIF(MKL_INCLUDE_DIRECTORIES)

# Check for mkl.h
FILE( WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl.cxx"
  "#include <iostream>\n #include <mkl.h>\n int main() { return 0; }\n" )
TRY_COMPILE(HAVE_MKL ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl.cxx
  CMAKE_FLAGS
  "-DINCLUDE_DIRECTORIES=${MKL_INCLUDE_DIRECTORIES} "
  LINK_LIBRARIES "${MKL_LIBRARIES}"
  OUTPUT_VARIABLE MKL_OUT)
IF( NOT HAVE_MKL )
  MESSAGE( "${MKL_OUT}" )
ENDIF( NOT HAVE_MKL )

# Check for mkl_vml_functions.h
FILE( WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_vml.cxx"
  "#include <iostream>\n #include <mkl_vml_functions.h>\n int main() { return 0; }\n" )
TRY_COMPILE(HAVE_MKL_VML ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_vml.cxx
  CMAKE_FLAGS
  "-DINCLUDE_DIRECTORIES=${MKL_INCLUDE_DIRECTORIES} "
  OUTPUT_VARIABLE MKL_OUT)

# Check for fftw3
FIND_PATH(MKL_FFTW3 "fftw3.h"
  HINTS ${MKL_INCLUDE_DIRECTORIES}
  PATH_SUFFIXES fftw)
IF(MKL_FFTW3)
  FILE( WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_fftw3.cxx"
    "#include <iostream>\n #include <fftw3.h>\n int main() { return 0; }\n" )
  TRY_COMPILE(HAVE_MKL_FFTW3 ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_fftw3.cxx
    CMAKE_FLAGS
    "-DINCLUDE_DIRECTORIES=${MKL_FFTW3} "
    LINK_LIBRARIES "${MKL_LIBRARIES}"
    OUTPUT_VARIABLE MKL_OUT)
ELSE(MKL_FFTW3)
  UNSET(HAVE_MKL_FFTW3)
ENDIF(MKL_FFTW3)

IF ( HAVE_MKL )
  SET( MKL_FOUND TRUE )
  MESSAGE(STATUS "MKL found: HAVE_MKL=${HAVE_MKL}, HAVE_MKL_VML=${HAVE_MKL_VML}, HAVE_MKL_FFTW3=${HAVE_MKL_FFTW3}")

  #Add BLAS_LAPACK header
  SET_TARGET_PROPERTIES(Math::BLAS_LAPACK PROPERTIES INTERFACE_COMPILE_DEFINITIONS "HAVE_MKL"
                                                     INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIRECTORIES}")

  IF( HAVE_MKL_VML )
    #create scalar_vector_functions target
    ADD_LIBRARY(Math::scalar_vector_functions INTERFACE IMPORTED)
    SET_TARGET_PROPERTIES(Math::scalar_vector_functions PROPERTIES INTERFACE_COMPILE_DEFINITIONS "HAVE_MKL;HAVE_MKL_VML"
                                                                   INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIRECTORIES}"
                                                                   INTERFACE_LINK_LIBRARIES "${MKL_LIBRARIES}")
  ENDIF( HAVE_MKL_VML )

  IF( HAVE_MKL_FFTW3 )
    #create FFTW3 target
    ADD_LIBRARY(Math::FFTW3 INTERFACE IMPORTED)
    SET_TARGET_PROPERTIES(Math::FFTW3 PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MKL_FFTW3}"
                                                 INTERFACE_COMPILE_DEFINITIONS "HAVE_MKL;HAVE_LIBFFTW"
                                                 INTERFACE_LINK_LIBRARIES "${MKL_LIBRARIES}")
  ENDIF( HAVE_MKL_FFTW3 )
ELSE( HAVE_MKL )
  SET( MKL_FOUND FALSE )
  MESSAGE(STATUS "MKL header files not found")
ENDIF( HAVE_MKL )
