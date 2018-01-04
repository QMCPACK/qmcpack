# Simple file to find MKL (if availible)
# This needs a lot of work to make it robust
INCLUDE( CheckCXXSourceCompiles )

# Check for mkl.h
FILE( WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl.cxx"
    "#include <iostream>\n #include <mkl.h>\n int main() { return 0; }\n" )
try_compile(HAVE_MKL ${CMAKE_BINARY_DIR}
      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl.cxx
      COMPILE_DEFINITIONS "-mkl" )

# Check for mkl_vml_functions.h
FILE( WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_vml.cxx"
    "#include <iostream>\n #include <mkl_vml_functions.h>\n int main() { return 0; }\n" )
try_compile(HAVE_MKL_VML ${CMAKE_BINARY_DIR}
      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_vml.cxx
      COMPILE_DEFINITIONS "-mkl" )

# Check for fftw3
FILE( WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_fftw3.cxx"
    "#include <iostream>\n #include <fftw/fftw3.h>\n int main() { return 0; }\n" )
try_compile(HAVE_MKL_FFTW3 ${CMAKE_BINARY_DIR}
      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_fftw3.cxx
      COMPILE_DEFINITIONS "-mkl" )

IF ( HAVE_MKL )
    SET( MKL_FOUND 1 )
    SET( MKL_FLAGS -mkl )
    SET( MKL_LIBRARIES )
    SET( MKL_LINKER_FLAGS -mkl )
    IF ( HAVE_MKL_FFTW3 )
        FILE(WRITE "${CMAKE_CURRENT_BINARY_DIR}/include/fftw3.h" "#include <fftw/fftw3.h>\n" )
        INCLUDE_DIRECTORIES( "${CMAKE_CURRENT_BINARY_DIR}/include" )
    ENDIF()
    MESSAGE(STATUS "MKL found: HAVE_MKL=${HAVE_MKL}, HAVE_MKL_VML=${HAVE_MKL_VML}, HAVE_MKL_FFTW3=${HAVE_MKL_FFTW3}")
ELSE()
    SET( MKL_FOUND 0 )
    SET( MKL_FLAGS )
    SET( MKL_LIBRARIES )
    SET( MKL_LINKER_FLAGS )
    MESSAGE(STATUS "MKL not found")
ENDIF()

