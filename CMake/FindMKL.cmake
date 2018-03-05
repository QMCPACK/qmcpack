# Simple file to find MKL (if availible)
# This needs a lot of work to make it robust
INCLUDE( CheckCXXSourceCompiles )
CMAKE_MINIMUM_REQUIRED(VERSION 3.2.0)
IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0056 NEW)
ENDIF(COMMAND cmake_policy)

# if MKL_ROOT is set, use that
if ( MKL_ROOT )
   message ("MKL_ROOT: ${MKL_ROOT}")
   set(MKL_FLAGS "-m64")
   set(MKL_LINK_DIRECTORIES "${MKL_ROOT}/lib")
   set(MKL_LINKER_FLAGS "-L${MKL_LINK_DIRECTORIES} -Wl,-rpath,${MKL_LINK_DIRECTORIES}")
   set(MKL_INCLUDE_DIRECTORIES "${MKL_ROOT}/include")
   set(MKL_FLAGS "-I${MKL_INCLUDE_DIRECTORIES}")
   set(MKL_LIBRARIES "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")
else ( MKL_ROOT )
# this takes away build control but is what the cmake does now
   string(TOLOWER CMAKE_CXX_COMPILER_ID LCASE_CMAKE_CXX_COMPILER_ID)
   if ( ${LCASE_CMAKE_CXX_COMPILER_ID} MATCHES "intel" )
     set(MKL_LIBRARIES "-mkl")
     set(MKL_FLAGS "-mkl")
   endif ( ${LCASE_CMAKE_CXX_COMPILER_ID} MATCHES "intel" )
endif ( MKL_ROOT )

set( org_CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" )
if ( CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin" AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang" )
  string(REPLACE "-fopenmp" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
endif ( CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin" AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang" )

 
# Check for mkl.h
FILE( WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl.cxx"
    "#include <iostream>\n #include <mkl.h>\n int main() { return 0; }\n" )
try_compile(HAVE_MKL ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl.cxx
  CMAKE_FLAGS ${MKL_FLAGS}
  "-DINCLUDE_DIRECTORIES=${MKL_INCLUDE_DIRECTORIES} "
  "-DLINK_DIRECTORIES=${MKL_LINK_DIRECTORIES}"
  LINK_LIBRARIES "${MKL_LIBRARIES}"
  OUTPUT_VARIABLE MKL_OUT)

if ( NOT HAVE_MKL )
  MESSAGE( "${MKL_OUT}" )
endif ( NOT HAVE_MKL )

# Check for mkl_vml_functions.h
FILE( WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_vml.cxx"
    "#include <iostream>\n #include <mkl_vml_functions.h>\n int main() { return 0; }\n" )
try_compile(HAVE_MKL_VML ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_vml.cxx
  CMAKE_FLAGS ${MKL_FLAGS}
  "-DINCLUDE_DIRECTORIES=${MKL_INCLUDE_DIRECTORIES} "
  "-DLINK_DIRECTORIES=${MKL_LINK_DIRECTORIES}"
  OUTPUT_VARIABLE MKL_OUT)

# Check for fftw3
FILE( WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_fftw3.cxx"
    "#include <iostream>\n #include <fftw/fftw3.h>\n int main() { return 0; }\n" )
try_compile(HAVE_MKL_FFTW3 ${CMAKE_BINARY_DIR}
      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_fftw3.cxx
      CMAKE_FLAGS ${MKL_FLAGS}
      "-DINCLUDE_DIRECTORIES=${MKL_INCLUDE_DIRECTORIES} "
      "-DLINK_DIRECTORIES=${MKL_LINK_DIRECTORIES}"
      LINK_LIBRARIES "${MKL_LIBRARIES}"
      OUTPUT_VARIABLE MKL_OUT)

set( CMAKE_CXX_FLAGS "${org_CMAKE_CXX_FLAGS}" )
    
IF ( HAVE_MKL )
    SET( MKL_FOUND 1 )
    SET( MKL_FLAGS ${MKL_FLAGS} )
    # if ( ${LCASE_CMAKE_CXX_COMPILER_ID} MATCHES "intel" )
    #    SET( MKL_LIBRARIES )
    #endif()
    IF ( HAVE_MKL_FFTW3 )
        FILE(WRITE "${CMAKE_CURRENT_BINARY_DIR}/include/fftw3.h" "#include <fftw/fftw3.h>\n" )
        INCLUDE_DIRECTORIES( "${CMAKE_CURRENT_BINARY_DIR}/include" )
    ENDIF()
    MESSAGE(STATUS "MKL found: HAVE_MKL=${HAVE_MKL}, HAVE_MKL_VML=${HAVE_MKL_VML}, HAVE_MKL_FFTW3=${HAVE_MKL_FFTW3}")
ELSE( HAVE_MKL )
    SET( MKL_FOUND 0 )
    SET( MKL_FLAGS )
    SET( MKL_LIBRARIES )
    SET( MKL_LINKER_FLAGS )
    MESSAGE(STATUS "MKL not found")
ENDIF( HAVE_MKL )


