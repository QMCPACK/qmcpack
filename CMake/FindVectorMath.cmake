#
# Find vectorized math libraries and enable appropriate flags
#
# Initial version only for MKL VML. Works for gcc+MKL case. libm and massv detection required.
#

SET( HAVE_VECTOR_MATH 0 )

IF ( HAVE_MKL_VML )
# We arrive here if MKL was detected earlier by FindMKL
  SET ( HAVE_VECTOR_MATH 1 )
  MESSAGE("Using MKL Vector Math functions")
ENDIF ()

IF ( NOT HAVE_VECTOR_MATH )
#MESSAGE(STATUS "Trying MKL VML")
# Check for mkl_vml_functions.h
  FILE( WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_vml.cxx"
      "#include <iostream>\n #include <mkl_vml_functions.h>\n int main() { return 0; }\n" )
  try_compile(HAVE_MKL_VML ${CMAKE_BINARY_DIR}
        ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_vml.cxx
        CMAKE_FLAGS "${CMAKE_CXX_FLAGS}" )
  IF (HAVE_MKL_VML)
    SET ( HAVE_VECTOR_MATH 1 )
    MESSAGE(STATUS "Using MKL Vector Math functions")
  ENDIF()
ENDIF()

IF ( NOT HAVE_VECTOR_MATH )
  MESSAGE(STATUS "No vector math library detected.")
ENDIF()


