#
# Find vectorized math libraries and enable appropriate flags
#
# Initial version only for MKL VML. Works for gcc+MKL case. libm and massv detection required.
#

SET( HAVE_VECTOR_MATH 0 )

IF ( HAVE_MKL_VML )
# We arrive here if MKL was detected earlier by FindMKL
  SET ( HAVE_VECTOR_MATH 1 )
  MESSAGE(STATUS "Using MKL Vector Math functions")
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
    # enable VML only when MKL libraries have been picked up
    IF (MKL_FOUND)
      SET ( HAVE_VECTOR_MATH 1 )
      MESSAGE(STATUS "Using MKL Vector Math functions (header file check passed)")
    ELSE(MKL_FOUND)
      MESSAGE(STATUS "Drop MKL Vector Math functions (header file check passed but libraries are not available)")
      SET ( HAVE_MKL_VML 0 )
    ENDIF(MKL_FOUND)
  ELSE ()
    IF (MKL_FOUND)
      MESSAGE(WARNING "mkl_vml_functions.h check failed but MKL libraries are available. At the risk of losing performance.")
    ENDIF()
  ENDIF()
ENDIF()

IF ( NOT HAVE_VECTOR_MATH )
  MESSAGE(STATUS "No usable vector math library detected.")
ENDIF()


