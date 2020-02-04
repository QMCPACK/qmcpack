INCLUDE(${CMAKE_ROOT}/Modules/CheckCXXSourceCompiles.cmake)

SET(CMAKE_REQUIRED_LIBRARIES Math::scalar_vector_functions)

SET( SINCOS_TEST_SRC
  "#include \"${SINCOS_INCLUDE}\"
int main(void) {
double input = 1.1;
double outsin, outcos;;
sincos(input,&outsin, &outcos);
}
")

CHECK_CXX_SOURCE_COMPILES("${SINCOS_TEST_SRC}" HAVE_SINCOS)

IF(HAVE_SINCOS)
  MESSAGE(STATUS "sincos found")
ELSE(HAVE_SINCOS)
  MESSAGE(STATUS "sincos not found")
ENDIF(HAVE_SINCOS)
