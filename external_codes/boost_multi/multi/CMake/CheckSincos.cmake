include(CheckCXXSourceCompiles)

set(CMAKE_REQUIRED_LIBRARIES Math::scalar_vector_functions)
set(SINCOS_TEST_SRC
    "#include \"${SINCOS_INCLUDE}\"
int main(void) {
double input = 1.1;
double outsin, outcos;;
sincos(input,&outsin, &outcos);
}
")

check_cxx_source_compiles("${SINCOS_TEST_SRC}" HAVE_SINCOS)
unset(CMAKE_REQUIRED_LIBRARIES)
