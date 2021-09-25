# Simple file to verify MKL VML
# Input: MKL_FOUND, MKL_INCLUDE, MKL_LIBRARIES
# Output: HAVE_MKL_VML

message(STATUS "Setting up Intel MKL Vector Math Library (VML)")

if(NOT MKL_FOUND)
  message(FATAL_ERROR "MKL was not found. Cannot proceed with Vector Math Library (VML) searching!")
endif()

include(CheckIncludeFileCXX)

# Check for mkl_vml_functions.h
set(CMAKE_REQUIRED_INCLUDES "${MKL_INCLUDE}")
check_include_file_cxx(mkl_vml_functions.h HAVE_MKL_VML)
unset(CMAKE_REQUIRED_INCLUDES)

if(HAVE_MKL_VML)
  target_compile_definitions(Math::scalar_vector_functions INTERFACE "HAVE_MKL;HAVE_MKL_VML")
  target_include_directories(Math::scalar_vector_functions INTERFACE "${MKL_INCLUDE}")
  target_link_libraries(Math::scalar_vector_functions INTERFACE "${MKL_LIBRARIES}")
  message(STATUS "MKL VML found")
else()
  message(FATAL_ERROR "MKL VML not found")
endif(HAVE_MKL_VML)
