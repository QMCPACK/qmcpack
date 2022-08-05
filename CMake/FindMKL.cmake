# Simple file to refine MKL search.
# It relies on FindLAPACK to locate MKL library files and set up linker options first.
include(CheckCXXSourceCompiles)

set(MKL_LIBRARIES ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})

message(STATUS "Looking for Intel MKL library header files")

# Finding and setting the MKL_INCLUDE based on MKL_ROOT, $ENV{MKLROOT}, $ENV{MKL_ROOT}, $ENV{MKL_HOME}
# Extremely Basic Support of common mkl module environment variables
find_path(
  MKL_INCLUDE "mkl.h"
  HINTS ${MKL_ROOT} $ENV{MKLROOT} $ENV{MKL_ROOT} $ENV{MKL_HOME}
  PATH_SUFFIXES include)
if(NOT MKL_INCLUDE)
  # Finding MKL headers in the system
  find_path(MKL_INCLUDE "mkl.h" PATH_SUFFIXES mkl)
endif()

if(MKL_INCLUDE)
  message(STATUS "MKL_INCLUDE: ${MKL_INCLUDE}")
else(MKL_INCLUDE)
  message(STATUS "mkl.h cannot be found")
  if(COMPILER MATCHES "Intel")
    message(
      FATAL_ERROR
        "Intel's standard compilervar.sh sets the env variable MKLROOT.\n"
        "If you are invoking icc without the customary environment\n"
        "you must set the the environment variable or pass cmake MKL_ROOT.")
  else(COMPILER MATCHES "Intel")
    message(FATAL_ERROR "Pass mkl root directory to cmake via MKL_ROOT.")
  endif(COMPILER MATCHES "Intel")
endif(MKL_INCLUDE)

# Check for mkl.h
file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl.cxx"
     "#include <iostream>\n #include <mkl.h>\n int main() { return 0; }\n")
try_compile(
  HAVE_MKL ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl.cxx
  CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${MKL_INCLUDE} "
  LINK_LIBRARIES "${MKL_LIBRARIES}"
  OUTPUT_VARIABLE MKL_OUT)
if(NOT HAVE_MKL)
  message("${MKL_OUT}")
endif(NOT HAVE_MKL)

# Check for fftw3
find_path(
  MKL_FFTW3 "fftw3.h"
  HINTS ${MKL_INCLUDE}
  PATH_SUFFIXES fftw)
if(MKL_FFTW3)
  file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_fftw3.cxx"
       "#include <iostream>\n #include <fftw3.h>\n int main() { return 0; }\n")
  try_compile(
    HAVE_MKL_FFTW3 ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_mkl_fftw3.cxx
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${MKL_FFTW3} "
    LINK_LIBRARIES "${MKL_LIBRARIES}"
    OUTPUT_VARIABLE MKL_OUT)
else(MKL_FFTW3)
  unset(HAVE_MKL_FFTW3)
endif(MKL_FFTW3)

if(HAVE_MKL)
  set(MKL_FOUND TRUE)
  message(STATUS "MKL found: HAVE_MKL=${HAVE_MKL}, HAVE_MKL_FFTW3=${HAVE_MKL_FFTW3}")

  #Add BLAS_LAPACK header
  add_library(MKL::BLAS_LAPACK INTERFACE IMPORTED)
  target_compile_definitions(MKL::BLAS_LAPACK INTERFACE "HAVE_MKL")
  target_include_directories(MKL::BLAS_LAPACK INTERFACE "${MKL_INCLUDE}")

  if(HAVE_MKL_FFTW3)
    add_library(MKL::FFTW3 INTERFACE IMPORTED)
    target_compile_definitions(MKL::FFTW3 INTERFACE "HAVE_MKL")
    target_include_directories(MKL::FFTW3 INTERFACE "${MKL_FFTW3}")
    target_link_libraries(MKL::FFTW3 INTERFACE "${MKL_LIBRARIES}")
  endif(HAVE_MKL_FFTW3)
else(HAVE_MKL)
  set(MKL_FOUND FALSE)
  message(STATUS "MKL header files not found")
endif(HAVE_MKL)

# check for mkl_sycl
if(HAVE_MKL AND ENABLE_SYCL)
  find_library(MKL_SYCL mkl_sycl
    HINTS ${MKL_ROOT} $ENV{MKLROOT} $ENV{MKL_ROOT} $ENV{MKL_HOME}
    PATH_SUFFIXES lib/intel64
    REQUIRED
  )

  if(MKL_SYCL)
    add_library(MKL::sycl INTERFACE IMPORTED)
    target_include_directories(MKL::sycl INTERFACE "${MKL_INCLUDE}")
    target_link_libraries(MKL::sycl INTERFACE ${MKL_SYCL})
  endif()
endif()
