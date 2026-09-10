find_path(MKL_SYCL_INCLUDE_DIR NAMES oneapi/mkl/spec.hpp
  HINTS
    "${MKL_ROOT}/include"
    "$ENV{MKLROOT}/include"
    "$ENV{MKL_ROOT}/include"
  PATH_SUFFIXES mkl
)

find_library(MKL_SYCL_LIB mkl_sycl
  HINTS $ENV{MKLROOT}
  PATH_SUFFIXES lib/intel64 lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKLsycl
  REQUIRED_VARS MKL_SYCL_LIB MKL_SYCL_INCLUDE_DIR
)

if(MKLsycl_FOUND AND NOT TARGET MKL::sycl)
  add_library(MKL::sycl INTERFACE IMPORTED)
  target_include_directories(MKL::sycl INTERFACE "${MKL_SYCL_INCLUDE_DIR}")
  target_link_libraries(MKL::sycl INTERFACE "${MKL_SYCL_LIB}")
endif()
