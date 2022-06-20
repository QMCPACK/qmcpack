# Output: HAVE_AMD_LIBM

message(STATUS "Setting up AMD LibM libraries")
message(WARNING "In limited benchmarks, AMD LibM slows down QMCPACK. "
                "Please check carefully before using it in production runs.")

include(CheckCXXSourceCompiles)

find_path(
  AMD_LIBM_INCLUDE_DIR amdlibm.h
  HINTS ${AMD_LIBM_ROOT} ${AOCL_ROOT}
  PATH_SUFFIXES include)

if(NOT AMD_LIBM_INCLUDE_DIR)
  message(FATAL_ERROR "AMD_LIBM_INCLUDE_DIR not set. Header file amdlib.h not found!")
else()
  message(STATUS "Header file amdlib.h found at  ${AMD_LIBM_INCLUDE_DIR}")
endif()

find_library(
  AMD_LIBM_LIBRARY amdlibm
  HINTS ${AMD_LIBM_ROOT} ${AOCL_ROOT}
  PATH_SUFFIXES lib)

if(NOT AMD_LIBM_LIBRARY)
  message(FATAL_ERROR "AMD_LIBM_LIBRARY not set. Library file amdlibm not found!")
else()
  message(STATUS "Library file amdlibm found as ${AMD_LIBM_LIBRARY}")
endif()

set(CMAKE_REQUIRED_INCLUDES "${AMD_LIBM_INCLUDE_DIR}")
set(CMAKE_REQUIRED_LIBRARIES "${AMD_LIBM_LIBRARY}")
check_cxx_source_compiles(
  "
#include <amdlibm.h>
int main() {
  double d_in(0), d_sin, d_cos;
  amd_sincos(d_in, &d_sin, &d_cos);
  float s_in(0), s_sin, s_cos;
  amd_sincosf(s_in, &s_sin, &s_cos);
  return 0;
}
"
  HAVE_AMD_LIBM)
unset(CMAKE_REQUIRED_INCLUDES)
unset(CMAKE_REQUIRED_LIBRARIES)

if(HAVE_AMD_LIBM)
  target_compile_definitions(Math::scalar_vector_functions INTERFACE "HAVE_AMD_LIBM")
  target_include_directories(Math::scalar_vector_functions INTERFACE "${AMD_LIBM_INCLUDE_DIR}")
  target_link_libraries(Math::scalar_vector_functions INTERFACE "${AMD_LIBM_LIBRARY}")
  message(STATUS "AMD LIBM found")
else()
  message(FATAL_ERROR "AMD LIBM not found")
endif()
