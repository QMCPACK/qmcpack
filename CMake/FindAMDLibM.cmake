message(STATUS "Looking for AMD LibM libraries")

find_path(
  AMD_LIBM_INCLUDE_DIR NAME amdlibm.h
  HINTS ${AMD_LIBM_ROOT} ${AOCL_ROOT}
  PATH_SUFFIXES include)

if(NOT AMD_LIBM_INCLUDE_DIR)
  message(FATAL_ERROR "AMD_LIBM_INCLUDE_DIR not set. Header file amdlib.h not found!")
endif()
message(STATUS "Header file amdlib.h found at  ${AMD_LIBM_INCLUDE_DIR}")

find_library(
  AMD_LIBM_LIBRARY NAME amdlibm
  HINTS ${AMD_LIBM_ROOT} ${AOCL_ROOT}
  PATH_SUFFIXES lib)

target_link_libraries(Math::scalar_vector_functions INTERFACE "${AMD_LIBM_LIBRARY}")

#if(HAVE_AMD_LIBM)
    target_include_directories(Math::scalar_vector_functions INTERFACE "${AMD_LIBM_INCLUDE_DIR}")
    target_compile_definitions(Math::scalar_vector_functions INTERFACE "HAVE_AMD_LIBM")
    set(SINCOS_INCLUDE amdlibm.h)
#endif()
