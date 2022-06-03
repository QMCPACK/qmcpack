set(CUDA_GNU_COMPATIBLE TRUE)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12.0)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS_EQUAL 11.7)
    set(CUDA_GNU_COMPLATIBLE FALSE)
  endif()
endif()

if(NOT CUDA_GNU_COMPATIBLE)
  message(FATAL_ERROR
          "CUDA version <=11.7 do not support GCC versions >= 12. Please change CMAKE_CXX_COMPILER"
          "to an older GCC version compatible with the CUDA version ${CMAKE_CUDA_COMPILER_VERSION}"
endif()

try_compile(PASS_MAPCUDA ${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR}/CMake/try_compile_sources/test_map.cu
	    CMAKE_FLAGS "${CMAKE_CUDA_FLAGS}" OUTPUT_VARIABLE COMPILE_OUTPUT)

if(NOT PASS_MAPCUDA)
  message(FATAL_ERROR "CUDA functionality failed!" "See compiler output at ${COMPILE_OUTPUT}")
endif()
