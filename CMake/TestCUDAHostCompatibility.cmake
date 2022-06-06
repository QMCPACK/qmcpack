try_compile(PASS_MAPCUDA ${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR}/CMake/try_compile_sources/test_map.cu
            CMAKE_FLAGS "${CMAKE_CUDA_FLAGS}" OUTPUT_VARIABLE COMPILE_OUTPUT)

if(NOT PASS_MAPCUDA)
  message(FATAL_ERROR "CUDA_HOST_COMPILER is not able to compile a simple test program!" "Compiler output: ${COMPILE_OUTPUT}")
endif()
