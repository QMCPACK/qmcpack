try_compile(PASS_MAPCUDA ${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR}/CMake/try_compile_sources/test_map.cu
            CMAKE_FLAGS "${CMAKE_CUDA_FLAGS}" OUTPUT_VARIABLE COMPILE_OUTPUT)

if(NOT PASS_MAPCUDA)
  message(FATAL_ERROR "Failed to compile host codes containing std::map by nvcc. Please check nvcc and host compiler compatibility. "
                      "You may need to explicitly set an nvcc compatible host compiler as CUDA_HOST_COMPILER. "
                      "Compiler output: ${COMPILE_OUTPUT}")
endif()
