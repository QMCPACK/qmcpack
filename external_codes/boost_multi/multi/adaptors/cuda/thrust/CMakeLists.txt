cmake_minimum_required(VERSION 3.11)
project(boost-multi-adaptor-cuda-thrust-test VERSION 0.1 LANGUAGES CXX CUDA)

message(STATUS "before: ${CUDA_ARCHITECTURES}")

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  # temporary workaround to detect CUDA arch, for between the deprecation of FindCUDA in CMake 3.10 and the availability of FindCudaToolkit in CMake 3.17 # https://stackoverflow.com/a/68223399/225186
  include(FindCUDA/select_compute_arch)
  CUDA_DETECT_INSTALLED_GPUS(INSTALLED_GPU_CCS_1)
  string(STRIP "${INSTALLED_GPU_CCS_1}" INSTALLED_GPU_CCS_2)
  string(REPLACE " " ";" INSTALLED_GPU_CCS_3 "${INSTALLED_GPU_CCS_2}")
  string(REPLACE "." "" CUDA_ARCH_LIST "${INSTALLED_GPU_CCS_3}")
  SET(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_LIST})
  set_property(GLOBAL CUDA_ARCHITECTURES "${CUDA_ARCH_LIST}")
endif()

message(STATUS "after: ${CUDA_ARCHITECTURES}")

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr -Xcudafe \"--display_error_number --diag_suppress=implicit_return_from_non_void_function --diag_suppress=class_and_member_name_conflict\"")

enable_testing()

find_program(MEMORYCHECK_COMMAND valgrind)
set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --error-exitcode=1")
include (CTest)

add_subdirectory(cuda/test)

add_subdirectory(test)

