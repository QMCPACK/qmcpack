cmake_minimum_required(VERSION 3.11)
project(boost-multi-adaptor-cuda-thrust-test VERSION 0.1 LANGUAGES CXX CUDA)

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr -Xcudafe \"--display_error_number --diag_suppress=implicit_return_from_non_void_function --diag_suppress=class_and_member_name_conflict\"")

enable_testing()

find_program(MEMORYCHECK_COMMAND valgrind)
set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --error-exitcode=1")
include (CTest)

add_subdirectory(cuda/test)

add_subdirectory(test)

