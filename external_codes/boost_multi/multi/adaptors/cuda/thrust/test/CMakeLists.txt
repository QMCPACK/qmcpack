# -*-indent-tabs-mode:nil;c-basic-offset:2;tab-width:4;autowrap:nil;-*-
cmake_minimum_required(VERSION 3.11)
project(boost-multi-adaptor-cuda-thrust-test VERSION 0.1 LANGUAGES CXX CUDA)

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr -Xcudafe \"--display_error_number --diag_suppress=implicit_return_from_non_void_function --diag_suppress=class_and_member_name_conflict\"")

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

enable_testing()

find_program(MEMORYCHECK_COMMAND valgrind)
set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --error-exitcode=1")
include (CTest)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  # temporary workaround to detect CUDA arch, for between the deprecation of FindCUDA in CMake 3.10 and the availability of FindCudaToolkit in CMake 3.17 # https://stackoverflow.com/a/68223399/225186
  include(FindCUDA/select_compute_arch)
  CUDA_DETECT_INSTALLED_GPUS(INSTALLED_GPU_CCS_1)
  string(STRIP "${INSTALLED_GPU_CCS_1}" INSTALLED_GPU_CCS_2)
  string(REPLACE " " ";" INSTALLED_GPU_CCS_3 "${INSTALLED_GPU_CCS_2}")
  string(REPLACE "." "" CUDA_ARCH_LIST "${INSTALLED_GPU_CCS_3}")
  SET(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_LIST})
  set_property(GLOBAL PROPERTY CUDA_ARCHITECTURES "${CUDA_ARCH_LIST}")
endif()

message(STATUS "after: ${CUDA_ARCHITECTURES}")

#file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

message(STATUS "before cu tests: ${CUDA_ARCH_LIST}")

if("30;35;50;52;60;61;70;70+PTX;75;75+PTX" STREQUAL "${CUDA_ARCH_LIST}" OR "" STREQUAL "${CUDA_ARCH_LIST}")
message(STATUS "gpu not detected")
set(TEST_SRCS
)
else()
message(STATUS "gpu detected")
set(TEST_SRCS
  array.cu
  vector.cu
)
endif()

find_package(Boost REQUIRED COMPONENTS unit_test_framework timer)

foreach(TEST_FILE ${TEST_SRCS})
  SET(TEST_EXE "${TEST_FILE}.x")
  add_executable            (${TEST_EXE} ${TEST_FILE})
  if(ENABLE_CUDA OR DEFINED CXXCUDA)
    set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
    target_compile_options  (${TEST_EXE} PRIVATE -std=c++17 --expt-relaxed-constexpr)
  endif()
#  target_compile_features   (${TEST_EXE} PUBLIC cxx_std_17) 
  target_compile_definitions(${TEST_EXE} PRIVATE "BOOST_PP_VARIADICS")
  target_compile_definitions(${TEST_EXE} PRIVATE ${Boost_DEFINITIONS})
  target_include_directories(${TEST_EXE} PRIVATE ${Boost_INCLUDE_DIRS})
  target_link_libraries     (${TEST_EXE} PRIVATE ${Boost_LIBRARIES})
  target_link_directories   (${TEST_EXE} PRIVATE ${Boost_LIBRARY_DIRS})
#  if(NOT ENABLE_CUDA)
#    target_compile_options  (${TEST_EXE} PRIVATE
#       -Werror -Wall -Wextra -fno-common
#       $<$<CXX_COMPILER_ID:GNU>:
#          -Wpedantic -Wformat-truncation -fstack-usage>#-Wconversion
#       $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
#          -Wpedantic -Wmove>
#       $<$<CXX_COMPILER_ID:Intel>:
#                                        -wd161 -diag-disable=remark -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wpointer-sign -Wreorder -Wno-return-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings -Werror -diag-error:3846
#       >
#       $<$<CXX_COMPILER_ID:MSVC>:
#          /W4>)
#  endif()
  add_test(NAME ${TEST_EXE} COMMAND ./${TEST_EXE})
endforeach()

