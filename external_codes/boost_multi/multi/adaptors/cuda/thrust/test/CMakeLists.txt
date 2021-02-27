cmake_minimum_required(VERSION 3.11)
project(boost-multi-adaptor-cuda-thrust-test VERSION 0.1 LANGUAGES CXX CUDA)

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr -Xcudafe \"--display_error_number --diag_suppress=implicit_return_from_non_void_function\"")

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

enable_testing()

find_program(MEMORYCHECK_COMMAND valgrind)
set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --error-exitcode=1")
include (CTest)

#file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
set(TEST_SRCS
  array.cu
  vector.cu
)

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

