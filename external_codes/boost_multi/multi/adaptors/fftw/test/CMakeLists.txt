# -*-indent-tabs-mode:nil;c-basic-offset:2;tab-width:4;autowrap:nil;-*-
#[=[Multi Test suite can be run like this:
  mkdir -p build
  cd build
  cmake .. [-DENABLE_CUDA=1]
  make -j
  ctest -j --output-on-error [-T memcheck]
  exit
#]=]
cmake_minimum_required(VERSION 3.11)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${qmcpack_SOURCE_DIR}/cmake")

set(CMAKE_VERBOSE_MAKEFILE ON)

project(boost-multi-adaptors-fftw-test VERSION 0.1 LANGUAGES CXX)

find_package(Boost REQUIRED COMPONENTS unit_test_framework timer)

find_package(FFTW REQUIRED COMPONENTS DOUBLE_LIB)
include_directories(${FFTW_INCLUDE_DIRS})
link_libraries(${FFTW_LIBRARIES})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(ENABLE_CUDA OR DEFINED CXXCUDA)
	enable_language(CUDA)
	set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr -Xcudafe \"--diag_suppress=implicit_return_from_non_void_function\"")
endif()

find_package(CUDA QUIET)

if (CUDA_FOUND)
    message("CUDA found")
    include_directories(${CUDA_INCLUDE_DIRS})
else()
    message("CUDA not found")
endif()

enable_testing()
list(APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure") # needs cmake 3.17
include(CTest)

#configure_file("config.hpp.in" ${CMAKE_BINARY_DIR}/config.hpp)

include_directories(${CMAKE_BINARY_DIR})

#file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
set(TEST_SRCS
  combinations.cpp
#  copy.cpp
  core.cpp
#  thrust.cpp
  transpose.cpp
  transpose_square.cpp
)

foreach(TEST_FILE ${TEST_SRCS})
  SET(TEST_EXE "${TEST_FILE}.x")
  add_executable            (${TEST_EXE} ${TEST_FILE})
  if(ENABLE_CUDA OR DEFINED CXXCUDA)
    set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
    target_compile_options  (${TEST_EXE} PRIVATE -std=c++17)
  endif()
#  target_compile_features   (${TEST_EXE} PUBLIC cxx_std_17) 
  target_compile_definitions(${TEST_EXE} PRIVATE "BOOST_PP_VARIADICS")
  target_compile_definitions(${TEST_EXE} PRIVATE ${Boost_DEFINITIONS})
  target_include_directories(${TEST_EXE} PRIVATE ${Boost_INCLUDE_DIRS})
  target_link_libraries     (${TEST_EXE} PRIVATE ${Boost_LIBRARIES})
  target_link_directories   (${TEST_EXE} PRIVATE ${Boost_LIBRARY_DIRS})
  if(NOT ENABLE_CUDA)
    target_compile_options  (${TEST_EXE} PRIVATE
       -Werror -Wall -Wextra -fno-common
       $<$<CXX_COMPILER_ID:GNU>:
          -Wpedantic -Wformat-truncation -fstack-usage>#-Wconversion
       $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
          -Wpedantic -Wmove>
       $<$<CXX_COMPILER_ID:Intel>:
                                        -wd161 -diag-disable=remark -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wpointer-sign -Wreorder -Wno-return-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings -Werror -diag-error:3846
       >
       $<$<CXX_COMPILER_ID:MSVC>:
          /W4>)
  endif()
  add_test(NAME ${TEST_EXE} COMMAND ./${TEST_EXE})
endforeach()

