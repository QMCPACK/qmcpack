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

if(EXISTS "CMakeLists.txt")
  message(FATAL_ERROR "You probably don't want to run CMake in a directory with CMakeLists.txt")
endif()

project(boost-multi-test VERSION 0.1 LANGUAGES CXX)

find_package(Boost REQUIRED)# COMPONENTS unit_test_framework)
link_libraries("-lboost_unit_test_framework")

if(NOT DEFINED CMAKE_CXX_STANDARD)
	set(CMAKE_CXX_STANDARD 14)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(ENABLE_CUDA OR DEFINED CXXCUDA)
  enable_language(CUDA)
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr")

  include(FindCUDA/select_compute_arch)
  CUDA_DETECT_INSTALLED_GPUS(INSTALLED_GPU_CCS_1)
  string(STRIP "${INSTALLED_GPU_CCS_1}" INSTALLED_GPU_CCS_2)
  string(REPLACE " " ";" INSTALLED_GPU_CCS_3 "${INSTALLED_GPU_CCS_2}")
  string(REPLACE "." "" CUDA_ARCH_LIST "${INSTALLED_GPU_CCS_3}")
  SET(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_LIST})

endif()

enable_testing()
include(CTest)

file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach(TEST_FILE ${TEST_SRCS})
  SET(TEST_EXE "${TEST_FILE}.x")
  add_executable            (${TEST_EXE} ${TEST_FILE})
  if(ENABLE_CUDA OR DEFINED CXXCUDA)
    set_property(TARGET ${TEST_EXE} PROPERTY "${CUDA_ARCH_LIST}")
    set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
  endif()
  target_compile_definitions(${TEST_EXE} PRIVATE "BOOST_PP_VARIADICS")
  target_compile_definitions(${TEST_EXE} PRIVATE ${Boost_DEFINITIONS})
  target_include_directories(${TEST_EXE} PRIVATE ${Boost_INCLUDE_DIRS})
  target_link_libraries     (${TEST_EXE} PRIVATE ${Boost_LIBRARIES})
  target_link_directories   (${TEST_EXE} PRIVATE ${Boost_LIBRARY_DIRS})
  if(NOT ENABLE_CUDA AND (NOT DART_COMPILER_NAME STREQUAL "nvcc"))
    target_compile_options  (${TEST_EXE} PRIVATE
       -Werror -Wall -Wshadow -Wsign-conversion -Wold-style-cast -Wcast-align -Wunused -Wnon-virtual-dtor -Woverloaded-virtual -Wcast-qual -Wconversion -Wnull-dereference -Wformat=2
       $<$<CXX_COMPILER_ID:GNU>:
                             -Wextra -Wpedantic -Wuninitialized -Wformat-truncation -Wdouble-promotion -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op>
       $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
                             -Wextra -Wpedantic -Wmove>
       $<$<CXX_COMPILER_ID:Intel>:
                             -Wextra                            -wd161 -diag-disable=remark -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wpointer-sign -Wreorder -Wno-return-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings -diag-error:3846>
       $<$<CXX_COMPILER_ID:PGI>: #PGI 20  doesn't recognize the -Wextra option
                                                         --diag_suppress=implicit_return_from_non_void_function>
       $<$<CXX_COMPILER_ID:MSVC>:
          /WX /W4>)
  endif()
  add_test(NAME ${TEST_EXE} COMMAND ./${TEST_EXE})
endforeach()

