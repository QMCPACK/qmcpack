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

project(
	boost-multi-test
	VERSION 0.1
	LANGUAGES CXX
)

find_package(Boost REQUIRED COMPONENTS serialization unit_test_framework)

if(NOT DEFINED CMAKE_CXX_STANDARD)
	set(CMAKE_CXX_STANDARD 14)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# if(ENABLE_CUDA OR DEFINED CXXCUDA) enable_language(CUDA) set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr") endif()

enable_testing()
include(CTest)

# file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

add_executable(serialization.cpp.x serialization.cpp)
target_link_libraries(serialization.cpp.x PRIVATE ${Boost_LIBRARIES} stdc++fs)
add_test(NAME serialization.cpp.x COMMAND serialization.cpp.x)

# foreach(TEST_FILE ${TEST_SRCS}) SET(TEST_EXE "${TEST_FILE}x") add_executable            (${TEST_EXE} ${TEST_FILE}) if(ENABLE_CUDA OR DEFINED CXXCUDA) set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA) endif()
# target_compile_definitions(${TEST_EXE} PRIVATE "BOOST_PP_VARIADICS") target_compile_definitions(${TEST_EXE} PRIVATE ${Boost_DEFINITIONS}) target_include_directories(${TEST_EXE} PRIVATE ${Boost_INCLUDE_DIRS}) target_link_libraries
# (${TEST_EXE} PRIVATE ${Boost_LIBRARIES}) target_link_directories   (${TEST_EXE} PRIVATE ${Boost_LIBRARY_DIRS}) if(NOT ENABLE_CUDA) target_compile_options  (${TEST_EXE} PRIVATE -Werror -Wall -Wextra -Wundef -fno-common -Wshadow
# #-Wfatal-errors $<$<CXX_COMPILER_ID:GNU>: -Wpedantic -Wdouble-promotion -Wconversion -Wformat-truncation -fstack-usage>
# $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
# -Wpedantic -Wdouble-promotion -Wmove> $<$<CXX_COMPILER_ID:Intel>: -wd161 -diag-disable=remark -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull
# -Wparentheses -Wpointer-sign -Wreorder -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings -Werror -diag-error:3846> $<$<CXX_COMPILER_ID:MSVC>: /W4>) endif()
# add_test(NAME ${TEST_EXE} COMMAND ./${TEST_EXE}) endforeach()
