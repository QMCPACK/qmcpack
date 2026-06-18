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

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

set(CMAKE_VERBOSE_MAKEFILE ON)

project(
	boost-multi-adaptors-blis-test
	VERSION 0.1
	LANGUAGES CXX
)

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()

include(CTest)

include_directories(${CMAKE_BINARY_DIR})

include_directories(../../../../../include)

set(TEST_SRCS matrix.cpp)

add_executable(totalview.x totalview.cpp)

target_compile_features(totalview.x PUBLIC cxx_std_17)

# target_compile_definitions(totalview.x PRIVATE "BOOST_PP_VARIADICS")
target_compile_definitions(totalview.x PRIVATE ${Boost_DEFINITIONS})
target_include_directories(totalview.x SYSTEM PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(totalview.x PRIVATE ${Boost_LIBRARIES})
target_link_directories(totalview.x PRIVATE ${Boost_LIBRARY_DIRS})

add_test(NAME totalview.x COMMAND $<TARGET_FILE:totalview.x>)
