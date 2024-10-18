# -*-indent-tabs-mode:nil;c-basic-offset:2;tab-width:4;autowrap:nil;-*-
#[=[Multi Test suite can be run like this:
  mkdir -p build
  cd build
  cmake .. [-DENABLE_CUDA=1]
  make -j
  ctest -j --output-on-error [-T memcheck]
  exit
#]=]
cmake_minimum_required(VERSION 3.18)

project(
	boost-multi-examples
	VERSION 0.1
	LANGUAGES CXX
)

include_directories(../include)

#find_package(multi)  # see https://gitlab.com/correaa/boost-multi#using-the-library-installation-and-tests
#...
#target_link_library(my_target PUBLIC multi)

find_package(Boost REQUIRED COMPONENTS serialization unit_test_framework timer)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()
include(CTest)

add_executable(serialization.cpp.x serialization.cpp)
target_link_libraries(serialization.cpp.x PRIVATE ${Boost_LIBRARIES} stdc++fs)
add_test(NAME serialization.cpp.x COMMAND serialization.cpp.x)

add_executable(gj_solve.cpp.x gj_solve.cpp)
target_link_libraries(gj_solve.cpp.x PRIVATE ${Boost_LIBRARIES} stdc++fs)
add_test(NAME gj_solve.cpp.x COMMAND gj_solve.cpp.x)
