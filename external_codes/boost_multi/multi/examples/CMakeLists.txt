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
	boost-multi-test
	VERSION 0.1
	LANGUAGES CXX
)

find_package(Boost REQUIRED COMPONENTS serialization unit_test_framework)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()
include(CTest)

add_executable(serialization.cpp.x serialization.cpp)
target_link_libraries(serialization.cpp.x PRIVATE ${Boost_LIBRARIES} stdc++fs)
add_test(NAME serialization.cpp.x COMMAND serialization.cpp.x)

