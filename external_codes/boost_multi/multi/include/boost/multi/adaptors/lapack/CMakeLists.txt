cmake_minimum_required(VERSION 3.18)  # 3.18 for LAPACK::LAPACK

project(
	boost-multi-adaptor-lapack
	VERSION 0.1
	LANGUAGES CXX
)

add_library(multi-lapack INTERFACE)

find_package(LAPACK REQUIRED)

target_link_libraries(multi-lapack INTERFACE multi LAPACK::LAPACK)

add_subdirectory(./test)
