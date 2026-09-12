# -*-indent-tabs-mode:nil;c-basic-offset:2;tab-width:4;autowrap:nil;-*-
cmake_minimum_required(VERSION 3.13)  # for reference Ubuntu 20.04 uses 3.16, 3.18 for BLAS::BLAS

include_directories(${CMAKE_BINARY_DIR})

add_subdirectory(test)
