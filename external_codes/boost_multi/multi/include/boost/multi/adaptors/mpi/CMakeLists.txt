cmake_minimum_required(VERSION 3.10)
# project(boost-multi-mpi VERSION 0.1)

find_package(MPI REQUIRED)

add_subdirectory(../mpi/test)
