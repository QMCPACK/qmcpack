cmake_minimum_required(VERSION 3.5)

find_package(MPI)

if(DISABLE_MPI)
  message(WARNING "MPI disabled, MPI adaptor will not be tested.")
else()
  if(MPI_FOUND)
    add_subdirectory(./mpi)
  else()
    message(WARNING "MPI not found, MPI adaptor will not be tested.")
  endif()
endif()
