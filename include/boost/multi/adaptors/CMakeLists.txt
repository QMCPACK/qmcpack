cmake_minimum_required(VERSION 3.10)
find_package(MPI)

if(DISABLE_MPI)
  message(WARNING "MPI disabled, MPI adaptor will not be tested.")
else()
  if(MPI_FOUND)
    add_subdirectory(./mpi)
  else()
    message(WARNING "MPI not found, MPI adaptor will not be tested. If you want this feature install MPI, for example please run:"
                    "\n    sudo apt install libopenmpi-dev  # or libmpich-dev, in Ubuntu"
                    "\n    sudo dnf install openmpi-devel   # or mpich-devel, in Fedora"
                    "\n    vcpkg install msmpi              # on Windows")
  endif()
endif()

# add_subdirectory(./cufft)

