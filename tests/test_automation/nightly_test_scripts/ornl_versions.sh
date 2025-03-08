#
# Versions should be consistent with setup script
#

# GCC
# Dates at https://gcc.gnu.org/releases.html
gcc_vnew=14.2.0 # Released 2024-08-01
gcc_vold=12.4.0 # Released 2024-06-20

gcc_vcuda=${gcc_vold}
gcc_vintel=${gcc_vold}
gcc_vnvhpc=${gcc_vold}
gcc_vllvmoffload=${gcc_vold}

# LLVM 
# Dates at https://releases.llvm.org/
llvm_vnew=19.1.7 # Released 2025-01-14
llvm_voffload=${llvm_vnew}
cuda_voffload=12.8.0 # CUDA version for offload builds

# HDF5
# Dates at https://portal.hdfgroup.org/display/support/Downloads
hdf5_vnew=1.14.5 # Released 2024-09-30
hdf5_vold=${hdf5_vnew}

# CMake 
# Dates at https://cmake.org/files/
cmake_vnew=3.31.5
cmake_vold=${cmake_vnew}

# OpenMPI
# Dates at https://www.open-mpi.org/software/ompi/v5.0/
ompi_vnew=5.0.6 # Released 2024-11-15
ompi_vold=${ompi_vold}

# Libxml2
libxml2_v=2.13.5 # Released 2024-11-12 See https://gitlab.gnome.org/GNOME/libxml2/-/releases

# FFTW
# Dates at http://www.fftw.org/release-notes.html
fftw_vnew=3.3.10 # Released 2021-09-15
fftw_vold=${fftw_vnew} # Released 2018-05-28

# BOOST
# Dates at https://www.boost.org/users/history/
boost_vnew=1.87.0 # Released 2024-12-12
boost_vold=1.79.0 # Released 2022-04-13

# Python
# Use a single version to reduce dependencies. Ideally the spack prefered version.
python_version=3.13.1

numpy_vnew=2.2.2
numpy_vold=2.2.2
#numpy_vold=1.26.4

