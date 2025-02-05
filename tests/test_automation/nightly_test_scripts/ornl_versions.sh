#
# Versions should be consistent with setup script
#

# GCC
# Dates at https://gcc.gnu.org/releases.html
gcc_vnew=14.2.0 # Released 2024-08-01
gcc_vold=12.4.0 # Released 2024-06-20

#gcc_vcuda=11.4.0 # https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#host-compiler-support-policy
gcc_vcuda=${gcc_vold}
#gcc_vintel=13.3.0 # Compiler for C++ library used by Intel compiler
#gcc_vnvhpc=13.3.0 # Use makelocalrc to configure NVHPC with this compiler 
#gcc_vllvmoffload=11.4.0 # Version for LLVM offload builds, should be compatible with CUDA version used
gcc_vintel=${gcc_vold}
gcc_vnvhpc=${gcc_vold}
gcc_vllvmoffload=${gcc_vold}

# LLVM 
# Dates at https://releases.llvm.org/
llvm_vnew=19.1.4 # Released 2024-11-19
llvm_voffload=${llvm_vnew}
cuda_voffload=12.4.0 # CUDA version for offload builds

# HDF5
# Dates at https://portal.hdfgroup.org/display/support/Downloads
hdf5_vnew=1.14.5 # Released 2024-09-30
hdf5_vold=${hdf5_vnew}

# CMake 
# Dates at https://cmake.org/files/
cmake_vnew=3.30.5
cmake_vold=${cmake_vnew}

# OpenMPI
# Dates at https://www.open-mpi.org/software/ompi/v5.0/
ompi_vnew=5.0.5 # Released 2024-07-22
ompi_vold=${ompi_vold}

# Libxml2
libxml2_v=2.13.4 # Released 2024-10 See https://gitlab.gnome.org/GNOME/libxml2/-/releases

# FFTW
# Dates at http://www.fftw.org/release-notes.html
fftw_vnew=3.3.10 # Released 2021-09-15
fftw_vold=${fftw_vnew} # Released 2018-05-28

# BOOST
# Dates at https://www.boost.org/users/history/
boost_vnew=1.86.0 # Released 2024-08-14
boost_vold=1.79.0 # Released 2022-04-13

# Python
# Use a single version to reduce dependencies. Ideally the spack prefered version.
python_version=3.13.0
#numpy_vnew=2.1.2
numpy_vnew=1.26.4
numpy_vold=1.26.4

