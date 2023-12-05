#
# Versions should be consistent with setup script
#

# GCC
# Dates at https://gcc.gnu.org/releases.html
gcc_vnew=13.2.0 # Released 2023-07-27
gcc_vold=11.4.0 # Released 2023-05-29

gcc_vcuda=11.4.0 # https://docs.nvidia.com/hpc-sdk/hpc-sdk-release-notes/index.html
gcc_vintel=11.4.0 # Compiler for C++ library used by Intel compiler
gcc_vnvhpc=11.4.0 # Use makelocalrc to configure NVHPC with this compiler 
gcc_vllvmoffload=9.5.0 # Version for LLVM offload builds, should be compatible with CUDA version used

# LLVM 
# Dates at https://releases.llvm.org/
llvm_vnew=16.0.6 # Released 2023-06-19
llvm_voffload=16.0.6
cuda_voffload=11.2.0 # CUDA version for offload builds

# HDF5
# Dates at https://portal.hdfgroup.org/display/support/Downloads
hdf5_vnew=1.14.2 # Released 2023-08-11
hdf5_vold=${hdf5_vnew}

# CMake 
# Dates at https://cmake.org/files/
cmake_vnew=3.27.4 # Relased 2023-08-23
cmake_vold=3.21.3 # Release 2021-09-20

# OpenMPI
# Dates at https://www.open-mpi.org/software/ompi/v4.1/
ompi_vnew=4.1.5 # Released 2023-02-23

# Libxml2
libxml2_v=2.10.3 # Released 2022-12? See https://gitlab.gnome.org/GNOME/libxml2/-/releases

# FFTW
# Dates at http://www.fftw.org/release-notes.html
fftw_vnew=3.3.10 # Released 2021-09-15
fftw_vold=3.3.8 # Released 2018-05-28

# BOOST
# Dates at https://www.boost.org/users/history/
boost_vnew=1.83.0 # Released 2023-08-11
boost_vold=1.77.0 # Released 2021-08-11

# Python
# Use a single version to reduce dependencies. Ideally the spack prefered version.
python_version=3.10.12




