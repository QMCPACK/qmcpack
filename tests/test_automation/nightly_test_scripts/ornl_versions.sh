#
# Versions should be consistent with setup script
#

# GCC
# Dates at https://gcc.gnu.org/releases.html
gcc_vnew=15.2.0 # Released 2025-08-08
gcc_vold=13.4.0 # Released 2025-06-06

# Verify vs cuda_voffload release notes e.g.
#  https://docs.nvidia.com/cuda/archive/12.9.0/cuda-installation-guide-linux/index.html#host-compiler-support-policy
#gcc_vllvmoffload=13.4.0
gcc_vllvmoffload=11.5.0 # Must match RHEL system supplied gcc

# LLVM 
# Dates at https://releases.llvm.org/
llvm_vnew=22.1.2 # Released 2026-03-25
llvm_voffload=21.1.8 # 22.1.1 not working as of 2026-04-03
cuda_voffload=12.6.0 # Same as system installed version

# HDF5
# Dates at https://portal.hdfgroup.org/display/support/Downloads
hdf5_vnew=1.14.6 # Released 2025-02-05
#hdf5_vnew=1.14.5
hdf5_vold=${hdf5_vnew}

# CMake 
# Dates at https://cmake.org/files/
#cmake_vnew=3.30.8 # Try older version for py-pyscf build
cmake_vnew=4.2.3
cmake_vold=3.27.9

# OpenMPI
# Dates at https://www.open-mpi.org/software/ompi/v5.0/
ompi_vnew=5.0.10 # Released 2026-02-23
ompi_vold=${ompi_vold}

# FFTW
# Dates at http://www.fftw.org/release-notes.html
fftw_vnew=3.3.10 # Released 2021-09-15
fftw_vold=${fftw_vnew} # Released 2018-05-28

# BOOST
# Dates at https://www.boost.org/users/history/
boost_vnew=1.90.0 # Released 2025-12-10
boost_vold=1.84.0 # Released 2023-12-06

# Python
# Use a single version to reduce dependencies. Ideally the spack prefered version.
python_version=3.14.3

#numpy_vnew=2.4.3
#numpy_vold=2.4.3
numpy_vnew=2.3.5 # PySCF < v2.12.1 is incompatible with >=2.4.0
numpy_vold=2.3.5 

