#
# Versions should be consistent with setup script
#

# GCC
# Dates at https://gcc.gnu.org/releases.html
#gcc_vnew=15.1.0 # Released 2025-04-25  # Too ambituous 2025-05-02
gcc_vnew=14.2.0 # Released 2024-08-01
gcc_vold=12.4.0 # Released 2024-06-20

gcc_vcuda=${gcc_vold}
gcc_vintel=${gcc_vold}
gcc_vnvhpc=${gcc_vold}
gcc_vllvmoffload=${gcc_vold}

# LLVM 
# Dates at https://releases.llvm.org/
llvm_vnew=20.1.4 # Released 2025-05-02
llvm_voffload=${llvm_vnew}
cuda_voffload=12.6.0 # https://releases.llvm.org/20.1.0/tools/clang/docs/ReleaseNotes.html#cuda-support #12.6 is limit for LLVM 20.1.0

# HDF5
# Dates at https://portal.hdfgroup.org/display/support/Downloads
#hdf5_vnew=1.14.6 # Released 2025-02-05
hdf5_vnew=1.14.5
hdf5_vold=${hdf5_vnew}

# CMake 
# Dates at https://cmake.org/files/
#cmake_vnew=3.30.8 # Try older version for py-pyscf build
cmake_vnew=3.31.6
cmake_vold=${cmake_vnew}

# OpenMPI
# Dates at https://www.open-mpi.org/software/ompi/v5.0/
ompi_vnew=5.0.6 # Released 2024-11-15
#ompi_vnew=5.0.7 # Released 2025-02-014
ompi_vold=${ompi_vold}

# Libxml2
libxml2_v=2.13.5 # Released 2024-11-12 See https://gitlab.gnome.org/GNOME/libxml2/-/releases

# FFTW
# Dates at http://www.fftw.org/release-notes.html
fftw_vnew=3.3.10 # Released 2021-09-15
fftw_vold=${fftw_vnew} # Released 2018-05-28

# BOOST
# Dates at https://www.boost.org/users/history/
boost_vnew=1.88.0 # Released 2025-04-10
boost_vold=1.82.0 # Released 2023-08-11

# Python
# Use a single version to reduce dependencies. Ideally the spack prefered version.
python_version=3.13.2

numpy_vnew=2.2.5
numpy_vold=2.2.5
#numpy_vold=1.26.4

