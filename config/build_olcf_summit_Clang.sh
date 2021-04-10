#!/bin/bash

echo "----------------------- WARNING ------------------------------------"
echo "This is **not** production ready and intended for development only!!"
echo "Use config/build_olcf_summit.sh for production on Summit."
echo "----------------------- WARNING ------------------------------------"

echo "Purging current module set"
module purge
echo "Loading QMCPACK dependency modules for summit"
module load gcc/8.1.1
module load spectrum-mpi
module load cmake
module load git
module load cuda
module load essl
module load netlib-lapack
module load hdf5
module load python/3.6.6-anaconda3-5.3.0
# private module until OLCF provides a new llvm build
if [[ ! -d /ccs/proj/mat151/opt/modules ]] ; then
  echo "Required module folder /ccs/proj/mat151/opt/modules not found!"
  exit 1
fi
module use /ccs/proj/mat151/opt/modules
module load llvm/master-latest

#the XL built fftw is buggy, use the gcc version
#module load fftw
export FFTW_HOME=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-6.4.0/fftw-3.3.8-5gcj2ic4el7acu3rqnfnh735jz2ez7j5
export BOOST_ROOT=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-6.4.0/boost-1.66.0-l3sghp3ggjzwi4vtvyb5yzsjm36npgrk

TYPE=Release
Compiler=Clang

for name in offload_cuda_real offload_cuda_real_MP offload_cuda_cplx offload_cuda_cplx_MP \
            cpu_real cpu_real_MP cpu_cplx cpu_cplx_MP
do

# not sure the about '-b;packed:16' for non-mpi tests.
CMAKE_FLAGS="-D CMAKE_BUILD_TYPE=$TYPE -D ENABLE_MASS=1 -D MASS_ROOT=/sw/summit/xl/16.1.1-5/xlmass/9.1.1 -D MPIEXEC_EXECUTABLE=`which jsrun` -D MPIEXEC_NUMPROC_FLAG='-n' -D MPIEXEC_PREFLAGS='-c;16;-g;1;-b;packed:16' -DMPIEXEC_NO_MPI_PREFLAGS='--smpiargs="-diasble_gpu_hooks";-g;1;-b;packed:16'"
if [[ $name == *"cplx"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D QMC_COMPLEX=1"
fi

if [[ $name == *"_MP"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D QMC_MIXED_PRECISION=1"
fi

if [[ $name == *"offload"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D ENABLE_OFFLOAD=ON -D USE_OBJECT_TARGET=ON -DOFFLOAD_ARCH=sm_70"
fi

if [[ $name == *"cuda"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D ENABLE_CUDA=1 -D CUDA_ARCH=sm_70 -D CUDA_HOST_COMPILER=`which gcc`"
fi

folder=build_summit_${Compiler}_${name}
echo "**********************************"
echo "$folder"
echo "$CMAKE_FLAGS"
echo "**********************************"
mkdir $folder
cd $folder
if [ ! -f CMakeCache.txt ] ; then
cmake $CMAKE_FLAGS -D CMAKE_C_COMPILER=mpicc -D CMAKE_CXX_COMPILER=mpicxx ..
cmake ..
fi
make -j16
cd ..

echo
done
