# Note: this script needs an OpenMP enabled OpenBLAS which is not avaible from
# modules provided on the machine. For this reason, this script is for archiving purpose.
# last update July 28th, 2021

module load PrgEnv-cray
module swap gcc gcc/8.1.0
module load hdf5/1.10.1

if [ `module list 2>&1 | grep "openblas/dynamic" | wc -l` -gt 0 ]; then
  echo please module unload openblas/dynamic which is not thread-safe for OpenMP
  exit
fi

module use /home/users/coe0097/opt/privatemodules
module load openblas-omp

export FFTW_HOME=/cm/shared/apps/fftw/openmpi/gcc/64/3.3.8
export BOOST_ROOT=/cm/shared/opt/boost/1.72.0

for build in V100_Clang_offload_real_MP V100_Clang_offload_cuda_real_MP \
	     V100_Cray_offload_real_MP V100_Cray_offload_cuda_real_MP \
             MI100_Cray_offload_real_MP MI100_Cray_offload_cplx_MP \
             MI60_Cray_offload_real_MP
do
echo
echo "###################################"
echo "Building $build"
echo "###################################"

if [[ $build == *"V100_Clang"* ]]; then
  module_hw=llvm/main-20210726
  module load cuda11.2
  CTEST_FLAGS="-DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DQMC_MPI=0"
  CTEST_FLAGS="$CTEST_FLAGS -DENABLE_OFFLOAD=ON -DOFFLOAD_ARCH=sm_70 -DUSE_OBJECT_TARGET=ON"
elif [[ $build == *"V100_Cray"* ]]; then
  module_hw=craype-accel-nvidia70
  module load cuda11.2
  CTEST_FLAGS="-DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment -DENABLE_OFFLOAD=ON"
elif [[ $build == *"MI100"* ]]; then
  module_hw=craype-accel-amd-gfx908
  module load rocm
  CTEST_FLAGS="-DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment -DENABLE_OFFLOAD=ON"
elif [[ $build == *"MI60"* ]]; then
  module_hw=craype-accel-amd-gfx906
  module load rocm
  CTEST_FLAGS="-DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment -DENABLE_OFFLOAD=ON"
fi

if [[ $build == *"cuda"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -DENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=70 -DCUDAToolkit_ROOT=$CUDA_ROOT -DCMAKE_CUDA_HOST_COMPILER=`which g++`"
fi

if [[ $build == *"cplx"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -DQMC_COMPLEX=ON"
fi

if [[ $build == *"_MP"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -DQMC_MIXED_PRECISION=ON"
fi

echo "CTEST_FLAGS $CTEST_FLAGS"
module load $module_hw
module list

folder=build_$build
mkdir $folder; cd $folder
cmake $CTEST_FLAGS -DBLA_VENDOR=OpenBLAS \
      -DQMC_DATA=/home/users/coe0097/opt/h5data \
      .. && make -j32
cd ..
module unload $module_hw
done
