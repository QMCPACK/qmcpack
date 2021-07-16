# Note: this script needs an OpenMP enabled OpenBLAS which is not avaible from
# modules provided on the machine. For this reason, it is for archiving purpose.
# last update July 3rd, 3021

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

echo
echo "###################################"
echo Building V100_Clang_offload_real_MP
echo "###################################"
module load llvm/main-20210112
folder=build_V100_Clang_offload_real_MP
mkdir $folder
cd $folder
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
      -DQMC_MPI=0 -DQMC_MIXED_PRECISION=1 -DENABLE_OFFLOAD=ON \
      -DENABLE_CUDA=OFF -DUSE_OBEJCT_TARGET=ON -DBLA_VENDOR=OpenBLAS \
      -DQMC_DATA=/home/users/coe0097/opt/h5data \
      .. && make -j32
cd ..
module unload llvm/main-20210112

for build in V100_Cray_offload_real_MP \
             MI100_Cray_offload_real_MP MI100_Cray_offload_cplx_MP \
             MI60_Cray_offload_real_MP
do
echo
echo "###################################"
echo "Building $build"
echo "###################################"

CTEST_FLAGS="-DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DENABLE_OFFLOAD=ON"

if [[ $build == *"V100"* ]]; then
  module_hw=craype-accel-nvidia70
elif [[ $build == *"MI100"* ]]; then
  module_hw=craype-accel-amd-gfx908
  module load rocm
elif [[ $build == *"MI60"* ]]; then
  module_hw=craype-accel-amd-gfx906
  module load rocm
fi

if [[ $sys == *"cuda"* ]]; then
  CTEST_FLAGS="-DENABLE_CUDA=ON -DCUDA_ARCH=sm_70"
fi

if [[ $sys == *"cplx"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -DQMC_COMPLEX=ON"
fi

if [[ $sys == *"_MP"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -DQMC_MIXED_PRECISION=ON"
fi

module load $module_hw
module list
folder=build_$build
mkdir $folder; cd $folder
cmake $CTEST_FLAGS -DBLA_VENDOR=OpenBLAS \
      -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment \
      -DQMC_DATA=/home/users/coe0097/opt/h5data \
      .. && make -j32
cd ..
module unload $module_hw
done
