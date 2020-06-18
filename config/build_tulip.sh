module load PrgEnv-cray
module load gcc
module load hdf5/1.10.1
module load openblas

export FFTW_HOME=/cm/shared/apps/fftw/openmpi/gcc/64/3.3.8
export BOOST_ROOT=/cm/shared/opt/boost/1.72.0

echo
echo "###################################"
echo Building V100_Cray_offload_real_MP
echo "###################################"
module load craype-accel-nvidia70
folder=build_V100_Cray_offload_real_MP
mkdir $folder
cd $folder
cmake -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC \
      -DQMC_MIXED_PRECISION=1 -DENABLE_OFFLOAD=ON \
      -DENABLE_CUDA=ON -DCUDA_ARCH=sm_70 -DCUDA_HOST_COMPILER=`which gcc` \
       -DCUDA_PROPAGATE_HOST_FLAGS=OFF -DCUDA_TOOLKIT_ROOT_DIR=$CUDA_HOME \
      -DENABLE_TIMERS=1 -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment \
      .. && make -j32
cd ..
module unload craype-accel-nvidia70

module load rocm-alt
module load craype-accel-amd-gfx906
echo
echo "###################################"
echo Building MI60_Cray_offload_real_MP
echo "###################################"
folder=build_MI60_Cray_offload_real_MP
mkdir $folder
cd $folder
cmake -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC \
      -DQMC_MIXED_PRECISION=1 -DENABLE_OFFLOAD=ON \
      -DENABLE_TIMERS=1 -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment \
      .. && make -j32
cd ..

echo
echo "###################################"
echo Building MI60_Cray_offload_cplx_MP
echo "###################################"
folder=build_MI60_Cray_offload_cplx_MP
mkdir $folder
cd $folder
cmake -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC \
      -DQMC_MIXED_PRECISION=1 -DENABLE_OFFLOAD=ON -DQMC_COMPLEX=1 \
      -DENABLE_TIMERS=1 -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment \
      .. && make -j32
cd ..
module unload craype-accel-amd-gfx906


echo
echo "###################################"
echo Building MI60_AOMP_offload_real_MP
echo "###################################"
echo Using clang++ from `which clang++`
folder=build_MI60_AOMP_offload_real_MP
mkdir $folder
cd $folder
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
      -DQMC_MPI=0 \
      -DCMAKE_C_FLAGS="-march=native" \
      -DCMAKE_CXX_FLAGS="-march=native -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx906" \
      -DQMC_MIXED_PRECISION=1 -DENABLE_OFFLOAD=ON -DOFFLOAD_TARGET="amdgcn-amd-amdhsa" \
      -DENABLE_TIMERS=1 \
      .. && make -j32
cd ..
