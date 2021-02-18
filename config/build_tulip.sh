module load PrgEnv-cray/1.0.6
module swap gcc gcc/8.1.0
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
      -DCUDA_TOOLKIT_ROOT_DIR=$CUDA_HOME -DBLA_VENDOR=OpenBLAS \
      -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment \
      .. && make -j32
cd ..
module unload craype-accel-nvidia70

module load rocm
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
      -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment -DBLA_VENDOR=OpenBLAS \
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
      -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment -DBLA_VENDOR=OpenBLAS \
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
      -DQMC_MPI=0 -DQMC_MIXED_PRECISION=1 \
      -DENABLE_OFFLOAD=ON -DOFFLOAD_TARGET=amdgcn-amd-amdhsa -DOFFLOAD_ARCH=gfx906 \
      .. && make -j32
cd ..
