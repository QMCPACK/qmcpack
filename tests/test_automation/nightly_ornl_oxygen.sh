#!/bin/bash

echo --- Script START `date`

export LONGCFG="-E long --timeout 1800"
#export LONGCFG="--timeout 7200"

# Directory in which to run tests. Should be an absolute path and fastest usable filesystem
test_path=/scratch/${USER}   # RAID FLASH on oxygen

test_dir=${test_path}/QMCPACK_CI_BUILDS_DO_NOT_REMOVE

export QMC_DATA=/data/pk7/QMC_DATA # Route to directory containing performance test files

if [ -e ${test_path} ]; then

if [ ! -e ${test_dir} ]; then
mkdir ${test_dir}
fi

if [ -e ${test_dir} ]; then
cd ${test_dir}

# Minimize load of GitHub by maintaining a local cloned git used for all builds
if [ ! -e qmcpack ]; then
echo --- Cloning QMCPACK git `date`
git clone https://github.com/QMCPACK/qmcpack.git --depth 1
else
cd qmcpack
echo --- Updating local QMCPACK git `date`
git pull
cd ..
fi

# Sanity check cmake config file present
if [ -e qmcpack/CMakeLists.txt ]; then

echo --- Starting test builds and tests

# Specify GPUs for testing. Enquire via "nvidia-smi -L"
# GPU 0: Tesla K40c (UUID: GPU-224da96d-fb1a-955e-b082-a0f2214877e3)
# GPU 2: Tesla K40c (UUID: GPU-229b6777-2095-3b48-ed2d-50a9f28f5118)
export OXYGEN_KEPLER=GPU-224da96d-fb1a-955e-b082-a0f2214877e3,GPU-229b6777-2095-3b48-ed2d-50a9f28f5118
export OXYGEN_VOLTA=GPU-6bf1c875-b5de-2486-fd0e-ed4bca724ba1
export CUDA_VISIBLE_DEVICES=$OXYGEN_KEPLER

#for sys in build_intel2019 build_intel2019_nompi  build_intel2019_soa build_intel2019_complex build_intel2019_complex_soa  build_intel2019_mixed build_intel2019_mixed_soa build_intel2019_complex_mixed build_intel2019_complex_mixed_soa build_intel2019_nompi_soa build_clang6 build_clang6_cuda build_volta_clang6_cuda build_clang6_cuda_complex build_clang6_cuda_full build_clang6_cuda_complex_full build_clang6 build_clang7 build_gcc8 build_gcc8_mkl build_gcc8_mkl_soa build_gcc8_complex build_gcc8_mkl_complex
for sys in build_intel2019 build_intel2019_soa build_intel2019_nompi build_intel2019_nompi_soa build_gcc8_nompi build_clang7_nompi build_clang6_cuda build_volta_clang6_cuda build_gcc8_mkl build_gcc8_mkl_soa build_gcc8_complex build_gcc8_mkl_complex build_clang7 build_pgi2018_nompi_mkl build_pgi2018_nompi_soa_mkl
do

echo --- START $sys `date`

cd ${test_dir}

if [ -e $sys ]; then
rm -r -f $sys
fi
mkdir $sys
cd $sys

export CUDAVER=10.0
export PATH=/opt/local/bin:/opt/local/sbin:/usr/local/cuda-${CUDAVER}/bin/:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin
export LD_LIBRARY_PATH=/usr/local/cuda-${CUDAVER}/lib64

module() { eval `/usr/bin/modulecmd bash $*`; }

export SPACK_ROOT=$HOME/apps/spack
. $SPACK_ROOT/share/spack/setup-env.sh

echo --- Spack list
spack find
echo --- Modules list
module list
echo --- End listings

#Cmake >=3.12.2 is needed for CUDA 10.0 (cublas_device errors result with earlier versions)
spack load cmake@3.12.2%gcc@4.8.5

case $sys in
    "build_gcc8")
        spack load gcc@8.2.0
        compilerversion=`gcc --version|grep ^gcc|sed 's/^.* //g'`
        if [ $compilerversion = "8.2.0" ]; then
        spack load openmpi%gcc@8.2.0
        export OMPI_CC=gcc
        export OMPI_CXX=g++
	export QMCPACK_TEST_SUBMIT_NAME=GCC8-Release
	export CTCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DBUILD_AFQMC=1 -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short.*diamondC_.x1x1_pp-.mc_.*' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
        spack unload openmpi%gcc@8.2.0
        else
        echo "Did not find expected gcc 8.2.0"
        fi
        spack unload gcc@8.2.0
	;;
    "build_gcc8_nompi")
        spack load gcc@8.2.0
        compilerversion=`gcc --version|grep ^gcc|sed 's/^.* //g'`
        if [ $compilerversion = "8.2.0" ]; then
	export QMCPACK_TEST_SUBMIT_NAME=GCC8-NoMPI-Release
	export CTCFG="-DQMC_MPI=0 -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DBUILD_AFQMC=0 -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
        else
        echo "Did not find expected gcc 8.2.0"
        fi
        spack unload gcc@8.2.0
	;;
    "build_clang6")
        spack load gcc@8.2.0
        spack load llvm@6.0.1
        compilerversion=`clang --version|grep ^clang|sed -e 's/^.* version//g' -e 's/(.*//g'`
        if [ $compilerversion = "6.0.1" ]; then
        spack load openmpi%gcc@8.2.0
        export OMPI_CC=clang
        export OMPI_CXX=clang++
	export QMCPACK_TEST_SUBMIT_NAME=CLANG6-Release
	export CTCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DBUILD_AFQMC=1 -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short.*diamondC_.x1x1_pp-.mc_.*' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
        spack unload openmpi%gcc@8.2.0
        spack unload llvm@6.0.1
        else
        echo "Did not find expected clang 6.0.1"
        fi
	spack unload gcc@8.2.0
	;;
    "build_clang6_cuda")
        spack load gcc@8.2.0
        spack load llvm@6.0.1
        compilerversion=`clang --version|grep ^clang|sed -e 's/^.* version//g' -e 's/(.*//g'`
        if [ $compilerversion = "6.0.1" ]; then
        spack load openmpi%gcc@8.2.0
        export OMPI_CC=clang
        export OMPI_CXX=clang++
	export QMCPACK_TEST_SUBMIT_NAME=CLANG6-CUDA-Release
	export OLD_CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES
	export CUDA_VISIBLE_DEVICES=$OXYGEN_KEPLER
	export CTCFG="-DQMC_CUDA=1 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DBUILD_AFQMC=1 -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
	export CUDA_VISIBLE_DEVICES=$OLD_CUDA_VISIBLE_DEVICES	
        spack unload openmpi%gcc@8.2.0
        spack unload llvm@6.0.1
        else
        echo "Did not find expected clang 6.0.1"
        fi
	;;
    "build_volta_clang6_cuda")
        spack load gcc@8.2.0
        spack load llvm@6.0.1
        compilerversion=`clang --version|grep ^clang|sed -e 's/^.* version//g' -e 's/(.*//g'`
        if [ $compilerversion = "6.0.1" ]; then
        spack load openmpi%gcc@8.2.0
        export OMPI_CC=clang
        export OMPI_CXX=clang++
	export QMCPACK_TEST_SUBMIT_NAME=Volta-CLANG6-CUDA-Release
	export OLD_CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES
	export CUDA_VISIBLE_DEVICES=$OXYGEN_VOLTA
	export CTCFG="-DQMC_CUDA=1 -DCUDA_ARCH="sm_70" -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DBUILD_AFQMC=1 -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
	export CUDA_VISIBLE_DEVICES=$OLD_CUDA_VISIBLE_DEVICES	
        spack unload openmpi%gcc@8.2.0
        spack unload llvm@6.0.1
        else
        echo "Did not find expected clang 6.0.1"
        fi
	spack unload gcc@8.2.0
	;;
    "build_clang6_cuda_soa")
        spack load gcc@8.2.0
        spack load llvm@6.0.1
        compilerversion=`clang --version|grep ^clang|sed -e 's/^.* version//g' -e 's/(.*//g'`
        if [ $compilerversion = "6.0.1" ]; then
        spack load openmpi%gcc@8.2.0
        export OMPI_CC=clang
        export OMPI_CXX=clang++
	export QMCPACK_TEST_SUBMIT_NAME=CLANG6-CUDA-SoA-Release
	export OLD_CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES
	export CUDA_VISIBLE_DEVICES=$OXYGEN_KEPLER
	ctest   -R 'deterministic|performance|short.*diamondC_.x1x1_pp-.mc_.*' -E 'long' --timeout 1800
	export CTCFG="-DQMC_CUDA=1 -DENABLE_SOA=1 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DBUILD_AFQMC=1 -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short.*diamondC_.x1x1_pp-.mc_.*' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
	export CUDA_VISIBLE_DEVICES=$OLD_CUDA_VISIBLE_DEVICES	
        spack unload openmpi%gcc@8.2.0
        spack unload llvm@6.0.1
        else
        echo "Did not find expected clang 6.0.1"
        fi
	spack unload gcc@8.2.0
	;;
    "build_clang6_cuda_complex")
        spack load gcc@8.2.0
        spack load llvm@6.0.1
        compilerversion=`clang --version|grep ^clang|sed -e 's/^.* version//g' -e 's/(.*//g'`
        if [ $compilerversion = "6.0.1" ]; then
        spack load openmpi%gcc@8.2.0
        export OMPI_CC=clang
        export OMPI_CXX=clang++
	export QMCPACK_TEST_SUBMIT_NAME=CLANG6-CUDA-Complex-Release
	export OLD_CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES
	export CUDA_VISIBLE_DEVICES=$OXYGEN_KEPLER
	export CTCFG="-DQMC_CUDA=1 -DQMC_COMPLEX=1 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DBUILD_AFQMC=1 -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short.*diamondC_.x1x1_pp-.mc_.*' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
	export CUDA_VISIBLE_DEVICES=$OLD_CUDA_VISIBLE_DEVICES	
        spack unload openmpi%gcc@8.2.0
        spack unload llvm@6.0.1
        else
        echo "Did not find expected clang 6.0.1"
        fi
	spack unload gcc@8.2.0
	;;
    "build_clang6_cuda_full")
        spack load gcc@8.2.0
        spack load llvm@6.0.1
        compilerversion=`clang --version|grep ^clang|sed -e 's/^.* version//g' -e 's/(.*//g'`
        if [ $compilerversion = "6.0.1" ]; then
        spack load openmpi%gcc@8.2.0
        export OMPI_CC=clang
        export OMPI_CXX=clang++
	export QMCPACK_TEST_SUBMIT_NAME=CLANG6-CUDA-Full-Release
	export OLD_CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES
	export CUDA_VISIBLE_DEVICES=$OXYGEN_KEPLER
	export CTCFG="-DQMC_CUDA=1 -DQMC_MIXED_PRECISION=0 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DBUILD_AFQMC=1 -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short.*diamondC_.x1x1_pp-.mc_.*' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
	export CUDA_VISIBLE_DEVICES=$OLD_CUDA_VISIBLE_DEVICES	
        spack unload openmpi%gcc@8.2.0
        spack unload llvm@6.0.1
        else
        echo "Did not find expected clang 6.0.1"
        fi
	spack unload gcc@8.2.0
	;;
    "build_clang6_cuda_complex_full")
        spack load gcc@8.2.0
        spack load llvm@6.0.1
        compilerversion=`clang --version|grep ^clang|sed -e 's/^.* version//g' -e 's/(.*//g'`
        if [ $compilerversion = "6.0.1" ]; then
        spack load openmpi%gcc@8.2.0
        export OMPI_CC=clang
        export OMPI_CXX=clang++
	export QMCPACK_TEST_SUBMIT_NAME=CLANG6-CUDA-Complex-Full-Release
	export OLD_CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES
	export CUDA_VISIBLE_DEVICES=$OXYGEN_KEPLER
	export CTCFG="-DQMC_CUDA=1 -DQMC_COMPLEX=1 -DQMC_MIXED_PRECISION=0 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DBUILD_AFQMC=1 -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short.*diamondC_.x1x1_pp-.mc_.*' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
	export CUDA_VISIBLE_DEVICES=$OLD_CUDA_VISIBLE_DEVICES	
        spack unload openmpi%gcc@8.2.0
        spack unload llvm@6.0.1
        else
        echo "Did not find expected clang 6.0.1"
        fi
	spack unload gcc@8.2.0
	;;
    "build_clang7")
        spack load llvm@7.0.0
        compilerversion=`clang --version|grep ^clang|sed -e 's/^.* version//g' -e 's/(.*//g'`
        if [ $compilerversion = "7.0.0" ]; then
        spack load openmpi%gcc@4.8.5
        export OMPI_CC=clang
        export OMPI_CXX=clang++
	export QMCPACK_TEST_SUBMIT_NAME=CLANG7-Release
	export CTCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DBUILD_AFQMC=1 -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short.*diamondC_.x1x1_pp-.mc_.*' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
	export CTCFG="-DCMAKE_CXX_COMPILER=mpicxx -DBUILD_AFQMC=1 -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV" 
        spack unload openmpi%gcc@4.8.5
        spack unload llvm@7.0.0
        else
        echo "Did not find expected clang 7.0.0"
        fi
	;;
    "build_clang7_nompi")
        spack load llvm@7.0.0
        compilerversion=`clang --version|grep ^clang|sed -e 's/^.* version//g' -e 's/(.*//g'`
        if [ $compilerversion = "7.0.0" ]; then
	export QMCPACK_TEST_SUBMIT_NAME=CLANG7-NoMPI-Release
	export CTCFG="-DQMC_MPI=0 DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DBUILD_AFQMC=0 -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short.*diamondC_.x1x1_pp-.mc_.*' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
        spack unload llvm@7.0.0
        else
        echo "Did not find expected clang 7.0.0"
        fi
	;;
    "build_gcc8_mkl")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2019/mkl/bin/mklvars.sh intel64
        spack load gcc@8.2.0
        compilerversion=`gcc --version|grep ^gcc|sed 's/^.* //g'`
        if [ $compilerversion = "8.2.0" ]; then
        spack load openmpi%gcc@8.2.0
        export OMPI_CC=gcc
        export OMPI_CXX=g++
	export QMCPACK_TEST_SUBMIT_NAME=GCC8-MKL-Release
	export CTCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DBLA_VENDOR=Intel10_64lp_seq -DCMAKE_PREFIX_PATH=$MKLROOT/lib -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
        spack unload openmpi%gcc@8.2.0
        else
        echo "Did not find expected gcc 8.2.0"
        fi
        spack unload gcc@8.2.0
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_gcc8_mkl_soa")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2019/mkl/bin/mklvars.sh intel64
        spack load gcc@8.2.0
        compilerversion=`gcc --version|grep ^gcc|sed 's/^.* //g'`
        if [ $compilerversion = "8.2.0" ]; then
        spack load openmpi%gcc@8.2.0
        export OMPI_CC=gcc
        export OMPI_CXX=g++
	export QMCPACK_TEST_SUBMIT_NAME=GCC8-MKL-SoA-Release
	export CTCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DENABLE_SOA=1 -DBLA_VENDOR=Intel10_64lp_seq -DCMAKE_PREFIX_PATH=$MKLROOT/lib -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
        spack unload openmpi%gcc@8.2.0
        else
        echo "Did not find expected gcc 8.2.0"
        fi
        spack unload gcc@8.2.0
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2019")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2019/bin/compilervars.sh intel64
        spack load gcc@8.2.0

	#For Intel2019 we also setup QE
	export QE_VERSION=6.3
        # QE version 6.x unpacks to qe-; Older versions 5.x uses espresso-
        export QE_PREFIX=qe-
	export QE_BIN=${test_dir}/${sys}_QE/${QE_PREFIX}${QE_VERSION}/bin
	echo --- QE_BIN set to ${QE_BIN}
	if [ ! -e ${QE_BIN}/pw.x ]; then
	    # Start from clean build if no executable present
	    if [ -e ${test_dir}/${sys}_QE ]; then
		rm -r -f ${test_dir}/${sys}_QE
	    fi
	    mkdir ${test_dir}/${sys}_QE
		
	    cd ${test_dir}/${sys}_QE
	    cp -p ../qmcpack/external_codes/quantum_espresso/*${QE_VERSION}* .
	    ./download_and_patch_qe${QE_VERSION}.sh
	    cd ${QE_PREFIX}${QE_VERSION}
            ./configure CC=mpiicc MPIF90=mpiifort F77=mpiifort --with-scalapack=intel --with-hdf5=/home/pk7/apps/hdf5-1.10.1-intel-mpi
            make pwall # No parallel build due to sometimes broken dependencies in QE build system
	    cd ..
	    cd ${test_dir}/${sys}
	else
	    echo -- Found existing QE ${QE_VERSION} executable
	fi
	export QMCPACK_TEST_SUBMIT_NAME=Intel2019-Release
	export CTCFG="-DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DQE_BIN=${QE_BIN} -DBUILD_AFQMC=1 -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel-mpi/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
        spack unload gcc@8.2.0
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2019_soa")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2019/bin/compilervars.sh intel64
        spack load gcc@8.2.0
	export QMCPACK_TEST_SUBMIT_NAME=Intel2019-SoA-Release
	export CTCFG="-DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DBUILD_AFQMC=1 -DENABLE_SOA=1 -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel-mpi/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
        spack unload gcc@8.2.0
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2019_nompi")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2019/bin/compilervars.sh intel64
        spack load gcc@8.2.0
	export QMCPACK_TEST_SUBMIT_NAME=Intel2019-NoMPI-Release
	export CTCFG="-DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DQMC_MPI=0 -DBUILD_AFQMC=0 -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
        spack unload gcc@8.2.0
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2019_nompi_soa")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2019/bin/compilervars.sh intel64
        spack load gcc@8.2.0
	export QMCPACK_TEST_SUBMIT_NAME=Intel2019-NoMPI-SoA-Release
	export CTCFG="-DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DQMC_MPI=0 -DENABLE_SOA=1 -DBUILD_AFQMC=0 -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel/ -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
        spack unload gcc@8.2.0
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_gcc8_complex")
        spack load gcc@8.2.0
        compilerversion=`gcc --version|grep ^gcc|sed 's/^.* //g'`
        if [ $compilerversion = "8.2.0" ]; then
        spack load openmpi%gcc@8.2.0
        export OMPI_CC=gcc
        export OMPI_CXX=g++
	export QMCPACK_TEST_SUBMIT_NAME=GCC8-Complex-Release
	export CTCFG="-DQMC_COMPLEX=1 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DBUILD_AFQMC=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short.*diamondC_.x1x1_pp-.mc_.*' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
        spack unload openmpi%gcc@8.2.0
        else
        echo "Did not find expected gcc 8.2.0"
        fi
        spack unload gcc@8.2.0
	;;
    "build_gcc8_mkl_complex")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2019/mkl/bin/mklvars.sh intel64
        spack load gcc@8.2.0
        compilerversion=`gcc --version|grep ^gcc|sed 's/^.* //g'`
        if [ $compilerversion = "8.2.0" ]; then
        spack load openmpi%gcc@8.2.0
        export OMPI_CC=gcc
        export OMPI_CXX=g++
	export QMCPACK_TEST_SUBMIT_NAME=GCC8-MKL-Complex-Release
	export CTCFG="-DQMC_COMPLEX=1 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DBLA_VENDOR=Intel10_64lp_seq -DCMAKE_PREFIX_PATH=$MKLROOT/lib -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short.*diamondC_.x1x1_pp-.mc_.*' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
        spack unload openmpi%gcc@8.2.0
        else
        echo "Did not find expected gcc 8.2.0"
        fi
        spack unload gcc@8.2.0
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2019_complex")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2019/bin/compilervars.sh intel64
        spack load gcc@8.2.0
	export QMCPACK_TEST_SUBMIT_NAME=Intel2019-Complex-Release
	export CTCFG="-DQMC_COMPLEX=1 -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DQE_BIN=${QE_BIN} -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel-mpi/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short.*diamondC_.x1x1_pp-.mc_.*' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
        spack unload gcc@8.2.0
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2019_complex_soa")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2019/bin/compilervars.sh intel64
        spack load gcc@8.2.0
	export QMCPACK_TEST_SUBMIT_NAME=Intel2019-Complex-SoA-Release
	export CTCFG="-DQMC_COMPLEX=1 -DENABLE_SOA=1 -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel-mpi/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short.*diamondC_.x1x1_pp-.mc_.*' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
        spack unload gcc@8.2.0
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2019_mixed")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2019/bin/compilervars.sh intel64
        spack load gcc@8.2.0
	export QMCPACK_TEST_SUBMIT_NAME=Intel2019-Mixed-Release
	export CTCFG="-DQMC_MIXED_PRECISION=1 -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel-mpi/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short.*diamondC_.x1x1_pp-.mc_.*' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
        spack unload gcc@8.2.0
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2019_mixed_soa")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2019/bin/compilervars.sh intel64
        spack load gcc@8.2.0
	export QMCPACK_TEST_SUBMIT_NAME=Intel2019-Mixed-SoA-Release
	export CTCFG="-DQMC_MIXED_PRECISION=1 -DENABLE_SOA=1 -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel-mpi/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
        spack unload gcc@8.2.0
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2019_complex_mixed")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2019/bin/compilervars.sh intel64
        spack load gcc@8.2.0
	export QMCPACK_TEST_SUBMIT_NAME=Intel2019-Complex-Mixed-Release
	export CTCFG="-DQMC_MIXED_PRECISION=1 -DQMC_COMPLEX=1 -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel-mpi/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short.*diamondC_.x1x1_pp-.mc_.*' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
        spack unload gcc@8.2.0
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2019_complex_mixed_soa")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2019/bin/compilervars.sh intel64
        spack load gcc@8.2.0
	export QMCPACK_TEST_SUBMIT_NAME=Intel2019-Complex-Mixed-SoA-Release
	export CTCFG="-DQMC_MIXED_PRECISION=1 -DQMC_COMPLEX=1 -DENABLE_SOA=1 -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel-mpi/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short.*diamondC_.x1x1_pp-.mc_.*' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
        spack unload gcc@8.2.0
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_pgi2018_nompi_mkl")
# On oxygen PGI is setup to use the gcc 7.2.0 files
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2019/mkl/bin/mklvars.sh intel64
        spack load gcc@7.2.0
        compilerversion=`gcc --version|grep ^gcc|sed 's/^.* //g'`
        if [ $compilerversion = "7.2.0" ]; then
#        spack load openmpi%gcc@8.2.0
#        export OMPI_CC=gcc
#        export OMPI_CXX=g++
	export QMCPACK_TEST_SUBMIT_NAME=PGI2018-NoMPI-MKL-Release
#	export CTCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DBLA_VENDOR=Intel10_64lp_seq -DCMAKE_PREFIX_PATH=$MKLROOT/lib -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
	export CTCFG="-DCMAKE_C_COMPILER=pgcc -DCMAKE_CXX_COMPILER=pgc++ -DQMC_MPI=0 -DBLA_VENDOR=Intel10_64lp_seq -DCMAKE_PREFIX_PATH=$MKLROOT/lib -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
#        spack unload openmpi%gcc@8.2.0
        else
        echo "Did not find expected gcc 7.2.0"
        fi
        spack unload gcc@7.2.0
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_pgi2018_nompi_soa_mkl")
# On oxygen PGI is setup to use the gcc 7.2.0 files
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2019/mkl/bin/mklvars.sh intel64
        spack load gcc@7.2.0
        compilerversion=`gcc --version|grep ^gcc|sed 's/^.* //g'`
        if [ $compilerversion = "7.2.0" ]; then
#        spack load openmpi%gcc@8.2.0
#        export OMPI_CC=gcc
#        export OMPI_CXX=g++
	export QMCPACK_TEST_SUBMIT_NAME=PGI2018-NoMPI-SoA-MKL-Release
#	export CTCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DBLA_VENDOR=Intel10_64lp_seq -DCMAKE_PREFIX_PATH=$MKLROOT/lib -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
	export CTCFG="-DCMAKE_C_COMPILER=pgcc -DCMAKE_CXX_COMPILER=pgc++ -DENABLE_SOA=1 -DQMC_MPI=0 -DBLA_VENDOR=Intel10_64lp_seq -DCMAKE_PREFIX_PATH=$MKLROOT/lib -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV"
 	ctest ${CTCFG} -R deterministic ${LONGCFG} --track Deterministic
        sed -i "s/Deterministic/Short/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'short' -LE unstable ${LONGCFG} --track Short
        sed -i "s/Short/Performance/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R 'performance' -LE unstable ${LONGCFG} --track Performance
        sed -i "s/Performance/Unstable/g" Testing/TAG  # Hack to change cached track
	ctest ${CTCFG} -R short -L unstable ${LONGCFG} --track Unstable
#        spack unload openmpi%gcc@8.2.0
        else
        echo "Did not find expected gcc 7.2.0"
        fi
        spack unload gcc@7.2.0
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    *)
	echo "ERROR: Unknown build type $sys"
	;;
esac

echo --- END $sys `date`
done

else
echo "ERROR: No CMakeLists. Bad git clone or update"
exit 1
fi

else
echo "ERROR: Unable to make test directory ${test_dir}"
exit 1
fi

else
echo "ERROR: No directory ${test_path}"
exit 1
fi
echo --- Script END `date`
