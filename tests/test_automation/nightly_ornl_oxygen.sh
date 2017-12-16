#!/bin/bash

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

# Specify GPU on oxygen for testing. Enquire via "nvidia-smi -L"
# GPU 0: Tesla K40c (UUID: GPU-224da96d-fb1a-955e-b082-a0f2214877e3)
export CUDA_VISIBLE_DEVICES=GPU-224da96d-fb1a-955e-b082-a0f2214877e3

#Caution: intel2015 also builds QE and sets QE_BIN directory. Should be run ahead of intel2015_complex, intel2015_cuda, intel2015_cuda_complex
#for sys in build_intel2017_nompi build_intel2017_nompi_soa build_intel2017 build_intel2017_soa build_intel2017_complex build_intel2017_complex_soa build_gcc_mkl build_gcc_cuda build_intel2017_mixed build_intel2017_mixed_soa build_intel2017_complex_mixed build_intel2017_complex_mixed_soa build_gcc_mkl_complex build_intel2015 build_intel2015_complex build_intel2015_cuda build_intel2015_cuda_complex build_gcc build_gcc_complex build_gcc_cuda_complex build_gcc_cuda_full build_gcc_cuda_complex_full build_gcc_mkl_soa build_gcc_cuda_soa 
for sys in build_intel2017_nompi build_intel2017_nompi_soa build_intel2017 build_intel2017_soa build_intel2017_complex build_intel2017_complex_soa build_gcc_mkl build_gcc_cuda build_gcc_cuda_complex build_intel2017_mixed build_intel2017_mixed_soa build_intel2017_complex_mixed build_intel2017_complex_mixed_soa build_gcc_mkl_complex build_intel2015 build_intel2015_complex build_gcc build_gcc_complex build_gcc_mkl_soa build_gcc_cuda_soa
do

echo --- Building for $sys `date`

cd ${test_dir}

if [ -e $sys ]; then
rm -r -f $sys
fi
mkdir $sys
cd $sys


export PATH=/opt/local/bin:/opt/local/sbin:/usr/local/cuda-8.0/bin/:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin
export LD_LIBRARY_PATH=/usr/local/cuda-8.0/lib64

case $sys in
    "build_gcc")
#	module() { eval `/usr/bin/modulecmd sh $*`; }
#	module load mpi
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export PATH=$HOME/apps/openmpi-2.0.2/bin/:$PATH
	export LD_LIBRARY_PATH=$HOME/apps/openmpi-2.0.2/lib/:$LD_LIBRARY_PATH
	export QMCPACK_TEST_SUBMIT_NAME=GCC-Release
	ctest -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DBUILD_AFQMC=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
#	module unload mpi
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	;;
    "build_gcc_mkl")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
#	module() { eval `/usr/bin/modulecmd sh $*`; }
#	module load mpi
	export PATH=$HOME/apps/openmpi-2.0.2/bin/:$PATH
	export LD_LIBRARY_PATH=$HOME/apps/openmpi-2.0.2/lib/:$LD_LIBRARY_PATH
	source /opt/intel2017/mkl/bin/mklvars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=GCC-MKL-Release
	ctest -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DBLA_VENDOR=Intel10_64lp_seq -DCMAKE_PREFIX_PATH=$MKLROOT/lib -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
#	module unload mpi
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_gcc_mkl_soa")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
#	module() { eval `/usr/bin/modulecmd sh $*`; }
#	module load mpi
	export PATH=$HOME/apps/openmpi-2.0.2/bin/:$PATH
	export LD_LIBRARY_PATH=$HOME/apps/openmpi-2.0.2/lib/:$LD_LIBRARY_PATH
	source /opt/intel2017/mkl/bin/mklvars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=GCC-MKL-SoA-Release
	ctest -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DENABLE_SOA=1 -DBLA_VENDOR=Intel10_64lp_seq -DCMAKE_PREFIX_PATH=$MKLROOT/lib -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
#	module unload mpi
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2017")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2017/bin/compilervars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=Intel2017-Release
	ctest -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DBUILD_AFQMC=1 -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel-mpi/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2017_soa")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2017/bin/compilervars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=Intel2017-SoA-Release
	ctest -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DBUILD_AFQMC=1 -DENABLE_SOA=1 -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel-mpi/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2017_nompi")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2017/bin/compilervars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=Intel2017-NoMPI-Release
	ctest -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DQMC_MPI=0 -DBUILD_AFQMC=0 -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2017_nompi_soa")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2017/bin/compilervars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=Intel2017-NoMPI-SoA-Release
	ctest -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DQMC_MPI=0 -DENABLE_SOA=1 -DBUILD_AFQMC=0 -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2015")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel/bin/compilervars.sh intel64
	source /opt/intel/impi_latest/bin64/mpivars.sh
	export QMCPACK_TEST_SUBMIT_NAME=Intel2015-Release
	#For Intel2015 we also setup QE
	export QE_VERSION=5.3.0
	export QE_BIN=${test_dir}/${sys}_QE/espresso-${QE_VERSION}/bin
	echo --- QE_BIN set to ${QE_BIN}
	if [ ! -e ${QE_BIN}/pw.x ]; then
	    # Start from clean build if no final executable present
	    if [ -e ${test_dir}/${sys}_QE ]; then
		rm -r -f ${test_dir}/${sys}_QE
	    fi
	    mkdir ${test_dir}/${sys}_QE
		
	    cd ${test_dir}/${sys}_QE
	    cp -p ../qmcpack/external_codes/quantum_espresso/*${QE_VERSION}* .
	    ./download_and_patch_qe${QE_VERSION}.sh
	    cd espresso-${QE_VERSION}
	    ./configure --with-hdf5 HDF5_DIR=/usr/local/ MPIF90=mpiifort F90=mpiifort
	    # Espresso incorrect assumes we are using OpenMPI when IntelMPI is enabled.
	    # Update make.sys so correct MKL BLACS linked
	    mv make.sys make.sys_orig
	    sed 's/mkl_blacs_openmpi_lp64/mkl_blacs_intelmpi_lp64/' make.sys_orig >make.sys
	    make -j 24 pwall
	    cd ..
	    
	    cd ${test_dir}/${sys}
	else
	    echo -- Found existing QE ${QE_VERSION} executable
	fi
	ctest -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DQE_BIN=${QE_BIN} -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_gcc_complex")
#	module() { eval `/usr/bin/modulecmd sh $*`; }
#	module load mpi
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export PATH=$HOME/apps/openmpi-2.0.2/bin/:$PATH
	export LD_LIBRARY_PATH=$HOME/apps/openmpi-2.0.2/lib/:$LD_LIBRARY_PATH
	export QMCPACK_TEST_SUBMIT_NAME=GCC-Complex-Release
	ctest -DQMC_COMPLEX=1 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
#	module unload mpi
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	;;
    "build_gcc_mkl_complex")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
#	module() { eval `/usr/bin/modulecmd sh $*`; }
#	module load mpi
	export PATH=$HOME/apps/openmpi-2.0.2/bin/:$PATH
	export LD_LIBRARY_PATH=$HOME/apps/openmpi-2.0.2/lib/:$LD_LIBRARY_PATH
	source /opt/intel2017/mkl/bin/mklvars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=GCC-MKL-Complex-Release
	ctest -DQMC_COMPLEX=1 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DBLA_VENDOR=Intel10_64lp_seq -DCMAKE_PREFIX_PATH=$MKLROOT/lib -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
#	module unload mpi
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2017_complex")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2017/bin/compilervars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=Intel2017-Complex-Release
	ctest -DQMC_COMPLEX=1 -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel-mpi/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2017_complex_soa")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2017/bin/compilervars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=Intel2017-Complex-SoA-Release
	ctest -DQMC_COMPLEX=1 -DENABLE_SOA=1 -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel-mpi/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2017_mixed")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2017/bin/compilervars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=Intel2017-Mixed-Release
	ctest -DQMC_MIXED_PRECISION=1 -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel-mpi/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2017_mixed_soa")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2017/bin/compilervars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=Intel2017-Mixed-SoA-Release
	ctest -DQMC_MIXED_PRECISION=1 -DENABLE_SOA=1 -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel-mpi/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2017_complex_mixed")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2017/bin/compilervars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=Intel2017-Complex-Mixed-Release
	ctest -DQMC_MIXED_PRECISION=1 -DQMC_COMPLEX=1 -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel-mpi/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2017_complex_mixed_soa")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel2017/bin/compilervars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=Intel2017-Complex-Mixed-SoA-Release
	ctest -DQMC_MIXED_PRECISION=1 -DQMC_COMPLEX=1 -DENABLE_SOA=1 -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DBUILD_LMYENGINE_INTERFACE=1 -DHDF5_ROOT=/home/pk7/apps/hdf5-1.10.1-intel-mpi/ -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2015_complex")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel/bin/compilervars.sh intel64
	source /opt/intel/impi_latest/bin64/mpivars.sh
	export QMCPACK_TEST_SUBMIT_NAME=Intel2015-Complex-Release
	ctest -DQMC_COMPLEX=1 -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DQE_BIN=${QE_BIN} -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_gcc_cuda")
#	module() { eval `/usr/bin/modulecmd sh $*`; }
#	module load mpi
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export PATH=$HOME/apps/openmpi-2.0.2/bin/:$PATH
	export LD_LIBRARY_PATH=$HOME/apps/openmpi-2.0.2/lib/:$LD_LIBRARY_PATH
	export QMCPACK_TEST_SUBMIT_NAME=GCC-CUDA-Release
	ctest -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_CUDA=1 -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
#	module unload mpi
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	;;
    "build_gcc_cuda_soa")
#	module() { eval `/usr/bin/modulecmd sh $*`; }
#	module load mpi
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export PATH=$HOME/apps/openmpi-2.0.2/bin/:$PATH
	export LD_LIBRARY_PATH=$HOME/apps/openmpi-2.0.2/lib/:$LD_LIBRARY_PATH
	export QMCPACK_TEST_SUBMIT_NAME=GCC-CUDA-SoA-Release
	ctest -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_CUDA=1 -DENABLE_SOA=1 -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
#	module unload mpi
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	;;
    "build_gcc_cuda_complex")
#	module() { eval `/usr/bin/modulecmd sh $*`; }
#	module load mpi
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export PATH=$HOME/apps/openmpi-2.0.2/bin/:$PATH
	export LD_LIBRARY_PATH=$HOME/apps/openmpi-2.0.2/lib/:$LD_LIBRARY_PATH
	export QMCPACK_TEST_SUBMIT_NAME=GCC-CUDA-Complex-Release
	ctest -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_CUDA=1 -DQMC_COMPLEX=1 -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
#	module unload mpi
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	;;
    "build_gcc_cuda_full")
#	module() { eval `/usr/bin/modulecmd sh $*`; }
#	module load mpi
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export PATH=$HOME/apps/openmpi-2.0.2/bin/:$PATH
	export LD_LIBRARY_PATH=$HOME/apps/openmpi-2.0.2/lib/:$LD_LIBRARY_PATH
	export QMCPACK_TEST_SUBMIT_NAME=GCC-CUDA-Full-Release
	ctest -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_CUDA=1 -DQMC_MIXED_PRECISION=0 -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
#	module unload mpi
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	;;
    "build_gcc_cuda_complex_full")
#	module() { eval `/usr/bin/modulecmd sh $*`; }
#	module load mpi
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export PATH=$HOME/apps/openmpi-2.0.2/bin/:$PATH
	export LD_LIBRARY_PATH=$HOME/apps/openmpi-2.0.2/lib/:$LD_LIBRARY_PATH
	export QMCPACK_TEST_SUBMIT_NAME=GCC-CUDA-Complex-Full-Release
	ctest -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_CUDA=1 -DQMC_COMPLEX=1 -DQMC_MIXED_PRECISION=0 -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
#	module unload mpi
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	;;
    "build_intel2015_cuda")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel/bin/compilervars.sh intel64
	source /opt/intel/impi_latest/bin64/mpivars.sh
	export QMCPACK_TEST_SUBMIT_NAME=Intel2015-CUDA-Release
	ctest -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DQMC_CUDA=1 -DQE_BIN=${QE_BIN} -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
	export PATH=$OLD_PATH
	export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
	export MANPATH=$OLD_MANPATH
	export NLSPATH=$OLD_NLSPATH
	export CPATH=$OLD_CPATH
	export LIBRARY_PATH=$OLD_LIBRARY_PATH
	export MKLROOT=$OLD_MKLROOT
	;;
    "build_intel2015_cuda_complex")
	export OLD_PATH=$PATH
	export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export OLD_MANPATH=$MANPATH
	export OLD_NLSPATH=$NLSPATH
	export OLD_CPATH=$CPATH
	export OLD_LIBRARY_PATH=$LIBRARY_PATH
	export OLD_MKLROOT=$MKLROOT
	source /opt/intel/bin/compilervars.sh intel64
	source /opt/intel/impi_latest/bin64/mpivars.sh
	export QMCPACK_TEST_SUBMIT_NAME=Intel2015-CUDA-Complex-Release
	ctest -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DQMC_CUDA=1 -DQMC_COMPLEX=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -VV -E 'long' --timeout 1800
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

