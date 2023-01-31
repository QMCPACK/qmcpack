#!/bin/bash

echo --- Script START `date`

localonly=no
#localonly=yes

# Type of job determined by filename of script

if [[ $0 == *"nightly"* ]]; then
    jobtype=nightly
else
    if [[ $0 == *"weekly"* ]]; then
	jobtype=weekly
    fi
fi
     
case "$jobtype" in
    nightly )
	echo --- Nightly job
	;;
    weekly )
	echo --- Weekly job
	;;
    * )
# If a new jobtype is added, add support in all similar case statements below    
	echo Unknown jobtype $jobtype
	exit 1
	;;
esac

if [[ $localonly == "yes" ]]; then
    echo --- Local CMake/Make/CTest only. No cdash drop.
fi

if [[ $jobtype == "weekly" ]]; then
    if [ -e `dirname "$0"`/ornl_update.sh ]; then
	echo --- Running compiler updates
	source `dirname "$0"`/ornl_update.sh
    fi
fi



if [ -e `dirname "$0"`/ornl_versions.sh ]; then
    source `dirname "$0"`/ornl_versions.sh
else
    echo Did not find version numbers script ornl_versions.sh
    exit 1
fi


plat=`lscpu|grep Vendor|sed 's/.*ID:[ ]*//'`
case "$plat" in
    GenuineIntel )
	ourplatform=Intel
	;;
    AuthenticAMD )
	ourplatform=AMD
	;;
    * )
	# ARM support should be trivial, but is not yet done
	echo Unknown platform
	exit 1
	;;
esac
echo --- Using $ourplatform architecture

ourhostname=`hostname|sed 's/\..*//g'`
echo --- Host is $ourhostname
case "$ourhostname" in
    sulfur )
	if [[ $jobtype == "nightly" ]]; then
	    buildsys="clangnewmpi gccnewnompi gccnewmpi clangnewmpi_complex gccnewnompi_complex gccnewmpi_complex clangnewmpi_mixed gccnewnompi_mixed gccnewmpi_mixed clangnewmpi_mixed_complex gccnewnompi_mixed_complex gccnewmpi_mixed_complex gccoldmpi gccoldmpi_legacycuda gccoldnompi_legacycuda gccoldnompi_legacycuda gccoldnompi_legacycuda_complex gccoldnompi_legacycuda_full gccoldnompi_legacycuda_full_complex"
	else
	    buildsys="clangnewmpi gccnewmpi clangnewmpi_complex clangnewmpi_mixed clangnewmpi_mixed_complex gccoldmpi_legacycuda gccoldmpi_legacycuda_full gccoldmpi_legacycuda_complex"
	fi
	export QMC_DATA=/scratch/${USER}/QMC_DATA_FULL # Route to directory containing performance test files
	;;
    nitrogen )
	if [[ $jobtype == "nightly" ]]; then
	    buildsys="gccoldnompi_legacycu2hip gccoldmpi_legacycu2hip gccoldnompi_legacycu2hip_complex gccoldnompi_legacycu2hip_full gccoldnompi_legacycu2hip_full_complex clangnewmpi gccnewnompi gccnewmpi clangnewmpi_complex gccnewnompi_complex gccnewmpi_complex clangnewmpi_mixed gccnewnompi_mixed gccnewmpi_mixed clangnewmpi_mixed_complex gccnewnompi_mixed_complex gccnewmpi_mixed_complex gccoldmpi gccoldmpi_legacycuda gccoldnompi_legacycuda gccoldnompi_legacycuda gccoldnompi_legacycuda_complex gccoldnompi_legacycuda_full gccoldnompi_legacycuda_full_complex"
	else
	    buildsys="clangnewmpi gccnewmpi clangnewmpi_complex clangnewmpi_mixed clangnewmpi_mixed_complex gccoldmpi_legacycuda gccoldmpi_legacycuda_full gccoldmpi_legacycuda_complex gccoldmpi_legacycu2hip gccoldmpi_legacycu2hip_full gccoldmpi_legacycu2hip_complex"
	fi
	export QMC_DATA=/scratch/${USER}/QMC_DATA_FULL # Route to directory containing performance test files
	export amdgpuarch=`rocminfo | awk '/gfx/ {print $2; exit;}'`
	;;
    nitrogen2 )
	if [[ $jobtype == "nightly" ]]; then
	    buildsys="gccoldnompi_legacycu2hip gccoldmpi_legacycu2hip gccoldnompi_legacycu2hip_complex gccoldnompi_legacycu2hip_full gccoldnompi_legacycu2hip_full_complex clangnewmpi gccnewnompi gccnewmpi gccoldmpi clangnewmpi_complex gccnewnompi_complex gccnewmpi_complex clangnewmpi_mixed gccnewnompi_mixed gccnewmpi_mixed clangnewmpi_mixed_complex gccnewnompi_mixed_complex gccnewmpi_mixed_complex"
	else
	    buildsys="gccoldmpi_legacycu2hip gccoldmpi_legacycu2hip_full gccoldmpi_legacycu2hip_complex clangnewmpi gccnewmpi clangnewmpi_complex clangnewmpi_mixed clangnewmpi_mixed_complex"
	fi
	export QMC_DATA=/scratch/${USER}/QMC_DATA_FULL # Route to directory containing performance test files
	export amdgpuarch=`rocminfo | awk '/gfx/ {print $2; exit;}'`
	;;
    * )
	echo Unknown host will use gccnew only
	buildsys="gccnewnompi"
	;;
esac

case "$jobtype" in
    weekly )
	export GLOBALTCFG="-j 48 --timeout 10800 -VV"
	export QMC_OPTIONS="-DQMC_PERFORMANCE_NIO_MAX_ATOMS=256"
	export LIMITEDTESTS=""
	export LESSLIMITEDTESTS=""
	;;
    nightly )
	export GLOBALTCFG="-j 48 --timeout 600 -VV"
	export QMC_OPTIONS="-DQMC_PERFORMANCE_NIO_MAX_ATOMS=64"
	export LIMITEDTESTS="--tests-regex deterministic -E long-"
	export LESSLIMITEDTESTS="-E long-"
	;;
    * )
	echo Unknown jobtype $jobtype
	exit 1
	;;
esac

# Directory in which to run tests. Should be an absolute path and fastest usable filesystem
test_path=/scratch/${USER}

test_dir=${test_path}/QMCPACK_CI_BUILDS
if [ ! -e ${test_dir} ]; then
    mkdir ${test_dir}
fi

export OMP_NUM_THREADS=16

# Intel2019.1 MPI configure setting to avoid MPI crash
# via https://software.intel.com/en-us/forums/intel-clusters-and-hpc-technology/topic/799716
#export FI_PROVIDER=sockets
export I_MPI_FABRICS=shm

# LLVM Offload bug workaround 2021-03-02
export LIBOMP_USE_HIDDEN_HELPER_TASK=OFF

module() { eval `/usr/bin/modulecmd bash $*`; }

export SPACK_ROOT=$HOME/apps/spack
export PATH=$SPACK_ROOT/bin:$PATH
. $SPACK_ROOT/share/spack/setup-env.sh


echo --- Spack list
spack find
echo --- Modules list
module list
echo --- End listings

spack load --first git

module list
if [ -e ${test_path} ]; then

if [ ! -e ${test_dir} ]; then
mkdir ${test_dir}
fi

if [ -e ${test_dir} ]; then
cd ${test_dir}

# Minimize load on GitHub by maintaining a local cloned git used for all builds
if [ ! -e qmcpack ]; then
echo --- Cloning QMCPACK git `date`
git clone https://github.com/QMCPACK/qmcpack.git --depth 1
else
cd qmcpack
echo --- Updating local QMCPACK git `date`
git pull
git status
cd ..
fi

# Sanity check cmake config file present
if [ -e qmcpack/CMakeLists.txt ]; then

export PYTHONPATH=${test_dir}/qmcpack/nexus/lib
echo --- PYTHONPATH=$PYTHONPATH

echo --- Starting test builds and tests
for sys in $buildsys
do

echo --- START build configuration $sys `date`
syscompilermpi=`echo $sys|sed 's/_.*//g'`
echo --- Compilermpi=$syscompilermpi
cd ${test_dir}

if [ -e build_$sys ]; then
rm -r -f build_$sys
fi
mkdir build_$sys
cd build_$sys

# Use subshell to allow compiler setup to contaminate the environment
# e.g. No "unload" capability is provided for Intel compiler scripts
(


ourenv=env${syscompilermpi}
echo --- Sourcing environment $ourenv
source $HOME/apps/spack/var/spack/environments/$ourenv/loads    
# Compiler sanity check:
which gcc
which clang
which mpicc
# Extra configuration for this build
# All base sw should be available via the environments

# Choose python on a per build type basis to minimize risk of contamination by e.g. older/newer HDF5 picked up via python modules
case "$sys" in
    gccnew*mpi*) echo $ourenv
		
		if [[ $sys != *"nompi"* ]]; then
		    # Make PySCF available
		    export PYTHONPATH=${test_dir}/qmcpack/utils/afqmctools/:$PYTHONPATH
		    export PYTHONPATH=${test_dir}/qmcpack/src/QMCTools/:$PYTHONPATH

		    # For debugging module availability etc. can check if afmctools are working here
		    #${test_dir}/qmcpack/utils/afqmctools/bin/pyscf_to_afqmc.py
		fi
		;;
    gccold*mpi*) echo $ourenv
	
		;;
    clangnew*mpi*) echo $ourenv
		;;
    
    *) echo "Problems: Unknown build environment"
	exit 1
;;
esac
module list


# Construct test name and configure flags
# Compiler and major version, MPI or not
if [[ $sys == *"gcc"* ]]; then
    if [[ $sys == *"gccdev"* ]]; then
	compilerversion=Dev
    else
	compilerversion=`gcc --version|grep ^gcc|sed 's/^.* //g'|sed 's/\..*//g'`
    fi
if [[ $sys == *"nompi"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=GCC${compilerversion}-NoMPI
CTCFG="-DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DQMC_MPI=0"
else
QMCPACK_TEST_SUBMIT_NAME=GCC${compilerversion}
CTCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_MPI=1"
export OMPI_CC=gcc
export OMPI_CXX=g++

if [[ $sys == *"gccnew"* ]]; then
# Add QE to any gccnew MPI builds
# Restrict to gccnew to avoid problems with mismatched libraries, mpi etc.
CTCFG="$CTCFG -DQE_BIN=${QE_BIN}" 
fi

if [[ $sys == *"gcclegacycuda"* ]]; then
# Add QE to any gcclegacycuda MPI builds
# GCC compilers will be mismatched from QE and QMCPACK builds
CTCFG="$CTCFG -DQE_BIN=${QE_BIN}" 
fi

fi
fi

#Clang/LLVM
if [[ $sys == *"clang"* ]]; then
    if [[ $sys == *"clangdev"* ]]; then
	compilerversion=Dev
    else
	compilerversion=`clang --version|grep ^clang|sed -e 's/^.* version //g' -e 's/(.*//g'|sed 's/\..*//g'`
    fi
    if [[ $sys == *"nompi"* ]]; then
	QMCPACK_TEST_SUBMIT_NAME=Clang${compilerversion}-NoMPI
	CTCFG="-DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DQMC_MPI=0"
    else
	QMCPACK_TEST_SUBMIT_NAME=Clang${compilerversion}
	CTCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_MPI=1"
	export OMPI_CC=clang
	export OMPI_CXX=clang++
    fi

# Clang OpenMP offload CUDA builds. Setup here due to clang specific arguments
    if [[ $sys == *"offloadcuda"* ]]; then
        QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Offload-CUDA
        CTCFG="$CTCFG -DCMAKE_CXX_FLAGS=-Wno-unknown-cuda-version"
		QMC_OPTIONS="${QMC_OPTIONS};-DENABLE_OFFLOAD=ON;-DUSE_OBJECT_TARGET=ON;-DENABLE_CUDA=ON;-DCMAKE_CUDA_ARCHITECTURES=sm_70;-DCMAKE_CUDA_HOST_COMPILER=`which gcc`"
		echo *** CHECK/FIX ME
		exit 1
    fi
fi

#AOMP (fork of Clang/LLVM)
if [[ $sys == *"aomp"* ]]; then
    compilerversion=`aompversion|sed 's/-.*//g'`
    if [[ $sys == *"nompi"* ]]; then
	QMCPACK_TEST_SUBMIT_NAME=AOMP${compilerversion}-Offload-NoMPI
	CTCFG="-DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DQMC_MPI=0"
	QMC_OPTIONS="${QMC_OPTIONS};-DENABLE_OFFLOAD=ON;-DOFFLOAD_TARGET=amdgcn-amd-amdhsa;-DOFFLOAD_ARCH=${amdgpuarch}"
    else
	QMCPACK_TEST_SUBMIT_NAME=AOMP${compilerversion}-Offload
	CTCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_MPI=1"
	QMC_OPTIONS="${QMC_OPTIONS};-DENABLE_OFFLOAD=ON;-DOFFLOAD_TARGET=amdgcn-amd-amdhsa;-DOFFLOAD_ARCH=${amdgpuarch}"
	export OMPI_CC=clang
	export OMPI_CXX=clang++
    fi
fi

# Intel
if [[ $sys == *"intel"* ]]; then

if [[ $sys == *"intel2020"* ]]; then
source /opt/intel2020/bin/compilervars.sh intel64
fi
if [[ $sys == *"intel2019"* ]]; then
source /opt/intel2019/bin/compilervars.sh intel64
fi
if [[ $sys == *"intel2018"* ]]; then
source /opt/intel2018/bin/compilervars.sh intel64
fi
# Assume the Intel dumpversion string has format AA.B.C.DDD
# AA is historically major year and a good enough unique identifier,
# but Intel 2020 packages currently (12/2019) come with the 19.1 compiler.
# Use 19.1 for this compiler only to avoid breaking test histories.
if [[ $sys == *"intel2020"* ]]; then
compilerversion=`icc -dumpversion|sed -e 's/\([0-9][0-9]\)\.\([0-9]\)\..*/\1.\2/g'` # AA.B
else
compilerversion=`icc -dumpversion|sed -e 's/\..*//g'` # AA year only
fi

if [[ $sys == *"nompi"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=Intel20${compilerversion}-NoMPI
CTCFG="-DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DQMC_MPI=0"
else
QMCPACK_TEST_SUBMIT_NAME=Intel20${compilerversion}
CTCFG="-DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DQMC_MPI=1"
fi
fi

# NVIDIA HPC SDK compiler
# TODO: Ensure consistent CUDA versions used for future nvhpc+cuda builds
if [[ $sys == *"nvhpc"* ]]; then
    export NVHPC=/opt/nvidia/hpc_sdk
if [ -e $NVHPC/Linux_x86_64/2021/compilers/bin/nvc++ ]; then
    echo --- Found nvc++ NVIDIA HPC SDK compiler. Adding to PATH
    export PATH=${NVHPC}/Linux_x86_64/2021/compilers/bin:${PATH}
else
    echo --- Did not find expected nvc++ compiler for nvhpc build. Error.
    exit 1
fi
compilerversion=`nvc -V|grep nvc|sed 's/^nvc //g'|sed 's/-.*//g'`
if [[ $sys == *"nompi"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=NVHPC${compilerversion}-NoMPI
CTCFG="-DCMAKE_C_COMPILER=nvc -DCMAKE_CXX_COMPILER=nvc++ -DQMC_MPI=0"
else
QMCPACK_TEST_SUBMIT_NAME=NVHPC${compilerversion}
CTCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_MPI=1"
export OMPI_CC=nvc
export OMPI_CXX=nvc++
fi
fi

# Legacy CUDA builds setup
# TODO: Ensure consistent CUDA versions for nvhpc+cuda, spack sourced compilers etc.
if [[ $sys == *"legacycuda"* ]]; then
#export CUDAVER=11.6
export CUDAVER=11.7 # HYPOTHESIS: Broken 2022-05-28
if [ -e /usr/local/cuda-${CUDAVER}/bin/nvcc ]; then
    echo --- Found nvcc from CUDA ${CUDAVER} . Adding to PATH
    export PATH=/usr/local/cuda-${CUDAVER}/bin:${PATH}
    export LD_LIBRARY_PATH=/usr/local/cuda-${CUDAVER}/lib64:${LD_LIBRARY_PATH}
else
    echo --- Did not find expected nvcc compiler for CUDA build. Error.
    exit 1
fi
# Specify GPUs for testing. Obtain device IDs via "nvidia-smi -L"
#export CUDA_VISIBLE_DEVICES=
QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Legacy-CUDA
CTCFG="$CTCFG -DQMC_CUDA=1"
fi

# Legacy CUDA2HIP builds setup
# TODO: Ensure consistent CUDA versions for nvhpc+cuda, spack sourced compilers etc.
if [[ $sys == *"legacycu2hip"* ]]; then
#    export ROCMVER=4.5
#    export PATH=$PATH:/opt/rocm-${ROCMVER}/bin:/opt/rocm-${ROCMVER}/opencl/bin
    export ROCM_PATH=/opt/rocm
    export PATH=${PATH}:${ROCM_PATH}/bin:${ROCM_PATH}/opencl/bin
    QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Legacy-CUDA2HIP
    QMC_OPTIONS="${QMC_OPTIONS};-DENABLE_CUDA=ON;-DQMC_CUDA2HIP=ON;-DCMAKE_HIP_ARCHITECTURES=${amdgpuarch}"
fi

# MKL
# MKLROOT set in sourced Intel mklvars.sh 
if [[ $sys == *"mkl"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-MKL
# MKL setup used by many builds for BLAS, LAPACK etc.
source /opt/intel2020/mkl/bin/mklvars.sh intel64
fi

# Complex
if [[ $sys == *"complex"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Complex
CTCFG="$CTCFG -DQMC_COMPLEX=1"
fi

# Mixed/Full precision
if [[ $sys == *"mixed"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Mixed
CTCFG="$CTCFG -DQMC_MIXED_PRECISION=1"
fi
if [[ $sys == *"full"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Full
CTCFG="$CTCFG -DQMC_MIXED_PRECISION=0"
fi

# Boilerplate for all tests
CTCFG="$CTCFG -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1"

# Selectively enable AFQMC
if [[ $sys == *"nompi"* ]]; then
    echo "AFQMC is disabled for this build without MPI."
    CTCFG="$CTCFG -DBUILD_AFQMC=0"
else
    if [[ $sys == *"cu2hip"* ]]; then
	echo "AFQMC is disabled for this HIP build."
	CTCFG="$CTCFG -DBUILD_AFQMC=0"
    else
	echo "AFQMC build option is enabled."
	CTCFG="$CTCFG -DBUILD_AFQMC=1"
    fi
fi

# Adjust which tests are run to control overall runtime
case "$sys" in
    *intel2020*|*gccnew*|*clangnew*|*nvhpc*|*gcc*legacycuda*|*gcc*cu2hip*) echo "Running full ("less limited") test set for $sys"
							     THETESTS=$LESSLIMITEDTESTS
							     ;;
    *) echo "Running limited test set for $sys"
       THETESTS=$LIMITEDTESTS
       ;;
esac
#THETESTS=$LIMITEDTESTS # for DEBUG. Remove for production.
echo $THETESTS

if [[ $sys == *"debug"* ]]; then
    export QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Debug
    CTCFG="-DCMAKE_BUILD_TYPE=Debug $CTCFG"
    ctestscriptarg=debug
else
    export QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Release
    ctestscriptarg=release
fi

echo $QMCPACK_TEST_SUBMIT_NAME
echo $CTCFG
if [[ $localonly == "yes" ]]; then
echo --- START cmake `date` 
cmake ${CTCFG} ${GLOBALTCFG} -DQMC_OPTIONS=${QMC_OPTIONS} ../qmcpack/ 
echo --- END cmake `date`
echo --- START make `date` 
make -j 96
echo --- END make `date`
echo --- START ctest `date` 
#Workaround CUDA concurrency problems
#case "$sys" in
#    *cuda*)
#	ctest ${GLOBALTCFG} ${THETESTS} -DN_CONCURRENT_TESTS=1
#	;;
#    *)
	ctest ${GLOBALTCFG} ${THETESTS}
#	;;
#esac
echo --- END ctest `date`
else
echo --- START ctest `date` 
echo ctest ${CTCFG} ${GLOBALTCFG} -DQMC_OPTIONS=${QMC_OPTIONS} -S $PWD/../qmcpack/CMake/ctest_script.cmake,$ctestscriptarg ${THETESTS}
#Workaround CUDA concurrency problems
#case "$sys" in
#    *cuda*)
#	ctest ${CTCFG} ${GLOBALTCFG} -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,$ctestscriptarg ${THETESTS} -DN_CONCURRENT_TESTS=1
#	;;
#    *)
	ctest ${CTCFG} ${GLOBALTCFG} -DQMC_OPTIONS=${QMC_OPTIONS} -S $PWD/../qmcpack/CMake/ctest_script.cmake,$ctestscriptarg ${THETESTS}
#	;;
#esac
echo --- END ctest `date`
fi

)

module purge

echo --- END build configuration $sys `date`
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
