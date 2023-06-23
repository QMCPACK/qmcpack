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
	    buildsys="clangnewmpi gccnewnompi gccnewmpi gccoldmpi clangoffloadnompi_offloadcuda clangoffloadmpi_offloadcuda clangoffloadmpi_offloadcuda_complex clangnewmpi_complex gccnewnompi_complex gccnewmpi_complex clangnewmpi_mixed gccnewnompi_mixed gccnewmpi_mixed clangnewmpi_mixed_complex gccnewnompi_mixed_complex gccnewmpi_mixed_complex"
	else
	    buildsys="clangnewmpi gccnewmpi clangoffloadmpi_offloadcuda clangnewmpi_complex clangnewmpi_mixed clangnewmpi_mixed_complex"
	fi
	export QMC_DATA=/scratch/${USER}/QMC_DATA_FULL # Route to directory containing performance test files
	;;
    nitrogen2 )
#	if [[ $jobtype == "nightly" ]]; then
#	    buildsys="amdclangnompi_offloadhip amdclangnompi gccnewnompi_legacycu2hip gccnewnompi gccnewmpi gccoldmpi gccnewmpi_legacycu2hip clangnewmpi clangnewmpi_complex gccnewnompi_complex gccnewmpi_complex clangnewmpi_mixed gccnewnompi_mixed gccnewmpi_mixed clangnewmpi_mixed_complex gccnewnompi_mixed_complex gccnewmpi_mixed_complex"
#	else
#	    buildsys="amdclangnompi_offloadhip amdclangnompi gccnewnompi_legacycu2hip gccnewmpi gccoldnompi clangnewmpi clangnewmpi_complex clangnewmpi_mixed clangnewmpi_mixed_complex"
#	fi
	if [[ $jobtype == "nightly" ]]; then
	    buildsys="amdclangnompi gccnewnompi clangnewmpi amdclangnompi_offloadhip"
	else
	    buildsys="amdclangnompi gccnewnompi clangnewmpi amdclangnompi_offloadhip"
	fi
	export QMC_DATA=/scratch/${USER}/QMC_DATA_FULL # Route to directory containing performance test files
#	export amdgpuarch=`/usr/bin/rocminfo | awk '/gfx/ {print $2; exit;}'` # Does not work (?)
	export amdgpuarch=gfx90a
	;;
    nitrogen )
#	if [[ $jobtype == "nightly" ]]; then
#	    buildsys="gccnewnompi gccnewmpi gccoldmpi clangnewmpi clangnewmpi_complex gccnewnompi_complex gccnewmpi_complex clangnewmpi_mixed gccnewnompi_mixed gccnewmpi_mixed clangnewmpi_mixed_complex gccnewnompi_mixed_complex gccnewmpi_mixed_complex"
#	else
#	    buildsys="clangnewmpi gccnewmpi clangnewmpi_complex clangnewmpi_mixed clangnewmpi_mixed_complex"
#	fi
	if [[ $jobtype == "nightly" ]]; then
	    buildsys="amdclangnompi gccnewnompi clangnewmpi amdclangnompi_offloadhip"
	else
	    buildsys="amdclangnompi gccnewnompi clangnewmpi amdclangnompi_offloadhip"
	fi
	export QMC_DATA=/scratch/${USER}/QMC_DATA_FULL # Route to directory containing performance test files
	export amdgpuarch=gfx906
	;;
    * )
	echo Unknown host will use gccnew only
	buildsys="gccnewnompi"
	;;
esac

case "$jobtype" in
    weekly )
	export GLOBALTCFG="-j 48 --timeout 7200 -VV"
#	export QMC_OPTIONS="-DQMC_PERFORMANCE_NIO_MAX_ATOMS=256"
	export QMC_OPTIONS="-DQMC_PERFORMANCE_NIO_MAX_ATOMS=256;-DQMC_PERFORMANCE_C_MOLECULE_MAX_ATOMS=64;-DQMC_PERFORMANCE_C_GRAPHITE_MAX_ATOMS=64"
	export LIMITEDTESTS=""
	export LESSLIMITEDTESTS=""
	;;
    nightly )
#	export GLOBALTCFG="-j 48 --timeout 900 -VV"
	export GLOBALTCFG="-j 48 --timeout 300 -VV"
#	export QMC_OPTIONS="-DQMC_PERFORMANCE_NIO_MAX_ATOMS=255;-DQMC_PERFORMANCE_C_MOLECULE_MAX_ATOMS=12;-DQMC_PERFORMANCE_C_GRAPHITE_MAX_ATOMS=16"
	export QMC_OPTIONS="-DQMC_PERFORMANCE_NIO_MAX_ATOMS=128;-DQMC_PERFORMANCE_C_MOLECULE_MAX_ATOMS=12;-DQMC_PERFORMANCE_C_GRAPHITE_MAX_ATOMS=16"
        export LIMITEDTESTS="-LE unstable --exclude-regex 'short-.*|long-.*|example.*'"
#	export LESSLIMITEDTESTS="--exclude-regex 'long-.*'"
	export LESSLIMITEDTESTS="-E long"
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

# LLVM Offload bug workaround
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
echo --- Activating environment $ourenv
spack env activate $ourenv
echo --- Sourcing environment $ourenv
if [ ! -e $HOME/apps/spack/var/spack/environments/$ourenv/loads ]; then
    echo Loads file missing for environment $ourenv
    exit 1
fi
source $HOME/apps/spack/var/spack/environments/$ourenv/loads    
# Compiler sanity check:
which gcc
which clang
which mpicc
# Extra configuration for this build
# All base sw should be available via the environments

# Ensure GNU C++ library available. Problem symptoms:
# $ bin/qmcpack
# bin/qmcpack: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.26' not found (required by bin/qmcpack)

echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`which gcc|sed 's/bin\/gcc/lib64/g'`
echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH

# Setup additional python paths when gccnew and therefore pyscf available. TO DO: test for module availability
case "$sys" in
    gccnew*mpi*) echo $ourenv
		
		if [[ $sys != *"nompi"* ]]; then
		    # Make PySCF workflows available
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
    clangoffload*mpi*) echo $ourenv
		;;
    amdclang*) echo $ourenv
	       export ROCM_PATH=/opt/rocm
	       export PATH=$PATH:$ROCM_PATH/bin
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
    if [[ $sys == *"amdclang"* ]]; then
	clangname=amdclang
    else
	clangname=clang
    fi
    
    if [[ $sys == *"clangdev"* ]]; then
	compilerversion=Dev
    else
	compilerversion=`$clangname --version|grep ^clang|sed -e 's/^.* version //g' -e 's/(.*//g'|sed 's/\..*//g'`
    fi
    if [[ $sys == *"nompi"* ]]; then
	QMCPACK_TEST_SUBMIT_NAME=Clang${compilerversion}-NoMPI
	CTCFG="-DCMAKE_C_COMPILER=$clangname -DCMAKE_CXX_COMPILER=$clangname++ -DQMC_MPI=0"
    else
	QMCPACK_TEST_SUBMIT_NAME=Clang${compilerversion}
	CTCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_MPI=1"
	export OMPI_CC=$clangname
	export OMPI_CXX=$clangname++
    fi
    if [[ $sys == *"amdclang"* ]]; then
	QMCPACK_TEST_SUBMIT_NAME=AMD${QMCPACK_TEST_SUBMIT_NAME}
    fi
    
# Clang OpenMP offload CUDA builds. Setup here due to clang specific arguments
    if [[ $sys == *"offloadcuda"* ]]; then
        QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Offload-CUDA
        CTCFG="$CTCFG -DCMAKE_CXX_FLAGS=-Wno-unknown-cuda-version"
#	QMC_OPTIONS="${QMC_OPTIONS};-DENABLE_OFFLOAD=ON;-DUSE_OBJECT_TARGET=ON;-DENABLE_CUDA=ON;-DCMAKE_CUDA_ARCHITECTURES=70;-DCMAKE_CUDA_HOST_COMPILER=`which gcc`"
	QMC_OPTIONS="${QMC_OPTIONS};-DENABLE_OFFLOAD=ON;-DUSE_OBJECT_TARGET=ON;-DENABLE_CUDA=ON;-DCMAKE_CUDA_ARCHITECTURES=70"
    fi
    if [[ $sys == *"offloadhip"* ]]; then
        QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Offload-CUDA2HIP
#	QMC_OPTIONS="${QMC_OPTIONS};-DENABLE_OFFLOAD=ON;-DENABLE_CUDA=ON;-DQMC_CUDA2HIP=ON;-DOFFLOAD_TARGET=amdgcn-amd-amdhsa;-DOFFLOAD_ARCH=$amdgpuarch"
	QMC_OPTIONS="${QMC_OPTIONS};-DENABLE_OFFLOAD=ON;-DENABLE_CUDA=ON;-DQMC_CUDA2HIP=ON;-DCMAKE_HIP_ARCHITECTURES=$amdgpuarch"
    fi
fi

#OLD#AOMP (fork of Clang/LLVM)
#OLDif [[ $sys == *"aomp"* ]]; then
#OLD    compilerversion=`aompversion|sed 's/-.*//g'`
#OLD    if [[ $sys == *"nompi"* ]]; then
#OLD	QMCPACK_TEST_SUBMIT_NAME=AOMP${compilerversion}-Offload-NoMPI
#OLD	CTCFG="-DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DQMC_MPI=0"
#OLD	QMC_OPTIONS="${QMC_OPTIONS};-DENABLE_OFFLOAD=ON;-DOFFLOAD_TARGET=amdgcn-amd-amdhsa;-DOFFLOAD_ARCH=${amdgpuarch}"
#OLD    else
#OLD	QMCPACK_TEST_SUBMIT_NAME=AOMP${compilerversion}-Offload
#OLD	CTCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_MPI=1"
#OLD	QMC_OPTIONS="${QMC_OPTIONS};-DENABLE_OFFLOAD=ON;-DOFFLOAD_TARGET=amdgcn-amd-amdhsa;-DOFFLOAD_ARCH=${amdgpuarch}"
#OLD	export OMPI_CC=clang
#OLD	export OMPI_CXX=clang++
#OLD    fi
#OLDfi

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

# General CUDA setup for offload and legacy cuda builds
# Use system installed CUDA since this will match the driver. May conflict will a different version spack installed cuda
# TODO: Ensure consistent CUDA versions for nvhpc+cuda, spack sourced compilers etc.

# ASSUME CORRECT CUDA VERSION ALREADY ON PATH , e.g. correct CUDA spack module loaded

#if [[ $sys == *"legacycuda"* ]]; then
#    if [ -e /usr/local/cuda/bin/nvcc ]; then
#        export CUDAVER=`cat /usr/local/cuda/version.json | python3 -c "import sys, json; print(json.load(sys.stdin)['cuda']['version'])"`
#        echo --- Found nvcc in /usr/local/cuda , apparent version $CUDAVER . Adding to PATH
#        export PATH=/usr/local/cuda/bin:${PATH}
#        export LD_LIBRARY_PATH=/usr/local/cuda/lib64:${LD_LIBRARY_PATH}
#    else
#        echo --- Did not find expected nvcc compiler for CUDA build. Error.
#        exit 1
#    fi
#else
#    if [[ $sys == *"offloadcuda"* ]]; then
#	echo --- FORCING CUDA 11.2 FOR OFFLOAD BUILD TO WORKAROUND https://github.com/llvm/llvm-project/issues/54633
#        if [ -e /usr/local/cuda-11.2/bin/nvcc ]; then
#            export CUDAVER=`cat /usr/local/cuda-11.2/version.json | python3 -c "import sys, json; print(json.load(sys.stdin)['cuda']['version'])"`
#            echo --- Found nvcc in /usr/local/cuda-11.2 , apparent version $CUDAVER . Adding to PATH
#            export PATH=/usr/local/cuda-11.2/bin:${PATH}
#            export LD_LIBRARY_PATH=/usr/local/cuda-11.2/lib64:${LD_LIBRARY_PATH}
#        else
#            echo --- Did not find expected nvcc compiler for CUDA build. Error.
#            exit 1
#        fi
#    fi
#fi

# Legacy CUDA builds setup
if [[ $sys == *"legacycuda"* ]]; then
# Specify GPUs for testing. Obtain device IDs via "nvidia-smi -L"
#export CUDA_VISIBLE_DEVICES=
QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Legacy-CUDA
CTCFG="$CTCFG -DQMC_CUDA=1"
fi

# Legacy CUDA2HIP builds setup
# TODO: Ensure consistent CUDA versions for nvhpc+cuda, spack sourced compilers etc.
if [[ $sys == *"legacycu2hip"* ]]; then
    export ROCM_PATH=/opt/rocm
    export PATH=${PATH}:${ROCM_PATH}/bin:${ROCM_PATH}/opencl/bin
    QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Legacy-CUDA2HIP
    QMC_OPTIONS="${QMC_OPTIONS};-DQMC_CUDA=ON;-DQMC_CUDA2HIP=ON;-DCMAKE_HIP_ARCHITECTURES=${amdgpuarch}"
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
    if [[ $sys == *"offload"* ]]; then
	echo "AFQMC is disabled for this offload build."
	CTCFG="$CTCFG -DBUILD_AFQMC=0"
    else
	echo "AFQMC build option is enabled."
	CTCFG="$CTCFG -DBUILD_AFQMC=1"
    fi
fi

# Adjust which tests are run to control overall runtime
case "$sys" in
    *intel2020*|*gccnew*|*clangnew*|*gcc*legacycuda*|*gcc*cu2hip*) echo "Running full ("less limited") test set for $sys"
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
# To workaround any concurrency problems put -DN_CONCURRENT_TESTS=1
ctest ${GLOBALTCFG} ${THETESTS}
echo --- END ctest `date`
else
echo --- START ctest `date` 
# To workaround any concurrency problems put -DN_CONCURRENT_TESTS=1
echo ctest ${CTCFG} ${GLOBALTCFG} -DQMC_OPTIONS=${QMC_OPTIONS} ${THETESTS} -S $PWD/../qmcpack/CMake/ctest_script.cmake,$ctestscriptarg
     ctest ${CTCFG} ${GLOBALTCFG} -DQMC_OPTIONS=${QMC_OPTIONS} ${THETESTS} -S $PWD/../qmcpack/CMake/ctest_script.cmake,$ctestscriptarg
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
