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
	    buildsys="gccnewnompi_debug_asan gccnewmpi_mkl clangnewmpi gccnewnompi gccnewmpi gccoldmpi clangnewmpi_complex gccnewnompi_complex gccnewmpi_complex clangnewmpi_mixed gccnewnompi_mixed gccnewmpi_mixed clangnewmpi_mixed_complex gccnewnompi_mixed_complex gccnewmpi_mixed_complex \
		clangoffloadmpi_offloadcuda \
		clangoffloadnompi_offloadcuda clangoffloadnompi_offloadcuda_debug \
		clangoffloadnompi_offloadcuda_complex clangoffloadnompi_offloadcuda_complex_debug \
		clangoffloadnompi_offloadcuda_mixed clangoffloadnompi_offloadcuda_mixed_debug \
		clangoffloadnompi_offloadcuda_mixed_complex clangoffloadnompi_offloadcuda_mixed_complex_debug"
	else
	    buildsys="gccnewmpi_mkl clangnewmpi gccnewmpi clangnewmpi_complex clangnewmpi_mixed clangnewmpi_mixed_complex clangoffloadmpi_offloadcuda"
	fi
	export QMC_DATA=/scratch/${USER}/QMC_DATA_FULL # Route to directory containing performance test files
	;;
    nitrogen2 )
	if [[ $jobtype == "nightly" ]]; then
	    buildsys="amdclangnompi_offloadhip_complex gccnewnompi gccnewnompi_aocl gccnewnompi_complex gccnewnompi_debug gccnewnompi_complex_debug gccnewnompi_mixed_debug gccnewnompi_mixed_complex_debug gccnewnompi_aocl_mixed_complex_debug gccnewmpi gccnewmpi_aocl clangnewmpi \
		amdclangnompi amdclangnompi_debug \
		amdclangnompi_offloadhip amdclangnompi_offloadhip_debug \
		amdclangnompi_offloadhip_complex_debug \
		amdclangnompi_offloadhip_mixed amdclangnompi_offloadhip_mixed_debug \
		amdclangnompi_offloadhip_mixed_complex amdclangnompi_offloadhip_mixed_complex_debug"
	else
	    buildsys="gccnewmpi gccnewmpi_aocl amdclangnompi gccnewnompi gccnewnompi_aocl clangnewmpi amdclangnompi_offloadhip"
	fi
	export QMC_DATA=/scratch/${USER}/QMC_DATA_FULL # Route to directory containing performance test files
	export amdgpuarch=`/usr/bin/rocminfo | awk '/gfx/ {print $2; exit 0;}'`
	;;
    nitrogen )
	if [[ $jobtype == "nightly" ]]; then
	    buildsys="gccnewmpi amdclangnompi gccnewnompi clangnewmpi amdclangnompi_offloadhip"
	else
	    buildsys="gccnewmpi amdclangnompi gccnewnompi clangnewmpi amdclangnompi_offloadhip"
	fi
	export QMC_DATA=/scratch/${USER}/QMC_DATA_FULL # Route to directory containing performance test files
	export amdgpuarch=`/usr/bin/rocminfo | awk '/gfx/ {print $2; exit 0;}'`
	;;
    * )
	echo Unknown host will use gccnew only
	buildsys="gccnewnompi"
	;;
esac

case "$jobtype" in
    weekly )
	export PARALLELCFG="-j 48"
	export QMC_OPTIONS="-DQMC_PERFORMANCE_NIO_MAX_ATOMS=256;-DQMC_PERFORMANCE_C_GRAPHITE_MAX_ATOMS=64"
	if [[ $sys == *"complex"* ]]; then
	    export QMC_OPTIONS="${QMC_OPTIONS};-DQMC_PERFORMANCE_C_MOLECULE_MAX_ATOMS=64"
	fi
	export LIMITEDTESTS="--exclude-regex long-"
	export LESSLIMITEDTESTS=""
	;;
    nightly )
	export PARALLELCFG="-j 48"
	export QMC_OPTIONS="-DQMC_PERFORMANCE_NIO_MAX_ATOMS=16;-DQMC_PERFORMANCE_C_GRAPHITE_MAX_ATOMS=16"
	if [[ $sys == *"complex"* ]]; then
	    export QMC_OPTIONS="${QMC_OPTIONS};-DQMC_PERFORMANCE_C_MOLECULE_MAX_ATOMS=12"
	fi
	export LIMITEDTESTS="--label-regex deterministic"
#	export LIMITEDTESTS="--exclude-regex 'long-|short-'"
	export LESSLIMITEDTESTS="--exclude-regex long-"
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

# Consider using a GCC toolset on Red Hat systems to use
# recent compilers with better architecture support.
# e.g. dnf install gcc-toolset-11
#if [ -e /opt/rh/gcc-toolset-12/enable ]; then
#    echo --- Using gcc-toolset-12 for newer compilers
#    source /opt/rh/gcc-toolset-12/enable 
#fi

export OMP_NUM_THREADS=16

# Intel2019.1 MPI configure setting to avoid MPI crash
# via https://software.intel.com/en-us/forums/intel-clusters-and-hpc-technology/topic/799716
#export FI_PROVIDER=sockets
export I_MPI_FABRICS=shm

# OpenMPI 5.x policies
export OMPI_MCA_rmaps_base_oversubscribe=true

# LLVM Offload bug workaround
export LIBOMP_USE_HIDDEN_HELPER_TASK=OFF

# Avoid use of staging buffer on AMD GPUs, see https://github.com/QMCPACK/qmcpack/pull/5339
export LIBOMPTARGET_AMDGPU_MAX_ASYNC_COPY_BYTES=0

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
# Ensure clang libraries both available and ahead of any system installed versions
		   echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH
		   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`which clang++|sed 's/bin\/clang++/lib\/x86_64-unknown-linux-gnu/g'`
		   echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH
		   ;;
    clangoffload*mpi*) echo $ourenv
# Ensure clang libraries both available and ahead of any system installed versions
		   echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH
		   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`which clang++|sed 's/bin\/clang++/lib\/x86_64-unknown-linux-gnu/g'`
		   echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH
		;;
    amdclang*) echo $ourenv
	       rocvlist=`find /opt/rocm-* -name "rocm-[1-9]*" -print|sort -V -r`
	       for rocp in /opt/rocm $rocvlist 
	       do
		   if [ -e $rocp/bin/rocminfo ]; then
		       echo Found rocminfo under $rocp
		       export ROCM_PATH=$rocp
		       break
		   fi
	       done
	       export PATH=$ROCM_PATH/bin:$ROCM_PATH/llvm/bin:$PATH
	       export LD_LIBRARY_PATH=$ROCM_PATH/lib:$LD_LIBRARYPATH
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
CMCFG="-DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DQMC_MPI=0"
else
QMCPACK_TEST_SUBMIT_NAME=GCC${compilerversion}-MPI
CMCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_MPI=1"
export QMC_OPTIONS="${QMC_OPTIONS};-DMPIEXEC_PREFLAGS=--bind-to\;none\;--oversubscribe"
export OMPI_CC=gcc
export OMPI_CXX=g++

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
	CMCFG="-DCMAKE_C_COMPILER=$clangname -DCMAKE_CXX_COMPILER=$clangname++ -DQMC_MPI=0"
    else
	QMCPACK_TEST_SUBMIT_NAME=Clang${compilerversion}-MPI
	CMCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_MPI=1"
	export QMC_OPTIONS="${QMC_OPTIONS};-DMPIEXEC_PREFLAGS=--bind-to\;none\;--oversubscribe"
	export OMPI_CC=$clangname
	export OMPI_CXX=$clangname++
    fi
    if [[ $sys == *"amdclang"* ]]; then
	QMCPACK_TEST_SUBMIT_NAME=AMD${QMCPACK_TEST_SUBMIT_NAME}
    fi
    
# Clang OpenMP offload CUDA builds. Setup here due to clang specific arguments
    if [[ $sys == *"offloadcuda"* ]]; then
        QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-NVGPU
        CMCFG="$CMCFG -DCMAKE_CXX_FLAGS=-Wno-unknown-cuda-version"
	QMC_OPTIONS="${QMC_OPTIONS};-DQMC_GPU_ARCHS=sm_70"
    fi
    if [[ $sys == *"offloadhip"* ]]; then
        QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-AMDGPU
	QMC_OPTIONS="${QMC_OPTIONS};-DQMC_GPU_ARCHS=$amdgpuarch"
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
CMCFG="-DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DQMC_MPI=0"
else
QMCPACK_TEST_SUBMIT_NAME=Intel20${compilerversion}-MPI
CMCFG="-DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DQMC_MPI=1" # Verify thread binding OK when Intel compiler added back to nightlies
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
CMCFG="-DCMAKE_C_COMPILER=nvc -DCMAKE_CXX_COMPILER=nvc++ -DQMC_MPI=0"
else
QMCPACK_TEST_SUBMIT_NAME=NVHPC${compilerversion}-MPI
CMCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_MPI=1"
export QMC_OPTIONS="${QMC_OPTIONS};-DMPIEXEC_PREFLAGS=--bind-to\;none\;--oversubscribe"
export OMPI_CC=nvc
export OMPI_CXX=nvc++
fi
fi

if [[ $sys == *"aocl"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-AOCL
export QMC_OPTIONS="${QMC_OPTIONS};-DBLA_VENDOR=AOCL"
else
    if [[ $sys == *"mkl"* ]]; then
	QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-MKL
	export QMC_OPTIONS="${QMC_OPTIONS};-DBLA_VENDOR=Intel10_64_dyn"
    else
	export QMC_OPTIONS="${QMC_OPTIONS};-DBLA_VENDOR=OpenBLAS"
    fi
fi

# Complex
if [[ $sys == *"complex"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Complex
CMCFG="$CMCFG -DQMC_COMPLEX=1"
fi

# Mixed/Full precision
if [[ $sys == *"mixed"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Mixed
CMCFG="$CMCFG -DQMC_MIXED_PRECISION=1"
fi
if [[ $sys == *"full"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Full
CMCFG="$CMCFG -DQMC_MIXED_PRECISION=0"
fi

# Boilerplate for all tests
CMCFG="$CMCFG -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1"

# Selectively enable AFQMC
if [[ $sys == *"nompi"* ]]; then
    echo "AFQMC is disabled for this build without MPI."
    CMCFG="$CMCFG -DBUILD_AFQMC=0"
else
    if [[ $sys == *"offload"* ]]; then
	echo "AFQMC is disabled for this offload build."
	CMCFG="$CMCFG -DBUILD_AFQMC=0"
    else
	echo "AFQMC build option is enabled."
	CMCFG="$CMCFG -DBUILD_AFQMC=1"
    fi
fi

# Adjust which tests are run to control overall runtime
case "$sys" in
    *offload*) echo "Running limited test set for $sys"
               THETESTS=${LIMITEDTESTS}
               ;;
    *) echo "Running full/lesslimited test set for $sys"
       THETESTS=${LESSLIMITEDTESTS}
       ;;
esac
echo THETESTS: ${THETESTS}

if [[ $sys == *"debug"* ]]; then
    if [[ $jobtype == *"nightly"* ]]; then
	    export TESTCFG="--timeout 3600 -VV"
	else
	    export TESTCFG="--timeout 10800 -VV"
	fi
    export QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Debug
    CMCFG="-DCMAKE_BUILD_TYPE=Debug $CMCFG"
    ctestscriptarg=debug
else
    if [[ $jobtype == *"nightly"* ]]; then
	    export TESTCFG="--timeout 900 -VV"
	else
	    export TESTCFG="--timeout 10800 -VV"
	fi
    export QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Release
    ctestscriptarg=release
fi

if [[ $sys == *"asan"* ]]; then
    export QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Asan
    export QMC_OPTIONS="${QMC_OPTIONS};-DENABLE_SANITIZER=asan"
fi

echo TEST SUBMIT NAME: $QMCPACK_TEST_SUBMIT_NAME
echo CMCFG: $CMCFG
echo PARALLELCFG: $PARALLELCFG
echo TESTCFG: $TESTCFG
echo QMC_OPTIONS: $QMC_OPTIONS
CMOPTIONS=`echo $QMC_OPTIONS|sed 's/;/ /g'`
echo CMOPTIONS: $CMOPTIONS
if [[ $localonly == "yes" ]]; then
echo --- START cmake `date` 
#cmake ${CMCFG} -DQMC_OPTIONS=${QMC_OPTIONS} ../qmcpack/ 
cmake ${CMCFG} ${CMOPTIONS} ../qmcpack/ 
echo --- END cmake `date`
echo --- START make `date` 
make -j ${PARALLELCFG}
echo --- END make `date`
echo --- START ctest `date` 
# To workaround any concurrency problems put -DN_CONCURRENT_TESTS=1
ctest ${PARALLELCFG} ${TESTCFG} ${THETESTS}
echo --- END ctest `date`
else
echo --- START ctest `date` 
# To workaround any concurrency problems put -DN_CONCURRENT_TESTS=1
echo ctest ${PARALLELCFG} ${CMCFG} ${TESTCFG} -DQMC_OPTIONS=${QMC_OPTIONS} ${THETESTS} -S $PWD/../qmcpack/CMake/ctest_script.cmake,$ctestscriptarg
     ctest ${PARALLELCFG} ${CMCFG} ${TESTCFG} -DQMC_OPTIONS=${QMC_OPTIONS} ${THETESTS} -S $PWD/../qmcpack/CMake/ctest_script.cmake,$ctestscriptarg
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
