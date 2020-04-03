#!/bin/bash
# Next line is specific to the Argonne CELS computing environment
# it is needed in order to set-up modules and load the compilers
# correctly, particular Intel and PGI.
. /etc/profile.d/z00_lmod.sh

# Uncomment below for VERY verbose output from BASH
set -x 

# Script is run as a nightly cronjob on naromero-desktop.cels.anl.gov.
# chronic from moreutils is needed to silence the script when everything
# runs correctly.
# https://packages.debian.org/unstable/utils/moreutils
# Example of the cronjob. It runs at 5 pm (local time)
# MAILTO=naromero@anl.gov
# 00 17 * * * chronic /home/naromero/spack_nightly_ubuntu.sh

# Spack Environment
export SPACK_ROOT=/nfs/gce/projects/naromero-workspace/spack
source $SPACK_ROOT/share/spack/setup-env.sh

# Spack Failures
# Contains failures, anytime we don't get an exit 0 from `spack install`
export SPACK_FAILS=$HOME/qmcpack_spack_failures.out

# Spack Conflicts
# Contains conflicts, anytime we don't get an exit 0 from `spack spec'
export SPACK_CONFLICTS=$HOME/qmcpack_spack_conflicts.out

# Test Variant
function test_variant {
    echo "####################"
    variant="$@"
    spack spec $variant
    RET_VAL=$?
    if [[ ${RET_VAL} -ne 0 ]]
    then
        known_conflict=1
	echo ${RET_VAL} ${known_conflict} $variant >> ${SPACK_CONFLICTS}

    else
        known_conflict=0
    fi

    if [[ ${known_conflict} -eq 0 ]]
    then
	echo "### Installing ###"
        spack install $variant
	RET_VAL=$?
        if [[ ${RET_VAL} -ne 0 ]]
        then
	    echo ${RET_VAL} ${known_conflict} $variant >> ${SPACK_FAILS}
        fi
    fi
}

# Here are all the variants that we test, plus dependencies
# If you are testing with CUDA, you need to verify that your
# CUDA version is compatible with your compiler
# https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html

# declare -a versions=("3.8.0" "3.7.0")
declare -a versions=("3.8.0")
declare -a compilers=("gcc@7.4.0" "intel@19.0.3.199" "pgi@19.7" "clang@9.0.0")
# declare -a compilers=("intel@19.0.3.199")
declare -a withqe=("+qe" "~qe")
declare -a withmpi=("+mpi" "~mpi")
declare -a spotypes=("+complex" "~complex")
declare -a withtimers=("+timers" "~timers")
declare -a withmixed=("+mixed" "~mixed")
declare -a withsoa=("+soa" "~soa")
declare -a withcuda=("+cuda" "~cuda")
declare -a blasproviders=("netlib-lapack" "intel-mkl" "openblas")
# declare -a blasproviders=("intel-mkl")
# cuda_arch value explicitly set to value of GPU card on naromero-desktop.cels.anl.gov
gpu_card=61
cuda_version=10.2.89
# hash for QE 6.4.1 patch
qe_hash=57cb1b06ee2653a87c3acc0dd4f09032fcf6ce6b8cbb9677ae9ceeb6a78f85e2

# Checkout develop and pull new version.
# Check return status along the way to make sure things are working.

cd $SPACK_ROOT
git checkout develop 
if [[ $? -ne 0 ]]
then
    echo "Not a git repository"
    exit 1
fi

git pull upstream develop
if [[ $? -ne 0 ]]
then
    echo "Not able to pull Spack upstream into local git repo"
    exit 1
fi


# Start from a clean slate each time. This is super slow
# but ensures that nothing breaks from a modification to
# the Spack core. This can take a while. 
# Note the 'spack uninstall -ay <spec>' is dangerous
# if there is <spec> cannot be found, everything is 
# uninstalled
spack uninstall -ay qmcpack
spack uninstall -ay quantum-espresso
spack clean

# File for logging conflicts and build failures
rm -f ${SPACK_FAILS}
touch ${SPACK_FAILS}
rm -f ${SPACK_CONFLICTS}
touch ${SPACK_CONFLICTS}

# test QMCPACK variants
# NOTES:
# - openmpi has an incompatibility with QE
#
# - Intel cannot compile m4
#
# - PGI cannot compile boost, pkgconf, perl, libpciaccess, cmake, findutils
#
# - Union of packages that should get compiled with GCC are: ^boost%gcc
#   ^pkgconf%gcc ^perl%gcc ^libpciaccess%gcc ^numactl%gcc ^cmake%gcc ^findutils%gcc ^m4%gcc
#
build_with_gcc='^boost%gcc@7.4.0 ^pkgconf%gcc@7.4.0 ^perl%gcc@7.4.0 ^libpciaccess%gcc@7.4.0 ^cmake%gcc@7.4.0 ^findutils%gcc@7.4.0 ^m4%gcc@7.4.0'
build_with_gcc_nompi='^boost%gcc@7.4.0 ^pkgconf%gcc@7.4.0 ^perl%gcc@7.4.0 ^cmake%gcc@7.4.0'

echo "##### Compiling QMCPACK variants without QE ######"
echo "##### 36 variants per compiler #####"
for version in ${versions[@]}; do
    for compiler in ${compilers[@]}; do
	for spotype in ${spotypes[@]}; do
	    for mixed in ${withmixed[@]}; do
		for soa in ${withsoa[@]}; do
		    for blas in ${blasproviders[@]}; do
			variant1='qmcpack~qe+mpi+timers~cuda'${spotype}${mixed}${soa}'@'${version}'%'${compiler}' ^'${blas}' ^mpich '${build_with_gcc}
			test_variant "${variant1}"

			# cuda version takes an extra arg
			variant2='qmcpack~qe+mpi+timers+cuda'${spotype}${mixed}${soa}'@'${version}'%'${compiler}' cuda_arch='${gpu_card}' ^cuda@'${cuda_version}' ^'${blas}' ^mpich '${build_with_gcc}
			# variant2='qmcpack~qe+mpi+timers+cuda'${spotype}${mixed}${soa}'@'${version}'%'${compiler}' ^cuda@'${cuda_version}' ^'${blas}' ^mpich '${build_with_gcc}
			test_variant ${variant2}
		    done
		done
	    done
	done
    done
done

# test QE and FFT variants
echo "##### Compiling QMCPACK variants with QE ######"
echo "##### 6 variants per compiler #####" 
echo "##### Test that QE patch is applied ####"
for version in ${versions[@]}; do
    for compiler in ${compilers[@]}; do
	for blas in ${blasproviders[@]}; do
	    variant3='qmcpack+qe+mpi~timers~cuda~complex~mixed~soa@'${version}'%'${compiler}' ^'${blas}' ^mpich '${build_with_gcc}
	    test_variant $variant3
	    
	    variant4='qmcpack+qe~mpi~phdf5~timers~cuda~complex~mixed~soa@'${version}'%'${compiler}' ^'${blas}
	    test_variant $variant4

	    # test that QMCPACK patch was REALLY applied
	    spack find quantum-espresso@6.4.1"%"${compiler} patches=${qe_hash}
	    if [[ $? -ne 0 ]]
	    then
		echo "QMCPACK patch was not applied to QE." >> ${SPACK_FAILS}
	    fi
       done
   done
done
