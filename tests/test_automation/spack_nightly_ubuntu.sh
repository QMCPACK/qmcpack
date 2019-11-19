#!/bin/bash
# Uncomment below for VERY verbose output from BASH
# set -x 

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
# Contains conflicts and failures, anytime we don't get an exit 0
export SPACK_FAILS=$HOME/qmcpack_spack_failures.out

# Here are all the variants that we test, plus dependencies
# declare -a versions=("3.8.0" "3.7.0")
declare -a versions=("3.8.0")
declare -a compilers=("gcc@8.3.0" "intel@19.0.3.199" "pgi@19.7" "clang@9.0.0")
# declare -a compilers=("intel@19.0.3.199")
declare -a withqe=("+qe" "~qe")
declare -a withmpi=("+mpi" "~mpi")
declare -a spotypes=("+complex" "~complex")
declare -a withtimers=("+timers" "~timers")
declare -a withmixed=("+mixed" "~mixed")
declare -a withsoa=("+soa" "~soa")
declare -a withcuda=("+cuda" "~cuda")
declare -a blasproviders=("netlib-lapack" "intel-mkl" "openblas")
# cuda_arch value explicitly set to value of GPU card on naromero-desktop.cels.anl.gov
gpu_card=61
cuda_version=10.1.243
# hash for QE 6.4.1 patch
qe_hash=57cb1b06ee2653a87c3acc0dd4f09032fcf6ce6b8cbb9677ae9ceeb6a78f85e2

# Checkout develop and pull new version.
# Check return status along the way to make sure things are working.

cd $SPACK_ROOT
git checkout develop 
if [ $? -ne 0 ]
then
    echo "Not a git repository"
    exit 1
fi

git pull upstream develop
if [ $? -ne 0 ]
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
build_with_gcc="^boost%gcc ^pkgconf%gcc ^perl%gcc ^libpciaccess%gcc ^cmake%gcc ^findutils%gcc ^m4%gcc"
build_with_gcc_nompi="^boost%gcc ^pkgconf%gcc ^perl%gcc ^cmake%gcc"

echo "##### Compiling QMCPACK variants without QE ######"
echo "##### 36 variants per compiler #####"
for version in ${versions[@]}; do
    for compiler in ${compilers[@]}; do
	for spotype in ${spotypes[@]}; do
	    for mixed in ${withmixed[@]}; do
		for soa in ${withsoa[@]}; do
		    for blas in ${blasproviders[@]}; do
			echo "####################"
			variant1="qmcpack~qe+mpi+timers~cuda${spotype}${mixed}${soa}"@"${version}"%"${compiler} ^${blas} ^mpich ${build_with_gcc}"
			echo $variant1
			spack install $variant1
			if [ $? -ne 0 ]
			then
			    echo $variant1 >> ${SPACK_FAILS}
			fi

			echo "####################"
			variant2="qmcpack~qe+mpi+timers+cuda${spotype}${mixed}${soa}"@"${version}"%"${compiler} cuda_arch=${gpu_card} ^cuda"@"${cuda_version} ^${blas} ^mpich ${build_with_gcc}"
			echo $variant2
			spack install $variant2
			if [ $? -ne 0 ]
			then
			    echo $variant2 >> ${SPACK_FAILS}
			fi
			
			echo "####################"
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
            echo "###############"
	    variant3="qmcpack+qe+mpi~timers~cuda~complex~mixed~soa"@"${version}"%"${compiler} ^${blas} ^mpich ${build_with_gcc}"
	    echo $variant3
	    spack install $variant3
	    if [ $? -ne 0 ]
	    then
		echo $variant3 >> ${SPACK_FAILS}
	    fi
	    
	    echo "###############"
	    variant4="qmcpack+qe~mpi~phdf5~timers~cuda~complex~mixed~soa"@"${version}"%"${compiler} ^${blas}"
	    echo $variant4
	    spack install $variant4
	    if [ $? -ne 0 ]
	    then
		echo $variant4 >> ${SPACK_FAILS}
	    fi	    

	    # test that QMCPACK patch was REALLY applied
	    spack find quantum-espresso@6.4.1"%"${compiler} patches=${qe_hash}
	    if $? -ne 0
	    then
		echo "QMCPACK patch was not applied to QE." >> ${SPACK_FAILS}
	    fi
       done
   done
done
