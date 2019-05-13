#!/bin/bash
# Begin LSF directives
#BSUB -P MAT189
#BSUB -J test
#BSUB -o tst.o%J
#BSUB -W 120
#BSUB -nnodes 1
#BSUB -alloc_flags smt1
# End LSF directives and begin shell commands

#!/bin/bash
#
# Setup for summit.olcf.ornl.gov
#
# Run the "short" nightlies
# 

export TEST_SITE_NAME=summit.olcf.ornl.gov
export N_PROCS_BUILD=32
export N_CONCURRENT_TESTS=1
export CC=mpicc
export CXX=mpicxx
export FFTW_HOME=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-6.4.0/fftw-3.3.8-5gcj2ic4el7acu3rqnfnh735jz2ez7j5
export BOOST_ROOT=/ccs/home/yeluo/opt/boost_1_69_0

#Must be an absolute path
place=/gpfs/alpine/mat189/scratch/yeluo/QMCPACK_CI_BUILDS_DO_NOT_REMOVE

#QE_BIN=/sandbox/opt/qe-stable/qe-6.3/bin
QMC_DATA=$place/h5data

#define and load compiler
compiler=gcc6.4

if [ ! -e $place ]; then
mkdir $place
fi

if [ -e $place ]; then
cd $place

echo --- Hostname --- $HOSTNAME
echo --- Checkout for $sys `date`

branch=develop
entry=qmcpack-${branch}

if [ ! -e $entry ]; then
echo --- Cloning QMCPACK git `date`
git clone --depth 1 https://github.com/QMCPACK/qmcpack.git $entry
else
echo --- Updating local QMCPACK git `date`
cd $entry
git pull
cd ..
fi

if [ -e $entry/CMakeLists.txt ]; then
cd $entry

git checkout $branch

for sys in Real-SoA #Real-Mixed-SoA Complex-SoA Complex-Mixed-SoA
do

folder=build_${compiler}_$sys

if [ -e $folder ]; then
rm -r $folder
fi
mkdir $folder
cd $folder

echo --- Building for $sys `date`

# create log file folder if not exist
mydate=`date +%y_%m_%d`
if [ ! -e $place/log/$entry/$mydate ];
then
  mkdir -p $place/log/$entry/$mydate
fi

CTEST_FLAGS="-DQMC_DATA=$QMC_DATA -DBLAS_essl_LIBRARY=$OLCF_ESSL_ROOT/lib64/libessl.so -DENABLE_TIMERS=1 -DQMC_OPTIONS=-DMPIEXEC=sh;-DMPIEXEC_NUMPROC_FLAG=$place/$entry/tests/scripts/jsrunhelper.sh"

if [[ $sys == *"Complex"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -D QMC_COMPLEX=1"
fi

if [[ $sys == *"-SoA"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -D ENABLE_SOA=1"
fi

if [[ $sys == *"-Mixed"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -D QMC_MIXED_PRECISION=1"
fi

export QMCPACK_TEST_SUBMIT_NAME=${compiler}-${sys}-Release

ctest $CTEST_FLAGS -S $PWD/../CMake/ctest_script.cmake,release --stop-time `date --date=now+117mins +%H:%M:%S` -VV -E 'long' --timeout 800 &> $place/log/$entry/$mydate/${QMCPACK_TEST_SUBMIT_NAME}.log

cd ..
echo --- Finished $sys `date`
done

else
echo  "ERROR: No CMakeLists. Bad git clone."
exit 1
fi

else
echo "ERROR: No directory $place"
exit 1
fi
