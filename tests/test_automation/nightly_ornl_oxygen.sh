#!/bin/bash

place=/scratch/pk7/QMCPACK_CI_BUILDS_DO_NOT_REMOVE

if [ -e /scratch/pk7 ]; then

if [ ! -e $place ]; then
mkdir $place
fi

if [ -e $place ]; then


for sys in build_gcc build_intel2016 build_intel2015 build_gcc_complex build_intel2016_complex build_intel2015_complex build_gcc_cuda build_intel2015_cuda
do

cd $place

if [ -e $sys ]; then
rm -r -f $sys
fi
mkdir $sys
cd $sys

echo --- Checkout for $sys `date`
svn checkout https://svn.qmcpack.org/svn/trunk
#svn checkout https://subversion.assembla.com/svn/qmcdev/trunk

if [ -e trunk/CMakeLists.txt ]; then
cd trunk
mkdir $sys
cd $sys
echo --- Building for $sys `date`

export PATH=/opt/local/bin:/opt/local/sbin:/usr/local/cuda/bin/:/usr/lib64/openmpi/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin

case $sys in
"build_gcc")
    module() { eval `/usr/bin/modulecmd sh $*`; }
    module load mpi
    export QMCPACK_TEST_SUBMIT_NAME=GCC-Release
    ctest -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -S $PWD/../CMake/ctest_script.cmake,release -R 'short|converter' -VV
    module unload mpi
    ;;
"build_intel2016")
    ORIGPATH=$PATH
    source /opt/intel2016/bin/compilervars.sh intel64
    source /opt/intel2016/impi/5.1.1.109/bin64/mpivars.sh
    export QMCPACK_TEST_SUBMIT_NAME=Intel2016-Release
    ctest -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -S $PWD/../CMake/ctest_script.cmake,release -R 'short|converter' -VV
    export PATH=$ORIGPATH
    ;;
"build_intel2015")
    ORIGPATH=$PATH
    source /opt/intel/bin/compilervars.sh intel64
    source /opt/intel/impi_latest/bin64/mpivars.sh
    export QMCPACK_TEST_SUBMIT_NAME=Intel2015-Release
    ctest -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -S $PWD/../CMake/ctest_script.cmake,release -R 'short|converter' -VV
    ;;
"build_gcc_complex")
    module() { eval `/usr/bin/modulecmd sh $*`; }
    module load mpi
    export QMCPACK_TEST_SUBMIT_NAME=GCC-Complex-Release
    ctest -DQMC_COMPLEX=1 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -S $PWD/../CMake/ctest_script.cmake,release -R 'short|converter' -VV
    module unload mpi
    ;;
"build_intel2016_complex")
    ORIGPATH=$PATH
    source /opt/intel2016/bin/compilervars.sh intel64
    source /opt/intel2016/impi/5.1.1.109/bin64/mpivars.sh
    export QMCPACK_TEST_SUBMIT_NAME=Intel2016-Complex-Release
    ctest -DQMC_COMPLEX=1 -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -S $PWD/../CMake/ctest_script.cmake,release -R 'short|converter' -VV
    export PATH=$ORIGPATH
    ;;
"build_intel2015_complex")
    ORIGPATH=$PATH
    source /opt/intel/bin/compilervars.sh intel64
    source /opt/intel/impi_latest/bin64/mpivars.sh
    export QMCPACK_TEST_SUBMIT_NAME=Intel2015-Complex-Release
    ctest -DQMC_COMPLEX=1 -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -S $PWD/../CMake/ctest_script.cmake,release -R 'short|converter' -VV
    export PATH=$ORIGPATH
    ;;
"build_gcc_cuda")
    module() { eval `/usr/bin/modulecmd sh $*`; }
    module load mpi
    export QMCPACK_TEST_SUBMIT_NAME=GCC-CUDA-Release
    ctest -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_CUDA=1 -S $PWD/../CMake/ctest_script.cmake,release -R 'short|converter' -VV
    module unload mpi
    ;;
"build_intel2015_cuda")
    ORIGPATH=$PATH
    source /opt/intel/bin/compilervars.sh intel64
    source /opt/intel/impi_latest/bin64/mpivars.sh
    export QMCPACK_TEST_SUBMIT_NAME=Intel2015-CUDA-Release
    ctest -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DQMC_CUDA=1 -S $PWD/../CMake/ctest_script.cmake,release -R 'short|converter' -VV
    export PATH=$ORIGPATH
    ;;
*)
    echo "ERROR: Unknown build type $sys"
    ;;
esac

else
echo  "ERROR: No CMakeLists. Bad svn checkout."
exit 1
fi

done

else
echo "ERROR: No directory $place"
exit 1
fi

fi
