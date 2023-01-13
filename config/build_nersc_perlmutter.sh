#!/usr/bin/env bash

################################################################
## * This script builds CPU versions of QMCPACK               ##
##   on NERSC Perlmutter.                                     ##
##                                                            ##
## * Execute this script from qmcpack directory as:           ##
##   ./config/build_nersc_perlmutter.sh                       ##
##                                                            ##
## * Last verified: 2023/01/13                                ##
################################################################

# Make LibXml2 available to QMCPACK on PATH
module load e4s
spack install libxml2
spack load libxml2
# Load the required modules
module load cpu
module load cray-hdf5-parallel
module load cray-fftw
module list

# On the above date, the loaded modules were:
#Currently Loaded Modules:
#  1) craype-x86-milan     4) xpmem/2.5.2-2.4_3.18__gd0f7936.shasta   7) cray-libsci/22.11.1.2  10) gcc/11.2.0              13) xalt/2.10.2    16) cpu/1.0
#  2) libfabric/1.15.2.0   5) PrgEnv-gnu/8.3.3                        8) cray-mpich/8.1.22      11) perftools-base/22.09.0  14) darshan/3.4.0  17) cray-hdf5-parallel/1.12.2.1
#  3) craype-network-ofi   6) cray-dsmml/0.2.2                        9) craype/2.7.19          12) cpe/22.11               15) e4s/22.05      18) cray-fftw/3.3.10.2

mkdir build_cpu_real
cd build_cpu_real
cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment -DQMC_COMPLEX=0 ..
nice make -j 16
ls -l bin/qmcpack
cd ..

mkdir build_cpu_cplx
cd build_cpu_cplx
cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment -DQMC_COMPLEX=1 ..
nice make -j 16
ls -l bin/qmcpack_complex
cd ..
