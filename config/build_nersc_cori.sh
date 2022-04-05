#!/usr/bin/env bash

################################################################
## * This script builds available configurations of QMCPACK   ##
##   on NERSC Cori.                                           ##
##                                                            ##
## * Execute this script in trunk/                            ##
##   ./config/build_nersc_cori_hsw.sh                         ##
##                                                            ##
## Last verified: Apr 19, 2021                                ##
################################################################



export CRAYPE_LINK_TYPE=dynamic

module unload cray-libsci
module load boost/1.70.0
module load cray-hdf5-parallel
module load cmake
module load gcc/7.3.0 # Make C++ 14 standard library available to the Intel compiler
module list

# module files resulting from module imports above:
#Currently Loaded Modulefiles:
#  1) modules/3.2.11.4                                 10) dmapp/7.1.1-7.0.1.1_4.72__g38cf134.ari           19) craype-haswell
#  2) altd/2.0                                         11) gni-headers/5.0.12.0-7.0.1.1_6.46__g3b1768f.ari  20) cray-mpich/7.7.10
#  3) darshan/3.2.1                                    12) xpmem/2.2.20-7.0.1.1_4.28__g0475745.ari          21) craype-hugepages2M
#  4) craype-network-aries                             13) job/2.2.4-7.0.1.1_3.55__g36b56f4.ari             22) boost/1.70.0
#  5) intel/19.0.3.199                                 14) dvs/2.12_2.2.167-7.0.1.1_17.11__ge473d3a2        23) cray-hdf5-parallel/1.10.5.2
#  6) craype/2.6.2                                     15) alps/6.6.58-7.0.1.1_6.30__g437d88db.ari          24) cmake/3.21.3
#  7) udreg/2.3.2-7.0.1.1_3.61__g8175d3d.ari           16) rca/2.2.20-7.0.1.1_4.74__g8e3fb5b.ari            25) gcc/7.3.0
#  8) ugni/6.0.14.0-7.0.1.1_7.63__ge78e5b0.ari         17) atp/2.1.3
#  9) pmi/5.0.14                                       18) PrgEnv-intel/6.0.5


# Haswell CPU Complex
mkdir build_nersc_cori_hsw_cmplx
cd build_nersc_cori_hsw_cmplx
cmake -DQMC_SYMLINK_TEST_FILES=0 -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment -DQMC_COMPLEX=1 ..
nice make -j 8
cd ..
appendage=
if [ -f "build_nersc_cori_hsw_cmplx/bin/qmcpack_complex" ]; then
    appendage=_complex
fi
ls -l build_nersc_cori_hsw_cmplx/bin/qmcpack${appendage}
ln -sf ./build_nersc_cori_hsw_cmplx/bin/qmcpack${appendage} ./qmcpack_nersc_cori_cpu_hsw_comp


# Haswell CPU Real
mkdir build_nersc_cori_hsw
cd build_nersc_cori_hsw
cmake -DQMC_SYMLINK_TEST_FILES=0 -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment -DQMC_COMPLEX=0 ..
nice make -j 8
ls -l bin/qmcpack
cd ..
ln -sf ./build_nersc_cori_hsw/bin/qmcpack ./qmcpack_nersc_cori_cpu_hsw

# Swap modules for KNL recipe
module swap craype-haswell craype-mic-knl
module list

# module files resulting from module imports above:
#Currently Loaded Modulefiles:
#  1) modules/3.2.11.4                                  6) craype/2.6.2                                     11) gni-headers/5.0.12.0-7.0.1.1_6.43__g3b1768f.ari  16) rca/2.2.20-7.0.1.1_4.65__g8e3fb5b.ari            21) craype-hugepages2M
#  2) altd/2.0                                          7) udreg/2.3.2-7.0.1.1_3.52__g8175d3d.ari           12) xpmem/2.2.20-7.0.1.1_4.23__g0475745.ari          17) atp/2.1.3                                        22) boost/1.70.0
#  3) darshan/3.2.1                                     8) ugni/6.0.14.0-7.0.1.1_7.54__ge78e5b0.ari         13) job/2.2.4-7.0.1.1_3.50__g36b56f4.ari             18) PrgEnv-intel/6.0.5                               23) cray-hdf5-parallel/1.10.5.2
#  4) craype-network-aries                              9) pmi/5.0.14                                       14) dvs/2.12_2.2.167-7.0.1.1_17.6__ge473d3a2         19) craype-mic-knl                                   24) cmake/3.14.4
#  5) intel/19.0.3.199                                 10) dmapp/7.1.1-7.0.1.1_4.64__g38cf134.ari           15) alps/6.6.58-7.0.1.1_6.22__g437d88db.ari          20) cray-mpich/7.7.10                                25) gcc/7.3.0

# KNL CPU Complex
mkdir build_nersc_cori_knl_cmplx
cd build_nersc_cori_knl_cmplx
cmake -DQMC_SYMLINK_TEST_FILES=0 -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment -DQMC_COMPLEX=1 ..
nice make -j 8
cd ..
appendage=
if [ -f "build_nersc_cori_knl_cmplx/bin/qmcpack_complex" ]; then
    appendage=_complex
fi
ls -l build_nersc_cori_knl_cmplx/bin/qmcpack${appendage}
ln -sf ./build_nersc_cori_knl_cmplx/bin/qmcpack${appendage} ./qmcpack_nersc_cori_cpu_knl_comp

# KNL CPU Real
mkdir build_nersc_cori_knl
cd build_nersc_cori_knl
cmake -DQMC_SYMLINK_TEST_FILES=0 -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment -DQMC_COMPLEX=0 ..
nice make -j 8
ls -l bin/qmcpack
cd ..
ln -sf ./build_nersc_cori_knl/bin/qmcpack ./qmcpack_nersc_cori_cpu_knl
