#!/bin/bash
echo "Loading QMCPACK dependency modules for archer2"
echo "https://docs.archer2.ac.uk/"
echo
#module swap PrgEnv-cray PrgEnv-gnu
module purge
module load PrgEnv-gnu
module load craype-network-ucx
module load cray-mpich-ucx
#module load openmpi/4.0.4
#module load openblas/0.3.12
#module load netlib-lapack
#module load netlib-scalapack
#module load cray-hdf5
module load cray-hdf5-parallel
module load cray-fftw
export FFTW_ROOT=$FFTW_DIR/..
#export FFTW_HOME=$FFTW_DIR/..
module load libxml2
module load cmake
module load boost
module load cray-python

