#!/bin/bash
echo "Loading QMCPACK dependency modules for archer2"
echo "https://docs.archer2.ac.uk/"
echo
module restore
module load PrgEnv-gnu
module load cray-hdf5-parallel
module load cray-fftw
export FFTW_ROOT=$FFTW_DIR/..
module load libxml2
module load cmake
module load boost
module load cray-python
module load craype-network-ucx
module load cray-mpich-ucx

