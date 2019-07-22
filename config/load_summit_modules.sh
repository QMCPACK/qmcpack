#!/bin/bash
echo "Loading QMCPACK dependency modules for summit"
module load gcc
module load spectrum-mpi
module load essl
module load netlib-lapack
module load netlib-scalapack
module load hdf5
module load fftw
export FFTW_ROOT=$OLCF_FFTW_ROOT
module load cmake
module load boost
module load cuda
module load python/2.7.15-anaconda2-5.3.0

