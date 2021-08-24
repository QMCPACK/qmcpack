#!/bin/bash
echo "Loading QMCPACK dependency modules for summit"
module load gcc/9.3.0
module load spectrum-mpi
module load essl
module load netlib-lapack
module load netlib-scalapack
module load hdf5
module load fftw
module load cmake
module load boost
module load cuda
module load python/3.8-anaconda3
