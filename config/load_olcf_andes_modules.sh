#!/bin/bash
echo "Loading QMCPACK dependency modules for andes"
echo "https://docs.olcf.ornl.gov/systems/andes_user_guide.html"
echo
module load gcc/9.3.0
#module load intel/19.0.3
module load openmpi/4.0.4
#module load essl
module load openblas/0.3.12
module load netlib-lapack
#module load netlib-scalapack
module load hdf5
module load fftw
export FFTW_ROOT=$OLCF_FFTW_ROOT
module load cmake/3.18.4
module load boost/1.74.0
#module load cuda
module load python/3.7-anaconda3

