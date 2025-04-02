#!/bin/bash
echo "Loading QMCPACK dependency modules for andes"
echo "https://docs.olcf.ornl.gov/systems/andes_user_guide.html"
echo

module load gcc/10.3.0
module load openmpi/4.1.2
module load openblas/0.3.19
module load netlib-lapack
module load hdf5
module load fftw
export FFTW_ROOT=$OLCF_FFTW_ROOT
module load cmake/3.24.4
module load boost/1.78.0
module load python/3.7-anaconda3
