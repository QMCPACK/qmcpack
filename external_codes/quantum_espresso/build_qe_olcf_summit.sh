#!/bin/bash

module purge
echo "Purging current module set"
. ../../config/load_summit_modules.sh

./download_and_patch_qe6.4.1.sh
cd qe-6.4.1

CXX=g++ CC=gcc FC=gfortran MPIF90=mpif90 F90=gfortran F77=gfortran ./configure --with-hdf5=$OLCF_HDF5_ROOT

make all
