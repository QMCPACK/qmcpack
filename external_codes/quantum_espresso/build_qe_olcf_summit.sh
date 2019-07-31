#!/bin/bash

module purge
echo "Purging current module set"
. ../../config/load_summit_modules.sh

./download_and_patch_qe6.4.1.sh
cd qe-6.4.1

export BLAS_LIBS="-L$OLCF_ESSL_ROOT/lib64 -lessl"
export LAPACK_LIBS="-L$OLCF_ESSL_ROOT/lib64 -lessl $OLCF_NETLIB_LAPACK_ROOT/lib64/liblapack.a"
export SCALAPACK_LIBS="-L$OLCF_NETLIB_SCALAPACK_ROOT/lib -lscalapack"

./configure --with-scalapack --with-hdf5=$OLCF_HDF5_ROOT CC=mpicc MPIF90=mpif90 F90=gfortran F77=gfortran

sed -i "/DFLAGS/s/__FFTW/__LINUX_ESSL/" make.inc

make all
