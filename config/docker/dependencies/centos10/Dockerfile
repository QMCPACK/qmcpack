FROM quay.io/centos/centos:10

RUN dnf install -y epel-release
RUN /usr/bin/crb enable
RUN dnf update -y
RUN dnf install -y make cmake gcc gcc-c++ gcc-gfortran openmpi-devel fftw-devel
RUN dnf install -y boost-devel libxml2-devel
RUN dnf install -y blas-devel lapack-devel
RUN dnf install -y hdf5-devel
RUN dnf install -y ninja-build
RUN dnf install -y rsync # Nexus hard requirement
RUN dnf install -y python3-numpy
RUN dnf install -y python3-h5py

## Add a user different from root for OpenMPI
RUN adduser -u 1000 -d /home/user -s /bin/bash user
USER user
WORKDIR /home/user

# To build QMCPACK:
#
# module load mpi/openmpi-x86_64
# cd /work/build
# cmake -GNinja -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpiCC -DQMC_DATA=/data ..
