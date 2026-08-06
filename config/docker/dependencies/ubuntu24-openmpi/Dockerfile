FROM ubuntu:24.04

RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get clean \
    && apt-get update \
    && apt-get upgrade -y apt-utils \
    && apt-get install -y gpg wget

RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
        gcc \
        g++ \
        gcc-9 \
        g++-9 \
        gcc-13 \
        g++-13 \
        clang \
        clang-format \
        clang-tidy \
        libomp-dev \
        cmake \
        curl \
        gcovr \
        git \
        hdf5-tools \
        libboost-all-dev \
        libfftw3-dev \
        libhdf5-openmpi-dev \
        libhdf5-serial-dev \
        libopenblas-openmp-dev \
        libopenmpi-dev \
        libxml2-dev \
        ninja-build \
        numdiff \
        python3 \
        rsync \
        software-properties-common \
        sudo \
        vim \
        wget

# Python packages for Nexus and testing
RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
        python3-cif2cell \
        python3-coverage \
        python3-h5py \
        python3-matplotlib \
        python3-numpy \
        python3-pandas \
        python3-pip \
        python3-pytest \
        python3-pytest-cov \
        python3-pytest-order \
        python3-scipy \
        python3-spglib

# Install gcovr from pip since the Ubuntu 24.04 apt version doesn't support files >9999 lines
# https://github.com/gcovr/gcovr/pull/883
RUN export DEBIAN_FRONTEND=noninteractive \
    && sudo pip3 install --break-system-packages gcovr==8.6

# must add a user different from root to run MPI executables
RUN useradd -ms /bin/bash user
# allow in sudoers to install packages
RUN adduser user sudo
RUN echo "user:user" | chpasswd

USER user
WORKDIR /home/user
