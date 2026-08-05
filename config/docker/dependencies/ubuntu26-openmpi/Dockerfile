FROM ubuntu:26.04

RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get clean \
    && apt-get update \
    && apt-get upgrade -y apt-utils \
    && apt-get install -y gpg wget

RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
        gcc \
        g++ \
        gcc-15 \
        g++-15 \
        gcc-16 \
        g++-16 \
        clang-22 \
        libomp-22-dev \
        liboffload-22-dev \
        clang-format-22 \
        clang-tidy-22 \
        gcovr \
        python3 \
        cmake \
        ninja-build \
        libboost-all-dev \
        git \
        libopenmpi-dev \
        libhdf5-openmpi-dev \
        libhdf5-serial-dev \
        hdf5-tools \
        libfftw3-dev \
        libopenblas-openmp-dev \
        libxml2-dev \
        sudo \
        curl \
        rsync \
        wget \
        software-properties-common \
        vim \
        numdiff

# Python packages for Nexus and testing
RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
        python3-numpy \
        python3-scipy \
        python3-matplotlib \
        python3-h5py \
        python3-pandas \
        python3-coverage \
        python3-pytest \
        python3-pytest-cov \
        python3-pytest-order \
        python3-cif2cell \
        python3-spglib \
        python3-pip

# Install gcovr from pip since the Ubuntu version doesn't support files >9999 lines
# https://github.com/gcovr/gcovr/pull/883
RUN export DEBIAN_FRONTEND=noninteractive &&\
    sudo pip3 install --break-system-packages gcovr==8.6

# must add a user different from root to run MPI executables
RUN useradd -ms /bin/bash user
# allow in sudoers to install packages
RUN adduser user sudo
RUN echo "user:user" | chpasswd

USER user
WORKDIR /home/user
