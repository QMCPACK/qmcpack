FROM ubuntu:impish-20210722
LABEL maintainer="williamfgc@yahoo.com"

RUN export DEBIAN_FRONTEND=noninteractive &&\
    apt-get update -y &&\
    apt-get upgrade -y apt-utils

# Dependencies
RUN export DEBIAN_FRONTEND=noninteractive &&\
    apt-get install gcc g++ \ 
    clang \
    clang-format \
    clang-tidy \
    gcovr \
    python3 \
    cmake \
    ninja-build \
    libboost-all-dev \
    git \
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
    -y

# Python packages for tests
RUN export DEBIAN_FRONTEND=noninteractive &&\
    apt-get install python3-numpy \
    python3-h5py \
    python3-pandas \
    -y

# must add a user different from root 
# to run MPI executables
RUN useradd -ms /bin/bash user
# allow in sudoers to install packages
RUN adduser user sudo
RUN echo "user:user" | chpasswd

USER user
WORKDIR /home/user
