FROM ubuntu:20.04
LABEL maintainer="williamfgc@yahoo.com"

RUN export DEBIAN_FRONTEND=noninteractive &&\
    apt-get clean &&\
    apt-get update -y &&\
    apt-get upgrade -y apt-utils &&\
    apt-get install -y gpg wget

# Dependencies
RUN wget https://apt.kitware.com/kitware-archive.sh &&\
    sh kitware-archive.sh

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
    numdiff \
    -y

# Python packages for tests
RUN export DEBIAN_FRONTEND=noninteractive &&\
    apt-get install python3-numpy \
    python3-h5py \
    python3-pandas \
    python3-pip \
    -y

RUN export DEBIAN_FRONTEND=noninteractive &&\
    pip3 install cif2cell

# must add a user different from root 
# to run MPI executables
RUN useradd -ms /bin/bash user
# allow in sudoers to install packages
RUN adduser user sudo
RUN echo "user:user" | chpasswd

USER user
WORKDIR /home/user
