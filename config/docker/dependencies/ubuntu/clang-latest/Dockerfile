FROM ubuntu:20.04
MAINTAINER William F Godoy williamfgc@yahoo.com

RUN export DEBIAN_FRONTEND=noninteractive &&\
    apt-get update -y &&\
    apt-get upgrade -y apt-utils

# Dependencies
RUN export DEBIAN_FRONTEND=noninteractive &&\
    apt-get install gcc g++ \ 
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

# add the latest clang development
RUN wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key| apt-key add - &&\
    apt-add-repository 'deb http://apt.llvm.org/focal/ llvm-toolchain-focal-12 main'
RUN apt-get update -y &&\     
    apt-get install clang-12 clang-tools-12 -y

# must add a user different from root 
# to run MPI executables
RUN useradd -ms /bin/bash user
# allow in sudoers to install packages
RUN adduser user sudo
RUN echo "user:user" | chpasswd

USER user
WORKDIR /home/user
