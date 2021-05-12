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
                    -y

# add the latest clang development
RUN bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)"           
RUN apt-get install clang-tools-12

# must add a user different from root 
# to run MPI executables
RUN useradd -ms /bin/bash user
# allow in sudoers to install packages
RUN adduser user sudo
RUN echo "user:user" | chpasswd

USER user
WORKDIR /home/user
