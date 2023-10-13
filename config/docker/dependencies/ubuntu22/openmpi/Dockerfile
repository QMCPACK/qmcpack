FROM ubuntu:22.04

RUN export DEBIAN_FRONTEND=noninteractive &&\
    apt-get clean &&\
    apt-get update -y &&\
    apt-get upgrade -y apt-utils &&\
    apt-get install -y gpg wget

# Dependencies
RUN wget https://apt.kitware.com/kitware-archive.sh &&\
    sh kitware-archive.sh

RUN export DEBIAN_FRONTEND=noninteractive &&\
    apt-get install gcc-9 g++-9 \ 
    clang-14 \
    clang-format-14 \
    clang-tidy-14 \
    libomp-14-dev \
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

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 100 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 100

# add clang-14 as clang
RUN update-alternatives --install /usr/bin/clang clang /usr/bin/clang-14 100 && \
    update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-14 100

# add clang-format and clang-tidy as well as libomp
RUN update-alternatives --install /usr/bin/clang-format clang-format /usr/bin/clang-format-14 100 && \
    update-alternatives --install /usr/bin/clang-tidy clang-tidy /usr/bin/clang-tidy-14 100 && \
    update-alternatives --install /usr/bin/clang-tidy-diff.py clang-tidy-diff.py /usr/bin/clang-tidy-diff-14.py 100

# must add a user different from root 
# to run MPI executables
RUN useradd -ms /bin/bash user
# allow in sudoers to install packages
RUN adduser user sudo
RUN echo "user:user" | chpasswd

USER user
WORKDIR /home/user
