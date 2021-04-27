.. _running_docker:

Running QMCPACK in Docker
=========================

This guide will briefly cover what it will take to quickly and effectively install QMCPack on your machine using a docker image.

Current Supported OS's
----------------------
* Linux

***************************
Getting Started with Docker
***************************

Firstly you will need to install the latest version of `Docker <https://www.docker.com/get-started>`_ on your machine

Once you have that installed and can run the following command without issue, we may proceed.


                    ``docker run hello-world``

Download an Appropriate Image
*****************************

QMCPack will be supporting a few flavors of linux, and will have a list of available images `here <http://>`_


For this tutorial we will be using one for `ubuntu 20.04` `williamfgc/qmcpack-ci:ubuntu20-openmpi <https://hub.docker.com/r/williamfgc/qmcpack-ci/tags?page=1&ordering=last_updated>`_

Use the following command to download your selected image:

    ``docker pull <your image name here>``

So I will run the command [`docker pull williamfgc/qmcpack-ci:ubuntu20-openmpi`]

SSH into a Docker Container
***************************

Run the following command and docker will spin up a container with using the image we just downloaded, giving us console access:

    ``docker run -u root -v <QMCPack Source Directory>:/home/user -it williamfgc/qmcpack-ci:ubuntu20-openmpi /bin/bash``

To explain the flags used(Note: The flags -i and -t are combined above):
    `-u` : For building we need to run as the root user so that docker has write permissions for the build

    `-v` : Replace `<QMCPack Source Directory>` with the direct path to your QMCPack directory, this maps it to our landing directory and gives docker access to the files

    `-i` : Specifies the image to use

    `-t` : Allocate a pseudo-tty, allows an instance of bash to pass commands to it

Build and Run Tests
*******************

***work in progress***

Run the following commands to build QMCPack, this may take a few minutes total::

    CC=clang CXX=clang++ cmake -GNinja \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx \
    -DENABLE_SANITIZER=ubsan \
    -DQMC_COMPLEX=1 \
    .

Then::

    ninja

and finally run some tests to confirm everything is in order::

    ctest -VV -R deterministic-unit_test_wavefunction_trialwf
