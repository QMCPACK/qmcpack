.. _running_docker:

Running QMCPACK on Docker Containers
====================================

This guide will briefly cover running QMCPACK on your machine using a docker container. Docker containers are a portable way to achieve reproducibility across all developers and users. Two main uses:

1. Debugging CI related issues running on Docker containers that are not reproducible in other environments.
2. Ease the learning curve by providing a ready-to-go environment in a single container (binary available in DockerHub) to new community members. 

Current Images
--------------

Docker containers are identified by `domain/image:tag` and stored using `DockerHub <https://hub.docker.com/>`_.
Currently available containers have pre-installed QMCPACK dependencies, see the Dockerfile file link for available dependencies on each image:

- `Linux containers <https://hub.docker.com/r/williamfgc/qmcpack-ci/tags>`_ 
   - williamfgc/qmcpack-ci:ubuntu20-openmpi: `Dockerfile <https://github.com/QMCPACK/qmcpack/blob/develop/config/docker/dependencies/ubuntu/openmpi/Dockerfile>`_


Running Docker Containers
-------------------------

1. **Install the Docker engine**: install the latest version of the `Docker <https://www.docker.com/get-started>`_ engine for your system. Please see the documentation for different Linux distros `here <https://docs.docker.com/engine/install/#server>`_. 

   After installation run the following command to verify the Docker engine is properly installed. **Note**: restart your system if necessary. 

   .. code-block:: bash
   
      docker run hello-world

2. **Pull an image** (optional, see 3): once Docker is properly installed and running on your system, use the following command to download a QMCPACK image and tag:

   .. code-block:: bash
   
      docker pull williamfgc/qmcpack-ci:ubuntu20-openmpi

3. **Run an image**: the `docker run` command will spin up a container with using the image we just downloaded from step 2. Alternatively, `docker run` will automatically fallback to pulling the image and tag from DockerHub (requires connection).

   For a quick and safe, non `sudo`, run:   

   .. code-block:: bash

      docker run -it williamfgc/qmcpack-ci:ubuntu20-openmpi /bin/bash

   The above will run the container in interactive mode dropping the default `user` to `/home/user` using the `bash` shell.

   **Run an image (for Development)** `docker run` has a few extra options that can be used to run QMCPACK: 

   .. code-block:: bash
    
      docker run -u root --env USER=`stat -c "%u" .`  -v <QMCPACK Source Directory>:/home/user -it williamfgc/qmcpack-ci:ubuntu20-openmpi /bin/bash


   Flags used by `docker run` (Note: The flags -i and -t are combined above):
    
    `-u` : For building we need to run as the root user so that docker has write permissions for the build (e.g. install additional packages, allocating shared volume permissions, ect.).

    `-v` : Replace `<QMCPACK Source Directory>` with the direct path to your QMCPACK directory, this maps it to our landing directory and gives docker access to the files

    `-i` : Specifies the image to use

    `-t` : Allocate a pseudo-tty, allows an instance of bash to pass commands to it
	
    `--env` : Adds an environment variable containing the current user's uid for use later

   As an example, if extra permissions are needed the container can be run with the `sudo` user (not recommended):

   .. code-block:: bash

      docker run -u root --env USER=`stat -c "%u" .`  -v path/to/QMCPACK:home/user -it williamfgc/qmcpack-ci:ubuntu20-openmpi /bin/bash


Build QMCPACK on Docker
-----------------------

The following steps just follow a regular QMCPACK build on any Linux environment

1. **Get QMCPACK**: use `https` as `ssh` requires extra authentication  

* Option 1 (fresh build):

   .. code-block:: bash

      git clone https://github.com/QMCPACK/qmcpack.git
      cd build

* Option 2 (for development):

    .. code-block:: bash

       sudo useradd -u $USER local && sudo -u local bash && exit
       cd build

    * Note: this gives the non-root docker `user` permissions to read/write to the directory supplied with the `-v` flag of `docker run` in the previous step, and then opens bash as non root `user` 


2. **Configure**:

   .. code-block:: bash

		  cmake -GNinja \
		   -DCMAKE_BUILD_TYPE=RelWithDebInfo \
		   -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
		   -DQMC_COMPLEX=0 \
		   ..

* Note: To reproduce the build in the Docker container used by GitHub Actions CI pipeline we provide an optimized build with debug symbols `-DCMAKE_BUILD_TYPE=RelWithDebInfo` , but users can select any other cmake build type(`Release` being default): 
            
            - `Debug`
            - `Release` 
            - `RelWithDebInfo`    

3. **Build**:

   .. code-block:: bash
    
      ninja

3. **Test**:

   .. code-block:: bash

      ctest -VV -R deterministic-unit_test_wavefunction_trialwf
      ctest -L deterministic


.. caution::

   OpenMPI strongly advises against running as a `root` user, see `docs <https://www.open-mpi.org/doc/v3.1/man1/mpirun.1.php#sect22>`_