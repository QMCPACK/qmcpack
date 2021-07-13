.. _github_actions:

=================================
Github Actions, QMCPACK, and You!
=================================

This guide will cover the purpose and usual interactions a QMCPACK contributor might experience in regards to the Github Actions.  For how Github Actions accomplishes these things and other implementation details please refer to the offical `Github Actions Docs <https://docs.github.com/en/actions/guides>`_ and our scripts, which at the time of writing this doc are located `here <https://github.com/QMCPACK/qmcpack/tree/develop/tests/test_automation/github-actions/ci>`_.


So, What Even is Github Actions?
================================


Good question! Github Actions is an event driven automation tool that allows us to automatically execute commands in response to QMCPACK repo related actions.  For example, merging a branch into master might then trigger our test scripts to run.

Neat! How We are Using It and Why?
==================================

One of the biggest hows and whys is that it saves time! There are certain jobs we always trigger upon merge, and having a designated person trigger and manage them takes a lot of time off their hands.

Currently we are using Github Actions to automatically handle a few different jobs. These jobs are then either ran on the Github provided build vm's or pushed to our supplied hardware.  Usually the jobs are only ran on our hardware when they require gpu's to run.

All Github Runner jobs currently use the williamfgc/qmcpack-ci:ubuntu20-openmpi docker image, if you would like to reproduce theses tests exactly using docker, please refer to running_docker.rst for QMCPACK docker setup.

Github Runners
--------------

linux (gcc-openmpi-real-coverage)
"""""""""""""""""""""""""""""""""
+-----------+--------------------+
|Compiler   |GCC                 |
+-----------+--------------------+
|MPI Wrapper|OpenMPI             |
+-----------+--------------------+
|Objective  |Test Coverage (Real)|
+-----------+--------------------+
|Duration   |~45 Minutes         |
+-----------+--------------------+

+-----------+---------------------------------------------------------------------------------------------------------------------------------------------+
|Build Flags|`cmake -GNinja -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx \-DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_GCOV=TRUE \-DQMC_COMPLEX=0 \.`   |
+-----------+---------------------------------------------------------------------------------------------------------------------------------------------+

+------------+--------------------------------------------------------------------------------------------------------------------------------------------+
|Test Command|`export OMPI_MCA_rmaps_base_oversubscribe=1 && export OMPI_MCA_hwloc_base_binding_policy=none && ctest --output-on-failure -L deterministic`|
+------------+--------------------------------------------------------------------------------------------------------------------------------------------+

linux (gcc-openmpi-complex-coverage)
""""""""""""""""""""""""""""""""""""
+-----------+-----------------------+
|Compiler   |GCC                    |
+-----------+-----------------------+
|MPI Wrapper|OpenMPI                |
+-----------+-----------------------+
|Objective  |Test Coverage (Complex)|
+-----------+-----------------------+
|Duration   |~50 Minutes            |
+-----------+-----------------------+

+-----------+---------------------------------------------------------------------------------------------------------------------------------------------+
|Build Flags|`cmake -GNinja -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx \-DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_GCOV=TRUE \-DQMC_COMPLEX=1 \.`   |
+-----------+---------------------------------------------------------------------------------------------------------------------------------------------+

+------------+--------------------------------------------------------------------------------------------------------------------------------------------+
|Test Command|`export OMPI_MCA_rmaps_base_oversubscribe=1 && export OMPI_MCA_hwloc_base_binding_policy=none && ctest --output-on-failure -L deterministic`|
+------------+--------------------------------------------------------------------------------------------------------------------------------------------+

linux (clang-openmpi-real-asan)
"""""""""""""""""""""""""""""""
+-----------+------------------------+
|Compiler   |Clang                   |
+-----------+------------------------+
|MPI Wrapper|OpenMPI                 |
+-----------+------------------------+
|Objective  |Address Sanitizer (Real)|
+-----------+------------------------+
|Duration   |~25 Minutes             |
+-----------+------------------------+

+-----------+---------------------------------------------------------------------------------------------------------------------------------------------+
|Build Flags|`CC=clang CXX=clang++ \cmake -GNinja -DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_SANITIZER=asan \-DQMC_MPI=0 \-DQMC_COMPLEX=0 \.`             |
+-----------+---------------------------------------------------------------------------------------------------------------------------------------------+

+------------+------------------------------------------------------------------------------------------------------------------------------------------------+
|Test Command|`export OMPI_MCA_rmaps_base_oversubscribe=1 && export OMPI_MCA_hwloc_base_binding_policy=none && ctest --output-on-failure -L unit -LE noasan`  |
+------------+------------------------------------------------------------------------------------------------------------------------------------------------+

linux (clang-openmpi-complex-asan)
""""""""""""""""""""""""""""""""""
+-----------+---------------------------+
|Compiler   |Clang                      |
+-----------+---------------------------+
|MPI Wrapper|OpenMPI                    |
+-----------+---------------------------+
|Objective  |Address Sanitizer (Complex)|
+-----------+---------------------------+
|Duration   |~30 Minutes                |
+-----------+---------------------------+

+-----------+---------------------------------------------------------------------------------------------------------------------------------------------+
|Build Flags|`CC=clang CXX=clang++ \cmake -GNinja -DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_SANITIZER=asan \-DQMC_MPI=0 \-DQMC_COMPLEX=1 \.`             |
+-----------+---------------------------------------------------------------------------------------------------------------------------------------------+

+------------+------------------------------------------------------------------------------------------------------------------------------------------------+
|Test Command|`export OMPI_MCA_rmaps_base_oversubscribe=1 && export OMPI_MCA_hwloc_base_binding_policy=none && ctest --output-on-failure -L unit -LE noasan`  |
+------------+------------------------------------------------------------------------------------------------------------------------------------------------+

linux (clang-openmpi-real-ubsan)
""""""""""""""""""""""""""""""""
+-----------+-----------------------------------+
|Compiler   |Clang                              |
+-----------+-----------------------------------+
|MPI Wrapper|OpenMPI                            |
+-----------+-----------------------------------+
|Objective  |Undefined Behavior Sanitizer (Real)|
+-----------+-----------------------------------+
|Duration   |~55 Minutes                        |
+-----------+-----------------------------------+

+-----------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|Build Flags|`CC=clang CXX=clang++ \cmake -GNinja -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx \-DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_SANITIZER=ubsan \-DQMC_COMPLEX=0 \.`|
+-----------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+

+------------+--------------------------------------------------------------------------------------------------------------------------------------------+
|Test Command|`export OMPI_MCA_rmaps_base_oversubscribe=1 && export OMPI_MCA_hwloc_base_binding_policy=none && ctest --output-on-failure -L deterministic`|
+------------+--------------------------------------------------------------------------------------------------------------------------------------------+

linux (clang-latest-openmp-offload)
"""""""""""""""""""""""""""""""""""
+-----------+-----------------------------------------+
|Compiler   |Clang                                    |
+-----------+-----------------------------------------+
|MPI Wrapper|OpenMPI                                  |
+-----------+-----------------------------------------+
|Objective  |Build for GPU Acceleration (Experimental)|
+-----------+-----------------------------------------+
|Duration   |~35 Minutes                              |
+-----------+-----------------------------------------+

+-----------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|Build Flags|`cmake -GNinja -DCMAKE_C_COMPILER=clang-12 -DCMAKE_CXX_COMPILER=clang++-12 \-DENABLE_OFFLOAD=ON -DOFFLOAD_TARGET=x86_64-pc-linux-gnu \-DUSE_OBJECT_TARGET=ON -DQMC_MPI=0 \.`|
+-----------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

+------------+------------------------------------------------------------------------------------------------------+
|Test Command|`export LD_LIBRARY_PATH=/usr/lib/llvm-12/lib/:${LD_LIBRARY_PATH} && ctest --output-on-failure -L unit`|
+------------+------------------------------------------------------------------------------------------------------+

Self-Hosted Runners
-------------------

gpu-cuda (gcc-real-gpu-cuda)
""""""""""""""""""""""""""""
+-----------+-----------------------+
|Compiler   |GCC                    |
+-----------+-----------------------+
|MPI Wrapper|OpenMPI?               |
+-----------+-----------------------+
|Objective  |Build for Nvidia (Real)|
+-----------+-----------------------+
|Duration   |~2 Minutes             |
+-----------+-----------------------+

+-----------+-----------------------------------------------------------------------+
|Build Flags|`cmake -GNinja -DQMC_CUDA=1 \-DQMC_MPI=0 \-DQMC_COMPLEX=0 \.`          |
+-----------+-----------------------------------------------------------------------+

+------------+-------------------------------------------------------------------------------------------------------------------------------------+
|Test Command|`export LD_LIBRARY_PATH=/usr/local/cuda/lib/:/usr/local/cuda/lib64/:${LD_LIBRARY_PATH} && ctest --output-on-failure -L deterministic`|
+------------+-------------------------------------------------------------------------------------------------------------------------------------+

gpu-cuda (gcc-complex-gpu-cuda)
"""""""""""""""""""""""""""""""
+-----------+--------------------------+
|Compiler   |GCC                       |
+-----------+--------------------------+
|MPI Wrapper|OpenMPI?                  |
+-----------+--------------------------+
|Objective  |Build for Nvidia (Complex)|
+-----------+--------------------------+
|Duration   |~2 Minutes                |
+-----------+--------------------------+

+-----------+-----------------------------------------------------------------------+
|Build Flags|`cmake -GNinja -DQMC_CUDA=1 \-DQMC_MPI=0 \-DQMC_COMPLEX=1 \.`          |
+-----------+-----------------------------------------------------------------------+

+------------+-------------------------------------------------------------------------------------------------------------------------------------+
|Test Command|`export LD_LIBRARY_PATH=/usr/local/cuda/lib/:/usr/local/cuda/lib64/:${LD_LIBRARY_PATH} && ctest --output-on-failure -L deterministic`|
+------------+-------------------------------------------------------------------------------------------------------------------------------------+


How does it Know What to Do?
============================

We define these jobs in the ci-github-actions.yaml files located in the .github/workflows directory.  Each of the jobs currently run through the yaml files utilizing steps defined `here <https://github.com/QMCPACK/qmcpack/tree/develop/tests/test_automation/github-actions/ci/run_steps.sh>`_.
The yaml allows for deviation in steps based on the job name(for instance the job needs to contain the keyword 'coverage' in order to trigger the Coverage step) and other boolean checks.

The currently defined steps are:

Checkout Action
---------------
Triggers `actions/checkout@v1` which is a predifed GithubAction for checking out the repo.

Configure
---------
Based on certain keywords in the job name, it will add job specific flags.

Build
-----
Simple, after configuration it just issues a build command.

Test
----
Runs tests appropriate to job name.(complex vs real, asan, ect)

Coverage
--------
Generate code coverate reports once all tests have reported.

Upload Coverage
---------------
Upload the generated code coverage to `CodeCov <https://codecov.io/gh/QMCPACK/qmcpack/tree/develop/src>`_ where the badges on our repo will then be updated.



TODO: Firgure out examples where contributors might wanna add their own jobs and stuff, and how exactly they're supposed to do that.
TODO: Maybe layout some standards to keep everything clean and managable?
TODO: Review process for contributions? (security and such?)
TODO: Are we going to cover the different external runners in this and how to access them in the CI?
