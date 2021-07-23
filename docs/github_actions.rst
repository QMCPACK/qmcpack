.. _github_actions:

============================
Github Actions CI on QMCPACK
============================

QMCPACK uses GitHub Actions as part of the suite of continuous integration (CI) checks before a pull request can be merged in the main `develop` branch. Github Actions is an event driven automation tool that allows us to automatically execute commands in response to QMCPACK repo related actions. For example, merging a branch into master might then trigger our test scripts to run.

This guide covers the purpose and usual interactions a QMCPACK contributor with Github Actions CI.  For how Github Actions accomplishes these things and other implementation details please refer to the offical `Github Actions Docs <https://docs.github.com/en/actions/guides>`_ and our scripts located `here <https://github.com/QMCPACK/qmcpack/tree/develop/tests/test_automation/github-actions/ci>`_.

Currently we are using Github Actions to automatically handle a few different jobs. These jobs are then either ran on the Github provided build vm's or pushed to our supplied hardware.  Usually the jobs are only ran on our hardware when they require gpu's to run.


Summary of Test Jobs
--------------------

The following is a summary of the jobs run in the CI process required for a PR:

+-------------------------------+--------+----------+---------------+------+----------+
| Job Name with                 | Runner | Compiler | Tests         | Time | Trigger  |
| Build Info                    | Host   |          | ctest -L      | min  | event    |
+-------------------------------+--------+----------+---------------+------+----------+
| gcc-openmpi-real-coverage*    | GitHub | gcc-9    | deterministic | 45   | PR/merge |
+-------------------------------+--------+----------+---------------+------+----------+
| gcc-openmpi-complex-coverage* | GitHub | gcc-9    | deterministic | 50   | PR/merge |
+-------------------------------+--------+----------+---------------+------+----------+
| clang-openmpi-real-asan       | GitHub | clang-10 | unit          | 25   | PR/merge |
+-------------------------------+--------+----------+---------------+------+----------+
| clang-openmpi-complex-asan    | GitHub | clang-10 | unit          | 25   | PR/merge |
+-------------------------------+--------+----------+---------------+------+----------+
| clang-openmpi-real-ubsan      | GitHub | clang-10 | unit          | 55   | PR/merge |
+-------------------------------+--------+----------+---------------+------+----------+
| clang-latest-openmp-offload   | GitHub | clang-12 | unit          | 35   | PR/merge |
+-------------------------------+--------+----------+---------------+------+----------+
| gcc-real-gpu-cuda             | sulfur | clang-11 | deterministic | 2    | manual   |
+-------------------------------+--------+----------+---------------+------+----------+
| gcc-complex-gpu-cuda          | sulfur | clang-11 | deterministic | 2    | manual   |
+-------------------------------+--------+----------+---------------+------+----------+

Jobs running on GitHub hosted runners are triggered automatically. An admin is required to run jobs on self-hosted runners (e.g. sulfur) for security reasons. In addition, jobs running on GitHub hosted runners run automatically in parallel and the time each job takes may vary depending on system utilization. For information on the underlying hardware see the GitHub Actions `docs on the topic <https://docs.github.com/en/actions/using-github-hosted-runners/about-github-hosted-runners>`_.  

All jobs Github Runner hosts currently use the `williamfgc/qmcpack-ci:ubuntu20-openmpi <https://hub.docker.com/r/williamfgc/qmcpack-ci>`_ docker image, if you would like to reproduce theses tests exactly using docker, please refer to `Running QMCPACK on Docker Containers <https://qmcpack.readthedocs.io/en/develop/running_docker.html>`_ section in the QMCPACK documentation.


.. note::

    Jobs with the *-coverage pattern upload a Coverage report to `Codecov <https://app.codecov.io/gh/QMCPACK/qmcpack>`_. There might be some delay between jobs sending indepenpendent results triggering a false positive coverage check, just wait until all jobs are done.  



linux (gcc-openmpi-real-coverage)
"""""""""""""""""""""""""""""""""
+-------------+--------------------------------------------------------------------------------------------------------------------------------------------+
| Compiler    | GCC                                                                                                                                        |
+-------------+--------------------------------------------------------------------------------------------------------------------------------------------+
|Build Command|`cmake -GNinja -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx \-DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_GCOV=TRUE \-DQMC_COMPLEX=0 \.`  |
+-------------+--------------------------------------------------------------------------------------------------------------------------------------------+
|Test Command |`export OMPI_MCA_rmaps_base_oversubscribe=1 && export OMPI_MCA_hwloc_base_binding_policy=none && ctest --output-on-failure -L deterministic`|
+-------------+--------------------------------------------------------------------------------------------------------------------------------------------+
| Objective   | Test Coverage (Real)                                                                                                                       |
+-------------+--------------------------------------------------------------------------------------------------------------------------------------------+
| Duration    | ~45 Minutes                                                                                                                                |
+-------------+--------------------------------------------------------------------------------------------------------------------------------------------+

linux (gcc-openmpi-complex-coverage)
""""""""""""""""""""""""""""""""""""
+-------------+--------------------------------------------------------------------------------------------------------------------------------------------+
| Compiler    | GCC                                                                                                                                        |
+-------------+--------------------------------------------------------------------------------------------------------------------------------------------+
|Build Command|`cmake -GNinja -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx \-DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_GCOV=TRUE \-DQMC_COMPLEX=1 \.`  |
+-------------+--------------------------------------------------------------------------------------------------------------------------------------------+
|Test Command |`export OMPI_MCA_rmaps_base_oversubscribe=1 && export OMPI_MCA_hwloc_base_binding_policy=none && ctest --output-on-failure -L deterministic`|
+-------------+--------------------------------------------------------------------------------------------------------------------------------------------+
| Objective   | Test Coverage (Complex)                                                                                                                    |
+-------------+--------------------------------------------------------------------------------------------------------------------------------------------+
| Duration    | ~50 Minutes                                                                                                                                |
+-------------+--------------------------------------------------------------------------------------------------------------------------------------------+

linux (clang-openmpi-real-asan)
"""""""""""""""""""""""""""""""
+-------------+------------------------------------------------------------------------------------------------------------------------------------------------+
| Compiler    | Clang                                                                                                                                          |
+-------------+------------------------------------------------------------------------------------------------------------------------------------------------+
|Build Command|`CC=clang CXX=clang++ \cmake -GNinja -DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_SANITIZER=asan \-DQMC_MPI=0 \-DQMC_COMPLEX=0 \.`                |
+-------------+------------------------------------------------------------------------------------------------------------------------------------------------+
|Test Command |`export OMPI_MCA_rmaps_base_oversubscribe=1 && export OMPI_MCA_hwloc_base_binding_policy=none && ctest --output-on-failure -L unit -LE noasan`  |
+-------------+------------------------------------------------------------------------------------------------------------------------------------------------+
| Objective   | Address Sanitizer (Real)                                                                                                                       |
+-------------+------------------------------------------------------------------------------------------------------------------------------------------------+
| Duration    | ~25 Minutes                                                                                                                                    |
+-------------+------------------------------------------------------------------------------------------------------------------------------------------------+

linux (clang-openmpi-complex-asan)
""""""""""""""""""""""""""""""""""
+-------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| Compiler    | Clang                                                                                                                                        |
+-------------+----------------------------------------------------------------------------------------------------------------------------------------------+
|Build Command|`CC=clang CXX=clang++ \cmake -GNinja -DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_SANITIZER=asan \-DQMC_MPI=0 \-DQMC_COMPLEX=1 \.`              |
+-------------+----------------------------------------------------------------------------------------------------------------------------------------------+
|Test Command |`export OMPI_MCA_rmaps_base_oversubscribe=1 && export OMPI_MCA_hwloc_base_binding_policy=none && ctest --output-on-failure -L unit -LE noasan`|
+-------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| Objective   | Address Sanitizer (Complex)                                                                                                                  |
+-------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| Duration    | ~30 Minutes                                                                                                                                  |
+-------------+----------------------------------------------------------------------------------------------------------------------------------------------+


linux (clang-openmpi-real-ubsan)
""""""""""""""""""""""""""""""""
+-------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Compiler    | Clang                                                                                                                                                                |
+-------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|Build Command|`CC=clang CXX=clang++ \cmake -GNinja -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx \-DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_SANITIZER=ubsan \-DQMC_COMPLEX=0 \.`|
+-------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|Test Command |`export OMPI_MCA_rmaps_base_oversubscribe=1 && export OMPI_MCA_hwloc_base_binding_policy=none && ctest --output-on-failure -L deterministic`                          |
+-------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Objective   | Undefined Behavior Sanitizer (Real)                                                                                                                                  |
+-------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Duration    | ~55 Minutes                                                                                                                                                          |
+-------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+

linux (clang-latest-openmp-offload)
"""""""""""""""""""""""""""""""""""
+-------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Compiler    | Clang                                                                                                                                                                      |
+-------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|Build Command|`cmake -GNinja -DCMAKE_C_COMPILER=clang-12 -DCMAKE_CXX_COMPILER=clang++-12 \-DENABLE_OFFLOAD=ON -DOFFLOAD_TARGET=x86_64-pc-linux-gnu \-DUSE_OBJECT_TARGET=ON -DQMC_MPI=0 \.`|
+-------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|Test Command |`export LD_LIBRARY_PATH=/usr/lib/llvm-12/lib/:${LD_LIBRARY_PATH} && ctest --output-on-failure -L unit`                                                                      |
+-------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Objective   | Build for GPU Acceleration (Experimental)                                                                                                                                  |
+-------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Duration    | ~35 Minutes                                                                                                                                                                |
+-------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Self-Hosted Runners
-------------------

gpu-cuda (gcc-real-gpu-cuda)
""""""""""""""""""""""""""""
+-------------+-------------------------------------------------------------------------------------------------------------------------------------+
| Compiler    | GCC                                                                                                                                 |
+-------------+-------------------------------------------------------------------------------------------------------------------------------------+
|Build Command|`cmake -GNinja -DQMC_CUDA=1 \-DQMC_MPI=0 \-DQMC_COMPLEX=0 \.`                                                                        |
+-------------+-------------------------------------------------------------------------------------------------------------------------------------+
|Test Command |`export LD_LIBRARY_PATH=/usr/local/cuda/lib/:/usr/local/cuda/lib64/:${LD_LIBRARY_PATH} && ctest --output-on-failure -L deterministic`|
+-------------+-------------------------------------------------------------------------------------------------------------------------------------+
| Objective   | Build for Nvidia (Real)                                                                                                             |
+-------------+-------------------------------------------------------------------------------------------------------------------------------------+
| Duration    | ~2 Minutes                                                                                                                          |
+-------------+-------------------------------------------------------------------------------------------------------------------------------------+


gpu-cuda (gcc-complex-gpu-cuda)
"""""""""""""""""""""""""""""""
+-------------+-------------------------------------------------------------------------------------------------------------------------------------+
| Compiler    | GCC                                                                                                                                 |
+-------------+-------------------------------------------------------------------------------------------------------------------------------------+
|Build Command|`cmake -GNinja -DQMC_CUDA=1 \-DQMC_MPI=0 \-DQMC_COMPLEX=1 \.`                                                                        |
+-------------+-------------------------------------------------------------------------------------------------------------------------------------+
|Test Command |`export LD_LIBRARY_PATH=/usr/local/cuda/lib/:/usr/local/cuda/lib64/:${LD_LIBRARY_PATH} && ctest --output-on-failure -L deterministic`|
+-------------+-------------------------------------------------------------------------------------------------------------------------------------+ 
| Objective   | Build for Nvidia (Complex)                                                                                                          |
+-------------+-------------------------------------------------------------------------------------------------------------------------------------+
| Duration    | ~2 Minutes                                                                                                                          |
+-------------+-------------------------------------------------------------------------------------------------------------------------------------+

Workflow Steps
==============

We define these jobs in the yaml files located in the .github/workflows directory.  Each of the jobs currently run through the yaml files utilizing steps defined in a `test/test_automation/github-actions/ci/run_step.sh <https://github.com/QMCPACK/qmcpack/tree/develop/tests/test_automation/github-actions/ci/run_steps.sh>`_ file.

This script applies workflow branching (if-else) based on the job name(for instance the job needs to contain the keyword 'coverage' in order to trigger the Coverage step) and other boolean checks.

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
Upload the generated code coverage to `CodeCov <https://codecov.io/gh/QMCPACK/qmcpack/tree/develop/src>`_ where the badges on our repo will then be updated. Only done by jobs with name `*-coverage`.
