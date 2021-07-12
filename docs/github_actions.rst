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

One of the biggest hows and whys is that it saves time! There are certain jobs we always trigger upon merge, and having a designated person trigger and manage them takes a lot of time off their hands.  Time they could be spending actually improving the project with. 

Currently we are using Github Actions to automatically handle the following jobs:

Github Runners (Free)
---------------------

linux (gcc-openmpi-real-coverage)
"""""""""""""""""""""""""""""""""
+-----------+--------------------+
|Compiler   |GCC                 |
+-----------+--------------------+
|MPI Wrapper|OpenMPI             |
+-----------+--------------------+
|Objective  |Test Coverage (Real)|
+-----------+--------------------+

linux (gcc-openmpi-complex-coverage)
""""""""""""""""""""""""""""""""""""
+-----------+-----------------------+
|Compiler   |GCC                    |
+-----------+-----------------------+
|MPI Wrapper|OpenMPI                |
+-----------+-----------------------+
|Objective  |Test Coverage (Complex)|
+-----------+-----------------------+

linux (clang-openmpi-real-asan)
"""""""""""""""""""""""""""""""
+-----------+------------------------+
|Compiler   |Clang                   |
+-----------+------------------------+
|MPI Wrapper|OpenMPI                 |
+-----------+------------------------+
|Objective  |Address Sanitizer (Real)|
+-----------+------------------------+

linux (clang-openmpi-complex-asan)
""""""""""""""""""""""""""""""""""
+-----------+---------------------------+
|Compiler   |Clang                      |
+-----------+---------------------------+
|MPI Wrapper|OpenMPI                    |
+-----------+---------------------------+
|Objective  |Address Sanitizer (Complex)|
+-----------+---------------------------+

linux (clang-openmpi-real-ubsan)
""""""""""""""""""""""""""""""""
+-----------+-----------------------------------+
|Compiler   |Clang                              |
+-----------+-----------------------------------+
|MPI Wrapper|OpenMPI                            |
+-----------+-----------------------------------+
|Objective  |Undefined Behavior Sanitizer (Real)|
+-----------+-----------------------------------+

linux (clang-latest-openmp-offload)
"""""""""""""""""""""""""""""""""""
+-----------+-----------------------------------------+
|Compiler   |Clang                                    |
+-----------+-----------------------------------------+
|MPI Wrapper|OpenMPI                                  |
+-----------+-----------------------------------------+
|Objective  |Build for GPU Acceleration (Experimental)|
+-----------+-----------------------------------------+

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

gpu-cuda (gcc-complex-gpu-cuda)
"""""""""""""""""""""""""""""""
+-----------+--------------------------+
|Compiler   |GCC                       |
+-----------+--------------------------+
|MPI Wrapper|OpenMPI?                  |
+-----------+--------------------------+
|Objective  |Build for Nvidia (Complex)|
+-----------+--------------------------+

How does it Know What to Do?
============================

We define these jobs in the ci-github-actions.yaml files located in the .github/workflows directory.  Each of the jobs currently run through the yaml files utilizing steps defined `here <https://github.com/QMCPACK/qmcpack/tree/develop/tests/test_automation/github-actions/ci/run_steps.sh>`_.
The yaml allows for deviation in steps based on the job name(for instance the job needs to contain the keyword 'coverage' in order to trigger) and other boolean checks.

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

This is Pretty Cool, How Do I Contribute?
=========================================

TODO: Firgure out examples where contributors might wanna add their own jobs and stuff, and how exactly they're supposed to do that.
TODO: Maybe layout some standards to keep everything clean and managable?
TODO: Review process for contributions? (security and such?)
TODO: Are we going to cover the different external runners in this and how to access them in the CI?
