.. _github_actions:


GitHub Actions CI on QMCPACK
============================

QMCPACK uses GitHub Actions as part of the suite of continuous integration (CI) checks before a pull request can be merged in the main `develop` branch. GitHub Actions is an event driven automation tool that allows us to automatically execute commands in response to QMCPACK repo related actions. For example, merging a branch into master might then trigger our test scripts to run.

This guide covers the purpose and usual interactions a QMCPACK contributor would have with GitHub Actions CI.  For more information on GitHub Actions please refer to the official `Github Actions Docs <https://docs.github.com/en/actions/guides>`_ and our scripts located `here <https://github.com/QMCPACK/qmcpack/tree/develop/tests/test_automation/github-actions/ci>`_.

Currently we are using GitHub Actions to automatically handle a few different jobs. These jobs are either run on the GitHub provided build VM's or are pushed to our hardware, e.g., for GPU tests.

Note: This is not necessarily the intended typical way for users to build QMCPACK, please refer to our getting started and other build documentation for that.

Jobs running on GitHub hosted runners are triggered automatically. Permission from an admin is required to run jobs on self-hosted runners for security reasons. In addition, jobs running on GitHub hosted runners run automatically in parallel and the time each job takes may vary depending on system utilization. For information on the underlying hardware see the GitHub Actions `docs on the topic <https://docs.github.com/en/actions/using-github-hosted-runners/about-github-hosted-runners>`_.  

All Linux jobs GitHub Runner hosts currently use the `ghcr.io/qmcpack/ubuntu22-openmpi:latest <https://github.com/orgs/QMCPACK/packages/container/package/ubuntu22-openmpi>`_ docker image, if you would like to reproduce theses tests exactly using docker, please refer to `Running QMCPACK on Docker Containers <https://qmcpack.readthedocs.io/en/develop/running_docker.html>`_ section in the QMCPACK documentation. The macOS job runs directly on the `macos-latest GitHub Actions VM runner <https://docs.github.com/en/actions/using-github-hosted-runners/about-github-hosted-runners#supported-runners-and-hardware-resources>`_


.. note::

    Jobs with the \*-coverage pattern upload a Coverage report to `Codecov <https://app.codecov.io/gh/QMCPACK/qmcpack>`. There might be some delay between jobs sending independent results triggering a false positive coverage check, just wait until all jobs are done.  


Workflow Steps
==============

We define these jobs in the yaml files located in the .github/workflows directory.  Each of the jobs currently runs through the yaml files utilizing steps defined in a `test/test_automation/github-actions/ci/run_step.sh <https://github.com/QMCPACK/qmcpack/tree/develop/tests/test_automation/github-actions/ci/run_step.sh>`_ file.

This script applies workflow branching (if-else) based on the job name(for instance the job needs to contain the keyword 'coverage' in order to trigger the Coverage step) and other boolean checks.

The currently defined steps are:

Checkout Action
---------------
Triggers `actions/checkout@v1` which is a predefined GitHub Action for checking out the repo.

Configure
---------
Based on certain keywords in the job name, it will add job-specific flags.

Build
-----
After configuration it issues a build command.

Test
----
Runs tests appropriate to job name.(complex vs real, asan, etc.)

Coverage
--------
Generate code coverage reports once all tests have reported.

Upload Coverage
---------------
Upload the generated code coverage to `CodeCov <https://codecov.io/gh/QMCPACK/qmcpack/tree/develop/src>`_ where the badges on our repo will then be updated. Only done by jobs with name `*-coverage`.


Static Analysis Workflow
========================

A manually triggered workflow on the GitHub Actions tab can generate the required checks using the `clang-tidy <https://clang.llvm.org/extra/clang-tidy/>`_ static analyzer. The current approach is to set checks in the `qmcpack/src/.clang-tidy` configuration file and run using `clang-tidy` v14 on GitHub Actions runners. The workflow is not part of CI, and it's currently used for reporting potential warnings on the GitHub Actions logs as they are addressed on the `develop` branch as part of refactoring efforts for code quality. 

To run the workflow:
- Go to the Actions tab
- Click on the `static` workflow on the left
- Click on `Run workflow` on the right
- Use workflow from `Branch:develop` and click on the `Run workflow` button

**Note:** the current `.clang-tidy` configuration file is compatible with clang v14 and runs on the `ghcr.io/qmcpack/ubuntu22-openmpi:latest` docker image. To run locally on a Linux system use: `docker run -it user ghcr.io/qmcpack/ubuntu22-openmpi:latest /bin/bash` or refer to the :ref:`running_docker` section.

To build locally enabling `clang-tidy`` static checks defined in `qmcpack/src/.clang-tidy` use the CMake `-DCMAKE_CXX_CLANG_TIDY` option as follows:

.. code-block:: bash

    cmake -GNinja \
          -DCMAKE_C_COMPILER=clang \
          -DCMAKE_CXX_COMPILER=clang++ \
          -DCMAKE_BUILD_TYPE=Debug \
          -DCMAKE_CXX_CLANG_TIDY='clang-tidy' \
          /path/to/qmcpack

