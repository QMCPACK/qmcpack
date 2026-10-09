.. _dev-env:

Setting up a Development Environment
====================================

If you are interested in contributing to Nexus, setting up a development environment is the first step. This is most easily done with the ``uv`` project manager, which will help you get both a supported version of Python and allow you to install all relevant packages for development.

Prerequisites
-------------

* You must have `uv installed <https://docs.astral.sh/uv/getting-started/installation/>`__ already.
* You must have ``git`` installed on your machine.
* You should already have a fork of the `QMCPACK repository <https://github.com/QMCPACK/qmcpack>`__.


Setup commands
--------------

.. code-block:: bash

    > git clone https://github.com/<your_username>/qmcpack.git
    > cd qmcpack
    > git remote add upstream https://github.com/QMCPACK/qmcpack
    > cd nexus
    > uv sync --dev --group="docs"
    > source .venv/bin/activate # Change to suit your shell, e.g. `source .venv/bin/activate.fish`

With ``uv``, your project is built and installed to a virtual environment in the working directory in an editable mode. This means that as you make changes to the source code, your environment will be automatically and instantly updated with those changes upon saving the file. In addition, the Nexus scripts (``qmca``, ``eshdf``, etc.) will be available directly on the command line.

It is highly recommended to not commit changes to your fork's ``develop`` branch, as this should be your link to the upstream QMCPACK repository. You can pull the latest changes to the upstream ``develop`` branch with this command:

.. code-block:: bash

    > git switch develop
    > git pull upstream develop

After this, you will need to push the changes from your local branch to your remote fork on GitHub:

.. code-block:: bash

    > git push


Working in a Branch
-------------------

If you want to contribute to Nexus's development, you should create a branch where you can make changes and not alter your ``develop`` branch. Since Nexus is just a part of the QMCPACK repository, you may find it beneficial to prefix your Nexus-specific branches with ``nxs-``, so they are not confused with any changes you may be making to QMCPACK. You can easily create a new branch and switch to it with

.. code-block:: bash

    > git branch nxs-branch
    > git switch nxs-branch

or you can use the more convenient

.. code-block:: bash

    > git switch -c nxs-branch

which will automatically create the new branch and switch to it.

As you make changes to code, be sure to regularly add your changes to your local branch with ``git add``, and then when you are ready to commit your changes you can use

.. code-block:: bash

    > git commit -m "Add a descriptive commit message!"

Be sure to check that all of the changes you would like to make are staged with ``git status``!

If you would like to push your changes to your remote fork on GitHub (which is required to open a PR), you must first set the upstream repository for that branch:

.. code-block:: bash

    > git push --set-upstream origin nxs-branch

Even if you have no commits, this will link your local branch to your remote repository, and if you do have commits, they will be included in the remote repository.

Before you open a pull request or push changes to an open PR, please run ``pytest`` to see if your changes have broken testing! Every time a change is pushed or a pull request is opened, the QMCPACK automated CI actions will run, and if your branch breaks Nexus's test suite the automated CI will almost certainly fail, which in the case of the QMCPACK CI tests can result in holding up runners that could be run on someone else's PR.

Additionally, if you are making changes to a branch that has an active PR, please commit all of your changes *before* pushing them to the remote. Every time you push you re-trigger the CI, and just three pushes can bind up dozens of CI runners.


Tips and Tricks
---------------

This section contains some tips/tricks that may help you contribute to Nexus.

Faster Testing with ``pytest-xdist``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you find that Nexus's test suite is taking too long to run for your liking, you can install the ``pytest-xdist`` package, which allows you to run several tests in parallel. For example, on an 8-core laptop, running the test suite can take 30-60 seconds normally, but with the command ``pytest -n 8`` (to use 8 processes), the total test time can be cut down to 15 seconds. This is not 8 times faster as one might assume, and this is due to several tests that create subprocesses, such as the executable scripts and the user example tests.
