.. _installation:

Nexus Installation
==================

Nexus aims to be very flexible in installation method. Options include use of a package manager and through manual cloning of the
repository and setting of an environment variable. Package manager-based installation eases installation of dependencies, including
Python itself, as well as interoperability with other Python-based software. The environment variable route is particularly suited
to development or (currently) if any scripts/"binaries" are needed. Administrative privileges are not required.

Note that while many supercomputer systems have an outdated system Python, a recent Python version is usually provided via a
loadable module: check the local documentation or discuss with your system administrators to find the simplest installation route.

.. contents::

Installation using ``pip``
--------------------------

This simple installation method can install Nexus and all of its dependencies provided that the underlying Python version is new
enough. A user-level only installation does not require any elevated privileges.


User-level installation (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To install Nexus only for yourself ("user level installation") we recommend creating a virtual environment. No privileges are
needed, making it the preferred route in most scenarios. Multiple distinct environments can be created, e.g. for different versions
of Nexus or for different projects needing distinct Python dependencies. Packages installed via this route do not affect the system
Python installation and are only available when the environment is activated.

.. code-block:: bash

    cd $HOME/somewhere # Choose an appropriate directory to create the environment
    python3 -m venv nexusenv # Give the environment a memorable and descriptive name
    source nexusenv/bin/activate
    pip install --upgrade pip # Optional unless pip is outdated
    pip install "nexus@git+https://github.com/QMCPACK/qmcpack.git@main#subdirectory=nexus"

This will install Nexus and the minimum required dependencies, currently only numpy. To take full advantage of Nexus's
capabilities, optional dependencies include ``scipy``, ``h5py``, ``matplotlib``, ``spglib``, ``cif2cell``,
``pydot``, and ``seekpath``. A complete install command would look like this

.. code-block:: bash

    pip install "nexus[full]@git+https://github.com/QMCPACK/qmcpack.git@main#subdirectory=nexus"

To subsequently use the environment:

.. code-block:: bash

    cd $HOME/somewhere # Wherever you created the environment
    source nexusenv/bin/activate

Note that any versions of Nexus found via the ``PYTHONPATH`` environment variable will take precedence over a ``pip`` installed version.

.. important::
    If you install Nexus into a virtual environment, the Nexus executables (e.g. ``qmca``, ``nxs-sim``, etc.) will only be accessible if you have activated that virtual environment.


System-wide installation (not recommended in general)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you wish to make Nexus available globally, such as within a container or on a personal laptop, a system-wide installation may be
appropriate. (The user-level only installation method described above is the preferred route in general and avoids any possibility
of interference with other installed Python packages.) You will generally need to use ``sudo`` for the installation commands. The
underlying Python needs to be sufficiently up-to-date; if not you must use the ``uv`` based approach (or an equivalent) described in
the next section. Do not modify the system Python.

The following will result in a system-wide installation: (prepending ``python3 -m`` if necessary)

.. code-block:: bash

    pip install "nexus@git+https://github.com/QMCPACK/qmcpack.git@main#subdirectory=nexus"

will install Nexus with ``numpy`` as the only dependency. This is the only required dependency. To take full advantage of Nexus's
capabilities, optional dependencies include ``scipy``, ``h5py``, ``matplotlib``, ``spglib``, ``cif2cell``,
``pydot``, and ``seekpath``. A complete install command would look like this

.. code-block:: bash

    pip install "nexus[full]@git+https://github.com/QMCPACK/qmcpack.git@main#subdirectory=nexus"

If you wish to install a subset of the optional dependencies, you can do so with ``pip`` in the same manner as shown in
:ref:`manual_install`.

Installation using the ``uv`` package manager
---------------------------------------------

``uv`` is a recently developed Python package and project manager. Under nearly all circumstances it can install Nexus and all of
its dependencies, with no other requirements. Importantly, besides being a drop-in and significantly faster replacement for ``pip``,
it can be installed by a standard user without any administrative privileges. Unlike ``pip``, it is also able to install specific
versions of Python, making it a solution for systems with outdated Python installations or simply a route to get the latest Python.
If you are not familiar with ``uv``, see the latest getting started and installation guide at https://docs.astral.sh/uv/ . ``uv``
can be installed with a single line. 

Standard practice is to use a virtual environment. Once ``uv`` is installed, you can create a virtual environment with the command

.. code-block:: bash

    uv venv .venv

which will create a virtual environment in the current directory inside a directory called ``.venv``. If you need or have a
preference to use a specific version of Python, see https://docs.astral.sh/uv/#python-versions .

After setting up the virtual environment you can simply prepend ``uv`` to the installation command:

.. code-block:: bash

    uv pip install "nexus@git+https://github.com/QMCPACK/qmcpack.git@main#subdirectory=nexus"

which will install Nexus to your current virtual environment. As described in the ``pip`` section above, this will only install the
minimum set of dependencies. You can also install the full set of optional dependencies via:

.. code-block:: bash

    uv pip install "nexus[full]@git+https://github.com/QMCPACK/qmcpack.git@main#subdirectory=nexus"

Note that when using ``uv`` as a package manager, your Python scripts can optionally use a modified `shebang
<https://en.wikipedia.org/wiki/Shebang_(Unix)>`__ that tells it to use ``uv`` to run it. For example

.. code-block:: rest

    #! /usr/bin/env -S uv run --script

Additionally, you can optionally use ``uv`` to set the dependencies of a script with the following command

.. code-block:: bash

    uv add --script <script_name>.py "nexus@git+https://github.com/QMCPACK/qmcpack.git@main#subdirectory=nexus"

This adds a piece of `Inline Script Metadata <https://packaging.python.org/en/latest/specifications/inline-script-metadata/>`__ to
the top of your file that specifies two main things, first the Python version, and second the dependencies of the script. For Nexus
this will look something like:

.. code-block:: python

    # /// script
    # requires-python = ">=3.14"
    # dependencies = [
    #     "nexus",
    # ]
    #
    # [tool.uv.sources]
    # nexus = { git = "https://github.com/QMCPACK/qmcpack.git", subdirectory = "nexus", rev "main" }
    # ///

This metadata ensures that ``uv`` can find the required dependencies for your project and will ensure that they are loaded at
runtime. This inline metadata will not exclude the script from being run directly with Python (e.g. ``python <script_name>.py``),
and so can be thought of simply as extra documentation. This can help when exchanging scripts with other users but is completely
optional.

.. tip::
    As with ``pip``, you can use ``uv`` to add any additional dependencies or packages you desire.

.. _manual_install:

Manual Installation of Nexus
----------------------------

.. warning::
    Although very simple, this method requires you maintain a consistent Python environment and environment variables. Be sure to
    consistently point Python at the version of Nexus that you intend.
    
Installing Nexus via this method requires downloading a QMCPACK release or cloning the `QMCPACK repository <https://github.com/QMCPACK/qmcpack>`__.

To make your Python installation aware of Nexus, simply set the ``PYTHONPATH`` environment variable to include the topmost nexus
directory.  For example, in bash this would look like:

.. code-block:: bash

    export PYTHONPATH=/your_download_path/nexus:$PYTHONPATH

Optionally add this to, *e.g.*, your ``.bashrc`` file to make Nexus available in future sessions. Be sure to specify the topmost
nexus directory, not the lower nexus directory ``/your_download_path/nexus/nexus``.

.. note::
    If your shell does not use ``export`` for environment variables, you will need to adapt the provided command for use with your shell. For example, the ``fish`` shell would need a command that looks like ``set -gx PYTHONPATH /your_download_path/nexus``.

If you want to use Nexus's command line tools, add them to your PATH:

.. code-block:: bash

    export PATH=/your_download_path/nexus/bin:$PATH


Installing Python dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When manually installing Nexus, you must also install the Python dependencies. The ``numpy`` module must be installed for Nexus to
function at a basic level. To realize the full range of functionality available, it is recommended that the ``scipy``,
``matplotlib``, ``h5py``, ``pydot``, ``spglib``, ``pycifrw``, ``cif2cell`` and ``seekpath`` modules be installed as well. In
supercomputing environments, most of these packages will *not* be available via system modules due to their specialized nature; you
will need to install them via ``pip``, ``uv``, or individual manual installation.

You can perform a user level installation via ``pip``. If using ``uv``, prepend ``uv``.

.. code-block:: bash

    pip3 install --user numpy
    pip3 install --user scipy
    pip3 install --user matplotlib
    pip3 install --user h5py
    pip3 install --user pydot
    pip3 install --user spglib
    pip3 install --user PyCifRW
    pip3 install --user cif2cell
    pip3 install --user seekpath

.. note::
    ``PyCifRW`` is a dependency of ``cif2cell``, thus if you install ``cif2cell`` you will already have ``PyCifRW`` installed as well.

While Nexus does not have strict version requirements, most recent dependency versions that have been tested and are known to work
can be found at ``qmcpack/nexus/pyproject.toml``. These specific library versions can be installed using the following command:

.. code-block:: bash

    pip3 install --user -r pyproject.toml

.. note::
    Additionally, for users of ``uv``, the project's ``uv.lock`` file can be used to install the exact versions of dependencies that Nexus is tested with using the command ``uv sync --all-extras --locked``.

If you are making a system-wide installation, on a Linux-based system, such as Ubuntu or Fedora, installation of these Python
modules can be accomplished using the distribution's package manager (e.g. ``apt`` for Debian-based systems like Ubuntu). For
example, here is how one may install some of the packages with ``apt``:

.. code-block:: bash

    sudo apt install python3-numpy
    sudo apt install python3-scipy python3-matplotlib python3-h5py
    sudo apt install python3-pydot
    sudo apt install python3-pip

Simple substitutions will be needed on distributions with different package managers (e.g. ``dnf``). Note that the most specialist
packages may not be available via this route.

The purpose of each library is described below:

``numpy`` 
    Needed throughout Nexus for array computation. Nexus will not function without NumPy.

``scipy`` 
    Used by the ``qmc-fit`` executable to perform least squares fits to QMCPACK DMC energies vs. timestep to perform extrapolation to zero timestep.

``matplotlib`` 
    Needed to view plots of QMCPACK data generated by the ``qmca`` and ``qmc-fit`` executables.

``h5py`` 
    Needed by the ``qdens`` tool to postprocess QMCPACK densities.

``pydot`` 
    Used to plot simulation workflow diagrams.

``spglib`` 
    Used to find crystalline primitive cells and perform k-point symmetrization.

``PyCifRW`` 
    Needed to read ``.cif`` crystal structure files.

``seekpath`` 
    Used to find high symmetry lines in reciprocal space for excited state calculations with QMCPACK.

Of course, to run full calculations, the simulation codes and converters involved must be installed as well. These include a patched
version of Quantum ESPRESSO (``pw.x``, ``pw2qmcpack.x``, optionally ``pw2casino.x``), QMCPACK (``qmcpack``, ``qmcpack_complex``,
``convert4qmc``, ``wfconvert``, ``ppconvert``), VASP, and/or GAMESS. Complete coverage of this task is beyond the scope of the
current document, but please see :ref:`install-code`.

Testing your Nexus installation
-------------------------------

Nexus's testing suite is designed to run with ``pytest >= 6.2.4``.
We additionally make use of the ``pytest-cov`` and ``pytest-order`` plugins, however these are both optional and will not cause test failure if they are not installed.

To run the tests with ``pytest`` (``pip install --user pytest``), enter the ``nexus`` directory and simply invoke the ``pytest`` command:

.. code-block:: bash

    > cd nexus/
    > pytest
    =========================== test session starts ============================
    platform linux -- Python 3.14.4, pytest-9.0.3, pluggy-1.6.0
    rootdir: /home/brock/Documents/github/qmcpack/nexus
    configfile: pyproject.toml
    plugins: order-1.3.0, cov-7.1.0
    collected 397 items

    test_versions.py .....                                               [  1%]
    test_required_dependencies.py .                                      [  1%]
    test_nexus_imports.py .                                              [  1%]
    test_testing.py ....                                                 [  2%]
    test_execute.py ..                                                   [  3%]
    test_memory.py ....                                                  [  4%]
    test_generic.py ...                                                  [  5%]
    test_developer.py ...                                                [  5%]
    test_unit_converter.py ...                                           [  6%]
    test_periodic_table.py ......                                        [  8%]
    test_numerics.py ...............                                     [ 11%]
    test_grid_functions.py ......................                        [ 17%]
    test_fileio.py .......                                               [ 19%]
    test_hdfreader.py ..                                                 [ 19%]
    test_xmlreader.py .......                                            [ 21%]
    test_structure.py ...................................                [ 30%]
    test_physical_system.py .......                                      [ 31%]
    test_basisset.py .....                                               [ 33%]
    test_pseudopotential.py ......                                       [ 34%]
    test_nexus_base.py .....                                             [ 36%]
    test_machines.py ......................                              [ 41%]
    test_simulation_module.py .........................................  [ 51%]
    test_bundle.py ..                                                    [ 52%]
    test_project_manager.py ...........                                  [ 55%]
    test_settings.py ..                                                  [ 55%]
    test_pwscf_input.py ...                                              [ 56%]
    test_pwscf_postprocessor_input.py .....                              [ 57%]
    test_gamess_input.py .......                                         [ 59%]
    test_pyscf_input.py ....                                             [ 60%]
    test_quantum_package_input.py ....                                   [ 61%]
    test_rmg_input.py ......                                             [ 62%]
    test_qmcpack_converter_input.py ..........                           [ 65%]
    test_qmcpack_input.py .............                                  [ 68%]
    test_vasp_analyzer.py ....                                           [ 69%]
    test_vasp_input.py .......                                           [ 71%]
    test_pwscf_analyzer.py ...                                           [ 72%]
    test_pwscf_postprocessor_analyzers.py ...                            [ 73%]
    test_gamess_analyzer.py ...                                          [ 73%]
    test_pyscf_analyzer.py ..                                            [ 74%]
    test_quantum_package_analyzer.py ..                                  [ 74%]
    test_rmg_analyzer.py ..                                              [ 75%]
    test_qmcpack_converter_analyzers.py ....                             [ 76%]
    test_qmcpack_analyzer.py ......                                      [ 77%]
    test_vasp_simulation.py .......                                      [ 79%]
    test_pwscf_simulation.py ......                                      [ 81%]
    test_gamess_simulation.py ......                                     [ 82%]
    test_pyscf_simulation.py .....                                       [ 83%]
    test_quantum_package_simulation.py ......                            [ 85%]
    test_rmg_simulation.py ..                                            [ 85%]
    test_pwscf_postprocessor_simulations.py ......                       [ 87%]
    test_qmcpack_converter_simulations.py ..................             [ 91%]
    test_qmcpack_simulation.py ......                                    [ 93%]
    test_observables.py ..                                               [ 93%]
    test_nxs_redo.py .                                                   [ 94%]
    test_nxs_sim.py .                                                    [ 94%]
    test_qmc_fit.py .                                                    [ 94%]
    test_qdens.py .                                                      [ 94%]
    test_qdens_radial.py .                                               [ 95%]
    test_qmca.py ...........                                             [ 97%]
    test_user_examples_alt.py ........                                   [100%]

    ==================== 397 passed, 43 warnings in 58.38s =====================

Some tests may be skipped depending on what dependencies you have available, or if they are marked to be skipped.
Additionally, you may see a number of warnings appear; some of these may be warnings about Nexus, but it is likely that the majority arise from a Nexus dependency.
In general it is safe to ignore these warnings as they likely do not affect the functionality of Nexus.

.. note::
    If you are planning on adding a new feature to Nexus, it is **strongly** encouraged to add a test for the new feature.
    Additionally, if you are planning on changing an existing Nexus feature, you must ensure that you either write a test for it, update the existing test, or ensure that your changes do not fail the existing test; your changes will not get merged otherwise!

If you have installed Nexus through ``pip`` or ``uv``, you can still use the ``nxs-test`` executable to test Nexus, however in the future this executable will be repurposed to simply call ``pytest``.

Developer Topics
----------------

Assessing Test Coverage
^^^^^^^^^^^^^^^^^^^^^^^

Code coverage can be assessed by using the ``pytest-cov`` plugin (``pip install --user pytest-cov``):

.. code-block:: bash

    > cd nexus
    > pytest --cov=nexus
    ...
    > coverage report

    Name                                  Stmts   Miss  Cover
    ---------------------------------------------------------
    nexus/__init__.py                       305     99    68%
    nexus/_bin.py                            25     25     0%
    nexus/basisset.py                       645    387    40%
    nexus/bin/nxs-redo                       90     31    66%
    nexus/bin/nxs-sim                       148     62    58%
    nexus/bin/qdens                         846    438    48%
    nexus/bin/qdens-radial                  282     84    70%
    nexus/bin/qmc-fit                       380    186    51%
    nexus/bin/qmca                          887    326    63%
    nexus/bundle.py                         191     68    64%
    nexus/debug.py                           12      6    50%
    nexus/developer.py                      265     88    67%
    nexus/execute.py                         14      2    86%
    nexus/fileio.py                        1019    407    60%
    nexus/gamess.py                         151     62    59%
    nexus/gamess_analyzer.py                306    131    57%
    nexus/gamess_input.py                   593    167    72%
    nexus/gaussian_process.py               943    943     0%
    nexus/generic.py                        842    184    78%
    nexus/grid_functions.py                1641    689    58%
    nexus/hdfreader.py                      205     59    71%
    nexus/machines.py                      2593    705    73%
    nexus/memory.py                          60      7    88%
    nexus/nexus_base.py                      76     10    87%
    nexus/nexus_version.py                    2      0   100%
    nexus/numerics.py                       904    461    49%
    nexus/numpy_extensions.py                 7      1    86%
    nexus/observables.py                    891    492    45%
    nexus/periodic_table.py                 301      4    99%
    nexus/physical_system.py                429     67    84%
    nexus/project_manager.py                233     37    84%
    nexus/pseudopotential.py               1675   1003    40%
    nexus/pwscf.py                          216     86    60%
    nexus/pwscf_analyzer.py                 653    283    57%
    nexus/pwscf_data_reader.py              130    120     8%
    nexus/pwscf_input.py                   1387    494    64%
    nexus/pwscf_postprocessors.py           517    113    78%
    nexus/pyscf_analyzer.py                   3      0   100%
    nexus/pyscf_input.py                    285     44    85%
    nexus/pyscf_sim.py                       64     14    78%
    nexus/qmcpack.py                       1044    789    24%
    nexus/qmcpack_analyzer.py               456    106    77%
    nexus/qmcpack_analyzer_base.py          326    136    58%
    nexus/qmcpack_converters.py             673    215    68%
    nexus/qmcpack_input.py                 4440   1906    57%
    nexus/qmcpack_method_analyzers.py       196     64    67%
    nexus/qmcpack_property_analyzers.py     204    100    51%
    nexus/qmcpack_quantity_analyzers.py    2095   1818    13%
    nexus/qmcpack_result_analyzers.py       289    144    50%
    nexus/quantum_package.py                253    141    44%
    nexus/quantum_package_analyzer.py         3      0   100%
    nexus/quantum_package_input.py          337    164    51%
    nexus/rmg.py                             33     10    70%
    nexus/rmg_analyzer.py                   293    265    10%
    nexus/rmg_input.py                      682    148    78%
    nexus/simulation.py                    1055    185    82%
    nexus/structure.py                     4121   2099    49%
    nexus/template_simulation.py             64     64     0%
    nexus/testing.py                        451     87    81%
    nexus/unit_converter.py                 120      2    98%
    nexus/utilities.py                       56     16    71%
    nexus/vasp.py                            93     15    84%
    nexus/vasp_analyzer.py                  547     73    87%
    nexus/vasp_input.py                     958    454    53%
    nexus/versions.py                       348     52    85%
    nexus/xmlreader.py                      298     53    82%
    ---------------------------------------------------------
    TOTAL                                 39651  17491    56%

The first column is the total number of statements, the second is the number not yet covered by the tests and the third is the percent covered. At the bottom is the sum total of all covered lines, all missed lines, and the average coverage percent.

To obtain an annotated view of the statements in the source that are not yet covered, run:

.. code-block:: bash

    > pytest --cov=nexus --cov-report=html

Open ``htmlcov/index.html`` in a browser to view the report. More information regarding the ``coverage`` tool can be found at https://coverage.readthedocs.io/en/latest/.

Creating Portable Nexus Builds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With recent changes to Nexus's structure, it has become possible to build binary distributions of Nexus in the form of a wheel file (extension ``.whl``).

Wheel files are compressed ZIP-format archives that are officially recommended by the `Python Packaging Authority (PyPA) <https://packaging.python.org/en/latest/specifications/binary-distribution-format/>`__ and can be built by anyone with the right tools.
To get started, you must have a tool capable of building a project; here we will be describing how to use ``uv`` for building distributions.
First, start by cloning the QMCPACK GitHub repository and navigating to the ``nexus`` directory. 
Next, though this step is not always necessary, it is recommended to run ``pytest`` to ensure you have a clean Nexus installation.
Finally, you can run the command ``uv build``, which will create a ``qmcpack/nexus/dist`` directory, and build two files to go in it.
The first is a tar archive of the source code, and the second is a ``.whl`` file.
This wheel file is then usable as an installation source, which you can access via the following command (if you are in a virtual environment):

.. code-block:: bash

    uv pip install qmcpack/nexus/dist/nexus-2.2.0-py3-none-any.whl

This will install Nexus into your virtual environment and make all of its libraries and executables available to you while you are in that virtual environment.

.. tip::
    This route for installing and distributing Nexus is especially useful if you are a system administrator as you can create consistent builds of Nexus at regular intervals and ensure that your users all have up-to-date versions available to them.
    It also allows for users to have access to their own version of Nexus if they choose to not use a global Python environment. This is because the wheels will install a copy of Nexus into their local/virtual environment, which they can edit as much as they please, without disturbing the global installation.
