.. _installation:

Nexus Installation
==================

Installation of Nexus can be accomplished in two ways, either through a package manager or through manual cloning of the QMCPACK repository and setting of an environment variable.

Each section below talks about the related upsides and downsides, however for those who are not actively developing Nexus or QMCPACK, it is recommended to use a package manager.

.. contents::

.. _ via a single download from ``qmcpack.org``

Installation with ``pip`` or ``uv``
-----------------------------------

1) ``pip``
^^^^^^^^^^

.. note::
    Installing via ``pip`` requires at least ``pip`` version 10.0.0 or above. Older versions can be upgraded with the command ``python3 -m pip install -U pip``.

If you are working in a system where you can manage your own Python environment (e.g. a workstation or personal laptop), installing Nexus via ``pip`` is simple. Running the command (prepending ``python3 -m`` if necessary)

.. code-block:: bash

    > pip install "nexus@git+https://github.com/QMCPACK/qmcpack.git@main#subdirectory=nexus"

will install Nexus with ``numpy`` as the only dependency. This is required for Nexus to function, however should the user want to take full advantage of Nexus's capabilities, we offer a set of optional dependencies. These include ``scipy``, ``h5py``, ``matplotlib``, ``spglib``, ``cif2cell``, ``pydot``, and ``seekpath``. A complete install command would look like this

.. code-block:: bash

    > pip install "nexus[full]@git+https://github.com/QMCPACK/qmcpack.git@main#subdirectory=nexus"

If you do not want to install all of the dependencies, you can do so with ``pip``, in the same manner as shown in :ref:`manual_install`.

.. caution::
    This method of installation is not recommended for those who wish to customize Nexus for a specific project as ``pip`` will default to installing Nexus to the global Python environment and thus any changes made to the source code there will affect all projects that use it.
    
    Additionally, if you have an existing Nexus installation that has modified your ``$PYTHONPATH`` environment variable then that existing installation will override a version installed by ``pip``

1) ``uv``
^^^^^^^^^

.. note::
    If you are not able to install software on the machine you are working in, please skip to :ref:`manual_install`.

If you are working in a system where you can not modify the global Python installation, the standard practice is to use a virtual environment. A common choice is ``uv``, which can be installed by following the instructions `here <https://docs.astral.sh/uv/getting-started/installation/>`__. With ``uv`` installed, you can create a virtual environment with the command

.. code-block:: bash

    > uv venv .venv

which will create a virtual environment in the current directory inside a directory called ``.venv``. Following this, you can simply prepend ``uv`` to the installation command,

.. code-block:: bash

    > uv pip install "nexus@git+https://github.com/QMCPACK/qmcpack.git@main#subdirectory=nexus"

which will install Nexus to your current virtual environment. Note that when using ``uv`` as a package manager, your Python scripts should use a modified `shebang <https://en.wikipedia.org/wiki/Shebang_(Unix)>`__ that tells it to use ``uv`` to run it. For example, on a local machine, the shebang would look like

.. code-block:: rest

    #! /usr/bin/env -S uv run --script

Additionally, you can use ``uv`` to automatically set the dependencies of a script with the following command

.. code-block:: bash

    uv add --script <script_name>.py "nexus@git+https://github.com/QMCPACK/qmcpack.git@main#subdirectory=nexus"

which ensures that, as long as the script is run with ``uv``, Nexus will be available. This does not, however, exclude the script from being run directly with Python (e.g. ``python <script_name>.py``)

.. tip::
    You can use ``uv`` to add any additional dependencies you desire, however any dependencies installed with Nexus (e.g. ``numpy``, and if you add ``[full]``, all of the optional dependencies) will already be installed if you have Nexus added to the script, so there is no need to add them again.

.. note::
    If you want to customize Nexus for a specific project but do not want to continually download QMCPACK, this option will create a version of Nexus that is local to each virtual environment, so you can modify it without fear of altering other virtual environments. Importantly however, this has two caveats, the first being that this installation method does indeed create duplicates of Nexus in each virtual environment (though this is a minimal side effect, as the entirety of Nexus is only a handful of megabytes in size), and secondly that there is no version tracking in a ``uv`` installation, so any changes made will not be tracked via ``git``.

.. _manual_install:

Manual Installation of Nexus
----------------------------

.. warning::
    This method of installing Nexus can lead to undefined behavior! For example, if you wish to update your Nexus installation, then you must make sure that your ``.bashrc`` file is updated to point at the correct path, otherwise you will still be using the previous version of Nexus. Installation via this route can also override existing installations of Nexus via ``pip`` or ``uv``.
    
    Additionally, due to the method `by which Python searches for modules <https://docs.python.org/3/library/sys_path_init.html>`__, you can get situations where another Python script in the same directory as your current script will override Nexus's modules. This problem is largely mitigated by `PR #5700 <https://github.com/QMCPACK/qmcpack/pull/5700>`__ which changed Nexus's import method to behave more like a Python package, however not all edge cases have been tested. 
    
    If you encounter unusual behavior, please open an issue at `the QMCPACK GitHub page <https://github.com/QMCPACK/qmcpack>`__.

To make your Python installation (must be Python 3.x, 2.x is no longer supported) aware of Nexus, simply set the ``PYTHONPATH`` environment variable.  For example, in bash this would look like:

.. code-block:: bash

    export PYTHONPATH=/your_download_path/nexus:$PYTHONPATH

Add this to, *e.g.*, your ``.bashrc`` file to make Nexus available in future sessions.

If you want to use Nexus's command line tools, add them to your path:

.. code-block:: bash

    export PATH=/your_download_path/nexus/bin:$PATH

Both of these environment variables can be set automatically by the ``install`` script packaged with Nexus.  To use the installer, instead of performing the manual installation above, simply type the following at the command line:

::

    /your_download_path/nexus/install

If you want the Nexus binaries to reside a location different than ``/your_download_path/nexus/bin``, simply provide this path to the installer:

::

    /your_download_path/nexus/install /some/other/location

Installing Python dependencies
------------------------------

In addition to the standard Python installation, the ``numpy`` module must be installed for Nexus to function at a basic level. To realize the full range of functionality available, it is recommended that the ``scipy``, ``matplotlib``, ``h5py``, ``pydot``, ``spglib``, ``pycifrw``, ``cif2cell`` and ``seekpath`` modules be installed as well. Many of these packages are already available in various supercomputing environments.

On a Linux systems, such as Ubuntu or Fedora, installation of these Python modules is easily accomplished by using your distribution's package manager (e.g. ``apt`` for Debian-based systems, ``dnf`` for Fedora, etc.). For example, here is how one may install some of the packages with ``apt``:

.. code-block:: bash

    sudo apt install python3-numpy
    sudo apt install python3-scipy python3-matplotlib python3-h5py
    sudo apt install python3-pydot
    sudo apt install python3-pip

Other Linux systems often come with their own package managers (such as ``dnf`` for Fedora), and can often be used as drop-in replacements for ``apt``.

To install the Python modules on other platforms (as well as those not listed with ``apt`` above on Debian systems), try ``pip`` or ``pip3``:

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

While Nexus does not have strict version requirements, most recent dependency versions that have been tested and are known to work can be found at ``qmcpack/nexus/requirements.txt``. These specific library versions can be installed using the following command:

.. code-block:: bash

    pip3 install --user -r requirements.txt

``qmcpack/nexus/requirements_minimal.txt`` can be used similarly but only contains a recently tested version of ``numpy``.

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

Of course, to run full calculations, the simulation codes and converters
involved must be installed as well. These include a modified version of
Quantum ESPRESSO (``pw.x``, ``pw2qmcpack.x``, optionally
``pw2casino.x``), QMCPACK (``qmcpack``, ``qmcpack_complex``,
``convert4qmc``, ``wfconvert``, ``ppconvert``), SQD (``sqd``, packaged
with QMCPACK), VASP, and/or GAMESS. Complete coverage of this task is
beyond the scope of the current document, but please see :ref:`install-code`.

Testing your Nexus installation
-------------------------------

Nexus is packaged with an extensive suite of tests that can be run with
either the ``nxs-test`` executable packaged with Nexus or with
``pytest``. If you have installed Nexus, the ``nxs-test`` tool should be
in your ``PATH``. Installation is successful if all tests pass:

::

    > nxs-test

    1/67 versions................................   Passed  0.01 sec
    2/67 required_dependencies...................   Passed  0.00 sec
    3/67 nexus_base..............................   Passed  0.00 sec
    4/67 nexus_imports...........................   Passed  0.02 sec
    5/67 testing.................................   Passed  0.03 sec
    6/67 execute.................................   Passed  0.00 sec
    7/67 memory..................................   Passed  0.00 sec
    8/67 generic.................................   Passed  0.01 sec
    9/67 developer...............................   Passed  0.00 sec
    10/67 unit_converter..........................   Passed  0.00 sec
    11/67 periodic_table..........................   Passed  0.00 sec
    12/67 numerics................................   Passed  0.02 sec
    13/67 grid_functions..........................   Passed  0.70 sec
    14/67 fileio..................................   Passed  0.02 sec
    15/67 hdfreader...............................   Passed  0.00 sec
    16/67 xmlreader...............................   Passed  0.00 sec
    17/67 structure...............................   Passed  0.50 sec
    18/67 physical_system.........................   Passed  0.03 sec
    19/67 basisset................................   Passed  0.02 sec
    20/67 pseudopotential.........................   Passed  0.29 sec
    21/67 machines................................   Passed  0.98 sec
    22/67 simulation_module.......................   Passed  0.34 sec
    23/67 bundle..................................   Passed  0.00 sec
    24/67 project_manager.........................   Passed  2.90 sec
    25/67 settings................................   Passed  0.00 sec
    26/67 vasp_input..............................   Passed  0.01 sec
    27/67 pwscf_input.............................   Passed  0.01 sec
    28/67 pwscf_postprocessor_input...............   Passed  0.00 sec
    29/67 gamess_input............................   Passed  0.00 sec
    30/67 pyscf_input.............................   Passed  0.01 sec
    31/67 quantum_package_input...................   Passed  0.06 sec
    32/67 rmg_input...............................   Passed  0.06 sec
    33/67 qmcpack_converter_input.................   Passed  0.00 sec
    34/67 qmcpack_input...........................   Passed  0.06 sec
    35/67 vasp_analyzer...........................   Passed  0.02 sec
    36/67 pwscf_analyzer..........................   Passed  0.16 sec
    37/67 pwscf_postprocessor_analyzers...........   Passed  0.00 sec
    38/67 gamess_analyzer.........................   Passed  0.00 sec
    39/67 pyscf_analyzer..........................   Passed  0.00 sec
    40/67 quantum_package_analyzer................   Passed  0.00 sec
    41/67 rmg_analyzer............................   Passed  0.00 sec
    42/67 qmcpack_converter_analyzers.............   Passed  0.00 sec
    43/67 qmcpack_analyzer........................   Passed  0.23 sec
    44/67 vasp_simulation.........................   Passed  0.02 sec
    45/67 pwscf_simulation........................   Passed  0.01 sec
    46/67 gamess_simulation.......................   Passed  0.00 sec
    47/67 pyscf_simulation........................   Passed  0.00 sec
    48/67 quantum_package_simulation..............   Passed  0.00 sec
    49/67 rmg_simulation..........................   Passed  0.00 sec
    50/67 pwscf_postprocessor_simulations.........   Passed  0.00 sec
    51/67 qmcpack_converter_simulations...........   Passed  0.01 sec
    52/67 qmcpack_simulation......................   Passed  0.20 sec
    53/67 observables.............................   Passed  0.00 sec
    54/67 nxs_redo................................   Passed  0.92 sec
    55/67 nxs_sim.................................   Passed  2.03 sec
    56/67 qmca....................................   Passed  2.90 sec
    57/67 qmc_fit.................................   Passed  0.00 sec
    58/67 qdens...................................   Passed  0.00 sec
    59/67 qdens_radial............................   Passed  0.00 sec
    60/67 example_gamess_H2O......................   Passed  1.01 sec
    61/67 example_pwscf_relax_Ge_T................   Passed  0.47 sec
    62/67 example_qmcpack_H2O.....................   Passed  0.54 sec
    63/67 example_qmcpack_LiH.....................   Passed  0.54 sec
    64/67 example_qmcpack_c20.....................   Passed  0.52 sec
    65/67 example_qmcpack_diamond.................   Passed  0.73 sec
    66/67 example_qmcpack_graphene................   Passed  0.59 sec
    67/67 example_qmcpack_oxygen_dimer............   Passed  0.50 sec

    100% tests passed, 0 tests failed out of 67

    Total test time = 17.54 sec

Only portions of Nexus consistent with your Python installed Python libraries will be tested.

To run the tests with ``pytest`` (``pip install --user pytest``), enter the unit test directory and simply invoke the ``pytest`` command:

.. code-block:: bash

    > cd nexus/tests/unit/
    > pytest
    =========================== test session starts ============================
    platform linux -- Python 3.14.2, pytest-9.0.2, pluggy-1.6.0
    rootdir: qmcpack/nexus
    configfile: pyproject.toml
    plugins: cov-7.0.0
    collected 393 items                                                        

    test_basisset.py .....                                               [  1%]
    test_bundle.py ..                                                    [  1%]
    test_developer.py ...                                                [  2%]
    test_execute.py ..                                                   [  3%]
    test_fileio.py .......                                               [  4%]
    test_gamess_analyzer.py ...                                          [  5%]
    test_gamess_input.py .......                                         [  7%]
    test_gamess_simulation.py ......                                     [  8%]
    test_generic.py ...                                                  [  9%]
    test_grid_functions.py ......................                        [ 15%]
    test_hdfreader.py ..                                                 [ 15%]
    test_machines.py ......................                              [ 21%]
    test_memory.py ....                                                  [ 22%]
    test_nexus_base.py .....                                             [ 23%]
    test_nexus_imports.py .                                              [ 23%]
    test_numerics.py ...............                                     [ 27%]
    test_nxs_redo.py .                                                   [ 27%]
    test_nxs_sim.py .                                                    [ 28%]
    test_observables.py ..                                               [ 28%]
    test_optional_dependencies.py .......                                [ 30%]
    test_periodic_table.py ...                                           [ 31%]
    test_physical_system.py .......                                      [ 33%]
    test_project_manager.py ...........                                  [ 35%]
    test_pseudopotential.py ......                                       [ 37%]
    test_pwscf_analyzer.py ...                                           [ 38%]
    test_pwscf_input.py ...                                              [ 38%]
    test_pwscf_postprocessor_analyzers.py ...                            [ 39%]
    test_pwscf_postprocessor_input.py .....                              [ 40%]
    test_pwscf_postprocessor_simulations.py ......                       [ 42%]
    test_pwscf_simulation.py ......                                      [ 44%]
    test_pyscf_analyzer.py ..                                            [ 44%]
    test_pyscf_input.py ....                                             [ 45%]
    test_pyscf_simulation.py .....                                       [ 46%]
    test_qdens.py .                                                      [ 47%]
    test_qdens_radial.py .                                               [ 47%]
    test_qmc_fit.py .                                                    [ 47%]
    test_qmca.py ...........                                             [ 50%]
    test_qmcpack_analyzer.py ......                                      [ 51%]
    test_qmcpack_converter_analyzers.py ....                             [ 52%]
    test_qmcpack_converter_input.py ..........                           [ 55%]
    test_qmcpack_converter_simulations.py ..................             [ 60%]
    test_qmcpack_input.py .............                                  [ 63%]
    test_qmcpack_simulation.py ......                                    [ 64%]
    test_quantum_package_analyzer.py ..                                  [ 65%]
    test_quantum_package_input.py ....                                   [ 66%]
    test_quantum_package_simulation.py ......                            [ 67%]
    test_required_dependencies.py .                                      [ 68%]
    test_rmg_analyzer.py ..                                              [ 68%]
    test_rmg_input.py ......                                             [ 70%]
    test_rmg_simulation.py ..                                            [ 70%]
    test_settings.py ..                                                  [ 71%]
    test_simulation_module.py .........................................  [ 81%]
    test_structure.py ...................................                [ 90%]
    test_testing.py ....                                                 [ 91%]
    test_unit_converter.py ...                                           [ 92%]
    test_vasp_analyzer.py ....                                           [ 93%]
    test_vasp_input.py .......                                           [ 95%]
    test_vasp_simulation.py .......                                      [ 96%]
    test_versions.py .....                                               [ 98%]
    test_xmlreader.py .......                                            [100%]

Assessing Test Coverage (Developer Topic)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Code coverage can be assessed by using the ``pytest-cov`` plugin (``pip install --user pytest-cov``):

.. code-block:: bash

    >cd nexus
    >pytest-cov --cov=nexus 
    ...
    >coverage report | grep nexus/nexus

    nexus/nexus/basisset.py                      631    375    41%
    nexus/nexus/bundle.py                        191     68    64%
    nexus/nexus/debug.py                          12      6    50%
    nexus/nexus/developer.py                     261     97    63%
    nexus/nexus/execute.py                        13      2    85%
    nexus/nexus/fileio.py                        957    373    61%
    nexus/nexus/gamess.py                        102     20    80%
    nexus/nexus/gamess_analyzer.py               305    149    51%
    nexus/nexus/gamess_input.py                  597    167    72%
    nexus/nexus/generic.py                       817    173    79%
    nexus/nexus/grid_functions.py               1192    435    64%
    nexus/nexus/hdfreader.py                     215     61    72%
    nexus/nexus/machines.py                     1887    463    75%
    nexus/nexus/memory.py                         60      7    88%
    nexus/nexus/nexus.py                         297    140    53%
    nexus/nexus/nexus_base.py                     74     11    85%
    nexus/nexus/numerics.py                      756    372    51%
    nexus/nexus/periodic_table.py               1505     24    98%
    nexus/nexus/physical_system.py               427     73    83%
    nexus/nexus/plotting.py                       22      7    68%
    nexus/nexus/project_manager.py               234     37    84%
    nexus/nexus/pseudopotential.py              1225    559    54%
    nexus/nexus/pwscf.py                         198     73    63%
    nexus/nexus/pwscf_analyzer.py                634    316    50%
    nexus/nexus/pwscf_data_reader.py             132    120     9%
    nexus/nexus/pwscf_input.py                  1261    563    55%
    nexus/nexus/pwscf_postprocessors.py          434     56    87%
    nexus/nexus/pyscf_analyzer.py                  3      0   100%
    nexus/nexus/pyscf_input.py                   181     26    86%
    nexus/nexus/pyscf_sim.py                      57      8    86%
    nexus/nexus/qmcpack.py                       344    146    58%
    nexus/nexus/qmcpack_analyzer.py              457    104    77%
    nexus/nexus/qmcpack_analyzer_base.py         327    137    58%
    nexus/nexus/qmcpack_converters.py            507     83    84%
    nexus/nexus/qmcpack_input.py                3605   1439    60%
    nexus/nexus/qmcpack_method_analyzers.py      198     64    68%
    nexus/nexus/qmcpack_property_analyzers.py    205     97    53%
    nexus/nexus/qmcpack_quantity_analyzers.py   2070   1789    14%
    nexus/nexus/qmcpack_result_analyzers.py      285    142    50%
    nexus/nexus/quantum_package.py               253    141    44%
    nexus/nexus/quantum_package_analyzer.py        3      0   100%
    nexus/nexus/quantum_package_input.py         338    164    51%
    nexus/nexus/simulation.py                   1019    169    83%
    nexus/nexus/structure.py                    3830   2055    46%
    nexus/nexus/superstring.py                   311    199    36%
    nexus/nexus/testing.py                       409     67    84%
    nexus/nexus/unit_converter.py                121      4    97%
    nexus/nexus/vasp.py                           94     15    84%
    nexus/nexus/vasp_analyzer.py                 548     73    87%
    nexus/nexus/vasp_input.py                    906    412    55%
    nexus/nexus/versions.py                      335     50    85%
    nexus/nexus/xmlreader.py                     260     54    79%

The first column is the total number of statements, the second is the number not yet covered by the tests and the third is the percent covered. By averaging the third column, you can tell roughly what percent of Nexus is covered by testing. From the testing shown above, code coverage is at ~67.5%.

To obtain an annotated view of the statements in the source that are not yet covered, run:

.. code-block:: bash

    > pytest --cov=nexus --cov-report html

Open ``htmlcov/index.html`` in a browser to view the report. More
information regarding the ``coverage`` tool can be found at
https://coverage.readthedocs.io/en/v4.5.x/.
