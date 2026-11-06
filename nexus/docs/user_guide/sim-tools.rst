.. _sim-tools:

Simulation Tools
================

Nexus provides two primary tools for working with a simulation, :py:mod:`nxs-redo` and :py:mod:`nxs-sim`.
These are used to back up and reset an already-completed simulation and modify a simulation's status, respectively.

``nxs-redo``
------------

``nxs-redo`` is used when you want to restart an existing simulation run but keep the old files. Consider the following directory structure:

.. code-block:: text

    nxs-redo-example
    ├── nexus_script.py
    ├── pseudopotentials
    │   ├── h_pbesol_v1.4.uspp.F.UPF
    │   └── o_pbesol_v1.2.uspp.F.UPF
    └── runs
        ├── water_relax
        │   ├── h_pbesol_v1.4.uspp.F.UPF
        │   ├── o_pbesol_v1.2.uspp.F.UPF
        │   ├── relax_output
        │   ├── sim_water_relax
        │   │   ├── input.p
        │   │   └── sim.p
        │   ├── water_relax.in
        │   ├── water_relax.struct.xsf
        │   └── water_relax.struct.xyz
        └── water_scf
            ├── h_pbesol_v1.4.uspp.F.UPF
            ├── o_pbesol_v1.2.uspp.F.UPF
            ├── scf_output
            ├── sim_water_scf
            │   ├── input.p
            │   └── sim.p
            ├── water_scf.in
            ├── water_scf.struct.xsf
            └── water_scf.struct.xyz

    9 directories, 17 files


With this, when we are in the ``nxs-redo-example`` directory, we can use ``nxs-redo`` to "reset" each simulation directory.

Here we will use the ``-n`` option to preview a dry run of ``nxs-redo``:

.. code-block:: bash

    > nxs-redo -n ./runs/*
    Backing up attempted sim runs...

      Backing up files in /tmp/nxs-redo-example/runs/water_relax
        Would have created backup of 7 files in:
          /tmp/nxs-redo-example/runs/water_relax/attempt1

      Backing up files in /tmp/nxs-redo-example/runs/water_scf
        Would have created backup of 7 files in:
          /tmp/nxs-redo-example/runs/water_scf/attempt1

    Backup Complete!

We can confirm that this is the desired behavior, and proceed now with a real run.

.. code-block:: bash

    > nxs-redo -v ./runs/*
    Backing up attempted sim runs...

      Backing up files in /tmp/nxs-redo-example/runs/water_relax
        Creating backup of 7 files in:
          /tmp/nxs-redo-example/runs/water_relax/attempt1
      Done!

      Backing up files in /tmp/nxs-redo-example/runs/water_scf
        Creating backup of 7 files in:
          /tmp/nxs-redo-example/runs/water_scf/attempt1
      Done!

    Backup Complete!
    > tree .
    .
    ├── nexus_script.py
    ├── pseudopotentials
    │   ├── h_pbesol_v1.4.uspp.F.UPF
    │   └── o_pbesol_v1.2.uspp.F.UPF
    └── runs
        ├── water_relax
        │   └── attempt1
        │       ├── h_pbesol_v1.4.uspp.F.UPF
        │       ├── o_pbesol_v1.2.uspp.F.UPF
        │       ├── relax_output
        │       ├── sim_water_relax
        │       │   ├── input.p
        │       │   └── sim.p
        │       ├── water_relax.in
        │       ├── water_relax.struct.xsf
        │       └── water_relax.struct.xyz
        └── water_scf
            └── attempt1
                ├── h_pbesol_v1.4.uspp.F.UPF
                ├── o_pbesol_v1.2.uspp.F.UPF
                ├── scf_output
                ├── sim_water_scf
                │   ├── input.p
                │   └── sim.p
                ├── water_scf.in
                ├── water_scf.struct.xsf
                └── water_scf.struct.xyz

    11 directories, 17 files

Now, we can safely re-run ``nexus_script.py`` without fear of touching the previous run.

.. tip::

    If you want more control over the files given to ``nxs-redo``, you can explicitly list the files you want to back up:

    .. code-block:: bash

        > nxs-redo -v ./runs/water_scf
    
    or you can provide a file with a list of directory names:

    .. code-block:: bash

        > echo "./runs/water_relax" > file.txt
        > nxs-redo -v file.txt


Controlling Simulation Status
-----------------------------

TODO
