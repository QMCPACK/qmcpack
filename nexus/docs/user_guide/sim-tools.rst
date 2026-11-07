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


``nxs-sim``
-----------

``nxs-sim`` has several subcommands useful for modifying the status of a saved simulation image or creating backups of a simulation image.

All examples here assume the following directory structure:

.. code-block:: text

    .
    └── runs
        ├── water_relax
        │   └── sim_water_relax
        │       └── sim.p
        └── water_scf
            └── sim_water_scf
                └── sim.p

    6 directories, 2 files

Subcommands
^^^^^^^^^^^

``show``
    Show the selected simulation status indicators.

    .. dropdown:: Example
        :color: info
        :icon: code-square

        .. code-block:: bash

            > nxs-sim show -v setup failed finished runs/*/sim_*/sim.p
            Loading sim images:
              ./runs/water_relax/sim_water_relax/sim.p
              ./runs/water_scf/sim_water_scf/sim.p
            Done!
            ./runs/water_relax/sim_water_relax/sim.p
              setup:    True
              failed:   True
              finished: True

            ./runs/water_scf/sim_water_scf/sim.p
              setup:    True
              failed:   True
              finished: True

``showall``
    Show all simulation properties.

    .. dropdown:: Example
        :color: info
        :icon: code-square

        .. code-block:: bash

            > nxs-sim showall -v runs/water_scf/sim_water_scf/sim.p
            Loading sim images:
              ./runs/water_scf/sim_water_scf/sim.p
            Done!
            ./runs/water_scf/sim_water_scf/sim.p
              analyzed        = True
              analyzer_image  = analyzer.p
              app_name        = pw.x
              app_props       = ['mpi', 'serial']
              block           = False
              block_subcascade = False
              errfile         = water_scf.err
              failed          = True
              files           = {...}
              finished        = True
              got_output      = False
              identifier      = water_scf
              image_dir       = sim_water_scf
              imlocdir        = ./runs/water_scf/sim_water_scf
              imremdir        = ./runs/water_scf/sim_water_scf
              imresdir        = ./runs/water_scf/sim_water_scf
              infile          = water_scf.in
              input_image     = input.p
              locdir          = ./runs/water_scf
              outfile         = water_scf.out
              outputs         = None
              path            = water_scf
              process_id      = None
              remdir          = ./runs/water_scf
              resdir          = ./runs/water_scf
              sent_files      = True
              setup           = True
              sim_image       = sim.p
              subcascade_finished = False
              submitted       = False
              timestamps      = {...}

``complete``
    Set simulation status flags to indicate it is completed.

    .. dropdown:: Example
        :color: info
        :icon: code-square

        .. code-block:: bash

            > nxs-sim complete -v runs/water_relax/sim_water_relax/sim.p
            Loading sim images:
              ./runs/water_relax/sim_water_relax/sim.p
            Done!
            ./runs/water_relax/sim_water_relax/sim.p
              setup:      False -> True
              sent_files: False -> True
              analyzed:   False -> True
              submitted:  False -> True
              finished:   False -> True
              got_output: False -> True

``reset``
    Reset all simulation status flags.

    .. dropdown:: Example
        :color: info
        :icon: code-square

        .. code-block:: bash

            > nxs-sim reset -v runs/water_relax/sim_water_relax/sim.p
            Loading sim images:
              ./runs/water_relax/sim_water_relax/sim.p
            Done!
            ./runs/water_relax/sim_water_relax/sim.p
              setup:      True -> False
              sent_files: True -> False
              submitted:  True -> False (and setting sim.process_id to None)
              finished:   True -> False
              failed:     True -> False
              got_output: True -> False
              analyzed:   True -> False

``set``
    Set the specified simulation status flags to ``True``.

    .. dropdown:: Example
        :color: info
        :icon: code-square

        .. code-block:: bash

            > nxs-sim set -v setup sent_files submitted runs/water_*/sim_*/sim.p
            Loading sim images:
              ./runs/water_relax/sim_water_relax/sim.p
              ./runs/water_scf/sim_water_scf/sim.p
            Done!
            ./runs/water_relax/sim_water_relax/sim.p
              setup:      False -> True
              sent_files: False -> True
              submitted:  False -> True
            ./runs/water_scf/sim_water_scf/sim.p
              setup:      False -> True
              sent_files: False -> True
              submitted:  False -> True

``unset``
    Set the specified simulation status flags to ``False``. The opposite of ``set``

    .. dropdown:: Example
        :color: info
        :icon: code-square

        .. code-block:: bash

            > nxs-sim unset -v setup sent_files submitted runs/water_*/sim_*/sim.p
            Loading sim images:
              ./runs/water_relax/sim_water_relax/sim.p
              ./runs/water_scf/sim_water_scf/sim.p
            Done!
            ./runs/water_relax/sim_water_relax/sim.p
              setup:      True -> False
              sent_files: True -> False
              submitted:  True -> False (and setting sim.process_id to None)
            ./runs/water_scf/sim_water_scf/sim.p
              setup:      True -> False
              sent_files: True -> False
              submitted:  True -> False (and setting sim.process_id to None)

``backup``
    Make a backup of the specified simulation pickle at ``sim_backup.p``.

    .. dropdown:: Example
        :color: info
        :icon: code-square

        .. code-block:: bash

            > nxs-sim backup -v runs/water_*/sim_*/sim.p
            Loading sim images:
              ./runs/water_relax/sim_water_relax/sim.p
              ./runs/water_scf/sim_water_scf/sim.p
            Done!
            Creating backup of ./runs/water_relax/sim_water_relax/sim.p
            Creating backup of ./runs/water_scf/sim_water_scf/sim.p
            > tree .
            .
            └── runs
                ├── water_relax
                │   └── sim_water_relax
                │       ├── backup_sim.p
                │       └── sim.p
                └── water_scf
                    └── sim_water_scf
                        ├── backup_sim.p
                        └── sim.p

            6 directories, 4 files

``restore``
    Take a backup of the specified simulation pickle at ``sim_backup.p``
    and copy it back to ``sim.p``.

    .. dropdown:: Example
        :color: info
        :icon: code-square

        .. code-block:: bash

            > tree .
            .
            └── runs
                ├── water_relax
                │   └── sim_water_relax
                │       ├── backup_sim.p
                │       └── sim.p
                └── water_scf
                    └── sim_water_scf
                        ├── backup_sim.p
                        └── sim.p

            6 directories, 4 files
            > nxs-sim restore -v runs/water_*/sim_*/sim.p
            Loading sim images:
              ./runs/water_relax/sim_water_relax/sim.p
              ./runs/water_scf/sim_water_scf/sim.p
            Done!
            Restoring backup of ./runs/water_relax/sim_water_relax/sim.p
            Restoring backup of ./runs/water_scf/sim_water_scf/sim.p
