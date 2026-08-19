.. _nexus_execution_behavior:

Nexus Execution Behavior
========================

This section describes how Nexus behaves during execution (e.g. after calling :py:func:`~.run_project`).
The steps that a simulation goes through are documented here, as well as how Nexus handles errors both during simulation execution and after.


.. _internal_nexus_errors:

Internal Nexus Errors
---------------------

If a simulation fails due to an exception in Nexus, a message will also be printed to the screen and a logfile with the suffix ``.nexus.log`` will be created in the simulation run directory that contains the error message.

.. dropdown:: For Developers
    :color: info
    :margin: auto auto 5 5
    :icon: code-square

    These errors occur most commonly in :py:meth:`~.ProjectManager.progress_cascades`, and are handled through the :py:class:`~.simulation.sim_err_handler` context manager.
