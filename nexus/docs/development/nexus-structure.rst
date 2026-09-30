.. _nexus-structure:

Nexus's Structure
=================

The layout of Nexus can be difficult to understand for potential contributors.
This document provides a summary of the overall layout of Nexus, starting with a brief mention of the class structure, then moving to essential calculation components, followed by the call tree for a normal Nexus run.


.. _nexus-class-layout:

Class Structure
---------------

Most Nexus classes inherit, directly or indirectly, from the base class :py:class:`~.DevBase`, which provides a dict-like interface for classes.
One such example is :py:class:`~.NexusCore`, which provides a base for many of the most critical Nexus classes:

* :py:class:`~.machines.Job`
* :py:class:`~.machines.Machine`
* :py:class:`~.simulation.SimulationInput`
* :py:class:`~.simulation.SimulationAnalyzer`
* :py:class:`~.simulation.SimulationImage`
* :py:class:`~.simulation.Simulation`
* :py:class:`~.project_manager.ProjectManager`
* :py:class:`~.project_manager.DynamicWorkflowManager`
* :py:class:`~nexus.Settings`

Some other classes that merit mention are

* :py:class:`~.pseudoset.PseudoSet`, which serves as the through-point for all pseudopotential handling in Nexus.
* :py:class:`~.hdfreader.HDFreader`, which enables parsing of HDF5 archives through the ``h5py`` package.
* :py:class:`~.xmlreader.XMLreader`, which enables parsing of XML files through the ``expat`` XML parser interface in Python.
* :py:class:`~.NexusConfig` is where all of Nexus's core configuration variables are stored.


.. _nexus-calculation:

Essential Calculation Components
--------------------------------

This section provides an overview on some of the critical components of any Nexus script.

User Settings
^^^^^^^^^^^^^

User settings are parsed through :py:class:`~nexus.Settings`, and are both assigned to its singleton instance ``settings`` and assigned to the relevant classes. Any valid configuration variable for a simulation class or for :py:class:`~.NexusConfig` are able to be passed to ``settings``. Some variables can be passed through the command line, in which case they will silently replace the values in the script's call to ``settings``.


Job Objects
^^^^^^^^^^^

Every Nexus simulation requires an instance of the :py:class:`~.Job` class, which controls many aspects of Nexus's execution. Some of the notable things this class handles are setting the executable command for the simulation, setting the number of nodes, cores, threads, walltime allocation, queue partition, and account. A complete listing of its parameters/attributes are in its documentation.

:py:class:`~.Job` objects work closely with the user's specified machine, which is usually either a :py:class:`~.Workstation` instance or an instance of a subclass of :py:class:`~.Supercomputer`. These classes work together to create job submission files and ensure that the job that a user requests is actually possible on their machine. The :py:class:`~.Machine` class (parent class of :py:class:`~.Workstation` and :py:class:`~.Supercomputer`) also defines the ``query_queue`` abstract method, which subclasses implement so that Nexus is aware of when a job stops running or leaves the queue.


.. _nexus-call-tree:

Nexus Call Tree
---------------

.. admonition:: Definitions

    * **Cascade**: An instance of a :py:class:`~.Simulation` subclass.

This is broken down for two kinds of runs, flat runs and nested runs. A flat run is where a user is running one or more independent simulations, and a nested run is with one or more dependent simulations.

Common to all runs (excluding the actively-developed dynamic processes) is the call to :py:func:`~nexus.run_project`. This function serves one purpose, to create a :py:class:`~.ProjectManager` instance, add simulations to it, and call :py:meth:`.ProjectManager.run_project`.

Flat Runs
^^^^^^^^^

The main way a user might create this sort of run is by placing their physical system/simulation generation in a ``for`` loop that iterates over structure files and creates a large number of independent simulations. For these runs, all of the simulations are managed by the :py:class:`~.ProjectManager`.

The "main loop" of Nexus, where the actual steps of a simulation are driven, is in :py:meth:`.ProjectManager.run_project`, and starts with a call to :py:meth:`.ProjectManager.init_cascades`, which screens for fake simulations (usually encountered in testing), checks for file collisions, propagates any blockages for dependent simulations (usually not encountered in flat runs), and gets any relevant dependencies (almost always encountered on restarts of Nexus workflows).

The project manager will then run through the following steps:

#. Log the start time of the Nexus run. (only if :py:attr:`.NexusConfig.monitor` is ``True``)
#. Query the queue (usually does nothing on the first pass)
#. Progress the simulations
#. Submit any jobs that are ready to go
#. Update the process IDs

If :py:attr:`.NexusConfig.monitor` is ``False``, Nexus stops and does not continue to monitor jobs. More commonly, however, is when :py:attr:`.NexusConfig.monitor` is ``True``, in which case the project manager will sleep for the duration of :py:attr:`.NexusConfig.sleep`, and repeat the steps above. For flat runs, this is the only spot that simulations are actually driven.

A good way to tell what has happened inside this loop is to look for any text that is between two of the polling messages.
