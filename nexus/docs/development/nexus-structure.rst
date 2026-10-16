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
    * **Parent**: A simulation which a child simulation depends on.
    * **Child**: A simulation that depends on a parent simulation.

This is broken down for two kinds of runs, flat runs and nested runs. A flat run is where a user is running one or more independent simulations, and a nested run is with one or more dependent simulations.

Common to all runs (excluding the actively-developed dynamic processes) is the call to :py:func:`~nexus.run_project`. This function serves one purpose, to create a :py:class:`~.ProjectManager` instance, add simulations to it, and call :py:meth:`.ProjectManager.run_project`.

Flat Runs
^^^^^^^^^

The main way a user might create this sort of run is by placing their physical system/simulation generation in a ``for`` loop that iterates over structure files and creates a large number of independent simulations. For these runs, all of the simulations are managed by the :py:class:`~.ProjectManager`.

The "main loop" of Nexus, where the actual steps of a simulation are driven, is in :py:meth:`.ProjectManager.run_project`, and starts with a call to :py:meth:`.ProjectManager.init_cascades`, which screens for fake simulations (usually encountered in testing), checks for file collisions, propagates any blockages for dependent simulations (usually not encountered in flat runs), and gets any relevant dependencies (almost always encountered on restarts of Nexus workflows).

The project manager will then run through the following steps:

1. Log the start time of the Nexus run. (only if :py:attr:`.NexusConfig.monitor` is ``True``)
2. Query the queue (usually does nothing on the first pass)
3. Progress the simulations
4. Submit any jobs that are ready to go
5. Update the process IDs

If :py:attr:`.NexusConfig.monitor` is ``False``, Nexus stops and does not continue to monitor jobs. More commonly, however, is when :py:attr:`.NexusConfig.monitor` is ``True``, in which case the project manager will sleep for the duration of :py:attr:`.NexusConfig.sleep`, and repeat the steps above. For flat runs, this is the only spot that simulations are actually driven.

A good way to tell what has happened inside this loop is to look for any text that is between two of the polling messages.

The 3rd step of the project manager's flow is, as noted, simulation progression. This is done via :py:meth:`.Simulation.progress`, which controls the flow of essentially all simulations in Nexus.

In flat runs, this function serves one purpose, which is to go through the various simulation stages, most of which are enumerated in :py:class:`~.nexus_base.SimStage`. Those which are not will always happen for non-blocked simulations, which will be marked below in **bold**.

**Create Directories**
""""""""""""""""""""""
:py:meth:`.Simulation.create_directories` will create both the simulation directory at :py:attr:`~.Simulation.locdir` and the simulation's image directory at :py:attr:`~.Simulation.imlocdir`.


**Get Dependencies**
""""""""""""""""""""
:py:meth:`.Simulation.get_dependencies` will search through the parent simulations that the child simulation depends on, then call the abstract method :py:meth:`.Simulation.get_result` for each relevant dependency. In flat runs, this is a no-op except for setting :py:attr:`.Simulation.got_dependencies` to ``True``.


**Write Inputs**
""""""""""""""""
:py:meth:`.Simulation.write_inputs` will call up the simulation's input class (an instance of :py:class:`.SimulationInput` defined by the subclass's :py:attr:`~.Simulation.input_type`) which will write the input file for the simulation. :py:meth:`.Simulation.write_inputs` will also attempt to write an XYZ and XSF file containing the system's structural information, but will never fail on these write attempts.


**Send Files**
""""""""""""""
:py:meth:`.Simulation.send_files` is often the last step before a simulation is truly setup and ready to be submitted. This function loops over :py:attr:`.Simulation.files`, searches for them in :py:attr:`.Simulation.locdir` and :py:attr:`.NexusConfig.file_locations`, then calls :py:meth:`.Simulation.copy_file` which will copy the file *only if it does not already exist in the correct location*.


**Submit Job**
""""""""""""""
:py:meth:`.Simulation.submit` serves a dual purpose. Its primary purpose is, as the name suggests, to submit jobs. This is done by calling :py:meth:`.Job.submit` for queued jobs, or :py:meth:`.Simulation.execute` for local workstation jobs.

If you are unfamiliar with Nexus, the :py:meth:`.Job.submit` command can appear confusing at first. Here is a breakdown of how it works:

1. :py:meth:`.Job.submit` gets the current machine.
2. :py:meth:`.Job.submit` calls the :py:meth:`.Machine.add_job`.
3. :py:meth:`.Machine.add_job` calls the subclass's implementation of :py:meth:`.Machine.write_job`.
4. For a supercomputer, this is :py:meth:`.Supercomputer.write_job`, which writes the job's submission file. For a workstation, this is :py:meth:`.Workstation.write_job`, which simply returns the output of :py:meth:`.Workstation.job_command`.
5. The :py:attr:`.Job.internal_id` is added to :py:attr:`.Machine.jobs` and also added to :py:attr:`.Machine.waiting`, which is the set of jobs waiting to be submitted.
   
   .. note::
        These jobs are not immediately submitted, even though :py:attr:`.Simulation.submitted` is set to ``True``

6. After the current call to :py:meth:`.ProjectManager.progress_cascades` is done, the project manager will call the relevant implementation of the :py:meth:`.Machine.submit_jobs` abstract method.

.. danger::
    Any errors that occur between the time that :py:meth:`.Simulation.submit` and :py:meth:`.Machine.submit_jobs` can result in a very confusing state, so it is critical that this flow is extremely robust and not interrupted at any point.

As mentioned earlier, :py:meth:`.Simulation.submit` serves two purposes. The first, described above, is submitting jobs, but the second is to check for job completion. This is done via a call to :py:meth:`.Simulation.check_status`, which goes through the following steps:

.. mermaid::

    flowchart TD
        A[check_status] --> B[pre_check_status]
        B --> C{"`self.job.finished`"}
        C -->|"`True`"| D["Output and error file exist?"]
        C -->|"`False`"| I
        D -->|"Both files exist"| F["Call `self.check_sim_status()`"]
        D -->|"One of the files is missing"| G["Log time of queue exit, see if job has reached timeout"]
        G -->|"Simulation has timed out"| H["Mark sim as failed"]
        H & F --> I{"Sim failed?"}
        I -->|"True"| J["Sim is finished"]
        I -->|"False"| K["Sim may be finished"]
        J & K -->|"Sim has failed"| L["Log failed timestamp"]
        J & K & L --> M["Finished or newly exited queue?"]
        M -->|"Finished"| N["Record finished timestamp, save image"]
        M -->|"Newly exited queue"| O["Save image"]
        N & O --> P["Return"]


**Get Output**
""""""""""""""
:py:meth:`.Simulation.get_output` is designed for those who wish to store the results of their simulations in an alternative directory. Whether or not this happens is dictated by the identity of :py:attr:`.NexusConfig.results`; if it is an empty string, the results are not moved, otherwise they are moved to :py:attr:`.Simulation.resdir`, and the simulation's image is moved to :py:attr:`.Simulation.imresdir`.


**Analyze Output**
""""""""""""""""""
Simulation analysis is a critical part of Nexus, but the actual :py:meth:`.Simulation.analyze` method is largely just a convenient way to incorporate the functionality of a subclass :py:meth:`.SimulationAnalyzer.analyze` into the main Nexus workflow. For a more in-depth description of this functionality, you will need to look at the analyzer subclasses themselves.






