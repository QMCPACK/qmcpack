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

* The :py:class:`~.pseudoset.PseudoSet` class, which serves as the through-point for all pseudopotential handling in Nexus.
* The :py:class:`~.hdfreader.HDFreader` class, which enables parsing of HDF5 archives through the ``h5py`` package.
* The :py:class:`~.xmlreader.XMLreader` class, which enables parsing of XML files through the ``expat`` XML parser interface in Python.
