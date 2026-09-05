.. _output-overview:

Output Overview
===============

QMCPACK writes several output files that report information about the simulation (e.g., the physical properties such as the energy), as well as information about the computational aspects of the simulation, checkpoints, and restarts.
The types of output files generated depend on the details of a calculation. The following list is not meant to be exhaustive but rather to highlight some salient features of the more common file types. Further details can be found in the description of the estimator of interest.

.. _scalardat-file:

The .scalar.dat file
--------------------

The most important output file is the ``scalar.dat`` file. This file contains the
output of block-averaged properties of the system such as the local
energy and other estimators. Each line corresponds to an average over
:math:`N_{walkers}*N_{steps}` samples. By default, the quantities
reported in the ``scalar.dat`` file include the following:

LocalEnergy
   The local energy.

LocalEnergy_sq
   The local energy squared.

LocalPotential
   The local potential energy.

Kinetic
   The kinetic energy.

ElecElec
   The electron-electron potential energy.

IonIon
   The ion-ion potential energy.

LocalECP
   The energy due to the pseudopotential/effective core potential.

NonLocalECP
   The nonlocal energy due to the pseudopotential/effective core
   potential.

MPC
   The modified periodic Coulomb potential energy.

BlockWeight
   The number of MC samples in the block.

BlockCPU
   The number of seconds to compute the block.

AcceptRatio
   The acceptance ratio.

QMCPACK includes a python utility, ``qmca``, that can be used to process these files. Details and examples are given in :ref:`analyzing`.

.. _optxml-file:

The .opt.xml file
-----------------

This file is generated after a VMC wavefunction optimization and contains the part of the input file that lists the optimized Jastrow factors.
Conveniently, this file is already formatted such that it can easily be incorporated into a DMC input file.

.. _qmc-file:

The .qmc.xml file
-----------------

This file contains information about the computational aspects of the simulation, for example, which parts of the code are being executed when. This file is generated only during an ensemble run in which QMCPACK runs multiple input files.

.. _dmc-file:

The .dmc.dat file
-----------------

This file contains information similar to the ``.scalar.dat`` file but also includes extra information about the details of a DMC calculation, for example, information about the walker population.

Index
   The block number.

LocalEnergy
   The local energy.

Variance
   The variance.

Weight
   The number of samples in the block.

NumOfWalkers
   The number of walkers times the number of steps.

AvgSentWalkers
   The average number of walkers sent. During a DMC simulation, walkers
   might be created or destroyed. At every step, QMCPACK will do some
   load balancing to ensure that the walkers are evenly distributed
   across nodes.

TrialEnergy
   The trial energy. See :ref:`dmc` for an explanation of
   trial energy.

DiffEff
   The diffusion efficiency.

LivingFraction
   The fraction of the walker population from the previous step that
   survived to the current step.

.. _bandinfo-file:

The .bandinfo.dat file
----------------------

This file contains information from the trial wavefunction about the band structure of the system,
including the available :math:`k`-points. This can
be helpful in constructing trial wavefunctions.

.. _checkpoint-files:

Checkpoint and restart files
----------------------------

The .cont.xml file
~~~~~~~~~~~~~~~~~~

This file enables continuation of the run.  It is mostly a copy of the input XML file with the series number incremented and the ``mcwalkerset`` element added to read the walkers from a config file.   The ``.cont.xml`` file is always created, but other files it depends on are  present only if checkpointing is enabled.

The .config.h5 file
~~~~~~~~~~~~~~~~~~~

This file contains stored walker configurations.

::
   
   * GROUP "root"
     * GROUP "state_0"
       * DATASET "block"
         * int
         * SCALAR  
       * DATASET "number_of_walkers"
         * size_t
         * SCALAR
       * DATASET "walker_partition"
         * int
         * ARRAY ( offsets )
       * DATASET "walker_weights"
         * double
         * ARRAY ( weights )
       * DATASET "walkers"
         * double
         * ARRAY ( configurations )
     * DATASET "version"
       * int
       * ARRAY ( major version number, minor version number )

The .random.h5 file
~~~~~~~~~~~~~~~~~~~

This file contains the state of the random number generator to allow restarts.
(Older versions used an XML file with a suffix of ``.random.xml``).
