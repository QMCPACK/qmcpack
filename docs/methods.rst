.. _qmcmethods:

Quantum Monte Carlo Methods
===========================

``qmc`` factory element:

  +-----------------+----------------------+
  | Parent elements | ``simulation, loop`` |
  +-----------------+----------------------+
  | type selector   | ``method`` attribute |
  +-----------------+----------------------+

  type options:

  +--------+-----------------------------------------------+
  | vmc    | Variational Monte Carlo                       |
  +--------+-----------------------------------------------+
  | linear | Wavefunction optimization with linear method  |
  +--------+-----------------------------------------------+
  | dmc    | Diffusion Monte Carlo                         |
  +--------+-----------------------------------------------+
  | rmc    | Reptation Monte Carlo                         |
  +--------+-----------------------------------------------+

  shared attributes:

  +----------------+--------------+--------------+-------------+---------------------------------+
  | **Name**       | **Datatype** | **Values**   | **Default** | **Description**                 |
  +================+==============+==============+=============+=================================+
  | ``method``     | text         | listed above | invalid     | QMC driver                      |
  +----------------+--------------+--------------+-------------+---------------------------------+
  | ``move``       | text         | pbyp, alle   | pbyp        | Method used to move electrons   |
  +----------------+--------------+--------------+-------------+---------------------------------+
  | ``gpu``        | text         | yes/no       | dep.        | Use the GPU                     |
  +----------------+--------------+--------------+-------------+---------------------------------+
  | ``trace``      | text         |              | no          | ???                             |
  +----------------+--------------+--------------+-------------+---------------------------------+
  | ``checkpoint`` | integer      | -1, 0, n     | -1          | Checkpoint frequency            |
  +----------------+--------------+--------------+-------------+---------------------------------+
  | ``record``     | integer      | n            | 0           | Save configuration ever n steps |
  +----------------+--------------+--------------+-------------+---------------------------------+
  | ``target``     | text         |              |             | ???                             |
  +----------------+--------------+--------------+-------------+---------------------------------+
  | ``completed``  | text         |              |             | ???                             |
  +----------------+--------------+--------------+-------------+---------------------------------+
  | ``append``     | text         | yes/no       | no          | ???                             |
  +----------------+--------------+--------------+-------------+---------------------------------+

Additional information:

-  ``move``: There are two ways to move electrons. The more used method
   is the particle-by-particle move. In this method, only one electron
   is moved for acceptance or rejection. The other method is the
   all-electron move; namely, all the electrons are moved once for
   testing acceptance or rejection.

-  ``gpu``: When the executable is compiled with CUDA, the target
   computing device can be chosen by this switch. With a regular
   CPU-only compilation, this option is not effective.

-  ``checkpoint``: This enables and disables checkpointing and
   specifying the frequency of output. Possible values are:

   - **[-1]** No checkpoint (default setting).

   - **[0]** Dump after the completion of a QMC section.

   - **[n]** Dump after every :math:`n` blocks.  Also dump at the end of the run.

The particle configurations are written to a ``.config.h5`` file.

.. code-block::
  :caption: The following is an example of running a simulation that can be restarted.
  :name: Listing 42

  <qmc method="dmc" move="pbyp"  checkpoint="0">
    <parameter name="timestep">         0.004  </parameter>
    <parameter name="blocks">           100   </parameter>
    <parameter name="steps">            400    </parameter>
  </qmc>

The checkpoint flag instructs QMCPACK to output walker configurations.  This also
works in VMC.  This outputs an h5 file with the name ``projectid.run-number.config.h5``.
Check that this file exists before attempting a restart.

To continue a run, specify the ``mcwalkerset`` element before your VMC/DMC block:

.. code-block::
  :caption: Restart (read walkers from previous run).
  :name: Listing 43

  <mcwalkerset fileroot="BH.s002" version="0 6" collected="yes"/>
   <qmc method="dmc" move="pbyp"  checkpoint="0">
     <parameter name="timestep">         0.004  </parameter>
     <parameter name="blocks">           100   </parameter>
     <parameter name="steps">            400    </parameter>
   </qmc>

``BH`` is the project id, and ``s002`` is the calculation number to read in the walkers from the previous run.

In the project id section, make sure that the series number is different from any existing ones to avoid overwriting them.

.. _vmc:

Variational Monte Carlo
-----------------------
