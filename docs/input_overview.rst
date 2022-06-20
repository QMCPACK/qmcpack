.. _input-overview:

Input file overview
===================

This chapter introduces XML as it is used in the QMCPACK input file.  The focus is on the XML file format itself and the general structure of the input file rather than an exhaustive discussion of all keywords and structure elements.

QMCPACK uses XML to represent structured data in its input file.  Instead of text blocks like

::

  begin project
    id     = vmc
    series = 0
  end project

  begin vmc
    move     = pbyp
    blocks   = 200
    steps    =  10
    timestep = 0.4
  end vmc

QMCPACK input looks like

::

  <project id="vmc" series="0">
  </project>

  <qmc method="vmc" move="pbyp">
     <parameter name="blocks"  >  200 </parameter>
     <parameter name="steps"   >   10 </parameter>
     <parameter name="timestep">  0.4 </parameter>
  </qmc>

XML elements start with ``<element_name>``, end with ``</element_name>}``, and can be nested within each other to denote substructure (the trial wavefunction is composed of a Slater determinant and a Jastrow factor, which are each further composed of :math:`...`).  ``id`` and ``series`` are attributes of the ``<project/>`` element.  XML attributes are generally used to represent simple values, like names, integers, or real values.  Similar functionality is also commonly provided by ``<parameter/>`` elements like those previously shown.

The overall structure of the input file reflects different aspects of the QMC simulation: the simulation cell, particles, trial wavefunction, Hamiltonian, and QMC run parameters.  A condensed version of the actual input file is shown as follows:

::

  <?xml version="1.0"?>
  <simulation>

  <project id="vmc" series="0">
    ...
  </project>

  <qmcsystem>

    <simulationcell>
      ...
    </simulationcell>

    <particleset name="e">
      ...
    </particleset>

    <particleset name="ion0">
      ...
    </particleset>

    <wavefunction name="psi0" ... >
      ...
      <determinantset>
        <slaterdeterminant>
          ..
        </slaterdeterminant>
      </determinantset>
      <jastrow type="One-Body" ... >
         ...
      </jastrow>
      <jastrow type="Two-Body" ... >
        ...
      </jastrow>
    </wavefunction>

    <hamiltonian name="h0" ... >
      <pairpot type="coulomb" name="ElecElec" ... />
      <pairpot type="coulomb" name="IonIon"   ... />
      <pairpot type="pseudo" name="PseudoPot" ... >
        ...
      </pairpot>
    </hamiltonian>

   </qmcsystem>

   <qmc method="vmc" move="pbyp">
     <parameter name="warmupSteps">   20 </parameter>
     <parameter name="blocks"     >  200 </parameter>
     <parameter name="steps"      >   10 </parameter>
     <parameter name="timestep"   >  0.4 </parameter>
   </qmc>

  </simulation>

The omitted portions ``...`` are more fine-grained inputs such as the axes of the simulation cell, the number of up and down electrons, positions of atomic species, external orbital files, starting Jastrow parameters, and external pseudopotential files.

Project
-------

The ``<project>`` tag uses the ``id`` and ``series`` attributes.
The value of ``id`` is the first part of the prefix for output file names.

Output file names also contain the series number, starting at the value given by the
``series`` tag.  After every ``<qmc>`` section, the series value will increment, giving each section a unique prefix.

For the input file shown previously, the output files will start with ``vmc.s000``, for example, ``vmc.s000.scalar.dat``.
If there were another ``<qmc>`` section in the input file, the corresponding output files would use the prefix ``vmc.s001``.

The ``<project>`` tag accepts additional control parameters (using the ``<parameters/>`` tag) that can set time limits and specify the driver version.

Time limits
~~~~~~~~~~~
Batched drivers check against ``max_seconds`` and make efforts to stop the execution cleanly at the end of a block before reaching the maximum time. Classic drivers can also take the now-deprecated ``maxcpusecs`` parameter for the same effect in the per driver XML section.

In addition, a file named ``id`` plus ``.STOP``, in this case ``vmc.STOP``, stops QMCPACK execution on the fly cleanly once being found in the working directory.


.. _driver-version-parameter:

Driver version
~~~~~~~~~~~~~~
The ``driver_version`` parameter selects between the new performance-portable batched drivers and the previous drivers (now referred to as the 'legacy drivers').
The values for this parameter are ``legacy`` or ``batch`` (alternately, ``batched``).


Random number initialization
----------------------------

The random number generator state is initialized from the ``random`` element using the ``seed`` attribute:

::

  <random seed="1000"/>

If the random element is not present, or the seed value is negative, the seed will be generated from the current time.

To initialize the many independent random number generators (one per thread and MPI process), the seed value is used (modulo 1024) as a starting index into a list of prime numbers.
Entries in this offset list of prime numbers are then used as the seed for the random generator on each thread and process.

If checkpointing is enabled, the random number state is written to an HDF file at the end of each block (suffix: ``.random.h5``).
This file will be read if the ``mcwalkerset`` tag is present to perform a restart.
For more information, see the ``checkpoint`` element in the QMC methods :ref:`qmcmethods` and :ref:`checkpoint-files` on checkpoint and restart files.
