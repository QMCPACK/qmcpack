.. _additional-tools:

Additional Tools
================

QMCPACK provides a set of lightweight executables that address certain
common problems in QMC workflow and analysis.  These range from conversion utilities between
different file formats and QMCPACK (e.g., ``ppconvert`` and ``convert4qmc``),
(qmc-extract-eshdf-kvectors) to postprocessing utilities (``trace-density`` and ``qmcfinitesize``) to many others.  In this section, we cover the use cases, syntax, and features of all additional tools provided with QMCPACK.

Initialization
--------------

qmc-get-supercell
~~~~~~~~~~~~~~~~~

Postprocessing
--------------

qmca
~~~~

``qmca`` is a versatile tool to analyze and plot the raw data from QMCPACK ``*.scalar.dat`` files.
It is a Python executable and part of the Nexus suite of tools.  It can be found in
``qmcpack/nexus/executables``. For details, see :ref:`qmca`.

qmc-fit
~~~~~~~

``qmc-fit`` is a curve fitting tool used to obtain statistical error bars on fitted parameters.
It is useful for DMC time step extrapolation.  For details, see :ref:`qmcfit`.

qdens
~~~~~

``qdens`` is a command line tool to produce density files from QMCPACK's ``stat.h5`` output files.  For details, see :ref:`qdens`.

qmcfinitesize
~~~~~~~~~~~~~

``qmcfinitesize`` is a utility to compute many-body finite-size corrections to the energy.  It
is a C++ executable that is built alongside the QMCPACK executable.  It can be found in
``build/bin``.

Converters
----------

.. _convert4qmc:

convert4qmc
~~~~~~~~~~~

``Convert4qmc`` allows conversion of orbitals and wavefunctions from
quantum chemistry output files to QMCPACK XML and HDF5 input files.
It is a small C++ executable that is built alongside the QMCPACK
executable and can be found in ``build/bin``.

To date, ``convert4qmc`` supports the following codes:
GAMESS :cite:`schmidt93`, PySCF :cite:`Sun2018` and QP2 :cite:`QP2` natively, and NWCHEM :cite:`NWCHEM`, TURBOMOLE :cite:`TURBOMOLE`, PSI4 :cite:`PSI4`, CFOUR 2.0beta :cite:`CFOUR`, ORCA 3.X - 4.X :cite:`ORCA`, DALTON2016 :cite:`DALTON2016`, MOLPRO :cite:`MOLPRO`, DIRAC :cite:`DIRAC`, RMG :cite:`RMG`, and QCHEM 4.X :cite:`QCHEM` through the molden2qmc converter (see :ref:`molden2qmc`).



General use
^^^^^^^^^^^

General use of ``convert4qmc`` can be prompted by running with no options:

::

  >convert4qmc

  Defaults : -gridtype log -first 1e-6 -last 100 -size 1001 -ci required -threshold 0.01 -TargetState 0 -prefix sample

   convert [-gaussian|gamess|-orbitals|-dirac|-rmg]
   filename
  [-nojastrow -hdf5 -prefix title -addCusp -production -NbImages NimageX NimageY NimageZ]
  [-psi_tag psi0 -ion_tag ion0 -gridtype log|log0|linear -first ri -last rf]
  [-size npts -ci file.out -threshold cimin -TargetState state_number
  -NaturalOrbitals NumToRead -optDetCoeffs]
  Defaults : -gridtype log -first 1e-6 -last 100 -size 1001 -ci required
  -threshold 0.01 -TargetState 0 -prefix sample
  When the input format is missing, the  extension of filename is used to determine
  the format
   *.Fchk -> gaussian; *.out -> gamess; *.h5 -> hdf5 format

As an example, to convert a GAMESS calculation using a single determinant, the following use is sufficient:

::

  convert4qmc -gamess MyGamessOutput.out

By default, the converter will generate multiple files:

  ``convert4qmc`` output:

    +-------------------------+---------------+-------------+----------------------------------------------------+
    | **output**              | **file type** | **default** | **description**                                    |
    +=========================+===============+=============+====================================================+
    | ``*.qmc.in-wfs.xml``    | XML           | default     | Main input file for QMCPACK                        |
    +-------------------------+---------------+-------------+----------------------------------------------------+
    | ``*.qmc.in-wfnoj.xml``  | XML           | default     | Main input file for QMCPACK                        |
    +-------------------------+---------------+-------------+----------------------------------------------------+
    | ``*.structure.xml``     | XML           | default     | File containing the structure of the system        |
    +-------------------------+---------------+-------------+----------------------------------------------------+
    | ``*.wfj.xml``           | XML           | default     | Wavefunction file with 1-, 2-, and 3-body Jastrows |
    +-------------------------+---------------+-------------+----------------------------------------------------+
    | ``*.wfnoj.xml``         | XML           | default     | Wavefunction file with no Jastrows                 |
    +-------------------------+---------------+-------------+----------------------------------------------------+
    | ``*.orbs.h5``           | HDF5          | with -hdf5  | HDF5 file containing all wavefunction data         |
    +-------------------------+---------------+-------------+----------------------------------------------------+

If no ``-prefix`` option is specified, the prefix is taken from
the input file name. For instance, if the GAMESS output file is
``Mysim``.out, the files generated by ``convert4qmc`` will use the
prefix ``Mysim`` and output files will be
``Mysim.qmc.in-wfs.xml``, ``Mysim.structure.xml``, and so on.

- Files ``.in-wfs.xml`` and ``.in-wfnoj.xml``

  These
  are the input files for QMCPACK.  The geometry and the
  wavefunction are stored in external files ``*.structure.xml``
  and ``*.wfj.xml`` (referenced from ``*.in-wfs.xml``) or
  ``*.qmc.wfnoj.xml`` (referenced from
  ``*.qmc.in-wfnoj.xml``). The Hamiltonian section is included,
  and the presence or lack of presence of an ECP is detected during the
  conversion. If use of an ECP is detected, a default ECP name is
  added (e.g., ``H.qmcpp.xml``), and it is the responsibility of
  the user to modify the ECP name to match the one used to generate
  the wavefunction.

  ::

      <?xml version="1.0"?>
    <simulation>
      <!--

    Example QMCPACK input file produced by convert4qmc

    It is recommend to start with only the initial VMC block and adjust
    parameters based on the measured energies, variance, and statistics.

    -->
      <!--Name and Series number of the project.-->
      <project id="gms" series="0"/>
      <!--Link to the location of the Atomic Coordinates and the location of
          the Wavefunction.-->
      <include href="gms.structure.xml"/>
      <include href="gms.wfnoj.xml"/>
      <!--Hamiltonian of the system. Default ECP filenames are assumed.-->
      <hamiltonian name="h0" type="generic" target="e">
        <pairpot name="ElecElec" type="coulomb" source="e" target="e"
                                                       physical="true"/>
        <pairpot name="IonIon" type="coulomb" source="ion0" target="ion0"/>
        <pairpot name="PseudoPot" type="pseudo" source="ion0" wavefunction="psi0"
                                                               format="xml">
          <pseudo elementType="H" href="H.qmcpp.xml"/>
          <pseudo elementType="Li" href="Li.qmcpp.xml"/>
        </pairpot>
      </hamiltonian>

    The ``qmc.in-wfnoj.xml`` file will have one VMC block with a
    minimum number of blocks to reproduce the HF/DFT energy used to
    generate the trial wavefunction.

    ::

        <qmc method="vmc" move="pbyp" checkpoint="-1">
          <estimator name="LocalEnergy" hdf5="no"/>
          <parameter name="warmupSteps">100</parameter>
          <parameter name="blocks">20</parameter>
          <parameter name="steps">50</parameter>
          <parameter name="substeps">8</parameter>
          <parameter name="timestep">0.5</parameter>
          <parameter name="usedrift">no</parameter>
        </qmc>
      </simulation>

  If the ``qmc.in-wfj.xml`` file is used, Jastrow optimization
  blocks followed by a VMC and DMC block are included. These blocks
  contain default values to allow the user to test the accuracy of a
  system; however, they need to be updated and optimized for each
  system. The initial values might only be suitable for a small molecule.

  ::

      <loop max="4">
        <qmc method="linear" move="pbyp" checkpoint="-1">
          <estimator name="LocalEnergy" hdf5="no"/>
          <parameter name="warmupSteps">100</parameter>
          <parameter name="blocks">20</parameter>
          <parameter name="timestep">0.5</parameter>
          <parameter name="walkers">1</parameter>
          <parameter name="samples">16000</parameter>
          <parameter name="substeps">4</parameter>
          <parameter name="usedrift">no</parameter>
          <parameter name="MinMethod">OneShiftOnly</parameter>
          <parameter name="minwalkers">0.0001</parameter>
        </qmc>
      </loop>
      <!--

    Example follow-up VMC optimization using more samples for greater accuracy:

    -->
      <loop max="10">
        <qmc method="linear" move="pbyp" checkpoint="-1">
          <estimator name="LocalEnergy" hdf5="no"/>
          <parameter name="warmupSteps">100</parameter>
          <parameter name="blocks">20</parameter>
          <parameter name="timestep">0.5</parameter>
          <parameter name="walkers">1</parameter>
          <parameter name="samples">64000</parameter>
          <parameter name="substeps">4</parameter>
          <parameter name="usedrift">no</parameter>
          <parameter name="MinMethod">OneShiftOnly</parameter>
          <parameter name="minwalkers">0.3</parameter>
        </qmc>
      </loop>
      <!--

    Production VMC and DMC:

    Examine the results of the optimization before running these blocks.
    For example, choose the best optimized jastrow from all obtained, put in the
    wavefunction file, and do not reoptimize.

    -->
      <qmc method="vmc" move="pbyp" checkpoint="-1">
        <estimator name="LocalEnergy" hdf5="no"/>
        <parameter name="warmupSteps">100</parameter>
        <parameter name="blocks">200</parameter>
        <parameter name="steps">50</parameter>
        <parameter name="substeps">8</parameter>
        <parameter name="timestep">0.5</parameter>
        <parameter name="usedrift">no</parameter>
        <!--Sample count should match targetwalker count for
          DMC. Will be obtained from all nodes.-->
        <parameter name="samples">16000</parameter>
      </qmc>
      <qmc method="dmc" move="pbyp" checkpoint="20">
        <estimator name="LocalEnergy" hdf5="no"/>
        <parameter name="targetwalkers">16000</parameter>
        <parameter name="reconfiguration">no</parameter>
        <parameter name="warmupSteps">100</parameter>
        <parameter name="timestep">0.005</parameter>
        <parameter name="steps">100</parameter>
        <parameter name="blocks">100</parameter>
        <parameter name="nonlocalmoves">yes</parameter>
      </qmc>
    </simulation>

- File ``.structure.xml``

  This file will be referenced from the main QMCPACK input. It contains the geometry of the system, position of the atoms, number of atoms, atomic types and charges, and number of electrons.

- Files ``.wfj.xml`` and ``.wfnoj.xml``

  These files contain the basis set detail, orbital coefficients, and
  the 1-, 2-, and 3-body Jastrow (in the case of ``.wfj.xml``). If the
  wavefunction is multideterminant, the expansion will be at the end of
  the file. We recommend using the option ``-hdf5`` when large molecules
  are studied to store the data more compactly in an HDF5 file.

- File ``.orbs.h5``
  This file is generated only if the option ``-hdf5`` is added as
  follows:

  ::

    convert4qmc -gamess MyGamessOutput.out -hdf5

  In this case, the ``.wfj.xml`` or ``.wfnoj.xml`` files will point to
  this HDF file. Information about the basis set, orbital coefficients,
  and the multideterminant expansion is put in this file and removed from
  the wavefunction files, making them smaller.

``convert4qmc`` input type:

  +-----------------+----------------------------------------------------------------------------+
  | **option name** | **description**                                                            |
  +=================+============================================================================+
  | ``-orbitals``   | Generic HDF5 input file. Mainly automatically generated from QP2, Pyscf and|
  |                 | all codes  in molden2qmc                                                   |
  +-----------------+----------------------------------------------------------------------------+
  | ``-gamess``     | Gamess code                                                                |
  +-----------------+----------------------------------------------------------------------------+
  | ``-gaussian``   | Gaussian code                                                              |
  +-----------------+----------------------------------------------------------------------------+
  | ``-dirac``      | get spinors from DIRAC code                                                |
  +-----------------+----------------------------------------------------------------------------+
  | ``-rmg``        | RMG code                                                                   |
  +-----------------+----------------------------------------------------------------------------+

Command line options
^^^^^^^^^^^^^^^^^^^^

  ``convert4qmc`` command line options:

    +-----------------+-----------+-------------+--------------------------------------------------------------+
    | **Option Name** | **Value** | **default** | **description**                                              |
    +=================+===========+=============+==============================================================+
    | ``-nojastrow``  | -         | -           | Force no Jastrow. ``qmc.in.wfj`` will not be generated       |
    +-----------------+-----------+-------------+--------------------------------------------------------------+
    | ``-hdf5``       | -         | -           | Force the wf to be in HDF5 format                            |
    +-----------------+-----------+-------------+--------------------------------------------------------------+
    | ``-prefix``     | string    | -           | All created files will have the name of the string           |
    +-----------------+-----------+-------------+--------------------------------------------------------------+
    | ``-multidet``   | string    | -           | HDF5 file containing a multideterminant expansion            |
    +-----------------+-----------+-------------+--------------------------------------------------------------+
    | ``-addCusp``    | -         | -           | Force to add orbital cusp correction (ONLY for all-electron) |
    +-----------------+-----------+-------------+--------------------------------------------------------------+
    | ``-production`` | -         | -           | Generates specific blocks in the input                       |
    +-----------------+-----------+-------------+--------------------------------------------------------------+
    | ``-psi_tag``    | string    | psi0        | Name of the electrons particles inside QMCPACK               |
    +-----------------+-----------+-------------+--------------------------------------------------------------+
    | ``-ion_tag``    | string    | ion0        | Name of the ion particles inside QMCPACK                     |
    +-----------------+-----------+-------------+--------------------------------------------------------------+

- ``-multidet``

  This option is to be used when a multideterminant expansion (mainly a CI expansion) is present in an HDF5 file. The trial wavefunction file will not display the full list of multideterminants and will add a path to the HDF5 file as follows (full example for the C2 molecule in qmcpack/tests/molecules/C2_pp).

  ::

    <?xml version="1.0"?>
    <qmcsystem>
      <wavefunction name="psi0" target="e">
        <determinantset type="MolecularOrbital" name="LCAOBSet" source="ion0" transform="yes" href="C2.h5">
          <sposet basisset="LCAOBSet" name="spo-up" size="58">
            <occupation mode="ground"/>
            <coefficient size="58" spindataset="0"/>
          </sposet>
          <sposet basisset="LCAOBSet" name="spo-dn" size="58">
            <occupation mode="ground"/>
            <coefficient size="58" spindataset="0"/>
          </sposet>
          <multideterminant optimize="no" spo_up="spo-up" spo_dn="spo-dn">
            <detlist size="202" type="DETS" nca="0" ncb="0" nea="4" neb="4" nstates="58" cutoff="1e-20" href="C2.h5"/>
          </multideterminant>
        </determinantset>
      </wavefunction>
    </qmcsystem>

  To generate such trial wavefunction, the converter has to be invoked as follows:

  ::

    > convert4qmc -orbitals C2.h5 -multidet C2.h5

- ``-nojastrow``

  This option generates only an input file, ``*.qmc.in.wfnoj.xml``, containing no Jastrow optimization blocks and references a wavefunction file, ``*.wfnoj.xml``, containing no Jastrow section.

- ``-hdf5``

  This option generates the ``*.orbs.h5`` HDF5 file containing the basis set and the orbital coefficients. If the wavefunction contains a multideterminant expansion from QP2, it will also be stored in this file. This option minimizes the size of the ``*.wfj.xml`` file, which points to the HDF file, as in the following example:

  ::

      <?xml version="1.0"?>
     <qmcsystem>
       <wavefunction name="psi0" target="e">
         <determinantset type="MolecularOrbital" name="LCAOBSet" source="ion0"
            transform="yes" href="test.orbs.h5">
           <slaterdeterminant>
             <determinant id="updet" size="39">
               <occupation mode="ground"/>
               <coefficient size="411" spindataset="0"/>
             </determinant>
             <determinant id="downdet" size="35">
               <occupation mode="ground"/>
               <coefficient size="411" spindataset="0"/>
             </determinant>
           </slaterdeterminant>
         </determinantset>
       </wavefunction>
     </qmcsystem>

  Jastrow functions will be included if the option "-nojastrow" was
  not specified. Note that when initially optimization a wavefunction, we recommend
  temporarily removing/disabling the 3-body Jastrow.

- **-prefix**

  Sets the prefix for all output generated by ``convert4qmc``.
  If not specified, ``convert4qmc`` will use the defaults for the
  following:

  -  **Gamess** If the Gamess output file is named “**Name**.out” or
     “**Name**.output,” all files generated by ``convert4qmc`` will carry
     **Name** as a prefix (i.e., **Name**.qmc.in.xml).

  -  **Generic HDF5 input** If a generic HDF5 file is named “**Name**.H5,” all files generated by
     ``convert4qmc`` will carry **Name** as a prefix (i.e.,
     **Name**.qmc.in.xml).

- **-addCusp**

  This option is very important for all-electron (AE) calculations. In
  this case, orbitals have to be corrected for the electron-nuclear
  cusp. The cusp correction scheme follows the algorithm described by Ma
  et al. :cite:`Ma2005` When this option is present, the
  wavefunction file has a new set of tags:

  ::

    qmcsystem>
     <wavefunction name="psi0" target="e">
       <determinantset type="MolecularOrbital" name="LCAOBSet" source="ion0"
         transform="yes" cuspCorrection="yes">
         <basisset name="LCAOBSet">

  The tag “cuspCorrection” in the ``wfj.xml`` (or ``wfnoj.xml``)
  wavefunction file will force correction of the orbitals at the
  beginning of the run.
  In the “orbitals“ section of the wavefunction file, a new tag
  “cuspInfo” will be added for orbitals spin-up and orbitals spin-down:

  ::

      <slaterdeterminant>
           <determinant id="updet" size="2"
               cuspInfo="../updet.cuspInfo.xml">
             <occupation mode="ground"/>
             <coefficient size="135" id="updetC">

     <determinant id="downdet" size="2"
              cuspInfo="../downdet.cuspInfo.xml">
             <occupation mode="ground"/>
             <coefficient size="135" id="downdetC">

  These tags will point to the files ``updet.cuspInfo.xml`` and
  ``downdet.cuspInfo.xml``. By default, the converter assumes that
  the files are located in the relative path
  ``../``. If the files are not
  present in the parent directory, QMCPACK will run the cusp correction
  algorithm to generate both files in the current run directory (not in ``../``).  If the files exist, then QMCPACK
  will apply the corrections to the orbitals.

  **Important notes:**

  The cusp correction implementations has been parallelized and performance improved.  However, since the correction needs
  to be applied for every ion and then for every orbital on that ion, this operation can be costly (slow) for large
  systems. We recommend saving and reusing the computed cusp correction files ``updet.cuspInfo.xml`` and
  ``downdet.cuspInfo.xml``, and transferring them between computer systems where relevant.

- **-psi_tag**

  QMCPACK builds the wavefunction as a named object. In the vast majority of cases, one wavefunction is simulated at a time, but there may be situations where we want to distinguish different parts of a wavefunction, or even use multiple wavefunctions. This option can change the name for these cases.

  ::

     <wavefunction name="psi0" target="e">

- **-ion_tag**

  Although similar to **-psi_tag**, this is used for the type of ions.

  ::

    <particleset name="ion0" size="2">

- **-production**

  Without this option, input files with standard optimization, VMC, and
  DMC blocks are generated. When the "-production" option is
  specified, an input file containing complex options that may be
  more suitable for large runs at HPC centers is generated. This option
  is for users who are already familiar with QMC and QMCPACK. We encourage feedback
  on the standard and production sample inputs.

The following options are specific to using MCSCF multideterminants from Gamess.

  ``convert4qmc`` MCSCF arguments:

  +----------------------+-----------+-------------+----------------------------------------------+
  | **Option Name**      | **Value** | **default** | **description**                              |
  +======================+===========+=============+==============================================+
  | ``-ci``              | String    | none        | Name of the file containing the CI expansion |
  +----------------------+-----------+-------------+----------------------------------------------+
  | ``-threshold``       | double    | 1e-20       | Cutoff of the weight of the determinants     |
  +----------------------+-----------+-------------+----------------------------------------------+
  | ``-TargetState``     | int       | none        | ?                                            |
  +----------------------+-----------+-------------+----------------------------------------------+
  | ``-NaturalOrbitals`` | int       | none        | ?                                            |
  +----------------------+-----------+-------------+----------------------------------------------+
  | ``-optDetCoeffs``    | -         | no          | Enables the optimization of CI coefficients  |
  +----------------------+-----------+-------------+----------------------------------------------+

-  keyword **-ci** Path/name of the file containing the CI expansion in
   a Gamess Format.

-  keyword **-threshold** The CI expansion contains coefficients
   (weights) for each determinant. This option sets the maximum
   coefficient to include in the QMC run. By default it is set to 1e-20
   (meaning all determinants in an expansion are taken into account). At
   the same time, if the threshold is set to a different value, for
   example :math:`1e-5`, any determinant with a weight
   :math:`|weight| < 1e-5` will be discarded and the determinant will
   not be considered.

-  keyword **-TargetState** ?

-  keyword **-NaturalOrbitals** ?

-  keyword **-optDetCoeffs** This flag enables optimization of the CI
   expansion coefficients. By default, optimization of the coefficients
   is disabled during wavefunction optimization runs.

Examples and more thorough descriptions of these options can be found in the lab section of this manual: :ref:`lab-advanced-molecules`.

Grid options
^^^^^^^^^^^^

These parameters control how the basis set is projected on a grid. The default parameters are chosen to be very efficient. Unless you have a very good reason, we do not recommend modifying them.

=============== =============== =========== ===========================

Tags
  **keyword**   **Value**       **default** **description**
  ``-gridtype`` log|log0|linear log         Grid type
  ``-first``    double          1e-6        First point of the grid
  ``-last``     double          100         Last point of the grid
  ``-size``     int             1001        Number of point in the grid
=============== =============== =========== ===========================

-  **-gridtype** Grid type can be logarithmic, logarithmic base 10, or
   linear

-  **-first** First value of the grid

-  **-last** Last value of the grid

-  **-size** Number of points in the grid between “first” and “last.”

Supported codes
^^^^^^^^^^^^^^^

- **PySCF**

  PySCF :cite:`Sun2018` is an all-purpose quantum chemistry
  code that can run calculations from simple Hartree-Fock to DFT, MCSCF,
  and CCSD, and for both isolated systems and periodic boundary
  conditions. PySCF can be downloaded from https://github.com/sunqm/pyscf.
  Many examples and tutorials can be found on the PySCF website, and all
  types of single determinants calculations are compatible with , thanks
  to active support from the authors of PySCF. A few additional steps are
  necessary to generate an output readable by ``convert4qmc``.

  This example shows how to run a Hartree-Fock calculation for the :math:`LiH`
  dimer molecule from PySCF and convert the wavefunction for QMCPACK.

  - **Python path**

    PySCF is a Python-based code. A Python module named **PyscfToQmcpack**
    containing the function **savetoqmcpack** is provided by and is located
    at ``qmcpack/src/QMCTools/PyscfToQmcpack.py``. To be accessible to the
    PySCF script, this path must be added to the PYTHONPATH environment
    variable. For the bash shell, this can be done as follows:

    ::

      export PYTHONPATH=/PATH_TO_QMCPACK/qmcpack/src/QMCTools:\$PYTHONPATH

  - **PySCF Input File**

    Copy and paste the following code in a file named LiH.py.

    ::

      #! /usr/bin/env python3
      from pyscf import gto, scf, df
      import numpy

      cell = gto.M(
         atom ='''
      Li 0.0 0.0 0.0
      H  0.0 0.0 3.0139239778''',
         basis ='cc-pv5z',
         unit="bohr",
         spin=0,
         verbose = 5,
         cart=False,
      )
      mf = scf.ROHF(cell)
      mf.kernel()

      ###SPECIFIC TO QMCPACK###
      title='LiH'
      from PyscfToQmcpack import savetoqmcpack

      savetoqmcpack(cell,mf,title)

    The arguments to the function **savetoqmcpack** are:

    -  **cell** This is the object returned from gto.M, containing the type
       of atoms, geometry, basisset, spin, etc.

    -  **mf** This is an object representing the PySCF level of theory, in
       this example, ROHF. This object contains the orbital coefficients of
       the calculations.

    -  **title** The name of the output file generated by PySCF. By default,
       the name of the generated file will be “default” if nothing is
       specified.

    |

    By adding the three lines below the “SPECIFIC TO QMCPACK” comment in the
    input file, the script will dump all the necessary data for QMCPACK into
    an HDF5 file using the value of “title” as an output name. PySCF is run
    as follows:

    ::

       >python LiH.py

    The generated HDF5 can be read by ``convert4qmc`` to generate the
    appropriate QMCPACK input files.

  - **Generating input files**

    As described in the previous section, generating input files for PySCF is as follows:

    ::

      > convert4qmc -pyscf LiH.h5

    The HDF5 file produced by “savetoqmcpack” contains the wavefunction in a
    form directly readable by QMCPACK. The wavefunction files from
    ``convert4qmc`` reference this HDF file as if the “-hdf5" option were
    specified (converting from PySCF implies the “-hdf5” option is always
    present).

Periodic boundary conditions with Gaussian orbitals from PySCF is fully supported for Gamma point and kpoints.

- **Quantum Package**

  QP2 :cite:`QP2` is a quantum chemistry code developed by the
  LCPQ laboratory in Toulouse, France, and Argonne National Laboratory for the PBC version.
  It can be downloaded from  https://github.com/QuantumPackage/qp2, and the tutorial within is
  quite extensive. The tutorial section of QP2 can guide you on how to
  install and run the code.

  After a QP2 calculation, the data needed for ``convert4qmc`` can be
  generated through

  ::

    qp_run save_for_qmcpack Myrun.ezfio 
    

  This command will generate an HDF5 file in the QMCPACK format named ``QP2QMCPACK.h5``
  ``convert4qmc`` can read this file and generate the ``*.structure.xml``, ``*.wfj.xml`` and other files needed to run QMCPACK. .  For example:

  ::

    convert4qmc -orbitals QP2QMCPACK.h5 -multidet QP2QMCPACK.h5 -prefix MySystem

  The main reason to use QP2 is to access the CIPSI algorithm to generate a
  multideterminant wavefunction. CIPSI is the preferred choice for
  generating a selected CI trial wavefunction for QMCPACK. An example on
  how to use QP2 for Hartree-Fock and selected CI can be found in
  :ref:`cipsi` of this manual. The converter code is actively
  maintained and codeveloped by both QMCPACK and QP2 developers.

- **Using -hdf5 tag**

  ::

    convert4qmc -gamess Myrun.out -hdf5

  This option is only used/useful with the gamess code as it is the only code not providing an HDF5 output
  The result will create QMCPACK input files but will also store all key data in the HDF5 format.

- **Mixing orbitals and multideterminants**


  Note that the ``QP2QMCPACK.h5`` combined with the tags ``-orbitals`` and
  ``-multidet`` allows the user to choose orbitals from a different code
  such as PYSCF and the multideterminant section from QP2. These two codes
  are fully compatible, and this route is also the only possible route for
  multideterminants for solids.

  ::

    convert4qmc -orbitals MyPyscfrun.h5 -multidet QP2QMCPACK.h5

- **GAMESS**

  QMCPACK can use the output of GAMESS :cite:`schmidt93` for any type of single determinant calculation (HF or DFT) or multideterminant (MCSCF) calculation. A description with an example can be found in the Advanced Molecular Calculations Lab (:ref:`lab-advanced-molecules`).

- **DIRAC**

  QMCPACK can use the output of DIRAC to run spin-orbit calculations using single-particle spinor wave functions for single-determinant calculations (DFT or closed-shell Dirac HF) or multideterminant complete open-shell configuration interaction (COSCI) wavefunctions. In the case of COSCI, the desired ground or excited state can be requested with ``-TargetState x``.

- **RMG**

  QMCPACK can use the HDF5 output of RMG DFT calculations. To generate this HDF5 output, set ``write_qmcpack_restart = "true"`` in the RMG input (file will be written to ``Waves/wave.out.h5``). ``convert4qmc`` will read the data from this HDF5 file and generate ``*.structure.xml``, ``*.wf{j,noj}.xml``, and ``*.qmc.in-wf{j,noj}.xml``. Pseudopotential files must be generated/moved manually by the user to ``X.qmcpp.xml``, where ``X`` is the appropriate element symbol (PP filename/path can be changed in the Hamiltonian section of ``*.qmc.in-wf{j,noj}.xml``).

  ::

    convert4qmc -rmg wave.out.h5

.. _pw2qmcpack:

pw2qmcpack.x
~~~~~~~~~~~~

``pw2qmcpack.x`` is an executable that converts PWSCF wavefunctions from the Quantum ESPRESSO (QE) package to QMCPACK readable
HDF5 format.  This utility is built alongside the QE postprocessing utilities. This utility is written in Fortran90 and is
distributed as a patch of the QE source code.  The patch, as well as automated QE download and patch scripts, can be found in
``qmcpack/external_codes/quantum_espresso``. Once built, we recommend also build QMCPACK with the QE_BIN option pointing to the
build pw.x and pw2qmcpack.x directory. This will enable workflow tests to be run.

pw2qmcpack can be used in serial in small systems and should be used in parallel with large systems for best performance. The K_POINT gamma optimization is not supported.

.. code-block::
  :caption: Sample ``pw2qmcpack.x`` input file ``p2q.in``
  :name: Listing 66

  &inputpp
    prefix     = 'bulk_silicon'
    outdir     = './'
    write_psir = .false.
  /

This example will cause ``pw2qmcpack.x`` to convert wavefunctions saved from
PWSCF with the prefix “bulk_silicon.” Perform the conversion via, for
example:

::

  mpirun -np 1 pw2qmcpack.x < p2q.in>& p2q.out

Because of the large plane-wave energy cutoffs in the pw.x calculation required by accurate PPs and the large system sizes of interest, one limitation of QE can be easily reached:
that ``wf_collect=.true.`` results in problems of writing and loading correct plane-wave coefficients on disks by pw.x because of the 32 bit integer limits. Thus, ``pw2qmcpack.x`` fails to convert the orbitals for QMCPACK. Since the release of QE v5.3.0, the converter has been fully parallelized to overcome this limitation completely.

By setting ``wf_collect=.false.`` (by default ``.false.`` in v6.1 and before and ``.true.`` since v6.2), pw.x does not collect the whole wavefunction into individual files for each k-point but instead writes one smaller file for each processor.
By running ``pw2qmcpack.x`` in the same parallel setup (MPI tasks and k-pools) as the last scf/nscf calculation with pw.x,
the orbitals distributed among processors will first be aggregated by the converter into individual temporal HDF5 files for each k-pool and then merged into the final file.
In large calculations, users should benefit from a significant reduction of time in writing the wavefunction by pw.x thanks to avoiding the wavefunction collection.

pw2qmcpack has been included in the test suite of QMCPACK (see instructions about how to activate the tests in :ref:`buildqe`).
There are tests labeled "no-collect" running the pw.x with the setting ``wf_collect=.false.``
The input files are stored at ``examples/solids/dft-inputs-polarized-no-collect``.
The scf, nscf, and pw2qmcpack runs are performed on 16, 12, and 12 MPI tasks with 16, 2, and 2 k-pools respectively.

convertpw4qmc
~~~~~~~~~~~~~

Convertpw4qmc is an executable that reads xml from a plane wave based DFT code and produces a QMCPACK readable
HDF5 format wavefunction.  For the moment, this supports both QBox and Quantum Epresso

In order to save the wavefunction from QBox so that convertpw4qmc can work on it, one needs to add a line to the
QBox input like

::

  save -text -serial basename.sample

after the end of a converged dft calculation.  This will write an ascii wavefunction file and will avoid
QBox's optimized parallel IO (which is not currently supported).

After the wavefunction file is written (basename.sample in this case) one can use convertpw4qmc as follows:

::

  convertpw4qmc basename.sample -o qmcpackWavefunction.h5

This reads the Qbox wavefunction and performs the Fourier transform before saving to a QMCPACK eshdf format wavefunction.  Currently multiple k-points are supported, but due to difficulties with the qbox wavefunction file format, the single particle orbitals do not have their proper energies associated with them.  This means that when tiling from a primitive cell to a supercell, the lowest n single particle orbitals from all necessary k-points will be used.  This can be problematic in the case of a metal and this feature should be used with EXTREME caution.

In the case of Quantum ESPRESSO, QE must be compiled with HDF support.  If this is the case, then an eshdf file can be generated by targeting the data-file-schema.xml file
generated in the output of Quantum ESPRESSO.  For example, if one is running a calculation with outdir = 'out' and prefix='Pt' then the converter can be invoked as:

::

  convertpw4qmc out/Pt.save/data-file-schema.xml -o qmcpackWavefunction.h5

Note that this method is insensitive to parallelization options given to Quantum ESPRESSO.  Additionally, it supports noncollinear magnetism and can be used to generate
wavefunctions suitable for qmcpack calculations with spin-orbit coupling.

.. _ppconvert:

ppconvert
~~~~~~~~~

``ppconvert`` is a utility to convert PPs between different commonly used formats. As with all operations on pseudopotentials,
great care should be exercised when using this tool. The tool is not yet considered to be fully robust and converted potentials
should be examined carefully. Please report any issues. Generally DFT-derived potentials should not be used with QMC. The main
intended use for the converter is to convert potentials generated for QMC calculations into formats acceptable to DFT and quantum
chemistry codes for trial wavefunction generation.

Currently it converts CASINO, FHI, UPF (generated by OPIUM), BFD, and GAMESS formats to several other formats including XML
(QMCPACK) and UPF (QE). See all the formats via ``ppconvert -h``.

For output formats requiring Kleinman-Bylander projectors, the atom will be solved with DFT if the projectors are not provided in
the input formats. This requires providing reference states and often needs extra tuning for heavy elements. To avoid ghost
states, the local channel can be changed via the ``--local_channel`` option. Ghost state considerations are similar to those of
DFT calculations but could be worse if ghost states were not considered during the original PP construction. To make the
self-consistent calculation converge, the density mixing parameter may need to be reduced via the ``--density_mix`` option. Note
that the reference state should include only the valence electrons. One reference state should be included for each channel in the
PP.

For example, for a sodium atom with a neon core, the reference state would be "1s(1)."
``--s_ref`` needs to include a 1s state, ``--p_ref`` needs to include a 2p state,
``--d_ref`` needs to include a 3d state, etc. If not specified, a corresponding state with zero occupation is added.
If the reference state is chosen as the neon core, setting empty reference states "" is technically correct.
In practice, reasonable reference states should be picked with care.
For PP with semi-core electrons in the valence, the reference state can be long.
For example, Ti PP has 12 valence electrons. When using the neutral atom state,
``--s_ref``, ``--p_ref``, and ``--d_ref`` are all set as "1s(2)2p(6)2s(2)3d(2)."
When using an ionized state, the three reference states are all set as "1s(2)2p(6)2s(2)" or "1s(2)2p(6)2s(2)3d(0)."

Unfortunately, if the generated UPF file is used in QE, the calculation may be incorrect because of the presence of "ghost"
states. Potentially these can be removed by adjusting the local channel (e.g., by setting ``--local_channel 1``, which chooses the
p channel as the local channel instead of d. For this reason, validation of UPF PPs is always required from the third row and is
strongly encouraged in general. For example, check that the expected ionization potential and electron affinities are obtained for
the atom and that dimer properties are consistent with those obtained by a quantum chemistry code or a plane-wave code that does
not use the Kleinman-Bylander projectors.

.. _molden2qmc:

molden2qmc
~~~~~~~~~~~

``molden2qmc`` is a tool used to convert molden files into an HDF5 file with the QMCPACK format.
Molden2qmc is a single program that can use multiple different quantum chemistry codes.
It is python code developed by Vladimir Konjkov originally for the CASINO code but then extended to QMCPACK.
This tool can be found at https://github.com/gjohnson3/molden2qmc.git.

Using molden2qmc
^^^^^^^^^^^^^^^^

General use of ``molden2qmc`` can be prompted by running ``molden2qmc.py`` and entering the corresponding quantum chemistry code number and the molden file name:

::

   number corresponding to the quantum chemistry code used to produce this MOLDEN file:
            0 -- TURBOMOLE
            1 -- PSI4
            2 -- CFOUR 2.0beta
            3 -- ORCA 3.X - 4.X
            4 -- DALTON2016
            5 -- MOLPRO
            6 -- NWCHEM
            7 -- QCHEM 4.X
            
Use the ``--qmcpack`` flag to create the file as an hdf5 file, suitable for QMCPACK.
Without the ``--qmcpack`` flag, the file will become a gwfn file for CASINO.            
Example: ``molden2qmc.py 5 n4.molden --qmcpack``.

Obtaining pseudopotentials
--------------------------

Pseudopotentiallibrary.org
~~~~~~~~~~~~~~~~~~~~~~~~~~

An open website collecting community developed and tested
pseudopotentials for QMC and other many-body calculations is being
developed at https://pseudopotentiallibrary.org. This site
includes potentials in QMCPACK format and an increasing range of
electronic structure and quantum chemistry codes. We recommend using
potentials from this site if available and suitable for your science
application.

.. _opium:

Opium
~~~~~

Opium is a pseudopotential generation code available from the website http://opium.sourceforge.net/.  Opium can generate pseudopotentials with either Hartree-Fock or DFT methods.  Once you have a useable pseudopotential param file (for example, Li.param), generate pseudopotentials for use in Quantum ESPRESSO with the upf format as follows:

.. code-block:
  :caption: Generate UPF-formatted pseudopotential with Opium
  :name: Listing 67

  opium Li.param Li.log all upf

This generates a UPF-formatted pseudopotential (``Li.upf``, in this case) for use in Quantum ESPRESSO.  The pseudopotential conversion tool ``ppconvert`` can then convert UPF to FSAtom xml format for use in QMCPACK:

.. code-block::
  :caption: Convert UPF-formatted pseudopotential to FSAtom xml format
  :name: Listing 68

  ppconvert --upf_pot Li.upf --xml Li.xml

.. _bfd:

Burkatzki-Filippi-Dolg
~~~~~~~~~~~~~~~~~~~~~~

Burkatzki *et al.* developed a set of energy-consistent pseudopotenitals
for use in QMC :cite:`Burkatzki07,Burkatzki08`, available at
http://www.burkatzki.com/pseudos/index.2.html. To convert for use in
QMCPACK, select a pseudopotential (choice of basis set is irrelevant to
conversion) in GAMESS format and copy the ending (pseudopotential) lines
beginning with(element symbol)-QMC GEN:

.. code-block::
  :caption: BFD Li pseudopotential in GAMESS format
  :name: Listing 69

  Li-QMC GEN 2 1
  3
  1.00000000 1 5.41040609
  5.41040609 3 2.70520138
  -4.60151975 2 2.07005488
  1
  7.09172172 2 1.34319829

Save these lines to a file (here, named ``Li.BFD.gamess``; the exact name may be anything as long as it is passed to ``ppconvert`` after --gamess_pot).  Then, convert using ``ppconvert`` with the following:

.. code-block::
  :caption: Convert GAMESS-formatted pseudopotential to FSAtom xml format
  :name: Listing 70

  ppconvert --gamess_pot Li.BFD.gamess --s_ref "2s(1)" --p_ref "2p(0)" --xml Li.BFD.xml

.. code-block::
  :caption: Convert GAMESS-formatted pseudopotential to Quantum ESPRESSO UPF format
  :name: Listing 71

  ppconvert --gamess_pot Li.BFD.gamess --s_ref "2s(1)" --p_ref "2p(0)" --log_grid --upf Li.BFD.upf

.. _CASINO:

CASINO
~~~~~~

The QMC code CASINO also makes available its pseudopotentials available at the website https://vallico.net/casinoqmc/pplib/. To use one in QMCPACK, select a pseudopotential and download its summary file (``summary.txt``), its tabulated form (``pp.data``), and (for ppconvert to construct the projectors to convert to Quantum ESPRESSO's UPF format) a CASINO atomic wavefunction for each angular momentum channel (``awfn.data_*``).  Then, to convert using ppconvert, issue the following command:

.. code-block::
  :caption: Convert CASINO-formatted pseudopotential to Quantum ESPRESSO UPF format
  :name: Listing 72

  ppconvert --casino_pot pp.data --casino_us awfn.data_s1_2S --casino_up awfn.data_p1_2P --casino_ud awfn.data_d1_2D --upf Li.TN-DF.upf

QMCPACK can directly read in the CASINO-formated pseudopotential (``pp.data``), but four parameters found in the pseudopotential summary file must be specified in the pseudo element (``l-local``, ``lmax``, ``nrule``, ``cutoff``)[see :ref:`nlpp` for details]:

.. code-block::
  :caption: XML syntax to use CASINO-formatted pseudopotentials in QMCPACK
  :name: Listing 73

  <pairpot type="pseudo" name="PseudoPot" source="ion0" wavefunction="psi0" format="xml">
     <pseudo elementType="Li" href="Li.pp.data" format="casino" l-local="s" lmax="2" nrule="2" cutoff="2.19"/>
     <pseudo elementType="H" href="H.pp.data" format="casino" l-local="s" lmax="2" nrule="2" cutoff="0.5"/>
  </pairpot>

.. _wftester:

wftester
~~~~~~~~

While not really a stand-alone application, wftester (short for “Wave
Function Tester") is a helpful tool for testing pre-existing and
experimental estimators and observables. It provides the user with
derived quantities from the Hamiltonian and wave function, but evaluated
at a small set of configurations.

The wftester is implemented as a QMCDriver, so one invokes QMCPACK in
the normal manner with a correct input XML, the difference being the
addition of an additional qmc input block. This is the main advantage of
this tool–it allows testing of realistic systems and realistic
combinations of observables. It can also be invoked before launching
into optimization, VMC, or DMC runs, as it is a valid <qmc> block.

As an example, the following code generates a random walker configuration and compares the trial wave function ratio computed in two different ways:

.. code-block::
  :caption: The following executes the wavefunction ratio test in "wftester"
  :name: Listing 74

  <qmc method="wftester">
    <parameter name="ratio">    yes    </parameter>
  </qmc>

Here's a summary of some of the tests provided:

-  Ratio Test. Invoked with

   ::

      <parameter name="ratio">yes</parameter>

   This computes the implemented wave function ratio associated with a
   single-particle move using two different methods.

-  Clone Test. Invoked with

   ::

      <parameter name="clone">yes</parameter>

   This checks the cloning of TrialWaveFunction, ParticleSet,
   Hamiltonian, and Walkers.

-  Elocal Test. Invoked with

   ::

      <parameter name="printEloc">yes</parameter>

   For an input electron configuration (can be random), print the value
   of TrialWaveFunction, LocalEnergy, and all local observables for this
   configuration.

-  Derivative Test. Invoked with

   ::

      <parameter name="ratio">deriv</parameter>}

   Computes electron gradients, laplacians, and wave function parameter
   derivatives using implemented calls and compares them to
   finite-difference results.

-  Ion Gradient Test. Invoked with

   ::

      <parameter name="source">ion0</parameter>

   Calls the implemented evaluateGradSource functions and compares them
   against finite-difference results.

-  “Basic Test". Invoked with

   ::

      <parameter name="basic">yes</parameter>

   Performs ratio, gradient, and laplacian tests against
   finite-difference and direct computation of wave function values.

The output of the various tests will be to standard out or "wftest.000" after successful execution of qmcpack.

.. bibliography:: /bibs/additional_tools.bib
