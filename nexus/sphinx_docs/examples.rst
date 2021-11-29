.. _examples:

Complete Examples
=================

.. note::
   Disclaimer: Please note that the examples given here do not generally qualify as
   production calculations because the supercell size, optimization
   process, DMC timestep and other key parameters may not be converged.
   Pseudopotentials are provided “as is” and should not be trusted without
   explicit validation.

Complete examples of calculations performed with Nexus are provided in
the following sections. These examples are intended to highlight basic
features of Nexus and act as templates for future calculations. If there
is an example you would like to contribute, or if you feel an example on
a particular topic is needed, please contact the developer at
krogeljt@ornl.gov to discuss the possibilities.

To perform the example calculations yourself, consult the ``examples``
directory in your Nexus installation:

.. code-block: rest

  /your_download_path/nexus/examples

The examples assume that you have working versions of ``pw.x``,
``pw2qmcpack.x``, ``qmcpack`` (real version), and ``qmcpack_complex``
(complex version) installed and in your ``PATH``. A brief description of
each example is given below.

Bulk Diamond VMC
   |
   | A representative bulk calculation. A simple workflow consisting of
     orbital generation with PWSCF, orbital conversion with pw2qmcpack,
     and a short VMC calculation with QMCPACK is performed.

Graphene Sheet DMC
   |
   | A representative slab calculation. The total DMC energy of a
     graphene “sheet” consisting of 8 atoms is computed. DFT is
     performed with PWSCF on the primitive cell followed by Jastrow
     optimization by QMCPACK and finally a supercell VMC+DMC calculation
     by QMCPACK.

C 20 Molecule DMC
   |
   | A representative molecular calculation. The total DMC energy of an
     ideal C 20 molecule is computed. DFT is performed with PWSCF on a
     periodic cell with some vacuum surrounding the molecule. QMCPACK
     optimization and VMC+DMC follow on the system with open boundary
     conditions.

   (Note that without the crystal field splitting afforded by the
   initial artificial periodicity, the Kohn-Sham HOMO would be
   degenerate, and so a production calculation would likely require more
   care in appropriately setting up the wavefunction.)

Oxygen Dimer DMC
   |
   | An example demonstrating automation of a simple parameter scan, in
     this case the interparticle spacing in an oxygen dimer. The reader
     will gain some experience modifying Nexus user scripts to produce
     automated workflows.

Bulk Diamond Excited States VMC
   |
   | A representative bulk excited states calculation. A simple workflow
     consisting of orbital generation with PWSCF, orbital conversion
     with pw2qmcpack, and a short VMC calculation with QMCPACK is
     performed.

.. _diamond-dmc:

Example 1: Bulk Diamond VMC
---------------------------

The files for this example are found in:

.. code-block:: rest

  /your_download_path/nexus/examples/qmcpack/diamond

By following the instructions contained in this section the reader can
execute a simple workflow with Nexus. The workflow presented here is
intended to illustrate basic Nexus usage. The example workflow has three
stages: (1) orbital generation in a primitive (2 atom) cell of diamond
with PWSCF, (2) conversion of the orbitals from the native PWSCF format
to the ESHDF format that QMCPACK reads, and (3) a minimal variational
Monte Carlo (VMC) run of a 16 atom supercell of diamond with QMCPACK.
The Nexus input script corresponding to this workflow is shown below.

The script is similar to the one discussed in :ref:`user-imports`, differing mainly in the name-by-name imports and
the description of the physical system. Instead of reading an external
VASP POSCAR file, the structure is specified in a format native to
Nexus. The physical system is specified by providing the unit system
(Angstrom in this case), the three vectors comprising the axes of the
simulation cell, and the names and positions (Cartesian coordinates) of
the atoms involved. The use of “``tiling=(2,2,2)``” communicates the
request that a :math:`2\times2\times2` supercell be constructed out of
the specified 2 atom primitive cell. The k-point grid applies to the
supercell and is comprised of a single k-point. The eight corresponding
primitive cell images of this k-point are determined automatically by
Nexus. The text “``C = 4``” specifies the number of valence electrons
for the carbon pseudopotential. The pseudopotential files used in this
example (``C.BFD.*``) have been adapted from an open access
pseudopotential database (see
http://www.burkatzki.com/pseudos/index.2.html) for use in PWSCF and
QMCPACK. Since the ``PhysicalSystem`` object, “``dia16``”, contains both
the supercell and its equivalent folded/primitive version, the PWSCF DFT
calculation will be performed in the primitive cell to save memory for
the subsequent VMC calculation of the full supercell performed with
QMCPACK.

.. code::

  #! /usr/bin/env python3

  # nexus imports
  from nexus import settings,Job,run_project
  from nexus import generate_physical_system
  from nexus import generate_pwscf
  from nexus import generate_pw2qmcpack
  from nexus import generate_qmcpack,vmc

  # general settings for nexus
  settings(
      pseudo_dir    = '../pseudopotentials',# directory with all pseudopotentials
      status_only   = 0,                    # only show status of runs
      generate_only = 0,                    # only make input files
      sleep         = 3,                    # check on runs every 3 seconds
      machine       = 'ws16'                # local machine is 16 core workstation
      )

  # generate diamond structure
  dia16 = generate_physical_system(
      units  = 'A',                      # Angstrom units
      axes   = [[1.785,1.785,0.   ],     # cell axes
                [0.   ,1.785,1.785],
                [1.785,0.   ,1.785]],
      elem   = ['C','C'],                # 2 C atoms
      pos    = [[0.    ,0.    ,0.    ],  # atomic positions
                [0.8925,0.8925,0.8925]],
      tiling = (2,2,2),                  # tile to 16 atom cell
      kgrid  = (1,1,1),                  # single supercell k-point
      kshift = (0,0,0),                  #  at gamma
      C      = 4                         # pseudo-C (4 val. elec.)
      )

  # scf run produces orbitals
  scf = generate_pwscf(
      identifier   = 'scf',           # identifier/file prefix
      path         = 'diamond/scf',   # directory for scf run
      job          = Job(cores=16,app='pw.x'),
      input_type   = 'generic',
      calculation  = 'scf',           # perform scf calculation
      input_dft    = 'lda',           # dft functional
      ecutwfc      = 200,             # planewave energy cutoff
      conv_thr     = 1e-8,            # scf convergence threshold
      nosym        = True,            # don't use symmetry
      wf_collect   = True,            # write orbitals
      system       = dia16,           # run diamond system
      pseudos      = ['C.BFD.upf'],   # pwscf PP for C
      )

  # convert orbitals for qmcpack
  conv = generate_pw2qmcpack(
      identifier   = 'conv',          # identifier/file prefix
      path         = 'diamond/scf',   # directory for conv job
      job          = Job(cores=1,app='pw2qmcpack.x'),
      write_psir   = False,           # output in k-space
      dependencies = (scf,'orbitals') # get orbitals from scf
      )

  # vmc run
  qmc = generate_qmcpack(
      identifier   = 'vmc',           # identifier/file prefix
      path         = 'diamond/vmc',   # directory for vmc run
      job          = Job(cores=16,threads=4,app='qmcpack'),
      input_type   = 'basic',
      system       = dia16,           # run diamond system
      pseudos      = ['C.BFD.xml'],   # qmcpack PP for C
      jastrows     = [],              # no jastrows, test run
      calculations = [
          vmc(                        # vmc inputs
              walkers     =   1,      #  one walker per core
              warmupsteps =  20,      #  20 steps before measurement
              blocks      = 200,      #  200 blocks
              steps       =  10,      #   of 10 MC steps each
              substeps    =   2,      #   2 substeps w/o measurement
              timestep    =  .4       #  0.4/Ha timestep
              )
          ],
      dependencies = (conv,'orbitals')# get orbitals from conv job
      )

  # nexus monitors all runs
  run_project(scf,conv,qmc)

To fully execute the usage example provided here, copies of PWSCF,
QMCPACK, and the orbital converter pw2qmcpack will need to be installed
on the local machine. The example assumes that the executables are in
the user’s ``PATH`` and are named ``pw.x``, ``qmcpack``, and
``pw2qmcpack.x``. See :ref:`install-code` for download and
installation instructions for these codes. A test of Nexus including the
generation of input files, but without actual job submission, can be
performed without installing these codes. However, Python itself and
NumPy are required to run Nexus (see :ref:`install-python`. The
example also assumes the local machine is a workstation with 16
available cores (“``ws16``”). If fewer than 16 cores are available,
*e.g.* 4, change the example files to reflect this:
``ws16``\ :math:`\rightarrow`\ ``ws4``,
``Job(cores=16,\ldots)``\ :math:`\rightarrow`\ ``Job(cores=4,\ldots)``.

In this example, we will run Nexus in three different modes:

#. status mode: print the status of each simulation and then exit
   (``status_only=1``).

#. generate mode: generate input files but do not execute workflows
   (``generate_only=1``).

#. execute mode: execute workflows by submitting jobs and monitoring
   simulation progress (``status_only=0, generate_only=0``).

Only the last mode requires executables for PWSCF and QMCPACK.

First, run Nexus in status mode. Enter the ``examples/qmcpack/diamond``
directory, open ``diamond.py`` with a text editor and set
“``status_only``\ =1”. Run the script by typing “``./diamond.py``” at
the command line and inspect the output. The output should be similar to
the text below (without the comments):

.. code:: rest

  Pseudopotentials
    reading pp:  ../pseudopotentials/C.BFD.upf  # dft PP found
    reading pp:  ../pseudopotentials/C.BFD.xml  # qmc PP found

  Project starting
  checking for file collisions        # files do not overlap
  loading cascade images              # load saved workflow state
    cascade 0 checking in             # only one workflow/cascade
  checking cascade dependencies       # match producers/consumers
    all simulation dependencies satisfied
  cascade status
    setup, sent_files, submitted, finished, got_output, analyzed
    000000  scf  ./runs/diamond/scf   # no work has been done yet
    000000  conv  ./runs/diamond/scf  # for any of the
    000000  vmc  ./runs/diamond/vmc   # three simulations
    setup, sent_files, submitted, finished, got_output, analyzed


The binary string “000000” indicates that none of the six stages of
simulation progression have been completed. These stages correpond to
the following actions/states: writing input files (“``setup``”), copying
pseudopotential files (“``sent_files``”), submitting simulation jobs for
execution (“``submitted``”), the completion of a simulation job
(“``finished``”), collecting output files (“``got_output``”), and
preprocessing output files for later analysis (“``analyzed``”) . In a
production setting, this mode is useful for checking the status of
current workflows/cascades prior to adding new ones. It is also useful
in general for detecting any problems with the Nexus input script
itself.

Next, run the example in generate mode. Set “``status_only`` =0” and
“``generate_only`` =1”, then run the example script again. Instead of
showing workflow status, Nexus will now perform a dry run of the
workflows by generating all of the run directories and input files. The
output should contain text similar to what is shown below:

.. code:: rest

  starting runs:                      # start submitting jobs
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  poll 0  memory 60.45 MB             # first poll cycle
    Entering ./runs/diamond/scf 0     # scf job
      writing input files  0 scf      # input file written
    Entering ./runs/diamond/scf 0
      sending required files  0 scf   # PP files copied
      submitting job  0 scf           # job is in virtual queue
    Entering ./runs/diamond/scf 0
      Would have executed:            # shows submission command
        export OMP_NUM_THREADS=1      # does not execute
        mpirun -np 16 pw.x -input scf.in

  poll 1  memory 60.72 MB
    Entering ./runs/diamond/scf 0
      copying results  0 scf          # output file copying stage
    Entering ./runs/diamond/scf 0
      analyzing  0 scf                # output analysis stage

  poll 2  memory 60.73 MB             # third poll cycle
    Entering ./runs/diamond/scf 1     # similar for conv job
      writing input files  1 conv
    Entering ./runs/diamond/scf 1
      sending required files  1 conv
      submitting job  1 conv
    Entering ./runs/diamond/scf 1
      Would have executed:
        export OMP_NUM_THREADS=1
        mpirun -np 1 pw2qmcpack.x<conv.in

  poll 3  memory 60.73 MB
    Entering ./runs/diamond/scf 1
      copying results  1 conv
    Entering ./runs/diamond/scf 1
      analyzing  1 conv

  poll 4  memory 60.73 MB             # fifth poll cycle
    Entering ./runs/diamond/vmc 2     # similar for vmc job
      writing input files  2 vmc
    Entering ./runs/diamond/vmc 2
      sending required files  2 vmc
      submitting job  2 vmc
    Entering ./runs/diamond/vmc 2
      Would have executed:
        export OMP_NUM_THREADS=4
        mpirun -np 4 qmcpack vmc.in.xml

  poll 5  memory 60.78 MB
    Entering ./runs/diamond/vmc 2
      copying results  2 vmc
    Entering ./runs/diamond/vmc 2
      analyzing  2 vmc

  Project finished                      # jobs finished

The output describes the progress of each simulation. The run submission
commands are also clearly shown as well as the amount of memory used by
Nexus. There should now be a “``runs``” directory containing the
generated input files with the following structure:

.. code:: rest

  runs/
  \-- diamond               # main diamond directory
      |-- scf               #  scf directory
      |-- |-- C.BFD.upf     #   pwscf PP file
      |-- |-- conv.in       #   conv job input file
      |-- |-- pwscf_output  #   pwscf output directory
      |-- |-- scf.in        #   scf job input file
      |-- |-- sim_conv      #   nexus directory for conv
      |-- |-- |-- input.p   #    stored input object
      |-- |-- \-- sim.p     #    simulation status file
      |-- \-- sim_scf       #   nexus directory for scf
      |--     |-- input.p   #    stored input object
      |--     \-- sim.p     #    simulation status file
      \-- vmc               #  vmc directory
          |-- C.BFD.xml     #   qmcpack PP file
          |-- sim_vmc       #   nexus directory for vmc
          |-- |-- input.p   #    stored input object
          |-- \-- sim.p     #    simulation status file
          \-- vmc.in.xml    #   vmc job input file

The “``sim.p``” files record the state of each simulation. Inspect the
input files generated by Nexus (``scf.in``, ``conv.in``, and
``vmc.in.xml``). Compare the files with the input provided to Nexus in
``diamond.py``.

Finally, run the example in execute mode. Remove the “``runs``” and
“``results``” directories, set “``status_only`` =0” and
“``generate_only`` =0”, and rerun the script. The output shown should
be similar to what was seen for generate mode, only now there may be
multiple workflow polls while a particular simulation is running. The
“``sleep``” keyword controls how often the polls occur (every 3 seconds
in this example). Note that the Nexus host process sleeps in between
polls so that a minimum of computational resources are occupied. Once
“``Project finished``” is displayed, all the simulation runs should be
complete. Confirm the success of the runs by checking the output files.
The text “``JOB DONE.``” should appear near the end of the PWSCF output
file ``scf.out``. QMCPACK has completed successfully if
“``Total Execution time``” appears near the end of the output in
``vmc.out``.

.. _graphene-dmc:

Example 2: Graphene Sheet DMC
-----------------------------

The files for this example are found in:

.. code:: rest

  /your_download_path/nexus/examples/qmcpack/graphene

Take a moment to study the “input file” script (``graphene.py``) and the
attendant comments (prefixed with #).

.. code:: rest

  #! /usr/bin/env python3

  from nexus import settings,Job,run_project
  from nexus import generate_physical_system
  from nexus import generate_pwscf
  from nexus import generate_pw2qmcpack
  from nexus import generate_qmcpack
  from nexus import loop,linear,vmc,dmc


  # general settings for nexus
  settings(
      pseudo_dir    = '../pseudopotentials',# directory with all pseudopotentials
      sleep         = 3,                    # check on runs every 'sleep' seconds
      generate_only = 0,                    # only make input files
      status_only   = 0,                    # only show status of runs
      machine       = 'ws16',               # local machine is 16 core workstation
      )



  # generate the graphene physical system
  graphene = generate_physical_system(
      lattice   = 'hexagonal',      # hexagonal cell shape
      cell      = 'primitive',      # primitive cell
      centering = 'P',              # primitive basis centering
      constants = (2.462,10.0),     # a,c constants
      units     = 'A',              # in Angstrom
      atoms     = ('C','C'),        # C primitive atoms
      basis     = [[ 0  , 0  , 0],  # basis vectors
                   [2./3,1./3, 0]],
      tiling    = (2,2,1),          # tiling of primitive cell
      kgrid     = (1,1,1),          # Monkhorst-Pack grid
      kshift    = (.5,.5,.5),       # and shift
      C         = 4                 # C has 4 valence electrons
      )


  # list of simulations in workflow
  sims = []

  # scf run produces charge density
  scf = generate_pwscf(
      # nexus inputs
      identifier   = 'scf',           # identifier/file prefix
      path         = 'graphene/scf',  # directory for scf run
      job          = Job(cores=16),   # run on 16 cores
      pseudos      = ['C.BFD.upf'],   # pwscf PP file
      system       = graphene,        # run graphene
      # input format selector
      input_type   = 'scf',           # scf, nscf, relax, or generic
      # pwscf input parameters
      input_dft    = 'lda',           # dft functional
      ecut         =  150 ,           # planewave energy cutoff (Ry)
      conv_thr     =  1e-6,           # scf convergence threshold (Ry)
      mixing_beta  =    .7,           # charge mixing factor
      kgrid        = (8,8,8),         # MP grid of primitive cell
      kshift       = (1,1,1),         #  to converge charge density
      wf_collect   = False,           # don't collect orbitals
      use_folded   = True             # use primitive rep of graphene
      )
  sims.append(scf)

  # nscf run to produce orbitals for jastrow optimization
  nscf_opt = generate_pwscf(
      # nexus inputs
      identifier   = 'nscf',          # identifier/file prefix
      path         = 'graphene/nscf_opt', # directory for nscf run
      job          = Job(cores=16),   # run on 16 cores
      pseudos      = ['C.BFD.upf'],   # pwscf PP file
      system       = graphene,        # run graphene
      # input format selector
      input_type   = 'nscf',          # scf, nscf, relax, or generic
      # pwscf input parameters
      input_dft    = 'lda',           # dft functional
      ecut         =  150 ,           # planewave energy cutoff (Ry)
      conv_thr     =  1e-6,           # scf convergence threshold (Ry)
      mixing_beta  =    .7,           # charge mixing factor
      nosym        = True,            # don't symmetrize k-points
      use_folded   = True,            # use primitive rep of graphene
      wf_collect   = True,            # write out orbitals
      kgrid        = (1,1,1),         # single k-point for opt
      kshift       = (0,0,0),         # gamma point
      # workflow dependencies
      dependencies = (scf,'charge_density')
      )
  sims.append(nscf_opt)

  # orbital conversion job for jastrow optimization
  p2q_opt = generate_pw2qmcpack(
      # nexus inputs
      identifier   = 'p2q',
      path         = 'graphene/nscf_opt',
      job          = Job(cores=1),
      # pw2qmcpack input parameters
      write_psir   = False,
      # workflow dependencies
      dependencies = (nscf_opt,'orbitals')
      )
  sims.append(p2q_opt)

  # Jastrow optimization
  opt = generate_qmcpack(
      # nexus inputs
      identifier   = 'opt',           # identifier/file prefix
      path         = 'graphene/opt',  # directory for opt run
      job          = Job(cores=16,app='qmcpack'),
      pseudos      = ['C.BFD.xml'],   # qmcpack PP file
      system       = graphene,        # run graphene
      # input format selector
      input_type   = 'basic',
      # qmcpack input parameters
      corrections  = [],
      jastrows     = [('J1','bspline',8),   # 1 body bspline jastrow
                      ('J2','bspline',8)],  # 2 body bspline jastrow
      calculations = [
          loop(max = 6,                        # No. of loop iterations
               qmc = linear(                   # linearized optimization method
                  energy               =  0.0, # cost function
                  unreweightedvariance =  1.0, #   is all unreweighted variance
                  reweightedvariance   =  0.0, #   no energy or r.w. var.
                  timestep             =  0.5, # vmc timestep (1/Ha)
                  warmupsteps          =  100, # MC steps before data collected
                  samples              = 16000,# samples used for cost function
                  stepsbetweensamples  =   10, # steps between uncorr. samples
                  blocks               =   10, # ignore this
                  minwalkers           =   0.1,#  and this
                  bigchange            =  15.0,#  and this
                  alloweddifference    =  1e-4 #  and this, for now
                  )
               )
          ],
      # workflow dependencies
      dependencies = (p2q_opt,'orbitals')
      )
  sims.append(opt)


  # nscf run to produce orbitals for final dmc
  nscf = generate_pwscf(
      # nexus inputs
      identifier   = 'nscf',          # identifier/file prefix
      path         = 'graphene/nscf', # directory for nscf run
      job          = Job(cores=16),   # run on 16 cores
      pseudos      = ['C.BFD.upf'],   # pwscf PP file
      system       = graphene,        # run graphene
      # input format selector
      input_type   = 'nscf',          # scf, nscf, relax, or generic
      # pwscf input parameters
      input_dft    = 'lda',           # dft functional
      ecut         =  150 ,           # planewave energy cutoff (Ry)
      conv_thr     =  1e-6,           # scf convergence threshold (Ry)
      mixing_beta  =    .7,           # charge mixing factor
      nosym        = True,            # don't symmetrize k-points
      use_folded   = True,            # use primitive rep of graphene
      wf_collect   = True,            # write out orbitals
      # workflow dependencies
      dependencies = (scf,'charge_density')
      )
  sims.append(nscf)

  # orbital conversion job for final dmc
  p2q = generate_pw2qmcpack(
      # nexus inputs
      identifier   = 'p2q',
      path         = 'graphene/nscf',
      job          = Job(cores=1),
      # pw2qmcpack input parameters
      write_psir   = False,
      # workflow dependencies
      dependencies = (nscf,'orbitals')
      )
  sims.append(p2q)

  # final dmc run
  qmc = generate_qmcpack(
      # nexus inputs
      identifier   = 'qmc',           # identifier/file prefix
      path         = 'graphene/qmc',  # directory for dmc run
      job          = Job(cores=16,app='qmcpack'),
      pseudos      = ['C.BFD.xml'],   # qmcpack PP file
      system       = graphene,        # run graphene
      # input format selector
      input_type   = 'basic',
      # qmcpack input parameters
      corrections  = [],              # no finite size corrections
      jastrows     = [],              # overwritten from opt
      calculations = [                # qmcpack input parameters for qmc
          vmc(                        # vmc parameters
              timestep      = 0.5,    # vmc timestep (1/Ha)
              warmupsteps   = 100,    # No. of MC steps before data is collected
              blocks        = 200,    # No. of data blocks recorded in scalar.dat
              steps         =  10,    # No. of steps per block
              substeps      =   3,    # MC steps taken w/o computing E_local
              samplesperthread = 40   # No. of dmc walkers per thread
              ),
          dmc(                        # dmc parameters
              timestep      = 0.01,   # dmc timestep (1/Ha)
              warmupsteps   =  50,    # No. of MC steps before data is collected
              blocks        = 400,    # No. of data blocks recorded in scalar.dat
              steps         =   5,    # No. of steps per block
              nonlocalmoves = True    # use Casula's T-moves
              ),                      #  (retains variational principle for NLPP's)
          ],
      # workflow dependencies
      dependencies = [(p2q,'orbitals'),
                      (opt,'jastrow')]
      )


  # nexus monitors all runs
  run_project(sims)


  # print out the total energy
  performed_runs = not settings.generate_only and not settings.status_only
  if performed_runs:
      # get the qmcpack analyzer object
      # it contains all of the statistically analyzed data from the run
      qa = qmc.load_analyzer_image()
      # get the local energy from dmc.dat
      le = qa.dmc[1].dmc.LocalEnergy  # dmc series 1, dmc.dat, local energy
      #  print the total energy for the 8 atom system
      print 'The DMC ground state energy for graphene is:'
      print '    {0} +/- {1} Ha'.format(le.mean,le.error)
  #end if

To run the example, navigate to the example directory and type

.. code:: rest

  ./graphene.py

or, alternatively,

.. code:: rest

  python ./graphene.py

You should see output like this (without the added # comments):

.. code:: rest

  Pseudopotentials   # reading pseudopotential files
      reading pp:  ../pseudopotentials/C.BFD.upf
      reading pp:  ../pseudopotentials/C.BFD.xml

  Project starting
    checking for file collisions  # ensure created files don't overlap
    loading cascade images        # load previous simulation state
      cascade 0 checking in
    checking cascade dependencies # ensure sim.'s have needed dep.'s
      all simulation dependencies satisfied

    starting runs:                # start submitting jobs
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    poll 0  memory 56.28 MB
      Entering ./runs/graphene/scf 0      # scf job
        writing input files  0 scf        # input file written
      Entering ./runs/graphene/scf 0
        sending required files  0 scf     # PP files copied
        submitting job  0 scf             # job is in virtual queue
      Entering ./runs/graphene/scf 0
        Executing:                        # job executed on workstation
          export OMP_NUM_THREADS=1
          mpirun -np 16 pw.x -input scf.in

    poll 1  memory 56.30 MB               # waiting for job to finish
    poll 2  memory 56.30 MB
    poll 3  memory 56.30 MB
    poll 4  memory 56.30 MB
      Entering ./runs/graphene/scf 0
        copying results  0 scf            # job is finished, copy results
      Entering ./runs/graphene/scf 0
        analyzing  0 scf                  # analyze output data

                                          # now do the same for
                                          # nscf job for Jastrow opt
                                          #   single k-point
                                          # nscf job for VMC/DMC
                                          #   multiple k-points

    poll 5  memory 56.31 MB
      Entering ./runs/graphene/nscf 1     # nscf dmc
        writing input files  1 nscf
        ...
      Entering ./runs/graphene/nscfopt 4  # nscf opt
        writing input files  4 nscf
        ...

                                          # now convert KS orbitals
                                          # to eshdf format
                                          # with pw2qmcpack.x
                                          # for nscf opt & nscf dmc

    poll 7  memory 56.32 MB
      Entering ./runs/graphene/nscf 2     # convert dmc orbitals
        sending required files  2 p2q
        ...
      Entering ./runs/graphene/nscfopt 4  # convert opt orbitals
        copying results  4 nscf
        ...

    poll 10  memory 56.32 MB
      Entering ./runs/graphene/opt 6      # submit jastrow opt
        writing input files  6 opt        # write input file
      Entering ./runs/graphene/opt 6
        sending required files  6 opt     # copy PP files
        submitting job  6 opt             # job is in virtual queue
      Entering ./runs/graphene/opt 6
        Executing:                        # run qmcpack
          export OMP_NUM_THREADS=1        # w/ complex arithmetic
          mpirun -np 16 qmcpack_complex opt.in.xml

    poll 11  memory 56.32 MB
    poll 12  memory 56.32 MB
    poll 13  memory 56.32 MB
    ...
    ...
    ...
    poll 793  memory 56.32 MB   # qmcpack opt finishes
    poll 794  memory 56.32 MB   # nearly an hour later
    poll 795  memory 56.32 MB
      Entering ./runs/graphene/opt 6
        copying results  6 opt            # copy output files
      Entering ./runs/graphene/opt 6
        analyzing  6 opt                  # analyze the results

    poll 796  memory 56.41 MB
      Entering ./runs/graphene/qmc 3      # submit dmc
        writing input files  3 qmc        # write input file
      Entering ./runs/graphene/qmc 3
        sending required files  3 qmc     # copy PP files
        submitting job  3 qmc             # job is in virtual queue
      Entering ./runs/graphene/qmc 3
        Executing:                        # run qmcpack
          export OMP_NUM_THREADS=1
          mpirun -np 16 qmcpack_complex qmc.in.xml

    poll 797  memory 57.31 MB
    poll 798  memory 57.31 MB
    poll 799  memory 57.31 MB
    ...
    ...
    ...
    poll 1041  memory 57.31 MB   # qmcpack dmc finishes
    poll 1042  memory 57.31 MB   # about 15 minutes later
    poll 1043  memory 57.31 MB
      Entering ./runs/graphene/qmc 3
        copying results  3 qmc            # copy output files
      Entering ./runs/graphene/qmc 3
        analyzing  3 qmc                  # analyze the results

  Project finished                        # all jobs are finished

  The DMC ground state energy for graphene is:
      -45.824960552 +/- 0.00498990689364 Ha    # one value from
                                               # qmcpack analyzer

If successful, you have just performed a start-to-finish DMC calculation.
The total energy quoted above probably will not match the one you
produce due to different compilation environments and the probabilistic
nature of DMC. They should not, however differ by three sigma.

Take some time to inspect the input files generated by Nexus and the
output files from PWSCF and QMCPACK. The runs were performed in
sub-directories of the ``runs`` directory. The order of execution of the
simulations is roughly ``scf``, ``nscf``, ``nscfopt``, ``opt``, then
``qmc``.

.. code:: rest

  runs
  └── graphene_test
      ├── nscf
      │   ├── nscf.in
      │   └── nscf.out
      ├── nscfopt
      │   ├── nscf.in
      │   └── nscf.out
      ├── opt
      │   ├── opt.in.xml
      │   └── opt.out
      ├── qmc
      │   ├── qmc.in.xml
      │   └── qmc.out
      └── scf
          ├── scf.in
          └── scf.out

The directories above contain all the files generated by the
simulations. Often one only wants to save the files with the most
important data, which are generally small. These are copied to the
``results`` directory which mirrors the structure of ``runs``.

.. code:: rest

  results
  └── runs
      └── graphene_test
          ├── nscf
          │   ├── nscf.in
          │   └── nscf.out
          ├── nscfopt
          │   ├── nscf.in
          │   └── nscf.out
          ├── opt
          │   ├── opt.in.xml
          │   └── opt.out
          ├── qmc
          │   ├── qmc.in.xml
          │   └── qmc.out
          └── scf
              ├── scf.in
              └── scf.out

Although this QMC run was performed at a single k-point, a
twist-averaged run could be performed simply by changing ``kgrid`` in
``generate_physical_system`` from ``(1,1,1)`` to ``(4,4,1)``, or
similar.

.. _c20-dmc:

Example 3: C 20 Molecule DMC
----------------------------

The files for this example are found in:

.. code:: rest

  /your_download_path/nexus/examples/qmcpack/c20

Take a moment to study the “input file” script (``c20_example.py``) and
the attendant comments (prefixed with #). The relevant differences from
the graphene example mostly involve how the structure is procured (it is
read from an XYZ file rather than being generated), the boundary
conditions (open BC’s, see ``bconds`` in the QMCPACK input parameters),
and the workflow involved.

::

  #! /usr/bin/env python3

  from nexus import settings,Job,run_project
  from nexus import Structure,PhysicalSystem
  from nexus import generate_pwscf
  from nexus import generate_pw2qmcpack
  from nexus import generate_qmcpack
  from nexus import loop,linear,vmc,dmc


  # general settings for nexus
  settings(
      pseudo_dir    = '../pseudopotentials',# directory with all pseudopotentials
      sleep         = 3,                    # check on runs every 'sleep' seconds
      generate_only = 0,                    # only make input files
      status_only   = 0,                    # only show status of runs
      machine       = 'ws16',               # local machine is 16 core workstation
      )


  #generate the C20 physical system
  # specify the xyz file
  structure_file = 'c20.cage.xyz'
  # make an empty structure object
  structure = Structure()
  # read in the xyz file
  structure.read_xyz(structure_file)
  # place a bounding box around the structure
  structure.bounding_box(
      box   = 'cubic',         # cube shaped cell
      scale = 1.5              # 50% extra space
      )
  # make it a gamma point cell
  structure.add_kmesh(
      kgrid      = (1,1,1),    # Monkhorst-Pack grid
      kshift     = (0,0,0)     # and shift
      )
  # add electronic information
  c20 = PhysicalSystem(
      structure = structure,   # C20 structure
      net_charge = 0,          # net charge in units of e
      net_spin   = 0,          # net spin in units of e-spin
      C          = 4           # C has 4 valence electrons
      )


  # list of simulations in workflow
  sims = []

  # scf run produces charge density
  scf = generate_pwscf(
      # nexus inputs
      identifier   = 'scf',           # identifier/file prefix
      path         = 'c20/scf',       # directory for scf run
      job          = Job(cores=16),   # run on 16 cores
      pseudos      = ['C.BFD.upf'],   # pwscf PP file
      system       = c20,             # run c20
      # input format selector
      input_type   = 'scf',           # scf, nscf, relax, or generic
      # pwscf input parameters
      input_dft    = 'lda',           # dft functional
      ecut         =  150 ,           # planewave energy cutoff (Ry)
      conv_thr     =  1e-6,           # scf convergence threshold (Ry)
      mixing_beta  =    .7,           # charge mixing factor
      nosym        = True,            # don't use symmetry
      wf_collect   = True,            # write out orbitals
      )
  sims.append(scf)

  # orbital conversion job for opt and dmc
  p2q = generate_pw2qmcpack(
      # nexus inputs
      identifier   = 'p2q',
      path         = 'c20/nscf',
      job          = Job(cores=1),
      # pw2qmcpack input parameters
      write_psir   = False,
      # workflow dependencies
      dependencies = (scf,'orbitals')
      )
  sims.append(p2q)


  # Jastrow optimization
  opt = generate_qmcpack(
      # nexus inputs
      identifier   = 'opt',           # identifier/file prefix
      path         = 'c20/opt',       # directory for opt run
      job          = Job(cores=16,app='qmcpack'),
      pseudos      = ['C.BFD.xml'],   # qmcpack PP file
      system       = c20,             # run c20
      # input format selector
      input_type   = 'basic',
      # qmcpack input parameters
      corrections  = [],
      jastrows     = [('J1','bspline',8,6),   # 1 body bspline jastrow
                      ('J2','bspline',8,8)],  # 2 body bspline jastrow
      calculations = [
          loop(max = 6,                        # No. of loop iterations
               qmc = linear(                   # linearized optimization method
                  energy               =  0.0, # cost function
                  unreweightedvariance =  1.0, #   is all unreweighted variance
                  reweightedvariance   =  0.0, #   no energy or r.w. var.
                  timestep             =  0.5, # vmc timestep (1/Ha)
                  warmupsteps          =  100, # MC steps before data collected
                  samples              = 16000,# samples used for cost function
                  stepsbetweensamples  =   10, # steps between uncorr. samples
                  blocks               =   10, # ignore this
                  minwalkers           =   0.1,#  and this
                  bigchange            =  15.0,#  and this
                  alloweddifference    =  1e-4 #  and this, for now
                  )
               )
          ],
      # workflow dependencies
      dependencies = (p2q,'orbitals')
      )
  sims.append(opt)


  # final dmc run
  qmc = generate_qmcpack(
      # nexus inputs
      identifier   = 'qmc',           # identifier/file prefix
      path         = 'c20/qmc',  # directory for dmc run
      job          = Job(cores=16,app='qmcpack'),
      pseudos      = ['C.BFD.xml'],   # qmcpack PP file
      system       = c20,             # run c20
      # input format selector
      input_type   = 'basic',
      # qmcpack input parameters
      corrections  = [],              # no finite size corrections
      jastrows     = [],              # overwritten from opt
      calculations = [                # qmcpack input parameters for qmc
          vmc(                        # vmc parameters
              timestep      = 0.5,    # vmc timestep (1/Ha)
              warmupsteps   = 100,    # No. of MC steps before data is collected
              blocks        = 200,    # No. of data blocks recorded in scalar.dat
              steps         =  10,    # No. of steps per block
              substeps      =   3,    # MC steps taken w/o computing E_local
              samplesperthread = 40   # No. of dmc walkers per thread
              ),
          dmc(                        # dmc parameters
              timestep      = 0.01,   # dmc timestep (1/Ha)
              warmupsteps   =  50,    # No. of MC steps before data is collected
              blocks        = 400,    # No. of data blocks recorded in scalar.dat
              steps         =   5,    # No. of steps per block
              nonlocalmoves = True    # use Casula's T-moves
              ),                      #  (retains variational principle for NLPP's)
          ],
      # workflow dependencies
      dependencies = [(p2q,'orbitals'),
                      (opt,'jastrow')]
      )



  # nexus monitors all runs
  run_project(sims)



  # print out the total energy
  performed_runs = not settings.generate_only and not settings.status_only
  if performed_runs:
      # get the qmcpack analyzer object
      # it contains all of the statistically analyzed data from the run
      qa = qmc.load_analyzer_image()
      # get the local energy from dmc.dat
      le = qa.dmc[1].dmc.LocalEnergy  # dmc series 1, dmc.dat, local energy
      #  print the total energy for the 20 atom system
      print 'The DMC ground state energy for C20 is:'
      print '    {0} +/- {1} Ha'.format(le.mean,le.error)
  #end if


To run the example, navigate to the example directory and type

.. code:: rest

  ./c20.py

or, alternatively,

.. code:: rest

  python ./c20.py

You should see output like this (without the added # comments):

.. code:: rest

  Pseudopotentials   # reading pseudopotential files
      reading pp:  ../pseudopotentials/C.BFD.upf
      reading pp:  ../pseudopotentials/C.BFD.xml

  Project starting
    checking for file collisions  # ensure created files don't overlap
    loading cascade images        # load previous simulation state
      cascade 0 checking in
    checking cascade dependencies # ensure sim.'s have needed dep.'s
      all simulation dependencies satisfied

    starting runs:                # start submitting jobs
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    poll 0  memory 56.21 MB
      Entering ./runs/c20/scf 0       # scf job
        writing input files  0 scf    # input file written
      Entering ./runs/c20/scf 0
        sending required files  0 scf # PP files copied
        submitting job  0 scf         # job is in the virtual queue
      Entering ./runs/c20/scf 0
        Executing:                    # job executed on workstation
          export OMP_NUM_THREADS=1
          mpirun -np 16 pw.x -input scf.in

    poll 1  memory 56.23 MB           # waiting for job to finish
    poll 2  memory 56.23 MB
    poll 3  memory 56.23 MB
    poll 4  memory 56.23 MB
    poll 5  memory 56.23 MB
    poll 6  memory 56.23 MB
    poll 7  memory 56.23 MB
    poll 8  memory 56.23 MB
      Entering ./runs/c20/scf 0
        copying results  0 scf        # job is finished, copy results
      Entering ./runs/c20/scf 0
        analyzing  0 scf              # analyze output data

    poll 9  memory 56.23 MB           # now convert KS orbitals
      Entering ./runs/c20/scf 1       # to eshdf format
        writing input files  1 p2q    # with pw2qmcpack.x
        ...

    poll 12  memory 56.23 MB
      Entering ./runs/c20/opt 3       # submit jastrow opt
        writing input files  3 opt    # write input file
      Entering ./runs/c20/opt 3
        sending required files  3 opt # copy PP files
        submitting job  3 opt         # job is in virtual queue
      Entering ./runs/c20/opt 3
        Executing:                    # run qmcpack
          export OMP_NUM_THREADS=1    # w/ real arithmetic
          mpirun -np 16 qmcpack opt.in.xml

    poll 13  memory 56.24 MB
    poll 14  memory 56.24 MB
    poll 15  memory 56.24 MB
    ...
    ...
    ...
    poll 204  memory 56.24 MB   # qmcpack opt finishes
    poll 205  memory 56.24 MB   # about 10 minutes later
    poll 206  memory 56.24 MB
      Entering ./runs/c20/opt 3
        copying results  3 opt        # copy output files
      Entering ./runs/c20/opt 3
        analyzing  3 opt              # analyze the results

    poll 207  memory 56.27 MB
      Entering ./runs/c20/qmc 2       # submit dmc
        writing input files  2 qmc    # write input file
      Entering ./runs/c20/qmc 2
        sending required files  2 qmc # copy PP files
        submitting job  2 qmc         # job is in virtual queue
      Entering ./runs/c20/qmc 2
        Executing:                    # run qmcpack
          export OMP_NUM_THREADS=1
          mpirun -np 16 qmcpack qmc.in.xml

    poll 208  memory 56.49 MB
    poll 209  memory 56.49 MB
    poll 210  memory 56.49 MB
    ...
    ...
    ...
    poll 598  memory 56.49 MB   # qmcpack dmc finishes
    poll 599  memory 56.49 MB   # about 20 minutes later
    poll 600  memory 56.49 MB
      Entering ./runs/c20/qmc 2
        copying results  2 qmc        # copy output files
      Entering ./runs/c20/qmc 2
        analyzing  2 qmc              # analyze the results

  Project finished                    # all jobs are finished

  The DMC ground state energy for C20 is:
      -112.890695404 +/- 0.0151688786226 Ha  # one value from
                                             # qmcpack analyzer

Again, the total energy quoted above probably will not match the one you
produce due to different compilation environments and the probabilistic
nature of QMC. The results should still be statistically comparable.

The directory trees generated by Nexus for C 20 have a similar structure
to the graphene example. Note the absence of the ``nscf`` runs. The
order of execution of the simulations is ``scf``, ``opt``, then ``qmc``.

.. code:: rest

  runs
  └── c20_test
      ├── opt
      │   ├── opt.in.xml
      │   └── opt.out
      ├── qmc
      │   ├── qmc.in.xml
      │   └── qmc.out
      └── scf
          ├── scf.in
          └── scf.out
  results
  └── runs
      └── c20_test
          ├── opt
          │   ├── opt.in.xml
          │   └── opt.out
          ├── qmc
          │   ├── qmc.in.xml
          │   └── qmc.out
          └── scf
              ├── scf.in
              └── scf.out

Example 4 Automated oxygen dimer binding curve
----------------------------------------------

The files for this example are found in:

.. code:: rest

  /your_download_path/nexus/examples/qmcpack/oxygen_dimer

Enter the ``examples/qmcpack/oxygen_dimer`` directory. Open
``oxygen_dimer.py`` with a text editor. The overall format is similar to
the example file shown in the prior sections. The header material,
including Nexus imports, settings, and the job parameters for QMC are
nearly identical.

Following the job parameters, inputs for the optimization method are
given. The keywords are identical to the parameters of QMCPACK’s XML
input file.

::

  linopt1 = linear(
      energy               = 0.0,
      unreweightedvariance = 1.0,
      reweightedvariance   = 0.0,
      timestep             = 0.4,
      samples              = 5000,
      warmupsteps          = 50,
      blocks               = 200,
      substeps             = 1,
      nonlocalpp           = True,
      usebuffer            = True,
      walkers              = 1,
      minwalkers           = 0.5,
      maxweight            = 1e9,
      usedrift             = True,
      minmethod            = 'quartic',
      beta                 = 0.025,
      exp0                 = -16,
      bigchange            = 15.0,
      alloweddifference    = 1e-4,
      stepsize             = 0.2,
      stabilizerscale      = 1.0,
      nstabilizers         = 3
      )

Requesting multiple loop's with different numbers of samples is more compact than in the native XML input file:

::

  linopt1 = ...

  linopt2 = linopt1.copy()
  linopt2.samples = 20000 # opt w/ 20000 samples

  linopt3 = linopt1.copy()
  linopt3.samples = 40000 # opt w/ 40000 samples

  opt_calcs = [loop(max=8,qmc=linopt1), # loops over opt's
               loop(max=6,qmc=linopt2),
               loop(max=4,qmc=linopt3)]

The VMC/DMC method inputs also mirror the XML:

::

  qmc_calcs = [
      vmc(
          walkers     =   1,
          warmupsteps =  30,
          blocks      =  20,
          steps       =  10,
          substeps    =   2,
          timestep    =  .4,
          samples     = 2048
          ),
      dmc(
          warmupsteps   = 100,
          blocks        = 400,
          steps         =  32,
          timestep      = 0.01,
          nonlocalmoves = True
          )
      ]

As in the prior examples, the oxygen dimer is generated with the ``generate_physical_system`` function:

::

  dimer = generate_physical_system(
      type       = 'dimer',
      dimer      = ('O','O'),
      separation = 1.2074*scale,
      Lbox       = 15.0,
      units      = 'A',
      net_spin   = 2,
      O          = 6
      )

Similar syntax can be used to generate crystal structures or to specify systems with arbitrary atomic configurations and simulation cells.  Notice that a "``scale``" variable has been introduced to stretch or compress the dimer.

Next, objects representing DFT calculations and orbital conversions are
constructed with the ``generate_pwscf`` and ``generate_pw2qmcpack``
functions.

::

  scf = generate_pwscf(
      identifier   = 'scf',
      ...
      )
  sims.append(scf)

  p2q = generate_pw2qmcpack(
      identifier   = 'p2q',
      ...
      dependencies = (scf,'orbitals')
      )
  sims.append(p2q)

Finally, objects representing QMCPACK simulations are constructed with the ``generate_qmcpack`` function:

Finally, objects representing QMCPACK simulations are constructed with
the ``generate_qmcpack`` function:

::

  opt = generate_qmcpack(
      identifier   = 'opt',
      ...
      jastrows     = [('J1','bspline',8,4.5),
                      ('J2','pade',0.5,0.5)],
      calculations = opt_calcs,
      dependencies = (p2q,'orbitals')
      )
  sims.append(opt)

  qmc = generate_qmcpack(
      identifier   = 'qmc',
      ...
      jastrows     = [],
      calculations = qmc_calcs,
      dependencies = [(p2q,'orbitals'),
                      (opt,'jastrow')]
      )
  sims.append(qmc)

Shared details such as the run directory, job, pseudopotentials, and
orbital file have been omitted (``...``). The “``opt``” run will
optimize a 1-body B-spline Jastrow with 8 knots having a cutoff of 4.5
Bohr and a 2-body Padé Jastrow with up-up and up-down “``B``” parameters
set to 0.5 1/Bohr. The Jastrow list for the DMC run is empty and a new
keyword is present: ``dependencies``. The usage of ``dependencies``
above indicates that the DMC run depends on the optimization run for the
Jastrow factor. Nexus will submit the “``opt``” run first and upon
completion it will scan the output, select the optimal set of
parameters, pass the Jastrow information to the “``qmc``” run and then
submit the DMC job. Independent job workflows are submitted in parallel
when permitted. No input files are written or job submissions made until
the “``run_project``” function is reached.

As written, ``oxygen_dimer.py`` will only perform calculations at the
equilibrium separation distance of 1.2074 Angstrom. Modify the file now
to perform DMC calculations across a range of separation distances with
each DMC run using the Jastrow factor optimized at the equilibrium
separation distance. The necessary Python ``for`` loop syntax should
look something like this:

::

  sims = []
  for scale in [1.00,0.90,0.95,1.05,1.10]:
      ...
      dimer = ...
      if scale==1.00:
          opt = ...
          ...
      #end if
      qmc = ...
      ...
  #end for
  run_project(sims)

Note that the text inside the ``for`` loop and the ``if`` block must be
indented by precisely four spaces. If you use Emacs, changes in
indentation can be performed easily with ``Cntrl-C >`` and ``Cntrl-C <``
after highlighting a block of text (other editors should have similar
functionality).

Change the “``status_only``” parameter in the “``settings``” function to
``1`` and type “./oxygen_dimer.py” at the command line. This will print
the status of all simulations:

.. code:: rest

  Project starting
    checking for file collisions
    loading cascade images
      cascade 0 checking in
    checking cascade dependencies
      all simulation dependencies satisfied
    cascade status
      setup, sent_files, submitted, finished, got_output, analyzed
      000000  scf  ./scale_1.0
      000000  scf  ./scale_0.9
      000000  scf  ./scale_0.95
      000000  scf  ./scale_1.05
      000000  scf  ./scale_1.1
      000000  p2q  ./scale_1.0
      000000  p2q  ./scale_0.9
      000000  p2q  ./scale_0.95
      000000  p2q  ./scale_1.05
      000000  p2q  ./scale_1.1
      000000  opt  ./scale_1.0
      000000  qmc  ./scale_1.0
      000000  qmc  ./scale_0.9
      000000  qmc  ./scale_0.95
      000000  qmc  ./scale_1.05
      000000  qmc  ./scale_1.1
      setup, sent_files, submitted, finished, got_output, analyzed

In this case, a single independent simulation “cascade” (workflow) has
been identified, containing one “``opt``” and five dependent “``qmc``”
runs. The six status flags (``setup``, ``sent_files``, ``submitted``,
``finished``, ``got_output``, ``analyzed``) each show ``0``, indicating
that no work has been done yet.

Now change “``status_only``” back to ``0``, set “``generate_only``” to
``1``, and run ``oxygen_dimer.py`` again. This will perform a dry-run of
all simulations. The dry-run should finish in about 20 seconds:

.. code:: rest

  Project starting
    checking for file collisions
    loading cascade images
      cascade 0 checking in
    checking cascade dependencies
      all simulation dependencies satisfied

    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    poll 0  memory 88.54 MB
      Entering ./scale_1.0 0
        writing input files  0 opt
      Entering ./scale_1.0 0
        sending required files  0 opt
        submitting job  0 opt
      Entering ./scale_1.0 1
        Would have executed:  qsub --mode script --env BG_SHAREDMEMSIZE=32 opt.qsub.in

    poll 1  memory 88.54 MB
      Entering ./scale_1.0 0
        copying results  0 opt
      Entering ./scale_1.0 0
        analyzing  0 opt

    poll 2  memory 88.87 MB
      Entering ./scale_1.0 1
        writing input files  1 qmc
      Entering ./scale_1.0 1
        sending required files  1 qmc
        submitting job  1 qmc
      ...
      Entering ./scale_1.0 2
        Would have executed:  qsub --mode script --env BG_SHAREDMEMSIZE=32 qmc.qsub.in
      ...

  Project finished

Nexus polls the simulation status every 3 seconds and sleeps in between.
The ``scale_`` directories should now contain several files:

.. code:: rest

  scale_1.0
  ├── O2.pwscf.h5
  ├── O.BFD.xml
  ├── opt.in.xml
  ├── opt.qsub.in
  ├── qmc.in.xml
  ├── qmc.qsub.in
  ├── sim_opt
  │   ├── analyzer.p
  │   ├── input.p
  │   └── sim.p
  └── sim_qmc
      ├── analyzer.p
      ├── input.p
      └── sim.p

Take a minute to inspect the generated input (``opt.in.xml``,
``qmc.in.xml``) and submission (``opt.qsub.in``, ``qmc.qsub.in``) files.
The pseudopotential file ``O.BFD.xml`` has been copied into each local
directory. Two additional directories have been created: ``sim_opt`` and
``sim_qmc``. The ``sim.p`` files in each directory contain the current
status of each simulation. If you run ``oxygen_dimer.py`` again, it
should not attempt to rerun any of the simulations:

.. code:: rest

  Project starting
    checking for file collisions
    loading cascade images
      cascade 0 checking in
      cascade 8 checking in
      cascade 2 checking in
      cascade 4 checking in
      cascade 6 checking in
    checking cascade dependencies
      all simulation dependencies satisfied

    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    poll 0  memory 60.10 MB
  Project finished

This way one can continue to add to the ``oxygen_dimer.py`` file (*e.g.*
adding more separation distances) without worrying about duplicate job
submissions.

Now actually submit the optimization and DMC jobs. Reset the state of
the simulations by removing the ``sim.p`` files
(“``rm ./scale*/sim*/sim.p``”), set “``generate_only``” to ``0``, and
rerun ``oxygen_dimer.py``. It should take about 15 minutes for all the
jobs to complete. You may wish to open another terminal to monitor the
progress of the individual jobs while the current terminal runs
``oxygen_dimer.py`` in the foreground.

After completion, try expanding to the full set of “``scale``”
values:``[0.90,0.925,0.95,0.975,1.00,1.025,1.05,1.075,1.10]``. Nexus
should only run the workflows that correspond to new values. In this
way, Nexus scripts can be expanded over time to include new aspects of
growing projects.

.. _excited:

Example 5: Bulk Diamond Excited States VMC
------------------------------------------

The files for this example are found in:

.. code:: rest

  /your_download_path/nexus/examples/qmcpack/rsqmc_misc/excited

Please study `Lab 5`_ in QMCPACK manual for an in-depth discussion of the
excited states calculations. The primitive cell for a structure is not
unique, therefore it must be standardized to generate reproducible
Brillouin zones (BZ). For this end, we use "SeeK-Path" python libraries
to generate standardized primitive cells and generate band structure
paths. Therefore, SeeK-Path libraries must be installed using "pip
install seekpath", prior to running scripts in this example. In this
example, input file script "vmc.py" is identical to "optical.py" in Lab5
of the QMCPACK manual.

We use neutral primitive cells at the wavefunction generation. However,
creation / annihilation operations are applied at VMC level to simulate
optical excitations. Compared to the ground state bulk calculations, a
tiling matrix that is commensurate with the wavevectors involved in the
excitation must be chosen. This process has been automatized in Nexus
using the "get_band_tiling" function. There are two VMC scripts in this
lab that generate the tiling matrix in different ways: ``vmc.py`` script 
uses a non-optimal tiling matrix from Lab 5 in QMCPACK, whereas 
``vmc-opt-tiling.py`` uses the "get_band_tiling" function. In this 
example, we will use ``vmc-opt-tiling.py``. Note, there is also an 
additional VMC script included ``vmc_excitation_alternatives.py`` which
does not use a tiling matrix but includes a variety of ways that
excitations can be specified with Nexus.

In `Lab 5 <https://qmcpack.readthedocs.io/en/develop/lab_excited.html>`_ of the QMCPACK manual we found that VBM is located at
:math:`\Gamma` and the CBM is located at :math:`\Delta` ([0.377, 0.,
0.377]), hence C-diamond is an indirect band gap material. However,
studying this transition at the given exact points would involve a very
large simulation cell. Therefore, we try to round :math:`\Delta` to the
nearest fraction with integer denominator. Tolerance of this rounding
procedure is controlled by "max_volfac" keyword in the "get_band_tiling"
function. "max_volfac" basically controls the maximum volume factor of
the supercell being searched, which inherently controls the highest
k-point grid density in one dimension.

| In the "qmc" object given below, pay attention to keyword "excitation
  = [’up’, ’0 3 1 4’]". Value for the excitation keyword can be given in
  list or tuple format. Here, first element of the list (’up’) denotes
  the spin channel used to perform the excitation. Whereas ’0 3 1 4’
  denotes the twist-index (ti) and band-index (bi) of annihilation and
  creation (an/cr) operations in :math:`ti_{an} bi_{an} ti_{cr} bi_{cr}`
  format. Another option is to use an indexing of the orbitals depending
  on their energetic ordering. In this example this would correspond to
  "excitation = [’up’, ’-11 +12’]". Band/twist index and energy indexes
  of the orbitals can be found in "einspline" files or they can be
  determined after parsing the "nscf.out" file using PwscfAnalyzer.
  In addition to these options, "excitation = ['up','lowest']" can also 
  be specified which will execute a homo-lumo excitation based on the
  energetic ordering of the orbitals. Nexus also allows singlet and
  triplet excitation types. Please refer to ``vmc_excitation_alternatives.py``
  for examples using the various excitation types.
  Examples are also provided in Lab 5 of the QMCPACK manual.

::

  #! /usr/bin/env python3

  from nexus import settings,job,run_project
  from nexus import generate_physical_system
  from nexus import generate_pwscf
  from nexus import generate_pw2qmcpack
  from nexus import generate_qmcpack,vmc
  from structure import *

  settings(
      pseudo_dir    = '../pseudopotentials',
      status_only   = 0,
      generate_only = 0,
      sleep         = 3,
      machine       = 'ws16'
      )

  #Input structure
  dia = generate_physical_system(
      units  = 'A',
      axes   = [[ 1.785,  1.785,  0.   ],
                [ 0.   ,  1.785,  1.785],
                [ 1.785,  0.   ,  1.785]],
      elem   = ['C','C'],
      pos    = [[ 0.    ,  0.    ,  0.    ],
                [ 0.8925,  0.8925,  0.8925]],
      C      = 4
      )

  # Standardized Primitive cell -- run rest of the calculations on this cell
  dia2_structure   = get_primitive_cell(structure=dia.structure)['structure']
  # get_band_tiling and get_primitiev_cell require "SeeK-path" python libraries

  # Returns commensurate optimum tiling matrix for the kpoints_rel wavevectors
  # max_volfac is default to 20
  tiling = get_band_tiling(structure   = dia2_structure,
                           kpoints_rel = [[0.000, 0.000, 0.000],
                                          [0.3768116, 0.,        0.3768116]])
  # Numerical value for the second wavevector can be different in different computers,
  # adjust accordingly

  # All wavevectors must be on the iBZ k_path! (get_kpath output in band.py)

  # K-points can also be defined using their labels
  # (for high symmetry k-points only). For example:
  # tiling = get_band_tiling(structure   = dia2_structure,
  #                         kpoints_label = ['GAMMA', 'X'])

  dia2 = generate_physical_system(
      structure    = dia2_structure,
      kgrid  = (1,1,1),
      kshift = (0,0,0), # Assumes we study transitions from Gamma.
      #For non-gamma tilings, use kshift appropriately
      tiling = tiling,
      C            = 4,
      )

  scf = generate_pwscf(
      identifier   = 'scf',
      path         = 'diamond/scf',
      job          = job(nodes=1, app='pw.x',hours=1),
      input_type   = 'generic',
      calculation  = 'scf',
      nspin        = 2,
      input_dft    = 'lda',
      ecutwfc      = 200,
      conv_thr     = 1e-8,
      nosym        = True,
      wf_collect   = True,
      system       = dia2,
      tot_magnetization = 0,
      pseudos      = ['C.BFD.upf'],
      )

  nscf = generate_pwscf(
      identifier   = 'nscf',
      path         ='diamond/nscf_opt_tiling',
      job          = job(nodes=1, app='pw.x',hours=1),
      input_type   = 'generic',
      calculation  = 'nscf',
      input_dft    = 'lda',
      ecutwfc      = 200,
      nspin        = 2,
      conv_thr     = 1e-8,
      nosym        = True,
      wf_collect   = True,
      system       = dia2,
      nbnd         = 8,      #a sensible nbnd value can be given
      verbosity    = 'high', #verbosity must be set to high
      pseudos      = ['C.BFD.upf'],
      dependencies = (scf, 'charge_density'),
  )

  conv = generate_pw2qmcpack(
      identifier   = 'conv',
      path         = 'diamond/nscf_opt_tiling',
      job          = job(cores=1,app='pw2qmcpack.x',hours=1),
      write_psir   = False,
      dependencies = (nscf,'orbitals'),
      )

  qmc = generate_qmcpack(
      det_format     = 'old',
      identifier     = 'vmc',
      path           = 'diamond/vmc_opt_tiling',
      job            = job(cores=16,threads=16,app='qmcpack',hours = 1),
      input_type     = 'basic',
      spin_polarized = True,
      system         = dia2,
      pseudos        = ['C.BFD.xml'],
      jastrows       = [],
      calculations   = [
          vmc(
              walkers     =  16,
              warmupsteps =  20,
              blocks      = 1000,
              steps       =  10,
              substeps    =   2,
              timestep    =  .4
              )
          ],
      dependencies = (conv,'orbitals'),
      )

  qmc_optical = generate_qmcpack(
      det_format     = 'old',
      identifier     = 'vmc',
      path           = 'diamond/vmc_opt_tiling_optical',
      job            = job(cores=16,threads=16,app='qmcpack',hours = 1),
      input_type     = 'basic',
      spin_polarized = True,
      system         = dia2,
      excitation     = ['up', '0 3 1 4'], #
      pseudos        = ['C.BFD.xml'],
      jastrows       = [],
      calculations   = [
          vmc(
              walkers     =  16,
              warmupsteps =  20,
              blocks      = 1000,
              steps       =  10,
              substeps    =   2,
              timestep    =  .4
              )
          ],
      dependencies = (conv,'orbitals'),
      )

  run_project(scf,nscf,conv,qmc)

