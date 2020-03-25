Nexus PySCF+QMCPACK Example 1: Water molecule (RHF)
===================================================

In this example, we show how to run a very simple PySCF restricted 
Hartree-Fock calculation with Nexus for an all-electron water molecule.

The Nexus script for this example is ``h2o_ae_hf.py``.

Similar to Nexus, the user interface to PySCF is in the form of pure 
Python scripts.  Because a PySCF script can come in a wide variety of 
custom forms (far exceeding that of a Quantum Espresso input file for 
example), Nexus takes the conservative approach of only adding portions 
to an existing user-written PySCF template script.

A simple template script that performs a restricted Hartree-Fock 
calculation is used for this example.  The contents of this script 
(``scf_template.py``) are shown below:

.. code-block:: python

    #! /usr/bin/env python3
    
    from pyscf import scf
    
    # Nexus expands this with Mole info
    $system
    
    mf = scf.RHF(mol)
    mf.kernel()

If you are familiar with PySCF, you will notice that information about 
the molecular structure and the gaussian basisset have been suppressed. 
Instead, ``$system`` appears in the template PySCF script.  This variable 
will be expanded by Nexus prior to the actual execution of the script. 

The Nexus script, ``h2o_ae_hf.py``, is shown below:

.. code-block:: python

    #! /usr/bin/env python3
    
    from nexus import settings,job,run_project,obj
    from nexus import generate_physical_system
    from nexus import generate_pyscf
    
    settings(
        results = '',
        sleep   = 3,
        machine = 'ws16',
        )
    
    system = generate_physical_system(
        structure = 'H2O.xyz',
        )
    
    scf = generate_pyscf(
        identifier = 'scf',               # log output goes to scf.out
        path       = 'h2o_ae_hf',         # directory to run in
        job        = job(serial=True),    # pyscf must run serially         
        template   = './scf_template.py', # pyscf template file
        system     = system,
        mole       = obj(                 # used to make Mole() inputs
            basis    = 'ccpvtz',
            symmetry = True,
            ),
        )
    
    run_project()

The overall flow of the script should be familiar from the Quantum Espresso 
examples presented earlier.  If it is not, please revisit these examples before 
proceeding.

One important difference is how the physical system is specified.  In the 
prior QE examples, the atomic positions, etc., were provided explicitly 
within the Nexus script.  In this case we instead use an external file, 
``H2O.xyz``:

.. code-block:: bash

    3
    
    O  0.000000  0.000000  0.000000 
    H  0.000000  0.757160  0.586260
    H  0.000000  0.757160 -0.586260

Additionally all use of MPI (via mpirun) in the job execution is suppressed 
by specifying ``job(serial=True)``.  PySCF can make use of threads, however, 
but this is unimportant for a simple molecule like H2O in open boundary 
conditions.

Finally, information regarding the basisset and use of symmetry are provided 
via the ``mole`` keyword to ``generate_pyscf``.  All other variables that 
can be assigned to the object outputted by the PySCF function ``gto.Mole()`` 
can be supplied there.

Let's run the example now:

.. code-block:: bash

    >./h2o_ae_hf.py 

    ...
    
    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 101.58 MB 
      ...
      Entering ./runs/h2o_ae_hf 0 
        Executing:  
          export OMP_NUM_THREADS=1
          python scf.py 
      ...
    Project finished

Next, let's look at the PySCF script produced by Nexus (see 
``./runs/h2o_ae_hf/scf.py``):

.. code-block:: python

    #! /usr/bin/env python3
    
    from pyscf import scf
    
    # Nexus expands this with Mole info
    
    ### generated system text ###
    from pyscf import gto as gto_loc
    mol = gto_loc.Mole()
    mol.atom     = '''
                   O    0.00000000   0.00000000   0.00000000
                   H    0.00000000   0.75716000   0.58626000
                   H    0.00000000   0.75716000  -0.58626000
                   '''
    mol.basis    = 'ccpvtz'
    mol.unit     = 'A'
    mol.charge   = 0
    mol.spin     = 0
    mol.symmetry = True
    mol.build()
    ### end generated system text ###

    
    mf = scf.RHF(mol)
    mf.kernel()

Information from Nexus' physical system object (from 
``generate_physical_system``) have been populated into ``mol``, including 
the distance units, net charge, net spin, atomic species, and atomic positions. 
The basis and symmetry information, provided separately as described above, 
have also been filled in.

For the PySCF RHF total energy for the all electron water molecule, you 
should get something very similar to the following:

.. code-block:: bash

    >cat runs/h2o_ae_hf/scf.out 

    converged SCF energy = -76.0302783714398

If you want to see what changes are required to run with pseudopotentials 
(BFD in this case) see the other Nexus script in this directory: 
``h2o_pp_hf.py``.

In the next example we will look at how to run PySCF in periodic boundary 
conditions by considering an RHF calculation of the diamond primitive cell. 

