Nexus PySCF+QMCPACK Example 2: Diamond primitive cell (RHF)
===========================================================

In this example, we show how to run a simple PySCF RHF pseudopotential 
calculation with Nexus for the primitive cell of diamond at the Gamma point.  

The Nexus script for this example is ``diamond_pp_hf_gamma.py``.

A simple template PySCF script that performs an RHF calculation is used 
for this example.  The contents of this script (``scf_template.py``) are 
shown below:

.. code-block:: python

    #! /usr/bin/env python3
    
    from pyscf.pbc import df, scf
    
    $system
    
    gdf = df.FFTDF(cell,kpts)
    gdf.auxbasis = 'weigend'
    
    mf = scf.KRHF(cell,kpts).density_fit()
    mf.exxdiv  = 'ewald'
    mf.with_df = gdf
    mf.kernel()

Similar to the last example, information about the atomic structure and 
gaussian basisset have been suppressed with the ``$system`` placeholder. 

The Nexus script, ``diamond_pp_hf_gamma.py``, is shown below:

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
        units    = 'A',
        axes     = '''1.785   1.785   0.000
                      0.000   1.785   1.785
                      1.785   0.000   1.785''',
        elem_pos = '''
                   C  0.0000  0.0000  0.0000
                   C  0.8925  0.8925  0.8925
                   ''',
        kgrid    = (1,1,1),
        kshift   = (0,0,0),
        C        = 4,
        )
    
    
    scf = generate_pyscf(
        identifier = 'scf',                      # log output goes to scf.out
        path       = 'diamond_pp_hf_gamma',      # directory to run in
        job        = job(serial=True,threads=16),# pyscf must run w/o mpi
        template   = './scf_template.py',        # pyscf template file
        system     = system,
        cell       = obj(                        # used to make Cell() inputs
            basis         = 'bfd-vdz',
            ecp           = 'bfd',
            drop_exponent = 0.1,
            verbose       = 5,
            ),
        )
    
    run_project()

Here, as before, the use of MPI in the job execution is suppressed, but now 
16 OpenMP threads have been requested.  This will help speed up the (very 
slow) PBC calculations.

Instead of ``mole``, the script uses the keyword ``cell``, mirroring the 
usage of ``gto.Cell()`` in PySCF for periodic systems.  Information 
about the pseudopotential and basisset has been placed there.

Now let's run the example.  This may take a few minutes:

.. code-block:: bash

    >./diamond_pp_hf_gamma.py 
    
    ...
  
    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 102.27 MB 
      ...
      Entering ./runs/diamond_pp_hf_gamma 0 
        Executing:  
          export OMP_NUM_THREADS=16
          python scf.py 
  
    elapsed time 3.0 s  memory 168.03 MB 
    ...
    elapsed time 527.8 s  memory 102.24 MB 
      Entering ./runs/diamond_pp_hf_gamma 0 
        copying results  0 scf 
      Entering ./runs/diamond_pp_hf_gamma 0 
        analyzing  0 scf 
  
    Project finished


Next, let's look at the PySCF script produced by Nexus (see 
``./runs/diamond_pp_hf_gamma/scf.py``):

.. code-block:: python

    #! /usr/bin/env python3
    
    from pyscf.pbc import df, scf
    
    ### generated system text ###
    from numpy import array
    from pyscf.pbc import gto as gto_loc
    cell = gto_loc.Cell()
    cell.a             = '''
                         1.78500000   1.78500000   0.00000000
                         0.00000000   1.78500000   1.78500000
                         1.78500000   0.00000000   1.78500000
                         '''
    cell.basis         = 'bfd-vdz'
    cell.dimension     = 3
    cell.ecp           = 'bfd'
    cell.unit          = 'A'
    cell.atom          = '''
                         C    0.00000000   0.00000000   0.00000000
                         C    0.89250000   0.89250000   0.89250000
                         '''
    cell.drop_exponent = 0.1
    cell.verbose       = 5
    cell.charge        = 0
    cell.spin          = 0
    cell.build()
    kpts = array([
        [0.0, 0.0, 0.0]])
    ### end generated system text ###
    
    gdf = df.FFTDF(cell,kpts)
    gdf.auxbasis = 'weigend'
    
    mf = scf.KRHF(cell,kpts).density_fit()
    mf.exxdiv  = 'ewald'
    mf.with_df = gdf
    mf.kernel()

Similar to the prior example, information regarding the atoms, basisset and 
pseudopotentials are populated into ``cell``.  An important addition is the 
``kpts`` array, which holds an explicit list of primitive cell k-points for 
the HF calculation (just Gamma in this case).  This will be more important 
in later examples where we will add QMC calculations for a supercell.  

For the PySCF RHF total energy for the diamond primitive cell, you 
should get something very similar to the following:

.. code-block:: bash

  >tail -n1 runs/diamond_pp_hf_gamma/scf.out
  
  converged SCF energy = -10.2472172154513

In the next example we will return to the water molecule, but now with 
the necessary additional steps to perform VMC with QMCPACK.

VMC calculations for diamond are covered in Example 4.


