Nexus PySCF+QMCPACK Example 4: Diamond supercell (RHF+VMC/DMC)
==============================================================

In this example, we expand on the prior PySCF RHF pseudopotential 
calculation for diamond in Example 2 to include the necessary steps 
to perform VMC/DMC for a tiled supercell.  Due to cost (essentially 
the slowness of PySCF), the example is restricted to a small 2x1x1 
tiled cell at the Gamma point.  

The Nexus script for this example is ``diamond_pp_hf_gamma.py``. 
The most important changes relative to Example 2 are shown below:
 
.. code-block:: python

    #! /usr/bin/env python3
    
    ...
    from nexus import generate_convert4qmc
    from nexus import generate_qmcpack
    
    settings(
        ...
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
        tiling   = (2,1,1),   # 2x1x1 supercell tiling
        kgrid    = (1,1,1),   # Gamma point supercell
        kshift   = (0,0,0),
        C        = 4,
        )
    system.change_units('B') # currently a bug in pyscf with A units
    
    scf = generate_pyscf(
        ...
        save_qmc   = True,                # save wfn data for qmcpack
        )
    
    c4q = generate_convert4qmc(
        ...
        no_jastrow   = True,
        hdf5         = True,              # use hdf5 format
        dependencies = (scf,'orbitals'),
        )
    
    opt = generate_qmcpack(
        ...
        J2              = True,           # two-body B-spline Jastrow
        qmc             = 'opt',          # oneshift energy opt
        minmethod       = 'oneshiftonly', 
        init_cycles     = 3,              # 3 cycles w/ 0.1 minwalkers
        init_minwalkers = 0.1,
        cycles          = 3,              # 3 cycles w/ 0.5 minwalkers
        samples         = 25600,          # VMC samples used in each cycle
        dependencies    = (c4q,'orbitals'),
        )
    
    qmc = generate_qmcpack(
        ...
        qmc          = 'vmc',             # vmc run, default inputs
        dependencies = [(c4q,'orbitals'),
                        (opt,'jastrow')],
        )
    
    qmc = generate_qmcpack(
        ...
        qmc          = 'dmc',             # dmc run, default inputs
        vmc_samples  = 800,               # 800 dmc walkers from vmc
        eq_dmc       = True,              # add dmc equil run
        dependencies = [(c4q,'orbitals'),
                        (opt,'jastrow')],
        )
    
    run_project()


The physical system has been changed to include the 2x1x1 supercell 
tiling (``tiling=(2,1,1)``).  Orbital conversion, Jastrow optimization, 
VMC and DMC steps have been added.  Block statements inserted in 
the script currently prevent the QMC steps from running.


First, let's run the RHF and orbital conversion steps:

.. code-block:: bash

    >./diamond_pp_hf_gamma.py 
    
    ...
  
    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 102.27 MB 
      ...
      Entering ./runs/diamond/scf 0 
        Executing:  
          export OMP_NUM_THREADS=16
          python scf.py 
    ...
    elapsed time 45.5 s  memory 143.36 MB 
      ...
      Entering ./runs/diamond/scf 1 
        Executing:  
          export OMP_NUM_THREADS=1
          mpirun -np 1 convert4qmc -pyscf scf.h5 -prefix c4q -hdf5 -nojastrow 
    elapsed time 48.5 s  memory 143.37 MB
    ...
    Project finished

Let's briefly check that these steps have completed succesfully:

.. code-block:: bash

    >./diamond_pp_hf_gamma.py --status_only
   
    ...
    
    cascade status 
      setup, sent_files, submitted, finished, got_output, analyzed, failed 
      111111  0  25060     scf     ./runs/diamond/scf  
      111111  0  25228     c4q     ./runs/diamond/scf  
      000000  0  ------    opt     ./runs/diamond/optJ2  
      000000  0  ------    vmc     ./runs/diamond/vmc  
      000000  0  ------    dmc     ./runs/diamond/dmc  
      setup, sent_files, submitted, finished, got_output, analyzed, failed 


The RHF energy has changed somewhat from the primitive cell example 
due to the introduction of a second primitive cell k-point to account 
for the 2x1x1 tiling:

.. code-block:: bash

    >grep 'SCF energy' runs/diamond/scf/scf.out 
    
    converged SCF energy = -10.6033677691633

Next, comment out all the block statements for the Jastrow 
optimization, VMC and DMC runs:

.. parsed-literal::
    opt = generate_qmcpack(
        **\#block        = True,**
        identifier   = 'opt',
        ...
        )
    
    qmc = generate_qmcpack(
        **\#block        = True,**
        identifier   = 'vmc',
        ...
        )
    
    qmc = generate_qmcpack(
        **\#block        = True,**
        identifier   = 'dmc',
        ...
        )


And then run the QMC portions:

.. code-block:: bash

    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 143.21 MB 
      ...
      Entering ./runs/diamond/optJ2 2 
        Executing:  
          export OMP_NUM_THREADS=4
          mpirun -np 4 qmcpack opt.in.xml 
    ...
    elapsed time 271.3 s  memory 145.63 MB 
      ...
      Entering ./runs/diamond/vmc 3 
        Executing:  
          export OMP_NUM_THREADS=4
          mpirun -np 4 qmcpack vmc.in.xml 
    ...
    elapsed time 347.6 s  memory 146.83 MB 
      ...
      Entering ./runs/diamond/dmc 4 
        Executing:  
          export OMP_NUM_THREADS=4
          mpirun -np 4 qmcpack dmc.in.xml 
    ...
    elapsed time 1070.4 s  memory 147.44 MB 
    ...
    Project finished

Check that the optimization is passable (the high variance here is 
largely due to the extremely short Jastrow cutoffs necessitated by 
the very small supercell):

.. code-block:: bash

    >qmca -q ev runs/diamond/optJ2/*scalar*
     
                                       LocalEnergy              Variance                ratio 
    runs/diamond/optJ2/opt  series 0  -21.209849 +/- 0.013044   3.741734 +/- 0.054848   0.1764 
    runs/diamond/optJ2/opt  series 1  -21.655938 +/- 0.007163   1.047075 +/- 0.018505   0.0484 
    runs/diamond/optJ2/opt  series 2  -21.677865 +/- 0.006881   1.160119 +/- 0.024401   0.0535 
    runs/diamond/optJ2/opt  series 3  -21.685325 +/- 0.008423   1.178533 +/- 0.023183   0.0543 
    runs/diamond/optJ2/opt  series 4  -21.679097 +/- 0.008870   1.209289 +/- 0.026446   0.0558 
    runs/diamond/optJ2/opt  series 5  -21.674319 +/- 0.007903   1.212488 +/- 0.021813   0.0559

Finally, let's look at the VMC and DMC energies and compare to 
the RHF energy found earlier:

.. code-block:: bash

    >grep 'SCF energy' runs/diamond/scf/scf.out 
    
    converged SCF energy = -10.6033677691633
    #                 x2 = -21.2067355383266

    >qmca -e 20 -q e runs/diamond/vmc/*scalar*
    
    runs/diamond/vmc/vmc  series 0  LocalEnergy =  -21.675239 +/- 0.003276 

    >qmca -e 20 -q e runs/diamond/dmc/*s002*scalar*
    
    runs/diamond/dmc/dmc  series 2  LocalEnergy =  -21.826760 +/- 0.004704

