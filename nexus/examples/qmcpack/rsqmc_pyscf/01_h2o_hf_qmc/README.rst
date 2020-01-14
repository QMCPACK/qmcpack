Nexus PySCF+QMCPACK Example 3: Water molecule (RHF+VMC/DMC)
===========================================================

In this example, we extend Example 1 to add the needed steps to run 
respectable VMC and DMC calculations with QMCPACK.  The Nexus script 
for this example is ``h2o_ae_hf_qmc.py``.  The steps performed by Nexus 
for this example are:

1. Restricted Hartree-Fock calculation with PySCF
2. Orbital conversion with convert4qmc
3. Nuclear cusp correction calculation with QMCPACK
4. Two- and three-body Jastrow optimization with QMCPACK
5. VMC with QMCPACK
6. DMC with QMCPACK


Important parts of the Nexus script, ``h2o_ae_hf_qmc.py``, that differ from 
Example 1 are shown below:

.. code-block:: python
    #! /usr/bin/env python3
    
    ...
    from nexus import generate_convert4qmc
    from nexus import generate_cusp_correction
    from nexus import generate_qmcpack
    
    ...
    
    # perform Hartree-Fock
    scf = generate_pyscf(
        ...
        save_qmc   = True,                # save wfn data for qmcpack
        )
    
    # convert orbitals to QMCPACK format
    c4q = generate_convert4qmc(
        ...
        hdf5         = True,              # use hdf5 format
        dependencies = (scf,'orbitals'),
        )
    
    # calculate cusp correction
    cc = generate_cusp_correction(
        ...
        dependencies = (c4q,'orbitals'),
        )
    
    # collect dependencies relating to orbitals
    orbdeps = [(c4q,'particles'), # pyscf changes particle positions
               (c4q,'orbitals' ),
               (cc ,'cuspcorr' )]
    
    # optimize 2-body Jastrow
    optJ2 = generate_qmcpack(
        ...
        J2                = True,         # 2-body B-spline Jastrow
        J1_rcut           = 4.0,          # 4 Bohr cutoff for J1
        J2_rcut           = 7.0,          # 7 Bohr cutoff for J2
        qmc               = 'opt',        # quartic variance optimization
        cycles            = 6,            # loop max of 6
        alloweddifference = 1e-3,         # increase allowed energy difference
        dependencies      = orbdeps,
        )
    
    # optimize 3-body Jastrow
    optJ3 = generate_qmcpack(
        ...
        J3                = True,         # 3-body polynomial Jastrow
        qmc               = 'opt',        # quartic variance optimization
        cycles            = 6,            # loop max of 6
        alloweddifference = 1e-3,         # increase allowed energy difference
        dependencies      = orbdeps+[(optJ2,'jastrow')],
        )
    
    # run VMC with QMCPACK
    qmc = generate_qmcpack(
        ...
        qmc          = 'vmc',             # vmc run
        dependencies = orbdeps+[(optJ3,'jastrow')],
        )
    
    # run DMC with QMCPACK
    qmc = generate_qmcpack(
        ...
        qmc          = 'dmc',             # dmc run
        eq_dmc       = True,              # add dmc equilibration
        dependencies = orbdeps+[(optJ3,'jastrow')],
        )
    
    run_project()

The addition of ``save_qmc=True`` to ``generate_pyscf`` will cause 
Nexus to add code to the end of the template PySCF file to write out 
orbital information in HDF5 format.  PySCF actually changes the 
atom coordinates (at least with symmetries on), so both particle 
and orbital information is collected from ``convert4qmc``.  

The nuclear cusp correction algorithm employed by QMCPACK takes some 
time for larger systems (no so much for this one) so it is included 
as a separate step as would be done for larger calculations.  The 
cusp information is shared by all subsequent QMCPACK runs.

You will also notice that the QMCPACK runs have ``block=True`` set. 
This prevents them from executing when the Nexus script is run.  In 
the rest of the example below, we will proceed stepwise, but if you 
prefer you can comment out all of the lines containing ``block`` 
and all steps will be run with a single execution of the script.

First, let's check the current run status.  Clearly visible are the 
run directories for the HF, orbital conversion, cusp correction, 
optimization, VMC, and DMC steps:

.. code-block:: bash

    >./h2o_ae_hf_qmc.py --status_only
    
    ...  
  
    cascade status 
      setup, sent_files, submitted, finished, got_output, analyzed, failed 
      000000  0  ------    scf     ./runs/H2O/hf  
      000000  0  ------    c4q     ./runs/H2O/hf  
      000000  0  ------    cusp    ./runs/H2O/cuspcorr  
      000000  0  ------    opt     ./runs/H2O/optJ2  
      000000  0  ------    opt     ./runs/H2O/optJ3  
      000000  0  ------    vmc     ./runs/H2O/vmc  
      000000  0  ------    dmc     ./runs/H2O/dmc  
      setup, sent_files, submitted, finished, got_output, analyzed, failed 


If you run the script as-is, then it will perform the Hartree-Fock, 
orbital conversion, and cusp correction steps:

.. code-block:: bash

    >./h2o_ae_hf_qmc.py
  
    ...
    
    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 102.26 MB 
      ...
      Entering ./runs/H2O/hf 0 
        Executing:  
          export OMP_NUM_THREADS=1
          python scf.py 
      ...
    elapsed time 6.1 s  memory 102.29 MB 
      ...
      Entering ./runs/H2O/hf 1 
        Executing:  
          export OMP_NUM_THREADS=1
          mpirun -np 1 convert4qmc -pyscf scf.h5 -prefix c4q -hdf5 -nojastrow 
      ...
    elapsed time 12.1 s  memory 102.31 MB 
      ...
      Entering ./runs/H2O/cuspcorr 2 
        Executing:  
          export OMP_NUM_THREADS=16
          mpirun -np 1 qmcpack cusp.in.xml 
    ...
    Project finished

Before proceeding to Jastrow optimization, check that these steps have 
completed successfully (indicated by ``0`` for the failure flag):

.. code-block:: bash

    >./h2o_ae_hf_qmc.py --status_only
  
    ...
    
    cascade status 
      setup, sent_files, submitted, finished, got_output, analyzed, failed 
      111111  0  4013      scf     ./runs/H2O/hf  
      111111  0  4028      c4q     ./runs/H2O/hf  
      111111  0  4049      cusp    ./runs/H2O/cuspcorr  
      000000  0  ------    opt     ./runs/H2O/optJ2  
      000000  0  ------    opt     ./runs/H2O/optJ3  
      000000  0  ------    vmc     ./runs/H2O/vmc  
      000000  0  ------    dmc     ./runs/H2O/dmc  
      setup, sent_files, submitted, finished, got_output, analyzed, failed 

Next, comment out the ``block`` variables for the optimization steps:

.. parsed-literal::

    # optimize 2-body Jastrow
    optJ2 = generate_qmcpack(
        **\#block             = True,**
        ...
        )
    
    # optimize 3-body Jastrow
    optJ3 = generate_qmcpack(
        **\#block             = True,**
        ...
        )

Then run the Jastrow optimization.  This will take a few minutes:

.. code-block:: bash

    >./h2o_ae_hf_qmc.py
    
    ...
  
    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 102.22 MB 
      ...
      Entering ./runs/H2O/optJ2 3 
        Executing:  
          export OMP_NUM_THREADS=1
          mpirun -np 16 qmcpack opt.in.xml 
  
    elapsed time 3.0 s  memory 838.70 MB 
    ...
    elapsed time 146.0 s  memory 104.71 MB 
      ...
      Entering ./runs/H2O/optJ3 4 
        Executing:  
          export OMP_NUM_THREADS=1
          mpirun -np 16 qmcpack opt.in.xml 
  
    elapsed time 149.1 s  memory 844.58 MB 
    ...
    elapsed time 638.6 s  memory 105.46 MB 
      Entering ./runs/H2O/optJ3 4 
        copying results  4 opt 
      Entering ./runs/H2O/optJ3 4 
        analyzing  4 opt 
  
    Project finished

When the optimization has finished, use ``qmca`` to check on the quality 
of the Jastrows.  It is generally harder to get a small variance/energy 
ratio for all-electron systems, but in this case a ratio of 0.03 Ha is 
reasonable:

.. code-block:: bash

    >qmca -q ev runs/H2O/optJ*/*scalar*
    
                                   LocalEnergy              Variance                ratio 
    runs/H2O/optJ2/opt  series 0  -76.082747 +/- 0.007535   5.664733 +/- 0.079907   0.0745 
    runs/H2O/optJ2/opt  series 1  -76.260574 +/- 0.009308   3.343816 +/- 0.057349   0.0438 
    runs/H2O/optJ2/opt  series 2  -76.272118 +/- 0.007465   3.446286 +/- 0.079691   0.0452 
    runs/H2O/optJ2/opt  series 3  -76.255950 +/- 0.007688   3.371955 +/- 0.062224   0.0442 
    runs/H2O/optJ2/opt  series 4  -76.269735 +/- 0.006397   3.401114 +/- 0.065665   0.0446 
    runs/H2O/optJ2/opt  series 5  -76.267095 +/- 0.006419   3.507579 +/- 0.103478   0.0460 
     
    runs/H2O/optJ3/opt  series 0  -76.264783 +/- 0.007077   3.503846 +/- 0.152511   0.0459 
    runs/H2O/optJ3/opt  series 1  -76.355497 +/- 0.004231   2.458338 +/- 0.215515   0.0322 
    runs/H2O/optJ3/opt  series 2  -76.355215 +/- 0.005806   2.134565 +/- 0.055624   0.0280 
    runs/H2O/optJ3/opt  series 3  -76.358826 +/- 0.011683   2.295939 +/- 0.135095   0.0301 
    runs/H2O/optJ3/opt  series 4  -76.349593 +/- 0.005223   2.265864 +/- 0.127677   0.0297 
    runs/H2O/optJ3/opt  series 5  -76.369746 +/- 0.005126   2.228445 +/- 0.069403   0.0292 

Next, let's perform VMC with the optimal three-body Jastrow selected by 
Nexus.  First, comment out the ``block`` statement for the VMC calculation:

.. parsed-literal::

    # run VMC with QMCPACK
    qmc = generate_qmcpack(
        **\#block        = True,**
        ...
        )

Then run VMC:

.. code-block:: bash

    >./h2o_ae_hf_qmc.py
    
    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 102.05 MB 
      ...
      Entering ./runs/H2O/vmc 5 
        Executing:  
          export OMP_NUM_THREADS=1
          mpirun -np 16 qmcpack vmc.in.xml 
  
    elapsed time 3.1 s  memory 591.75 MB 
    ...
    elapsed time 112.4 s  memory 105.53 MB 
      ...
    Project finished

Similarly for DMC, comment out the ``block`` statement:

.. parsed-literal::

    # run DMC with QMCPACK
    qmc = generate_qmcpack(
        **\#block        = True,**
        ...
        )

Then run DMC:

.. code-block:: bash

    >./h2o_ae_hf_qmc.py

    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 102.24 MB 
      ...
      Entering ./runs/H2O/dmc 6 
        Executing:  
          export OMP_NUM_THREADS=1
          mpirun -np 16 qmcpack dmc.in.xml 
  
    elapsed time 3.1 s  memory 654.33 MB 
    ...
    elapsed time 167.0 s  memory 105.77 MB 
      ...
    Project finished

Let's take a moment to review the energy gains going from Hartree-Fock 
to VMC and DMC:

.. code-block:: bash

    >grep 'SCF energy' runs/H2O/hf/scf.out 
    
    converged SCF energy = -76.0302783714398
    
    >qmca -e 20 -q e runs/H2O/vmc/*scalar*
    
    runs/H2O/vmc/vmc  series 0  LocalEnergy =  -76.359209 +/- 0.003101 
    
    >qmca -e 20 -q e runs/H2O/dmc/*s002*scalar*
    
    runs/H2O/dmc/dmc  series 2  LocalEnergy =  -76.405950 +/- 0.001801 

Overall VMC with an optimal J3 lowers the energy by about 329(3) mHa 
over RHF, while DMC lowers it by an additional 47(4) mHa (w/o timestep 
extrapolation).

In the next example, we return to the diamond system and extend the 
example to include supercell tiling and VMC/DMC.
