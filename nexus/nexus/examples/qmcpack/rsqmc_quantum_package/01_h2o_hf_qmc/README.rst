Nexus QP+QMCPACK Example 3: Water molecule (RHF+VMC/DMC)
===========================================================

In this example, we extend Example 1 to add the needed steps to run 
respectable VMC and DMC calculations with QMCPACK.  The Nexus script 
for this example is ``h2o_ae_hf_qmc.py``.  The steps performed by Nexus 
for this example are:

1. Restricted Hartree-Fock calculation with QP
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
    scf = generate_quantum_package(
        ...
        save_for_qmcpack = True,          # write wfn h5 file for qmcpack
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
    
    # optimize 2-body Jastrow
    optJ2 = generate_qmcpack(
        ...
        J2              = True,          # 2-body B-spline Jastrow
        J2_rcut         = 8.0,           # 8 Bohr cutoff for J2
        qmc             = 'opt',         # oneshift energy optimization
        minmethod       = 'oneshiftonly',
        init_cycles     = 3,             # 3 cycles w/ 0.1 minwalkers
        init_minwalkers = 0.1,
        cycles          = 3,             # 3 cycles w/ 0.5 minwalkers
        samples         = 25600,         # small set of samples
        dependencies    = [(c4q,'orbitals'),
                           (cc,'cuspcorr')],

        )
    
    # optimize 3-body Jastrow
    optJ3 = generate_qmcpack(
        ...
        J3                = True,         # 3-body polynomial Jastrow
        qmc               = 'opt',        # oneshift energy optimization
        minmethod       = 'oneshiftonly',
        init_cycles     = 3,              # 3 cycles w/ 0.1 minwalkers
        init_minwalkers = 0.1,
        cycles          = 3,              # 3 cycles w/ 0.5 minwalkers
        samples         = 51200,          # moderate set of samples
        dependencies    = [(c4q,'orbitals'),
                           (cc,'cuspcorr'),
                           (optJ2,'jastrow')],
        )
    
    # run VMC with QMCPACK
    qmc = generate_qmcpack(
        ...
        qmc          = 'vmc',             # vmc run
        dependencies = [(c4q,'orbitals'),
                        (cc,'cuspcorr'),
                        (optJ3,'jastrow')],
        )
    
    # run DMC with QMCPACK
    qmc = generate_qmcpack(
        ...
        qmc          = 'dmc',             # dmc run
        eq_dmc       = True,              # add dmc equilibration
        dependencies = [(c4q,'orbitals'),
                        (cc,'cuspcorr'),
                        (optJ3,'jastrow')],
        )
    
    run_project()

The addition of ``save_for_qmcpack=True`` to ``generate_quantum_package`` 
will cause Nexus to add an execution of ``qp_run save_for_qmcpack`` 
immediatly following the SCF run.

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
      000000  0  ------    hf      ./runs/H2O/hf  
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

    >source /home/ubuntu/apps/qp2/quantum_package.rc

    >./h2o_ae_hf_qmc.py
  
    ...
    
    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 102.26 MB 
      ...
      Entering ./runs/H2O/hf 0 
        Executing:  
          export OMP_NUM_THREADS=16
          
          source /home/ubuntu/apps/qp2/quantum_package.rc
          
          mpirun -np 1 qp_run scf h2o.ezfio >hf.out 2>hf.err
          
          qp_run save_for_qmcpack h2o.ezfio >hf_savewf.out 2>hf_savewf.err
    ...  
    elapsed time 12.3 s  memory 102.31 MB 
      ...
      Entering ./runs/H2O/hf 1 
        Executing:  
          export OMP_NUM_THREADS=1
          mpirun -np 1 convert4qmc -QP hf_savewf.out -prefix c4q -hdf5 
    ...
    elapsed time 18.3 s  memory 102.32 MB 
      ...
      Entering ./runs/H2O/cuspcorr 2 
        Executing:  
          export OMP_NUM_THREADS=16
          mpirun -np 1 qmcpack cusp.in.xml 
    ...
    elapsed time 24.4 s  memory 102.65 MB 
    ...
    Project finished


Before proceeding to Jastrow optimization, check that these steps have 
completed successfully (indicated by ``0`` for the failure flag):

.. code-block:: bash

    >./h2o_ae_hf_qmc.py --status_only
  
    ...
    
    cascade status 
      setup, sent_files, submitted, finished, got_output, analyzed, failed 
      111111  0  3122      hf      ./runs/H2O/hf  
      111111  0  3755      c4q     ./runs/H2O/hf  
      111111  0  3774      cusp    ./runs/H2O/cuspcorr  
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
    elapsed time 75.9 s  memory 104.71 MB 
      ...
      Entering ./runs/H2O/optJ3 4 
        Executing:  
          export OMP_NUM_THREADS=1
          mpirun -np 16 qmcpack opt.in.xml 
  
    elapsed time 79.0 s  memory 844.58 MB 
    ...
    elapsed time 203.7 s  memory 105.46 MB 
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
    runs/H2O/optJ2/opt  series 0  -76.062943 +/- 0.029738   5.257077 +/- 0.111715   0.0691 
    runs/H2O/optJ2/opt  series 1  -76.311701 +/- 0.020065   3.675307 +/- 0.253241   0.0482 
    runs/H2O/optJ2/opt  series 2  -76.279908 +/- 0.018269   4.365413 +/- 0.278379   0.0572 
    runs/H2O/optJ2/opt  series 3  -76.288324 +/- 0.019347   4.171734 +/- 0.143760   0.0547 
    runs/H2O/optJ2/opt  series 4  -76.320356 +/- 0.021946   4.348944 +/- 0.172051   0.0570 
    runs/H2O/optJ2/opt  series 5  -76.283239 +/- 0.026230   4.062358 +/- 0.126687   0.0533 
     
    runs/H2O/optJ3/opt  series 0  -76.301726 +/- 0.017285   4.305849 +/- 0.212146   0.0564 
    runs/H2O/optJ3/opt  series 1  -76.343321 +/- 0.021465   2.127170 +/- 0.146499   0.0279 
    runs/H2O/optJ3/opt  series 2  -76.350092 +/- 0.012522   2.506410 +/- 0.242378   0.0328 
    runs/H2O/optJ3/opt  series 3  -76.348510 +/- 0.008705   2.127781 +/- 0.100740   0.0279 
    runs/H2O/optJ3/opt  series 4  -76.357728 +/- 0.011500   2.131485 +/- 0.087828   0.0279 
    runs/H2O/optJ3/opt  series 5  -76.346919 +/- 0.011912   2.195001 +/- 0.083254   0.0288 

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

    >grep SCF runs/H2O/hf/hf.out

    * SCF energy                               -76.03027837147563   
    
    >qmca -e 20 -q e runs/H2O/vmc/*scalar*
    
    runs/H2O/vmc/vmc  series 0  LocalEnergy =  -76.354799 +/- 0.003241
    
    >qmca -e 20 -q e runs/H2O/dmc/*s002*scalar*
    
    runs/H2O/dmc/dmc  series 2  LocalEnergy =  -76.407511 +/- 0.001400

Overall VMC with an optimal J3 lowers the energy by about 325(3) mHa 
over RHF, while DMC lowers it by an additional 53(4) mHa (w/o timestep 
extrapolation).

In the next example, we return to the oxygen dimer system and extend the 
selected-CI example to include VMC and DMC with QMCPACK.
