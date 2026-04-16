Nexus QP+QMCPACK Example 4: Oxygen dimer (selected-CI+VMC/DMC)
==============================================================

In this example, we extend Example 2 to add the needed steps to run 
respectable VMC and DMC calculations with QMCPACK.  The Nexus script 
for this example is ``o2_selci_vmc_dmc.py``.  The steps performed by 
Nexus for this example are:

1. Hartree-Fock calculation with QP
2. Selected-CI based on HF with QP
3. Generation of selected-CI natural orbitals with QP
4. Selected-CI based on the natural orbitals with QP
5. Wavefunction conversion with convert4qmc
6. Nuclear cusp correction calculation with QMCPACK
7. Two- and three-body Jastrow optimization with QMCPACK
8. VMC with QMCPACK
9. DMC with QMCPACK

Important parts of the Nexus script, ``o2_selci_vmc_dmc.py``, that differ 
from Example 2 are shown below:

.. code-block:: python
    #! /usr/bin/env python3
    
    ...
    from nexus import generate_convert4qmc
    from nexus import generate_cusp_correction
    from nexus import generate_qmcpack
    
    ...
    
    # run Hartree-Fock
    scf = generate_quantum_package(
        ...
        )
    
    # initial selected CI run
    fci0 = generate_quantum_package(
        ...
        )
    
    # final selected CI based on natural orbitals
    fci = generate_quantum_package(
        ...
        save_for_qmcpack = True,          # write wfn h5 file for qmcpack
        )
    
    # convert orbitals to QMCPACK format
    c4q = generate_convert4qmc(
        ...
        hdf5         = True,              # use hdf5 format
        dependencies = (fci,'orbitals'),
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
        vmc_samples  = 1024,              # dmc walkers from vmc
        dependencies = [(c4q,'orbitals'),
                        (cc,'cuspcorr'),
                        (optJ3,'jastrow')],
        )
    
    run_project()

The addition of ``save_for_qmcpack=True`` to ``generate_quantum_package`` 
will cause Nexus to add an execution of ``qp_run save_for_qmcpack`` 
immediatly following the selected-CI run.

The nuclear cusp correction algorithm employed by QMCPACK takes some 
time for multideterminant systems so it is included as a separate step 
here. This saves several minutes from the overall calculation process 
as the cusp information is shared by all subsequent QMCPACK runs.

You will also notice that the QMCPACK runs have ``block=True`` set. 
This prevents them from executing when the Nexus script is run.  In 
the rest of the example below, we will proceed stepwise, but if you 
prefer you can comment out all of the lines containing ``block`` 
and all steps will be run with a single execution of the script.

First, let's check the current run status.  Clearly visible are the 
run directories for the HF, selected-CI, wavefunction conversion, 
cusp correction, optimization, VMC, and DMC steps:

.. code-block:: bash

    >./o2_selci_vmc_dmc.py --status_only
    
    ...  
  
    cascade status 
      setup, sent_files, submitted, finished, got_output, analyzed, failed 
      000000  0  ------    scf     ./runs/O_dimer/selci  
      000000  0  ------    fci0    ./runs/O_dimer/selci  
      000000  0  ------    fci     ./runs/O_dimer/selci  
      000000  0  ------    c4q     ./runs/O_dimer/selci  
      000000  0  ------    cusp    ./runs/O_dimer/cuspcorr  
      000000  0  ------    opt     ./runs/O_dimer/optJ2  
      000000  0  ------    opt     ./runs/O_dimer/optJ3  
      000000  0  ------    vmc     ./runs/O_dimer/vmc  
      000000  0  ------    dmc     ./runs/O_dimer/dmc  
      setup, sent_files, submitted, finished, got_output, analyzed, failed 


If you run the script as-is, then it will perform the Hartree-Fock, 
selected-CI, wavefunction conversion, and cusp correction steps:

.. code-block:: bash

    >source /home/ubuntu/apps/qp2/quantum_package.rc

    >./o2_selci_vmc_dmc.py
  
    ...
    
    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 102.22 MB 
      ..
      Entering ./runs/O_dimer/selci 0 
        Executing:  
          export OMP_NUM_THREADS=16
          
          source /home/ubuntu/apps/qp2/quantum_package.rc
          
          mpirun -np 1 qp_run scf fci.ezfio >scf.out 2>scf.err
          
          echo "Write" > fci.ezfio/mo_two_e_ints/io_mo_two_e_integrals
          qp_run four_idx_transform fci.ezfio >scf_fit.out 2>scf_fit.err
    ...
    elapsed time 12.3 s  memory 102.40 MB 
      ...
      Entering ./runs/O_dimer/selci 1 
        Executing:  
          export OMP_NUM_THREADS=16
          
          source /home/ubuntu/apps/qp2/quantum_package.rc
          
          mpirun -np 1 qp_run fci fci.ezfio >fci0.out 2>fci0.err
          
          qp_run save_natorb fci.ezfio >fci0_natorb.out 2>fci0_natorb.err
          
          echo "Write" > fci.ezfio/mo_two_e_ints/io_mo_two_e_integrals
          qp_run four_idx_transform fci.ezfio >fci0_fit.out 2>fci0_fit.err
    ...
    elapsed time 33.5 s  memory 102.41 MB 
      ...
      Entering ./runs/O_dimer/selci 2 
        Executing:  
          export OMP_NUM_THREADS=16
          
          source /home/ubuntu/apps/qp2/quantum_package.rc
          
          mpirun -np 1 qp_run fci fci.ezfio >fci.out 2>fci.err
          
          qp_run save_for_qmcpack fci.ezfio >fci_savewf.out 2>fci_savewf.err
    ...
    elapsed time 48.8 s  memory 102.41 MB 
      ...
      Entering ./runs/O_dimer/selci 3 
        Executing:  
          export OMP_NUM_THREADS=1
          mpirun -np 1 convert4qmc -QP fci_savewf.out -prefix c4q -hdf5 
    ...
    elapsed time 54.9 s  memory 102.41 MB 
      ...
      Entering ./runs/O_dimer/cuspcorr 4 
        Executing:  
          export OMP_NUM_THREADS=16
          mpirun -np 1 qmcpack cusp.in.xml 
    ...
    elapsed time 235.9 s  memory 102.72 MB 
    ...
    Project finished


Before proceeding to Jastrow optimization, check that these steps have 
completed successfully (indicated by ``0`` for the failure flag):

.. code-block:: bash

    >./o2_selci_vmc_dmc.py --status_only
  
    ...
    
    cascade status 
      setup, sent_files, submitted, finished, got_output, analyzed, failed 
      111111  0  5263      scf     ./runs/O_dimer/selci  
      111111  0  5572      fci0    ./runs/O_dimer/selci  
      111111  0  5947      fci     ./runs/O_dimer/selci  
      111111  0  6257      c4q     ./runs/O_dimer/selci  
      111111  0  6276      cusp    ./runs/O_dimer/cuspcorr  
      000000  0  ------    opt     ./runs/O_dimer/optJ2  
      000000  0  ------    opt     ./runs/O_dimer/optJ3  
      000000  0  ------    vmc     ./runs/O_dimer/vmc  
      000000  0  ------    dmc     ./runs/O_dimer/dmc  
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

Then run the Jastrow optimization.  This will take several minutes:

.. code-block:: bash

    >./o2_selci_vmc_dmc.py
    
    ...
  
    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 102.22 MB 
      ...
      Entering ./runs/O_dimer/optJ2 3 
        Executing:  
          export OMP_NUM_THREADS=1
          mpirun -np 16 qmcpack opt.in.xml 
  
    elapsed time 3.0 s  memory 838.70 MB 
    ...
    elapsed time 392.8 s  memory 104.71 MB 
      ...
      Entering ./runs/O_dimer/optJ3 4 
        Executing:  
          export OMP_NUM_THREADS=1
          mpirun -np 16 qmcpack opt.in.xml 
    ...
    elapsed time 1209.6 s  memory 105.46 MB 
      Entering ./runs/O_dimer/optJ3 4 
        copying results  4 opt 
      Entering ./runs/O_dimer/optJ3 4 
        analyzing  4 opt 
  
    Project finished

When the optimization has finished, use ``qmca`` to check on the quality 
of the Jastrows.  It is generally harder to get a small variance/energy 
ratio for all-electron systems, but in this case a ratio of 0.05 Ha is 
reasonable:

.. code-block:: bash

    >qmca -q ev runs/O_dimer/optJ*/*scalar*
     
                                       LocalEnergy               Variance                 ratio 
    runs/O_dimer/optJ2/opt  series 0  -150.018615 +/- 0.040825   14.594573 +/- 2.613162   0.0973 
    runs/O_dimer/optJ2/opt  series 1  -150.280100 +/- 0.039179    9.031973 +/- 0.485841   0.0601 
    runs/O_dimer/optJ2/opt  series 2  -150.201156 +/- 0.046328    9.124909 +/- 0.607106   0.0608 
    runs/O_dimer/optJ2/opt  series 3  -150.210220 +/- 0.034927    8.957077 +/- 0.476612   0.0596 
    runs/O_dimer/optJ2/opt  series 4  -150.181311 +/- 0.036444   10.004254 +/- 0.614452   0.0666 
    runs/O_dimer/optJ2/opt  series 5  -150.171909 +/- 0.034095    8.589145 +/- 0.323540   0.0572 
     
    runs/O_dimer/optJ3/opt  series 0  -150.209635 +/- 0.029633    9.904128 +/- 0.645692   0.0659 
    runs/O_dimer/optJ3/opt  series 1  -150.272138 +/- 0.026971    7.670001 +/- 0.511987   0.0510 
    runs/O_dimer/optJ3/opt  series 2  -150.234821 +/- 0.024377    7.263859 +/- 0.330053   0.0484 
    runs/O_dimer/optJ3/opt  series 3  -150.235615 +/- 0.022895    7.370637 +/- 0.434385   0.0491 
    runs/O_dimer/optJ3/opt  series 4  -150.253519 +/- 0.018772    7.380013 +/- 0.366188   0.0491 
    runs/O_dimer/optJ3/opt  series 5  -150.268180 +/- 0.019884    7.065689 +/- 0.266124   0.0470 


Next, let's perform VMC with the optimal three-body Jastrow selected by 
Nexus.  First, comment out the ``block`` statement for the VMC calculation:

.. parsed-literal::

    # run VMC with QMCPACK
    qmc = generate_qmcpack(
        **\#block        = True,**
        ...
        )

Then run VMC (be patient, this takes a while):

.. code-block:: bash

    >./o2_selci_vmc_dmc.py
    
    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 102.05 MB 
      ...
      Entering ./runs/O_dimer/vmc 5 
        Executing:  
          export OMP_NUM_THREADS=1
          mpirun -np 16 qmcpack vmc.in.xml 
  
    elapsed time 3.1 s  memory 591.75 MB 
    ...
    elapsed time 2029.0 s  memory 105.53 MB 
      ...
    Project finished

Similarly for DMC, comment out the ``block`` statement:

.. parsed-literal::

    # run DMC with QMCPACK
    qmc = generate_qmcpack(
        **\#block        = True,**
        ...
        )

Then run DMC (be patient, this takes a while):

.. code-block:: bash

    >./o2_selci_vmc_dmc.py

    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 102.24 MB 
      ...
      Entering ./runs/O_dimer/dmc 6 
        Executing:  
          export OMP_NUM_THREADS=1
          mpirun -np 16 qmcpack dmc.in.xml 
  
    elapsed time 3.1 s  memory 654.33 MB 
    ...
    elapsed time 2759.3 s  memory 105.77 MB 
      ...
    Project finished

Let's take a moment to review the energy gains going from Hartree-Fock 
to variational selected-CI to selected-CI+PT2 to VMC and DMC:

.. code-block:: bash

    >grep SCF runs/O_dimer/hf/scf.out 
    
    * SCF energy          -149.6199872983760    
    
    >grep 'E               =' runs/O_dimer/selci/fci0.out | tail -n1
    
     E                =   -149.96976111218555     
    
    >grep 'E               =' runs/O_dimer/selci/fci.out | tail -n1
    
     E                =   -149.98213334936918     

    >grep 'E+PT2            =' runs/O_dimer/selci/fci0.out | tail -n1
    
     E+PT2            =   -150.02457005565802       +/-    1.7052281470379021E-004
    
    >grep 'E+PT2            =' runs/O_dimer/selci/fci.out | tail -n1
    
     E+PT2            =   -150.02759522657587       +/-    7.0329356682808259E-005


    >qmca -e 20 -q e runs/O_dimer/vmc/*scalar*
    
    runs/O_dimer/vmc/vmc  series 0  LocalEnergy =  -150.240194 +/- 0.006105 

    >qmca -e 20 -q e runs/O_dimer/dmc/*s002*scalar*
    
    runs/O_dimer/dmc/dmc  series 2  LocalEnergy =  -150.331795 +/- 0.003518 

These energies are lower than RHF by:

.. code-block:: bash
    
    CIPSI@HF      350    mHa
    CIPSI@NO      362    mHa
    CIPSI@HF+PT2  405    mHa
    CIPSI@NO+PT2  408    mHa
    VMC@CIPSI@NO  620(6) mHa
    DMC@CIPSI@NO  712(4) mHa
    
As before, to take this example to full production, converge the total 
energies with respect to the number of determinants.  In general, the 
DMC energy should converge faster than the CIPSI energy.
