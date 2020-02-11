Nexus QP+QMCPACK Example 2: Oxygen dimer (selected-CI)
======================================================

In this example, we show how to run selected-CI (CIPSI) calculations 
with Quantum Package via Nexus for an all-electron spin-polarized 
oxygen dimer.  The steps performed by Nexus in this example are:

1. Hartree-Fock calculation with QP
2. Selected-CI based on HF with QP
3. Generation of selected-CI natural orbitals with QP
4. Selected-CI based on the natural orbitals with QP

The Nexus script for this example is ``o2_selci.py``.  The contents of 
the script are shown below:


.. code-block:: python

    #! /usr/bin/env python3
    
    from nexus import settings,job,run_project
    from nexus import generate_physical_system
    from nexus import generate_quantum_package
    
    # note: you must source the QP config file before running this script
    #   source /home/ubuntu/apps/qp2/quantum_package.rc
    
    settings(
        results       = '',
        status_only   = 0,
        generate_only = 0,
        sleep         = 3,
        machine       = 'ws16',
        qprc          = '/home/ubuntu/apps/qp2/quantum_package.rc',
        )
    
    # define run details
    qp_job = job(cores=16,threads=16)
    
    # read in structure for oxygen dimer
    dimer = generate_physical_system(
        structure = './O2.xyz',
        net_spin  = 2,           # nup-ndown=2
        )
    
    # path, job, system details are shared across runs
    qp_shared = dict(
        path   = 'O_dimer_selected_CI',
        job    = qp_job,
        system = dimer,
        prefix = 'fci', # single shared ezfio, rsync if different
        )
    
    # run Hartree-Fock
    scf = generate_quantum_package(
        identifier            = 'scf',
        run_type              = 'scf',
        ao_basis              = 'aug-cc-pvdz',
        io_ao_two_e_integrals = 'Write', # write 2e integrals
        four_idx_transform    = True,    # compute 2e integrals
        **qp_shared
        )
    
    # initial selected CI run
    fci0 = generate_quantum_package(
        identifier         = 'fci0',
        run_type           = 'fci',
        n_det_max          = 5000, # max determinant count
        save_natorb        = True, # write natural orbitals
        four_idx_transform = True, # compute 2e integrals
        dependencies       = (scf,'other'),
        **qp_shared
        )
    
    # final selected CI based on natural orbitals
    fci = generate_quantum_package(
        identifier    = 'fci',
        run_type      = 'fci',
        n_det_max     = 5000,
        dependencies  = (fci0,'other'),
        **qp_shared
        )
    
    run_project()

In the script above, as before, we are reading the structure from an xyz 
file.  The net spin of the system (number of up minus number of down 
electrons) is set to ``2`` in ``generate_physical_system``.  Most of the 
QP runs share some details, so we store these in a Python ``dict`` called 
``qp_shared`` which is reused in each QP run (see ``**qp_shared``). In 
the HF run, we request that the two electron integrals be computed and 
stored, which is important for the following selected CI run (which is 
denoted confusingly as ``fci`` in QP).  For the selected-CI we use 
a determinant count threshold of 5000.  The selected-CI calculation will 
terminate once this threshold is exceeded.  The request to compute the 
natural orbitals is found where ``save_natorb=True`` and the two electron 
integrals are again recomputed.  Finally, selected-CI is performed based 
on the natural orbitals.

Now let's run the Nexus script and see how these steps are executed:

.. code-block:: bash

    >source /home/ubuntu/apps/qp2/quantum_package.rc

    >./o2_selci.py

    ...
    
    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 101.89 MB 
      ...
      Entering ./runs/O_dimer_selected_CI 0 
        Executing:  
          export OMP_NUM_THREADS=16
          
          source /home/ubuntu/apps/qp2/quantum_package.rc
          
          mpirun -np 1 qp_run scf fci.ezfio >scf.out 2>scf.err
          
          echo "Write" > fci.ezfio/mo_two_e_ints/io_mo_two_e_integrals
          qp_run four_idx_transform fci.ezfio >scf_fit.out 2>scf_fit.err
    ...
    elapsed time 12.4 s  memory 101.98 MB 
      ...
      Entering ./runs/O_dimer_selected_CI 1 
        Executing:  
          export OMP_NUM_THREADS=16
          
          source /home/ubuntu/apps/qp2/quantum_package.rc
          
          mpirun -np 1 qp_run fci fci.ezfio >fci0.out 2>fci0.err
          
          qp_run save_natorb fci.ezfio >fci0_natorb.out 2>fci0_natorb.err
          
          echo "Write" > fci.ezfio/mo_two_e_ints/io_mo_two_e_integrals
          qp_run four_idx_transform fci.ezfio >fci0_fit.out 2>fci0_fit.err
    ...
    elapsed time 33.6 s  memory 101.98 MB 
      ...
      Entering ./runs/O_dimer_selected_CI 2 
        Executing:  
          export OMP_NUM_THREADS=16
          
          source /home/ubuntu/apps/qp2/quantum_package.rc
          
          mpirun -np 1 qp_run fci fci.ezfio 
    ...
    elapsed time 45.8 s  memory 101.98 MB 
    Project finished

You should obtain variational energies similar to the following for HF, 
CIPSI\@HF and CIPSI\@NO:

.. code-block:: bash

    >grep SCF runs/O_dimer_selected_CI/scf.out 
    
    * SCF energy         -149.6199872983760    
    
    >grep 'E               =' runs/O_dimer_selected_CI/fci0.out | tail -n1
    
     E               =   -149.96976111218555     
    
    >grep 'E               =' runs/O_dimer_selected_CI/fci.out | tail -n1
    
     E               =   -149.98213334936918     

With PT2 corrections, better (but non-variational) estimates of the ground 
state energy within this basis can be obtained:

.. code-block:: bash

    >grep 'E+PT2            =' runs/O_dimer_selected_CI/fci0.out | tail -n1
    
     E+PT2            =   -150.02457005565802       +/-    1.7052281470379021E-004
    
    >grep 'E+PT2            =' runs/O_dimer_selected_CI/fci.out | tail -n1
    
     E+PT2            =   -150.02759522657587       +/-    7.0329356682808259E-005

The only significant change to this example that is needed to obtain 
production level results is to perform a series of calculations with 
increasing maximum determinant counts until convergence is reached.

In the next example, we return to the water molecule and add the necessary 
steps to perform VMC and DMC with QMCPACK.

