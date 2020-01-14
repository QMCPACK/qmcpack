Nexus QP+QMCPACK Example 1: Water molecule (RHF)
================================================

In this example, we show how to run a simple restricted Hartree-Fock 
calculation with Quantum Package via Nexus for an all-electron water 
molecule.

Quantum Package takes a different philosophy from many electronic 
structure codes.  QP doesn't have a single file that serves as input. 
Instead, QP uses a directory tree with single text files storing single 
input variable values, as well as general program state.  Also unusual, 
QP presents an assemblage of executables to interact with the input 
and drive stages of the calculation. Nexus wraps several of these steps 
together, and in the end the Nexus interface to QP resembles its interfaces 
to other Gaussian-based codes like GAMESS or PySCF.

The Nexus script for this example is ``h2o_ae_hf.py``.  The contents of 
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
    
    scf_job = job(cores=16,threads=16)
    
    system = generate_physical_system(
        structure = 'H2O.xyz',
        )
    
    scf = generate_quantum_package(
        identifier   = 'hf',        # log output goes to hf.out
        path         = 'h2o_ae_hf', # directory to run in
        job          = scf_job,
        system       = system,
        prefix       = 'h2o',       # create/use h2o.ezfio
        run_type     = 'scf',       # qprun scf h2o.ezfio
        ao_basis     = 'cc-pvtz',   # use cc-pvtz basis
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

The QP specific inputs are ``prefix``, which specifies the name of the 
``ezfio`` directory (the input "file"), ``run_type``, which in this case 
selects a Hartree-Fock run ('scf'), and ``ao_basis``, which selects the 
Gaussian basis set.

Before running any Nexus script that involves QP (or before running QP in 
general) you need to ensure that the QP executables are all visible to the 
executing shell.  This is done by sourcing the QP configuration file:

.. code-block:: bash

    >source /home/ubuntu/apps/qp2/quantum_package.rc

Now lets run the Nexus script:

.. code-block:: bash

    >./h2o_ae_hf.py 

    ...
    
    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 101.90 MB
      Entering ./runs/h2o_ae_hf 0
        writing input files  0 hf
      ...
      Entering ./runs/h2o_ae_hf 0
        Executing:
          export OMP_NUM_THREADS=16
  
          source /home/ubuntu/apps/qp2/quantum_package.rc
  
          mpirun -np 1 qp_run scf h2o.ezfio
  
    elapsed time 3.2 s  memory 142.65 MB
    ...
    Project finished

The creation of the ``ezfio`` directory (via ``qp_create_ezfio`` and 
``qp_edit -c``) was performed during the "writing input files" step shown 
above.  You will also notice that the configuration file is sourced again 
by Nexus.  This is done both here (local workstation) and within batch job 
submission scripts on supercomputers because a new shell is entered. 

Because the ``ezfio`` directory/file contents are in constant flux and 
represent only current program state, Nexus stores changes it makes to the 
QP input in text files as a record of these actions, partially fulfilling 
the role a traditional input file would.  This "input file" for the HF run 
is shown below:

.. code-block:: bash

    >cat runs/h2o_ae_hf/hf.in 
  
    ao_basis
      ao_basis        = cc-pvtz
    end ao_basis
    determinants
      n_det_max       = 5000
    end determinants
    electrons
      elec_alpha_num  = 5
      elec_beta_num   = 5
    end electrons
    run_control
      four_idx_transform = False
      postprocess     = []
      prefix          = h2o
      run_type        = scf
      save_for_qmcpack = False
      save_natorb     = False
      sleep           = 30
    end run_control

In this case, the direct inputs are ``ao_basis``, ``elec_alpha_num``, and 
``elec_beta_num``.  The electron counts have been inferred from the Nexus 
physical system object.  The number of determinants is a default value and 
is not active in the present SCF case.  The contents of ``run_control`` 
relate to Nexus' actions with ``qp_run`` and related commands, in this case 
noting the ``ezfio`` file prefix and the run type as "scf".

For the QP RHF total energy for the all electron water molecule, you 
should get something very similar to the following:

.. code-block:: bash

    >grep SCF runs/h2o_ae_hf/hf.out 
    
    * SCF energy                                        -76.03027837147572  

In the next example we will move beyond Hartree-Fock to perform selected-CI
(CIPSI) calculations with QP and Nexus for a spin polarized oxygen dimer.

