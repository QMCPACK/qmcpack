Nexus QE+QMCPACK Example 3: Diamond cubic cell (twist averaged VMC)
===================================================================

In this example, we expand on the previous example to include twist 
averaging.  The NSCF calculation is updated to generate orbitals for 
a 2x2x2 supercell twist grid.  Using a Jastrow factor optimized only 
at the Gamma point, VMC calculations are performed for each twist.

The Nexus script used for this example is ``diamond_lda_vmc_twistavg.py``.
Important differences between this script and the Gamma-point only 
one from the prior example are shown below:

.. parsed-literal::

    system = generate_physical_system(
        ...
        **kgrid  = (2,2,2), \# 2x2x2 supercell twist grid
        kshift = (0,0,0),**
        ...
        )
    
    scf = generate_pwscf(
        ...
        )
    
    nscf = generate_pwscf(
        **path = 'diamond/nscf_twist', \# new nscf run**
        ...
        )
    
    conv = generate_pw2qmcpack(
        **path = 'diamond/nscf_twist', \# new orb conversion**
        ...
        )
    
    opt = generate_qmcpack(
        ...
        )
    
    qmc = generate_qmcpack(
        **path = 'diamond/vmc_twist', \# new vmc run
        \# run with 8 mpi tasks, one per twist
        job  = job(cores=16,threads=2,app='qmcpack'),**
        ...
        )

The NSCF run will be performed in a new directory.  When NSCF is 
performed in a separate directory, as in this case, Nexus copies the 
SCF output into the new NSCF directory via rsync.  In this way, the 
converged charge density from a single SCF run can be re-used to make 
orbitals for a number of different QMC supercells.

Another important difference is that the VMC job has been modified 
so that one mpi task will be assigned to each of the eight twists in 
the batched VMC run.

Before running this example, let's copy in the completed runs from 
the prior example:

.. code-block:: bash

    >rsync -av ../02_qe_diamond_dft_vmc/runs ./

Now confirm that the prior runs are complete from Nexus' point of view:

.. code-block:: bash

    >./diamond_lda_vmc_twistavg.py --status_only
  
    ...
  
    cascade status 
      setup, sent_files, submitted, finished, got_output, analyzed, failed 
      111111  0  17009     scf     ./runs/diamond/scf  
      000000  0  ------    nscf    ./runs/diamond/nscf_twist  
      000000  0  ------    conv    ./runs/diamond/nscf_twist  
      111111  0  17362     opt     ./runs/diamond/optJ2  
      000000  0  ------    vmc     ./runs/diamond/vmc_twist  
      setup, sent_files, submitted, finished, got_output, analyzed, failed 

Next, go ahead and run the NSCF and VMC steps combined:

.. code-block:: bash

    >./diamond_lda_vmc_twistavg.py 
  
    ...
    
    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 102.56 MB 
      ...
      Entering ./runs/diamond/nscf_twist 1 
        Executing:  
          export OMP_NUM_THREADS=1
          mpirun -np 16 pw.x -input nscf.in 
    ...
    elapsed time 6.1 s  memory 102.79 MB 
      ...
      Entering ./runs/diamond/nscf_twist 2 
        Executing:  
          export OMP_NUM_THREADS=1
          mpirun -np 16 pw2qmcpack.x<conv.in 
    ...
    elapsed time 12.2 s  memory 102.79 MB 
      ...
      Entering ./runs/diamond/vmc_twist 4 
        Executing:  
          export OMP_NUM_THREADS=2
          mpirun -np 8 qmcpack vmc.in 
  
    elapsed time 15.3 s  memory 554.23 MB 
    ...
  
    elapsed time 51.8 s  memory 105.50 MB 
      Entering ./runs/diamond/vmc_twist 4 
        copying results  4 vmc 
      Entering ./runs/diamond/vmc_twist 4 
        analyzing  4 vmc 
  
    Project finished

The input file to QMCPACK (``vmc.in``) is in the form for an ensemble 
or batched run, i.e. it is simply a list of other input files in the 
same directory:

.. code-block:: bash

    >cat runs/diamond/vmc_twist/vmc.in

    vmc.g000.twistnum_0.in.xml
    vmc.g001.twistnum_1.in.xml
    vmc.g002.twistnum_2.in.xml
    vmc.g003.twistnum_3.in.xml
    vmc.g004.twistnum_4.in.xml
    vmc.g005.twistnum_5.in.xml
    vmc.g006.twistnum_6.in.xml
    vmc.g007.twistnum_7.in.xml

These files differ from each other only by the twist index that is 
selected for each MPI group:

.. parsed-literal::

    >diff runs/diamond/vmc_twist/vmc.g000.twistnum_0.in.xml runs/diamond/vmc_twist/vmc.g007.twistnum_7.in.xml 
    47c47
    <  <sposet_builder type="bspline" href="../nscf_twist/pwscf_output/pwscf.pwscf.h5" tilematrix="1 -1 1 1 1 -1 -1 1 1" **twistnum="0"** source="ion0" version="0.10" meshfactor="1.0" precision="float" truncate="no">
    ---
    >  <sposet_builder type="bspline" href="../nscf_twist/pwscf_output/pwscf.pwscf.h5" tilematrix="1 -1 1 1 1 -1 -1 1 1" **twistnum="7"** source="ion0" version="0.10" meshfactor="1.0" precision="float" truncate="no">

The VMC energy for each twist can be estimated with ``qmca``:

.. code-block:: bash

    >qmca -e 40 -q e runs/diamond/vmc_twist/*scalar*
    
    runs/diamond/vmc_twist/vmc.g000  series 0  LocalEnergy  =  -45.075486 +/- 0.008989 
    runs/diamond/vmc_twist/vmc.g001  series 0  LocalEnergy  =  -45.860851 +/- 0.008472 
    runs/diamond/vmc_twist/vmc.g002  series 0  LocalEnergy  =  -45.880189 +/- 0.009115 
    runs/diamond/vmc_twist/vmc.g003  series 0  LocalEnergy  =  -45.880318 +/- 0.009575 
    runs/diamond/vmc_twist/vmc.g004  series 0  LocalEnergy  =  -45.874489 +/- 0.008558 
    runs/diamond/vmc_twist/vmc.g005  series 0  LocalEnergy  =  -45.884590 +/- 0.008889 
    runs/diamond/vmc_twist/vmc.g006  series 0  LocalEnergy  =  -45.870434 +/- 0.009012 
    runs/diamond/vmc_twist/vmc.g007  series 0  LocalEnergy  =  -45.302198 +/- 0.009218 

Finally, we can also obtain the twist averaged energy:

.. code-block:: bash

    >qmca -a -e 40 -q e runs/diamond/vmc_twist/*scalar*
    
    avg  series 0  LocalEnergy  =  -45.703569 +/- 0.003261 

In the next example we will return to the Gamma point case, but instead 
perform timestep extrapolation with DMC.

 
