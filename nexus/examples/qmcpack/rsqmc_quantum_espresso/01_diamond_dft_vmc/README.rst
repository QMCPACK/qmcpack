Nexus QE+QMCPACK Example 2: Diamond cubic cell (DFT+VMC)
========================================================

In this example, we expand on the simple SCF run of the diamond primitive 
cell (2 atoms) to include the steps necessary to perform VMC in a tiled 
cubic supercell (8 atoms).

The calculation steps performed in this example are:

1. SCF run with QE to obtain the charge density.
2. NSCF run with QE to obtain select orbitals to form the 8 atom supercell.
3. Orbital conversion step with pw2qmcpack to convert the NSCF orbitals into the HDF5 format that QMCPACK reads.
4. Jastrow optimization with QMCPACK for the 8 atom supercell via the linear method.
5. VMC run with QMCPACK for the 8 atom supercell.

Apart from these steps, a main difference from the prior example is the 
tiling of the cell.  The tiling matrix from the fcc primitive cell to the 
simple cubic supercell is:

::

   1 -1  1
   1  1 -1
  -1  1  1

The supercell k-point (twist angle) we will use is the Gamma point. Nexus 
maintains a simultaneous representation of the folded primitive cell and its 
tiled supercell counterpart.  Eight k-points in the primitive cell BZ are 
equivalent to the supercell Gamma point.  The primitive cell, and its set of 
8 k-points, is used automatically in the NSCF calculation, while the supercell 
(with its single k-point/twist) is automatically used in the QMC calculations. 
The use of twist averaged boundary conditions will is covered in the next 
example.

The Nexus script used in this example, ``diamond_lda_vmc.py``, is shown below. 
Important differences from the prior example are bolded.

.. parsed-literal:: 

    #! /usr/bin/env python3
    
    from nexus import settings,job,run_project
    from nexus import generate_physical_system
    from nexus import generate_pwscf
    from nexus import **generate_pw2qmcpack**  **\# used for pw2qmcpack** 
    from nexus import **generate_qmcpack**     **\# used for qmcpack**
    
    settings(
        pseudo_dir = '../../pseudopotentials',
        results    = '',
        sleep      = 3,
        machine    = 'ws16',
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
        **tiling   = [[ 1, -1,  1],       \# tiling matrix
                    [ 1,  1, -1],
                    [-1,  1,  1]],
        kgrid    = (1,1,1),             \# supercell k-point grid
        kshift   = (0,0,0),**
        C        = 4,
        )
    
    **\# SCF run with QE to get charge density**
    scf = generate_pwscf(
        identifier   = 'scf',
        path         = 'diamond/scf',
        job          = job(cores=16,app='pw.x'),
        input_type   = 'generic',
        calculation  = 'scf',
        input_dft    = 'lda', 
        ecutwfc      = 200,   
        conv_thr     = 1e-8, 
        system       = system,
        pseudos      = ['C.BFD.upf'],
        kgrid        = (4,4,4),
        kshift       = (0,0,0),
        )
    
    **\# NSCF run with QE to get orbitals**
    nscf = generate_pwscf(
        identifier   = 'nscf',
        **path         = 'diamond/nscf',  \# nscf directory**
        job          = job(cores=16,app='pw.x'),
        input_type   = 'generic',
        **calculation  = 'nscf',          \# nscf calculation**
        input_dft    = 'lda', 
        ecutwfc      = 200,   
        conv_thr     = 1e-8, 
        system       = system,
        pseudos      = ['C.BFD.upf'],
        **nosym        = True,            \# no symmetry for k-points
        dependencies = (scf,'charge_density'), \# depends on scf**
        )
    
    **\# Orbital conversion with pw2qmcpack**
    conv = generate_pw2qmcpack(
        identifier   = 'conv',
        **path         = 'diamond/nscf',  \# run in nscf directory**
        job          = job(cores=16,app='pw2qmcpack.x'),
        write_psir   = False,
        **dependencies = (nscf,'orbitals'), \# depends on nscf**
        )
    
    **\# Jastrow optimization with QMCPACK**
    opt = generate_qmcpack(
        **block        = True,            \# don't run for now**
        identifier   = 'opt',
        path         = 'diamond/optJ2',
        job          = job(cores=16,threads=4,app='qmcpack'),
        input_type   = 'basic',
        system       = system,
        pseudos      = ['C.BFD.xml'],
        **J2           = True,            \# use two-body B-spline Jastrow
        qmc          = 'opt',           \# optimization run, variance opt
        cycles       = 6,               \# loop max=6
        samples      = 51200,           \# VMC samples used in each cycle
        dependencies = (conv,'orbitals'), \# depends on conversion**
        )
    
    **\# VMC with QMCPACK**
    qmc = generate_qmcpack(
        **block        = True,            \# don't run for now**
        identifier   = 'vmc',
        path         = 'diamond/vmc',
        job          = job(cores=16,threads=4,app='qmcpack'),
        input_type   = 'basic',
        system       = system,
        pseudos      = ['C.BFD.xml'],
        J2           = True,
        **qmc          = 'vmc',           \# vmc run, default inputs
        dependencies = [(conv,'orbitals'), \# depends on conversion
                        (opt,'jastrow')],  \# and on optimization**
        )
    
    run_project()


As before, let's check check the status before running.  Now all five 
simulation steps are shown in the workflow:

.. code-block:: bash
    >./diamond_lda_vmc.py --status_only
  
    ...  
  
    cascade status 
      setup, sent_files, submitted, finished, got_output, analyzed, failed 
      000000  0  ------    scf     ./runs/diamond/scf  
      000000  0  ------    nscf    ./runs/diamond/nscf  
      000000  0  ------    conv    ./runs/diamond/nscf  
      000000  0  ------    opt     ./runs/diamond/optJ2  
      000000  0  ------    vmc     ./runs/diamond/vmc  
      setup, sent_files, submitted, finished, got_output, analyzed, failed 


Next, let's run the steps needed to generate the orbitals for the 8 atom 
cubic cell for QMCPACK.  The ``block=True`` statements in the script will 
prevent the QMCPACK optimization and VMC steps from running. 

.. code-block:: bash
    >./diamond_lda_vmc.py 
  
    ...
    
    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 102.71 MB 
      ...
      Entering ./runs/diamond/scf 0 
        Executing:  
          export OMP_NUM_THREADS=1
          mpirun -np 16 pw.x -input scf.in 
    ...
    elapsed time 6.1 s  memory 102.79 MB 
      ...
      Entering ./runs/diamond/nscf 1 
        Executing:  
          export OMP_NUM_THREADS=1
          mpirun -np 16 pw.x -input nscf.in 
    ...
    elapsed time 12.1 s  memory 102.80 MB 
      ...
      Entering ./runs/diamond/nscf 2 
        Executing:  
          export OMP_NUM_THREADS=1
          mpirun -np 16 pw2qmcpack.x<conv.in 
    ...
  
    Project finished

The HDF5 file containing the orbitals should now be available:

.. parsed-literal::

    >ls runs/diamond/nscf/pwscf_output/

    nexus_sync_record  pwscf.wfc11  pwscf.wfc2  pwscf.wfc8
    pwscf.ptcl.xml     pwscf.wfc12  pwscf.wfc3  pwscf.wfc9
    **pwscf.pwscf.h5**     pwscf.wfc13  pwscf.wfc4  pwscf.wfs.xml
    pwscf.save         pwscf.wfc14  pwscf.wfc5  pwscf.xml
    pwscf.wfc1         pwscf.wfc15  pwscf.wfc6
    pwscf.wfc10        pwscf.wfc16  pwscf.wfc7

With the orbitals successfully converted, we can now proceed to optimize 
the Jastrow factor.  Edit ``diamond_lda_vmc.py`` by commenting out the 
``block`` input for the optimization step:

.. parsed-literal::

    opt = generate_qmcpack(
        **\#block        = True,**
        ...
        )

Now rerun the script to perform the optimization step.  This step may 
take several minutes.

.. code-block:: bash

    >./diamond_lda_vmc.py 
    
    ...
      
    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 102.91 MB 
      ...
      Entering ./runs/diamond/optJ2 3 
        Executing:  
          export OMP_NUM_THREADS=4
          mpirun -np 4 qmcpack opt.in.xml 
  
    elapsed time 3.0 s  memory 360.59 MB 
    ...
    elapsed time 601.5 s  memory 103.06 MB 
      Entering ./runs/diamond/optJ2 3 
        copying results  3 opt 
      Entering ./runs/diamond/optJ2 3 
        analyzing  3 opt 
    
    Project finished


The generated QMCPACK input file contains an appropriately tiled orbital set 
as well as a one- and two-body B-spline Jastrow factor.  For the Jastrow, 
the cutoffs are set by default to match the Wigner radius of the supercell 
and one B-spline parameter (knot) has been introduced for every half Bohr 
(see ``./runs/diamond/optJ2/opt.in.xml``):

.. code-block:: xml

    <wavefunction name="psi0" target="e">
       <sposet_builder type="bspline" href="../nscf/pwscf_output/pwscf.pwscf.h5" tilematrix="1 -1 1 1 1 -1 -1 1 1" twistnum="0" source="ion0" version="0.10" meshfactor="1.0" precision="float" truncate="no">
          <sposet type="bspline" name="spo_ud" size="16" spindataset="0"/>
       </sposet_builder>
       <determinantset>
          <slaterdeterminant>
             <determinant id="updet" group="u" sposet="spo_ud" size="16"/>
             <determinant id="downdet" group="d" sposet="spo_ud" size="16"/>
          </slaterdeterminant>
       </determinantset>
       <jastrow type="One-Body" name="J1" function="bspline" source="ion0" print="yes">
          <correlation elementType="C" size="7" rcut="3.37316115" cusp="0.0">
             <coefficients id="eC" type="Array">                  
               0 0 0 0 0 0 0
             </coefficients>
          </correlation>
       </jastrow>
       <jastrow type="Two-Body" name="J2" function="bspline" print="yes">
          <correlation speciesA="u" speciesB="u" size="7" rcut="3.37316115">
             <coefficients id="uu" type="Array">                  
               0 0 0 0 0 0 0
             </coefficients>
          </correlation>
          <correlation speciesA="u" speciesB="d" size="7" rcut="3.37316115">
             <coefficients id="ud" type="Array">                  
               0 0 0 0 0 0 0
             </coefficients>
          </correlation>
       </jastrow>
    </wavefunction>

Default inputs for the quartic variant of the linear optimizer have also 
been populated:

.. code-block:: xml

    <loop max="6">
       <qmc method="linear" move="pbyp" checkpoint="-1">
          <cost name="energy"              >    0.0                </cost>
          <cost name="unreweightedvariance">    1.0                </cost>
          <cost name="reweightedvariance"  >    0.0                </cost>
          <parameter name="warmupSteps"         >    300                </parameter>
          <parameter name="blocks"              >    100                </parameter>
          <parameter name="steps"               >    1                  </parameter>
          <parameter name="subSteps"            >    10                 </parameter>
          <parameter name="timestep"            >    0.3                </parameter>
          <parameter name="useDrift"            >    no                 </parameter>
          <parameter name="samples"             >    51200              </parameter>
          <parameter name="MinMethod"           >    quartic            </parameter>
          <parameter name="minwalkers"          >    0.3                </parameter>
          <parameter name="nonlocalpp"          >    yes                </parameter>
          <parameter name="useBuffer"           >    yes                </parameter>
          <parameter name="alloweddifference"   >    0.0001             </parameter>
          <parameter name="exp0"                >    -6                 </parameter>
          <parameter name="bigchange"           >    10.0               </parameter>
          <parameter name="stepsize"            >    0.15               </parameter>
          <parameter name="nstabilizers"        >    1                  </parameter>
       </qmc>
    </loop>

Check that the optimization has completed successfully by using the ``qmca`` 
tool.  In this case, a variance to energy ratio of 0.025 Ha is acceptable.

.. code-block:: bash

    >qmca -q ev runs/diamond/optJ2/*scalar*
     
                                       LocalEnergy              Variance                ratio 
    runs/diamond/optJ2/opt  series 0  -44.042987 +/- 0.014884   7.043654 +/- 0.077766   0.1599 
    runs/diamond/optJ2/opt  series 1  -45.099480 +/- 0.006295   1.111803 +/- 0.026875   0.0247 
    runs/diamond/optJ2/opt  series 2  -45.088210 +/- 0.005014   1.090274 +/- 0.012517   0.0242 
    runs/diamond/optJ2/opt  series 3  -45.087752 +/- 0.005664   1.095030 +/- 0.021429   0.0243 
    runs/diamond/optJ2/opt  series 4  -45.094245 +/- 0.006613   1.104652 +/- 0.018217   0.0245 
    runs/diamond/optJ2/opt  series 5  -45.105876 +/- 0.009184   1.117030 +/- 0.022910   0.0248 

With optimization completed successfully, let's proceed with VMC. Edit 
``diamond_lda_vmc.py`` by commenting out the ``block`` input for the 
VMC step:

.. parsed-literal::

    qmc = generate_qmcpack(
        **\#block        = True,**
        ...
        )

Now rerun the script to perform the VMC run:  

.. code-block:: bash
  
    >./diamond_lda_vmc.py 
    
    ...
    
    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 102.59 MB 
      ...
      Entering ./runs/diamond/vmc 4 
        Executing:  
          export OMP_NUM_THREADS=4
          mpirun -np 4 qmcpack vmc.in.xml 
  
    elapsed time 3.1 s  memory 363.95 MB 
    ...
    elapsed time 51.6 s  memory 105.23 MB 
      Entering ./runs/diamond/vmc 4 
        copying results  4 vmc 
      Entering ./runs/diamond/vmc 4 
        analyzing  4 vmc 
  
    Project finished

In this case, Nexus will automatically select the best optimized Jastrow 
factor (according to the target cost function) from among the six candidates 
generated during the optimization run (see ``./runs/diamond/vmc/vmc.in.xml``):

.. code-block:: xml
    <jastrow type="One-Body" name="J1" function="bspline" source="ion0" print="yes">
       <correlation elementType="C" size="7" rcut="3.37316115" cusp="0.0">
          <coefficients id="eC" type="Array">                  
            -0.2966053085 -0.238597639 -0.1863207071 -0.1314790098 -0.07964389729 -0.03769253253 -0.01051452959
          </coefficients>
       </correlation>
    </jastrow>
    <jastrow type="Two-Body" name="J2" function="bspline" print="yes">
       <correlation speciesA="u" speciesB="u" size="7" rcut="3.37316115">
          <coefficients id="uu" type="Array">                  
            0.2929092277 0.2146386942 0.1438214568 0.09359039597 0.05663886553 0.02874308877 0.01356022399
          </coefficients>
       </correlation>
       <correlation speciesA="u" speciesB="d" size="7" rcut="3.37316115">
          <coefficients id="ud" type="Array">                  
            0.4658691438 0.3166737337 0.1999575811 0.1193902802 0.06713398775 0.03323672696 0.01554003667
          </coefficients>
       </correlation>
    </jastrow>


A default VMC XML block has also been populated in the input:

.. code-block:: xml
    <qmc method="vmc" move="pbyp" checkpoint="-1">
       <parameter name="walkers"             >    1               </parameter>
       <parameter name="warmupSteps"         >    50              </parameter>
       <parameter name="blocks"              >    800             </parameter>
       <parameter name="steps"               >    10              </parameter>
       <parameter name="subSteps"            >    3               </parameter>
       <parameter name="timestep"            >    0.3             </parameter>
    </qmc>

Finally, let's look at the the VMC result for the total energy.  As expected, the 
result does not differ greatly from the optimization results, but has a smaller 
statistical error bar:

.. code-block:: bash

    >qmca -e 40 -q e runs/diamond/vmc/*scalar*
    
    runs/diamond/vmc/vmc  series 0  LocalEnergy           =  -45.093478 +/- 0.003094 

In the next example, we will build on these results to show how to 
perform twist-averaged VMC calculations with Nexus.

