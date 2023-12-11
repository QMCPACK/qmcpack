Nexus QE+QMCPACK Example 4: Diamond cubic cell (DMC timestep extrapolation)
===========================================================================

In this example, we expand on Example 2 to include DMC timestep extrapolation.

The Nexus script used for this example is ``diamond_lda_dmc_extrap.py``.
Important differences between this script and the VMC only one from Example 2
 are shown below:

.. parsed-literal::

    qmc = generate_qmcpack(
        identifier   = 'dmc',
        path         = 'diamond/dmc',
        job          = job(cores=16,threads=4,app='qmcpack'),
        input_type   = 'basic',
        system       = system,
        pseudos      = ['C.BFD.xml'],
        J2           = True,
        **qmc          = 'dmc',  \# dmc run inputs
        vmc_samples  = 400,    \# dmc target pop, sampled from vmc
        eq_dmc       = True,   \# add an equilibration run
        timestep     = 0.02,   \# starting timestep
        ntimesteps   = 3,      \# use 3 timesteps, reducing each time**
        dependencies = [(conv,'orbitals'),
                        (opt,'jastrow')],
        )

For demonstration purposes, we use a relatively small DMC walker 
population of 400 (use 1000 or greater for production runs).  
Three successively halved timesteps are used (0.02, 0.01, and 0.005 Ha^-1) 
and the number of DMC steps is doubled relative to the prior one for 
each (Note: If one stricly wished to keep the error bar constant with 
decreasing timestep then 4x more steps would be used instead).

Before running this example, let's copy in the completed runs from 
Example 2:

.. code-block:: bash

    >rsync -av ../02_qe_diamond_dft_vmc/runs ./

Now confirm that the prior runs are complete from Nexus' point of view:

.. code-block:: bash

    >./diamond_lda_dmc_textrap.py --status_only
  
    ...  
  
    cascade status 
      setup, sent_files, submitted, finished, got_output, analyzed, failed 
      111111  0  17009     scf     ./runs/diamond/scf  
      111111  0  17082     nscf    ./runs/diamond/nscf  
      111111  0  17154     conv    ./runs/diamond/nscf  
      111111  0  17362     opt     ./runs/diamond/optJ2  
      000000  0  ------    dmc     ./runs/diamond/dmc  
      setup, sent_files, submitted, finished, got_output, analyzed, failed 

Next we will run the DMC.  This step will take several minutes.

.. code-block:: bash

    >./diamond_lda_dmc_textrap.py 
    
    ...
  
    starting runs:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    elapsed time 0.0 s  memory 102.76 MB 
      Entering ./runs/diamond/dmc 4 
        writing input files  4 dmc 
      Entering ./runs/diamond/dmc 4 
        sending required files  4 dmc 
        submitting job  4 dmc 
      Entering ./runs/diamond/dmc 4 
        Executing:  
          export OMP_NUM_THREADS=4
          mpirun -np 4 qmcpack dmc.in.xml 
  
    elapsed time 3.1 s  memory 366.20 MB 
    ...
    elapsed time 1737.0 s  memory 105.34 MB 
      Entering ./runs/diamond/dmc 4 
        copying results  4 dmc 
      Entering ./runs/diamond/dmc 4 
        analyzing  4 dmc 
  
    Project finished

The generated DMC input file has the DMC timestep XML blocks as described 
above, including warmup VMC and DMC sampling (see ``./runs/diamond/dmc/dmc.in.xml``):

.. code-block:: xml
    <qmc method="vmc" move="pbyp" checkpoint="-1">
       <parameter name="walkers"             >    1               </parameter>
       <parameter name="warmupSteps"         >    30              </parameter>
       <parameter name="blocks"              >    40              </parameter>
       <parameter name="steps"               >    10              </parameter>
       <parameter name="subSteps"            >    3               </parameter>
       <parameter name="timestep"            >    0.3             </parameter>
       <parameter name="samples"             >    400             </parameter>
    </qmc>
    <qmc method="dmc" move="pbyp" checkpoint="-1">
       <parameter name="warmupSteps"         >    20              </parameter>
       <parameter name="blocks"              >    20              </parameter>
       <parameter name="steps"               >    5               </parameter>
       <parameter name="timestep"            >    0.02            </parameter>
    </qmc>
    <qmc method="dmc" move="pbyp" checkpoint="-1">
       <parameter name="warmupSteps"         >    20              </parameter>
       <parameter name="blocks"              >    200             </parameter>
       <parameter name="steps"               >    10              </parameter>
       <parameter name="timestep"            >    0.02            </parameter>
    </qmc>
    <qmc method="dmc" move="pbyp" checkpoint="-1">
       <parameter name="warmupSteps"         >    40              </parameter>
       <parameter name="blocks"              >    200             </parameter>
       <parameter name="steps"               >    20              </parameter>
       <parameter name="timestep"            >    0.01            </parameter>
    </qmc>
    <qmc method="dmc" move="pbyp" checkpoint="-1">
       <parameter name="warmupSteps"         >    80              </parameter>
       <parameter name="blocks"              >    200             </parameter>
       <parameter name="steps"               >    40              </parameter>
       <parameter name="timestep"            >    0.005           </parameter>
    </qmc>

Total energies for each DMC series can be obtained with ``qmca``:

.. code-block:: bash

    >qmca -e 15 -q e runs/diamond/dmc/*scalar*
     
    runs/diamond/dmc/dmc  series 0  LocalEnergy           =  -45.059154 +/- 0.014558 
    runs/diamond/dmc/dmc  series 1  LocalEnergy           =  -45.292932 +/- 0.010916 
    runs/diamond/dmc/dmc  series 2  LocalEnergy           =  -45.272315 +/- 0.005772 
    runs/diamond/dmc/dmc  series 3  LocalEnergy           =  -45.268339 +/- 0.006606 
    runs/diamond/dmc/dmc  series 4  LocalEnergy           =  -45.276465 +/- 0.005578 

In this case it is clear that longer runs are desirable to distinguish 
the different timesteps.  Nevertheless we can use this data to illustrate 
the use of the ``qmc-fit`` timestep extrapolation tool.  The tool uses the 
jack-knife method to obtain an error bar for the DMC energy as it is 
extrapolated to zero timestep.  It can be used in the following way:

.. code-block:: bash

    >qmc-fit ts -e 20 -t '0.02 0.01 0.005' runs/diamond/dmc/*s002*scalar.dat  runs/diamond/dmc/*s003*scalar.dat  runs/diamond/dmc/*s004*scalar.dat 
    
    fit function  : linear
    fitted formula: (-45.2794 +/- 0.0055) + (0.48 +/- 0.42)*t
    intercept     : -45.2794 +/- 0.0055  Ha

Here the "intercept" is the zero timestep value.


