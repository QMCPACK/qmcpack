.. _qmcpack-input-mod:

Modifying Existing QMCPACK Input Files with Nexus Tools
=======================================================

Below are a few examples on using the ``QmcpackInput.modify`` function to 
modify existing QMCPACK input files, e.g. those made by the converters, 
previously generated with Nexus, or created by hand.


Reading and Writing a QMCPACK Input File with Nexus
---------------------------------------------------

.. code-block:: python

  # import
  from nexus import read_input

  # read 
  #   input can be included in workflows as
  #   generate_qmcpack(input=input,...)  
  input = read_input('qmc.in.xml',format='qmcpack')

  # print
  #   the output is shown in the next section
  print(input.write())

  # modify
  #   see examples in subsequent sections
  input.modify(...)

  # write
  #   a new file is created with the modified input
  input.write('qmc_mod.in.xml')



Example QMCPACK Input File with No Changes
------------------------------------------

Here is an example QMCPACK input file we will modify in the subsequent examples:

.. code-block:: xml

  <?xml version="1.0"?>
  <simulation>
     <project id="MoS2_2H" series="0">
        <parameter name="driver_version"      >    batch           </parameter>
     </project>
     <qmcsystem>
        <simulationcell>
           <parameter name="lattice">
                    5.15349857       -2.97537378        0.00000000
                    0.00000000        5.95074757        0.00000000
                    0.00000000        0.00000000       38.67891432
           </parameter>
           <parameter name="bconds">
              p p p
           </parameter>
           <parameter name="LR_dim_cutoff"       >    15                 </parameter>
        </simulationcell>
        <particleset name="ion0" size="3">
           <attrib name="ionid" datatype="stringArray">
              Mo S S
           </attrib>
           <attrib name="position" datatype="posArray">
                    1.71783286        2.97537378       19.33945716
                    3.43566571        0.00000000       16.28377002
                    3.43566571        0.00000000       22.39514430
           </attrib>
           <group name="MO">
              <parameter name="charge"              >    14                    </parameter>
              <parameter name="valence"             >    14                    </parameter>
              <parameter name="atomicnumber"        >    42                    </parameter>
           </group>
           <group name="S">
              <parameter name="charge"              >    6                     </parameter>
              <parameter name="valence"             >    6                     </parameter>
              <parameter name="atomicnumber"        >    16                    </parameter>
           </group>
        </particleset>
        <particleset name="e" random="yes" randomsrc="ion0">
           <group name="u" size="13">
              <parameter name="charge"              >    -1                    </parameter>
           </group>
           <group name="d" size="13">
              <parameter name="charge"              >    -1                    </parameter>
           </group>
        </particleset>
        <wavefunction name="psi0" target="e">
           <determinantset type="MolecularOrbital" href="MoS2_2H.h5" twist="0 0 0" 
                           source="ion0" transform="yes" name="LCAOBSet" PBCimages="8 8 8">
              <sposet basisset="LCAOBSet" name="spo-up" size="243">
                 <occupation mode="ground">                             </occupation>
                 <coefficient size="243" spindataset="0">                                               </coefficient>
              </sposet>
              <sposet basisset="LCAOBSet" name="spo-dn" size="243">
                 <occupation mode="ground">                              </occupation>
                 <coefficient size="243" spindataset="0">                                               </coefficient>
              </sposet>
              <multideterminant optimize="no" spo_up="spo-up" spo_dn="spo-dn">
                 <detlist size="1657026" type="DETS" nca="0" ncb="0" nea="13" neb="13" 
                          nstates="243" cutoff="1e-20" ext_level="0" href="MoS2_2H.det.h5"/>
              </multideterminant>
           </determinantset>
           <jastrow type="Two-Body" name="J2" function="Bspline" print="yes">
              <correlation speciesA="u" speciesB="u" size="10">
                 <coefficients id="uu" type="Array">                  
  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                 </coefficients>
              </correlation>
              <correlation speciesA="u" speciesB="d" size="10">
                 <coefficients id="ud" type="Array">                  
  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                 </coefficients>
              </correlation>
           </jastrow>
           <jastrow type="One-Body" name="J1" function="Bspline" source="ion0" print="yes">
              <correlation elementType="MO" size="10" cusp="0">
                 <coefficients id="eMo" type="Array">                  
  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                 </coefficients>
              </correlation>
              <correlation elementType="S" size="10" cusp="0">
                 <coefficients id="eS" type="Array">                  
  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                 </coefficients>
              </correlation>
           </jastrow>
           <jastrow type="eeI" name="J3" function="polynomial" print="yes" source="ion0">
              <correlation ispecies="MO" especies="u" isize="3" esize="3">
                 <coefficients id="uuMo" type="Array" optimize="yes">                                               </coefficients>
              </correlation>
              <correlation ispecies="MO" especies1="u" especies2="d" isize="3" esize="3">
                 <coefficients id="udMo" type="Array" optimize="yes">                                               </coefficients>
              </correlation>
              <correlation ispecies="S" especies="u" isize="3" esize="3">
                 <coefficients id="uuS" type="Array" optimize="yes">                                               </coefficients>
              </correlation>
              <correlation ispecies="S" especies1="u" especies2="d" isize="3" esize="3">
                 <coefficients id="udS" type="Array" optimize="yes">                                               </coefficients>
              </correlation>
           </jastrow>
        </wavefunction>
     </qmcsystem>
     <hamiltonian name="h0" type="generic" target="e">
        <pairpot type="coulomb" name="ElecElec" source="e" target="e" physical="yes"/>
        <pairpot type="coulomb" name="IonIon" source="ion0" target="ion0"/>
        <pairpot type="pseudo" name="PseudoPot" source="ion0" wavefunction="psi0" format="xml">
           <pseudo elementType="MO" href="Mo.ccECP.AREP.xml"/>
           <pseudo elementType="S" href="S.ccECP.xml"/>
        </pairpot>
     </hamiltonian>
     <qmc method="vmc" move="pbyp">
        <parameter name="walkers_per_rank"    >    1680            </parameter>
        <parameter name="warmupSteps"         >    10              </parameter>
        <parameter name="blocks"              >    10              </parameter>
        <parameter name="steps"               >    5               </parameter>
        <parameter name="subSteps"            >    2               </parameter>
        <parameter name="timestep"            >    0.1             </parameter>
        <parameter name="useDrift"            >    no              </parameter>
     </qmc>
     <loop max="4">
        <qmc method="linear" move="pbyp" checkpoint="-1">
           <parameter name="walkers_per_rank"    >    1680               </parameter>
           <parameter name="warmupSteps"         >    10                 </parameter>
           <parameter name="blocks"              >    20                 </parameter>
           <parameter name="subSteps"            >    5                  </parameter>
           <parameter name="timestep"            >    0.5                </parameter>
           <parameter name="useDrift"            >    no                 </parameter>
           <parameter name="MinMethod"           >    OneShiftOnly         </parameter>
           <parameter name="minwalkers"          >    0.1                </parameter>
           <estimator name="LocalEnergy" hdf5="no"/>
        </qmc>
     </loop>
     <loop max="10">
        <qmc method="linear" move="pbyp" checkpoint="-1">
           <parameter name="walkers_per_rank"    >    1680               </parameter>
           <parameter name="warmupSteps"         >    10                 </parameter>
           <parameter name="blocks"              >    40                 </parameter>
           <parameter name="subSteps"            >    5                  </parameter>
           <parameter name="timestep"            >    0.5                </parameter>
           <parameter name="useDrift"            >    no                 </parameter>
           <parameter name="MinMethod"           >    OneShiftOnly         </parameter>
           <parameter name="minwalkers"          >    0.5                </parameter>
           <estimator name="LocalEnergy" hdf5="no"/>
        </qmc>
     </loop>
     <qmc method="vmc" move="pbyp" checkpoint="-1">
        <parameter name="walkers_per_rank"    >    1680            </parameter>
        <parameter name="warmupSteps"         >    10              </parameter>
        <parameter name="blocks"              >    20              </parameter>
        <parameter name="subSteps"            >    30              </parameter>
        <parameter name="timestep"            >    0.1             </parameter>
        <parameter name="useDrift"            >    no              </parameter>
        <estimator name="LocalEnergy" hdf5="no"/>
     </qmc>
     <qmc method="dmc" move="pbyp" checkpoint="20">
        <parameter name="walkers_per_rank"    >    1680            </parameter>
        <parameter name="warmupSteps"         >    80              </parameter>
        <parameter name="blocks"              >    1000            </parameter>
        <parameter name="steps"               >    30              </parameter>
        <parameter name="timestep"            >    0.001           </parameter>
        <parameter name="nonlocalmoves"       >    v3              </parameter>
        <estimator name="LocalEnergy" hdf5="no"/>
     </qmc>
  </simulation>





Modification Example
--------------------

This example touches on a subset of the features made available through 
the ``modify`` function.  

As a set of simultaneous changes, we demonstrate 1) changing the orbital 
file path in ``<determinantset/>``, 2) changing the multideterminant filepath 
in ``<multideterminant``, 3) updating the multideterminant ``cutoff``, 
4) generating a new Jastrow factor, 5) disabling optimization of the 
multideterminant while enabling it for the Jastrow, 5) providing new pseudopotential 
files 6) generating new ``<qmc/>`` input sections as well as how to provide your 
own custom ones.

For a full listing of the input parameters and behaviors, see the docstring:

::

  >python3
  >>>
  >>>from nexus import QmcpackInput
  >>>help(QmcpackInput.modify)

Particularly note there the variety of ways to generate optimization 
calculations using QMCPACK's supported optimizers as well as the 
generation options for VMC/DMC calculations.

Here is the code that creates the changes listed above:

.. code-block:: python 

  # import
  from nexus import read_input

  # read 
  input = read_input('qmc.in.xml',format='qmcpack')

  # modify
  input.modify(
      # orbital file
      orbitals_h5     = 'orbs.h5',
      # jastrow generation and optimization
      J2              = True,
      jastrow_opt     = True,
      # multideterminant file, cutoff, optimization
      multidet_h5     = 'multi.h5',
      multidet_cutoff = 1e-6,
      multidet_opt    = False,
      # new pseudopotential files
      pseudo_files    = dict(Mo='Mo.BFD.xml',
                             S ='S.BFD.xml'),
      # generate standard DMC calculation
      #   including equilibration w/ larger timestep
      qmc             = 'dmc',
      eq_dmc          = True,
      )

  # print
  print(input.write())

Here is the input with the modifications:

.. code-block:: xml

  <?xml version="1.0"?>
  <simulation>
     <project id="MoS2_2H" series="0">
        <parameter name="driver_version"      >    batch           </parameter>
     </project>
     <qmcsystem>
        <simulationcell>
           <parameter name="lattice">
                    5.15349857       -2.97537378        0.00000000
                    0.00000000        5.95074757        0.00000000
                    0.00000000        0.00000000       38.67891432
           </parameter>
           <parameter name="bconds">
              p p p
           </parameter>
           <parameter name="LR_dim_cutoff"       >    15                 </parameter>
        </simulationcell>
        <particleset name="ion0" size="3">
           <attrib name="ionid" datatype="stringArray">
              Mo S S
           </attrib>
           <attrib name="position" datatype="posArray">
                    1.71783286        2.97537378       19.33945716
                    3.43566571        0.00000000       16.28377002
                    3.43566571        0.00000000       22.39514430
           </attrib>
           <group name="MO">
              <parameter name="charge"              >    14                    </parameter>
              <parameter name="valence"             >    14                    </parameter>
              <parameter name="atomicnumber"        >    42                    </parameter>
           </group>
           <group name="S">
              <parameter name="charge"              >    6                     </parameter>
              <parameter name="valence"             >    6                     </parameter>
              <parameter name="atomicnumber"        >    16                    </parameter>
           </group>
        </particleset>
        <particleset name="e" random="yes" randomsrc="ion0">
           <group name="u" size="13">
              <parameter name="charge"              >    -1                    </parameter>
           </group>
           <group name="d" size="13">
              <parameter name="charge"              >    -1                    </parameter>
           </group>
        </particleset>
        <wavefunction name="psi0" target="e">
           <determinantset type="MolecularOrbital" href="orbs.h5" twist="0 0 0" source="ion0" transform="yes" name="LCAOBSet" PBCimages="8 8 8">
              <sposet basisset="LCAOBSet" name="spo-up" size="243">
                 <occupation mode="ground">                                               </occupation>
                 <coefficient size="243" spindataset="0">                                               </coefficient>
              </sposet>
              <sposet basisset="LCAOBSet" name="spo-dn" size="243">
                 <occupation mode="ground">                                               </occupation>
                 <coefficient size="243" spindataset="0">                                               </coefficient>
              </sposet>
              <multideterminant optimize="no" spo_up="spo-up" spo_dn="spo-dn">
                 <detlist size="1657026" type="DETS" nca="0" ncb="0" nea="13" neb="13" nstates="243" cutoff="1e-06" ext_level="0" href="multi.h5" optimize="False"/>
              </multideterminant>
           </determinantset>
           <jastrow type="One-Body" name="J1" function="bspline" source="ion0" print="yes">
              <correlation elementType="MO" size="6" rcut="2.97537378312769" cusp="0.0">
                 <coefficients id="eMo" type="Array" optimize="yes">                  
  0 0 0 0 0 0
                 </coefficients>
              </correlation>
              <correlation elementType="S" size="6" rcut="2.97537378312769" cusp="0.0">
                 <coefficients id="eS" type="Array" optimize="yes">                  
  0 0 0 0 0 0
                 </coefficients>
              </correlation>
           </jastrow>
           <jastrow type="Two-Body" name="J2" function="bspline" print="yes">
              <correlation speciesA="u" speciesB="u" size="6" rcut="2.97537378312769">
                 <coefficients id="uu" type="Array" optimize="yes">                  
  0 0 0 0 0 0
                 </coefficients>
              </correlation>
              <correlation speciesA="u" speciesB="d" size="6" rcut="2.97537378312769">
                 <coefficients id="ud" type="Array" optimize="yes">                  
  0 0 0 0 0 0
                 </coefficients>
              </correlation>
           </jastrow>
        </wavefunction>
     </qmcsystem>
     <hamiltonian name="h0" type="generic" target="e">
        <pairpot type="coulomb" name="ElecElec" source="e" target="e" physical="yes"/>
        <pairpot type="coulomb" name="IonIon" source="ion0" target="ion0"/>
        <pairpot type="pseudo" name="PseudoPot" source="ion0" wavefunction="psi0" format="xml">
           <pseudo elementType="MO" href="Mo.BFD.xml"/>
           <pseudo elementType="S" href="S.BFD.xml"/>
        </pairpot>
     </hamiltonian>
     <qmc method="vmc" move="pbyp">
        <parameter name="total_walkers"       >    2048            </parameter>
        <parameter name="warmupSteps"         >    30              </parameter>
        <parameter name="blocks"              >    40              </parameter>
        <parameter name="steps"               >    10              </parameter>
        <parameter name="subSteps"            >    3               </parameter>
        <parameter name="timestep"            >    0.3             </parameter>
        <parameter name="useDrift"            >    no              </parameter>
     </qmc>
     <qmc method="dmc" move="pbyp">
        <parameter name="total_walkers"       >    2048            </parameter>
        <parameter name="warmupSteps"         >    20              </parameter>
        <parameter name="blocks"              >    20              </parameter>
        <parameter name="steps"               >    5               </parameter>
        <parameter name="timestep"            >    0.02            </parameter>
     </qmc>
     <qmc method="dmc" move="pbyp">
        <parameter name="total_walkers"       >    2048            </parameter>
        <parameter name="warmupSteps"         >    20              </parameter>
        <parameter name="blocks"              >    200             </parameter>
        <parameter name="steps"               >    10              </parameter>
        <parameter name="timestep"            >    0.01            </parameter>
     </qmc>
  </simulation>

