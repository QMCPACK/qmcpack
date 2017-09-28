#! /usr/bin/env python

from nexus import settings
from nexus import Structure,PhysicalSystem
from nexus import generate_pwscf,Job
from nexus import run_project



# set global parameters of nexus
settings(
    pseudo_dir    = '../pseudopotentials',# directory with pseudopotentials
    #generate_only   = False,
    # Complicated setting only so examples can be run in test harness.
    # For real runs, use the plain setting of 'generate_only' above.
    generate_only   = globals().get('override_generate_only_setting',False),

    status_only   = 0,                    # only show run status, T/F
    machine       = 'ws16'                # local machine is 16 core workstation
    )



# describe the physical system
T_structure = Structure()             # empty structure
T_structure.read_xyz('./Ge_T_16.xyz') # read in Ge T interstitial structure

T_structure.reset_axes([              # specify cell axes (in Angstrom)
        [ 5.66,  5.66,  0.  ],
        [ 0.  ,  5.66,  5.66],
        [ 5.66,  0.  ,  5.66]
        ])

T_system = PhysicalSystem(            # make the physical system
    structure = T_structure,          # out of the T interstitial structure
    Ge        = 4                     # pseudo-Ge has 4 valence electrons
    )


# specify MP k-point grids for successive relaxations
supercell_kgrids = [(1,1,1),  #   1 k-point
                    (2,2,2),  #   8 k-points
                    (4,4,4),  #  64 k-points
                    (6,6,6)]  # 216 k-points



# describe the relaxation calculations
# and link them together into a simulation cascade
relaxations = []                        # list of relax simulation objects
for kgrid in supercell_kgrids:          # loop over supercell kgrids
    relax = generate_pwscf(             # make each relax simulation
        identifier = 'relax',               # file prefix
                                            # run directory
        path       = 'relax/kgrid_{0}{1}{2}'.format(*kgrid),
        job        = Job(cores=16),         # will run with mpirun -np 16
        input_type = 'relax',               # this is a relax calculation
        input_dft  = 'pbe',                 # PBE functional
        ecut       = 50,                    # 50 Ry planewave cutoff
        conv_thr   = 1e-6,                  # convergence threshold
        kgrid      = kgrid,                 # supercell k-point grid
        kshift     = (1,1,1),               # grid centered at supercell L point
        pseudos    = ['Ge.pbe-kjpaw.UPF'],  # PBE pseudopotential
        system     = T_system               # the interstitial system
        )
                                        # link together the simulation cascade
                                        #   current relax gets structure from previous
    if len(relaxations)>0:              #   if it exists
        relax.depends(relaxations[-1],'structure')
    #end if
    relaxations.append(relax)           # add relax simulation to the list
#end for



# perform the simulations
run_project(relaxations)



# analyze the results
performed_runs = not settings.generate_only and not settings.status_only
if performed_runs:
    print
    print 'Relaxation results:'
    print '-------------------'
    print '    kgrid     starting force   max force    # of cycles'
    for ik in range(len(supercell_kgrids)):
        kgrid = supercell_kgrids[ik]
        relax = relaxations[ik]
        pa = relax.load_analyzer_image()
        start_force = pa.tot_forces[0]
        max_force   = pa.tot_forces.max()
        ncycles     = len(pa.tot_forces)
        print '  {0:10}  {1:10}     {2:10}  {3:8}'.format(kgrid,start_force,max_force,ncycles)
    #end for
    print
    print
    print 'The final structure is:'
    print 
    print pa.structures.list()[-1].positions
#end if


