#! /usr/bin/env python3

from nexus import generate_physical_system, generate_pwscf, job, run_project, settings

# set global parameters of nexus
settings(
    pseudo_dir    = '../pseudopotentials',# directory with pseudopotentials
    generate_only = 0,                    # only write input files, T/F
    status_only   = 0,                    # only show run status, T/F
    machine       = 'ws16'                # local machine is 16 core workstation
    )


# describe the physical system
T_system = generate_physical_system(     # make the physical system
    structure = './Ge_T_16.xyz',         # out of the T interstitial structure
    axes      = [[ 5.66,  5.66,  0.  ],  # specify cell axes (in Angstrom)
                 [ 0.  ,  5.66,  5.66],
                 [ 5.66,  0.  ,  5.66]],
    Ge        = 4,                       # pseudo-Ge has 4 valence electrons
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
        path       = f'relax/kgrid_{kgrid[0]}{kgrid[1]}{kgrid[2]}',
        job        = job(cores=16),         # will run with mpirun -np 16
        input_type = 'relax',               # this is a relax calculation
        input_dft  = 'pbe',                 # PBE functional
        ecut       = 50,                    # 50 Ry planewave cutoff
        conv_thr   = 1e-6,                  # convergence threshold
        kgrid      = kgrid,                 # supercell k-point grid
        kshift     = (1,1,1),               # grid centered at supercell L point
        pseudos    = ['Ge.pbe-kjpaw.UPF'],  # PBE pseudopotential
        system     = T_system,              # the interstitial system
        )
                                        # link together the simulation cascade
                                        #   current relax gets structure from previous
    if len(relaxations)>0:              #   if it exists
        relax.depends(relaxations[-1],'structure')
    #end if
    relaxations.append(relax)           # add relax simulation to the list
#end for


# perform the simulations
run_project()


# analyze the results
performed_runs = not settings.generate_only and not settings.status_only
if performed_runs:
    print()
    print('Relaxation results:')
    print('-------------------')
    print('    kgrid     starting force   max force    # of cycles')
    for ik in range(len(supercell_kgrids)):
        kgrid = supercell_kgrids[ik]
        relax = relaxations[ik]
        pa = relax.load_analyzer_image()
        forces      = pa.forces('Ry/B')
        total_force = (forces**2).sum(axis=(1,2))**0.5
        start_force = total_force[0]
        max_force   = total_force.max()
        ncycles     = len(total_force)
        kgrid_text  = str(kgrid)
        print(f'  {kgrid_text:10}  {start_force:10}     {max_force:10}  {ncycles:8}')
    #end for
    print()
    print()
    print('The final structure is:')
    print()
    print(pa.relaxed_structure('B').pos)
#end if
