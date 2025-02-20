#! /usr/bin/env python3

import os
cores = os.cpu_count() // 2

# Paths to executables.  *** EDIT ME ***
rmg_bin="rmg-cpu"
qmc_bin="qmcpack_complex"

# Import Nexus stuff
from nexus import settings,job,run_project,obj
from nexus import generate_physical_system
from nexus import ppset
from nexus import generate_qmcpack
from nexus import generate_rmg

# machine settings
settings(
    pseudo_dir = '../../qmcpack/pseudopotentials',         # Pseudopotential directory
    generate_only = 0,                                     # only write input files, T/F
    results    = '',                                       # Don't store results separately
    sleep      = 5,                                        # Workflow polling frequency (sec)
    machine    = 'ws'+str(cores),                          # Executing on simple workstation
    account     = 'PSFMat_2',
    )

scf_job  = job(cores=cores, app=rmg_bin)           # How to run RMG
opt_job = job(cores=cores, app=qmc_bin)


ppset(
    label   = 'BFD',
    qmcpack = ['C.BFD.xml'],
    rmg     = ['C.BFD.xml'],

    )

for struct in ['Diamond']:
    primcell = generate_physical_system(
            units    = 'A',
            structure = 'diamond',
            cell      = 'prim',
            symm_kgrid=True,
            C=4,
            )
    #Number of Supercells
    tiled = primcell.structure.tile_opt(2)

   
    scf=generate_rmg(
            path                   = 'RMG/scf',
            system                 = primcell,
            virtual_frac           = 1.0,
            wf_grid_spacing        = 0.15,
            # nexus inputs
            identifier             = 'scf',
            job                    = job(cores=1,app='rmg-cpu'),
            pseudos                = 'BFD',
            input_type             = 'generic',
            # control options
            calculation_mode       = 'Quench Electrons',
            compressed_infile      = False,
            compressed_outfile     = False,
            description            = 'diamond',
            energy_convergence_criterion = 1.0e-09,
            max_scf_steps          = 100,
            #start_mode             = 'Restart From File',
            write_data_period      = 10,
            # cell parameter options
            atomic_coordinate_type = 'Cell Relative',
            kpoint_is_shift        = (0,0,0),
            kpoint_mesh            = (3,2,1),
            potential_grid_refinement = 2,
            # pseudopotential related options
            localize_localpp       = False,
            localize_projectors    = False,
            # kohn sham solver options
            kohn_sham_mucycles     = 3,
            kohn_sham_solver       = 'davidson',
            # orbital occupation options
            occupations_type       = 'Fixed',
            # charge density mixing options
            charge_density_mixing  = 0.5,
            charge_mixing_type     = 'Broyden',
            potential_acceleration_constant_step = 1.0,
            # diagonalization options
            subdiag_driver         = 'lapack',
            exchange_correlation_type = 'pbe',
            # miscellaneous options
            kpoint_distribution    = 1,
            #qmcpack output
            write_qmcpack_restart = True,
            )
    opt1=1
    #Creating a 3x3x3 twist grid with a 2 Supercell 
    for grid in range(3,4):
        supercell = generate_physical_system(
                structure  = tiled.copy(),
                kgrid      = (grid,grid,grid),
                kshift     = (1,1,1),
                symm_kgrid = True,
                C         = 4,
                )
          
        nscf = generate_rmg(
                path                   = 'RMG/nscf-'+str(grid)+'x'+str(grid)+'x'+str(grid),
                system                 = supercell,
                virtual_frac           = 1.0,
                wf_grid_spacing        = 0.15,
                # nexus inputs
                identifier             = 'nscf',
                job                    = job(cores=1,app='rmg-cpu'),
                 pseudos                = 'BFD',
                input_type             = 'generic',
                # control options
                calculation_mode       = 'NSCF',
                compressed_infile      = False,
                compressed_outfile     = False,
                description            = 'diamond',
                # cell parameter options
                atomic_coordinate_type = 'Cell Relative',
                kohn_sham_mucycles     = 3,
                kohn_sham_solver       = 'davidson',
                # diagonalization options
                subdiag_driver         = 'lapack',
                exchange_correlation_type = 'pbe',
                # miscellaneous options
                kpoint_distribution    = 1,
                kpoint_is_shift        = (1,1,1),
                #Kpoint_mesh needs to have a negative number!!!!!
                kpoint_mesh            = (-3,2,1),
                #qmcpack output
                write_qmcpack_restart = True,
                dependencies = (scf,'wavefunctions'),
                )


        if(opt1==1):
            optJ12  = generate_qmcpack(
                    #block           = True,
                    identifier      = 'opt',
                    path            = 'RMG/optJ12',                      # Run directory
                    job             = opt_job,
                    input_type      = 'basic',
                    system          = supercell,                   # System to calculate
                    pseudos         = 'BFD',
                    twistnum        = 0,
                    meshfactor      = 1.0,  # you choose this and hybrid params
                    J1              = True,         # Add a 1-body B-spline Jastrow
                    J2              = True,         # Add a 2-body B-spline Jastrow
                    #J1_rcut         = 6.842691799768850,
                    qmc             = 'opt',        # Do a wavefunction optimization
                    minmethod       = 'oneshift',   # Optimization algorithm (assumes energy minimization)
                    init_cycles     = 2,            # First 4 iterations allow large parameter changes
                    cycles          = 5,           # 8 subsequent iterations with smaller parameter changes
                    warmupsteps     = 8,            # First 8 steps are not recorded
                    blocks          = 100,          # Number of blocks to write in the .scalar.dat file
                    timestep        = 0.1,          # MC step size (nothing to do with time for VMC)
                    init_minwalkers = 0.01,         # Smaller values -> bigger parameter change
                    minwalkers      = 0.5,          # 
                    samples         = 80000,        # VMC samples per iteration
                    use_nonlocalpp_deriv = False,
                    dependencies    = (nscf,'orbitals'),
                    )


            optJ123 = generate_qmcpack(
                    #block           = True,
                    identifier      = 'opt',
                    path            = 'RMG/optJ123',                      # Run directory
                    job             = opt_job,
                    input_type      = 'basic',
                    system          = supercell,                   # System to calculate
                    twistnum        = 0,
                    J3              = True,         # Add a 2-body B-spline Jastrow
                    J3_rcut         = 3.00,         # Cutoff for J3
                    pseudos         = 'BFD',
                    qmc             = 'opt',        # Do a wavefunction optimization
                    minmethod       = 'oneshift',   # Optimization algorithm (assumes energy minimization)
                    init_cycles     = 2,            # First 4 iterations allow large parameter changes
                    init_minwalkers = 0.01,         # Smaller value -> bigger parameter change
                    cycles          = 5,           # Subsequent iterations with smaller parameter changes
                    minwalkers      = 0.5,          # Larger value -> smaller parameter change
                    warmupsteps     = 4,            # First steps are not recorded
                    blocks          = 100,          # Number of blocks to write in the .scalar.dat file
                    timestep        = 0.1,          # MC step size (nothing to do with time for VMC)
                    samples         = 160000,        # VMC samples per iteration
                    use_nonlocalpp_deriv = False,   # Don't include nonlocal pseudo derivatives in optimization (this is nice to have for deep semicore states but expensive!)
                    dependencies    = [(nscf,'orbitals'),
                                     (optJ12, 'jastrow')],
                    )
            opt1=2
       
        # run DMC with 1,2 and 3 Body Jastrow function 
        qmc = generate_qmcpack(
                #block        = True,
                identifier   = 'dmc',
                path         = 'RMG/dmc-'+str(grid)+'x'+str(grid)+'x'+str(grid),                      # Run directory
                job          = job(cores=cores),                  # Submit with the number of cores available
                system       = supercell,                   # System to calculate
                pseudos      = 'BFD',
                jastrows     = [],                                                                                      
                qmc          = 'dmc',                             # dmc run
                vmc_samples  = 16000,                             # Number of Samples (selected from a VMC step)
                vmc_warmupsteps = 100,                            # Number of Equilibration steps
                warmupsteps  = 100,                              # Number of Equilibration steps
                vmc_blocks   = 20,                               # Number of VMC blocks (To generate the DMC samples) 
                vmc_steps    = 20,                                # Number of VMC steps (To generate DMC samples)
                vmc_timestep = 0.1,                               # VMC Timestep (To Generate DMC samples)
                timestep     = 0.01,                            # DMC timestep
                steps           = 40,                             # start with small number for large timesteps [autocorrelation]
                blocks          = 1000,                            # Number of DMC blocks
                nonlocalmoves     = 'v3',
                dependencies    = [(nscf,'orbitals'),
                                 (optJ123, 'jastrow')],
                )
   
   
run_project()
   
   
   
   
   
   
