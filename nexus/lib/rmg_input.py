
import numpy as np

from generic import obj
from developer import DevBase,error
from unit_converter import convert
from pseudopotential import pp_elem_label
from structure import generate_structure
from simulation import SimulationInput


class RmgInputSettings(DevBase):
    enforce_min_value = True
    enforce_max_value = True
    enforce_allowed   = True
    check_on_write    = True
#end class RmgInputSettings


# raw input spec below
#   taken directly from
#     https://github.com/RMGDFT/rmgdft/wiki/Input-File-Options
#   changes made from website values
#     write_data_period
#       Max value: 50     -> 500
#     pseudopotential
#       Key type : string -> formatted
#     Hubbard_U
#       Key type : string -> formatted

raw_input_spec = '''
Control options

    Key name:     a_length
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      0.000000e+00
    Description:  First lattice constant. 

    Key name:     b_length
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      0.000000e+00
    Description:  Second lattice constant. 

    Key name:     c_length
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      0.000000e+00
    Description:  Third lattice constant. 

    Key name:     calculation_mode
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Quench Electrons"
    Allowed:      "Exx Only" "NEB Relax" "Band Structure Only" "Psi Plot" "Plot" 
                  "Constant Pressure And Energy" "TDDFT" "Dimer Relax" "Constant 
                  Temperature And Energy" "Constant Volume And Energy" "Relax 
                  Structure" "Quench Electrons"   
    Description:  Type of calculation to perform. 

    Key name:     cell_relax
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  flag to control unit cell relaxation 

    Key name:     coalesce_factor
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    1
    Max value:    16
    Default:      4
    Description:  Grid coalescing factor. 

    Key name:     coalesce_states
    Required:     no
    Key type:     boolean
    Expert:       Yes
    Experimental: No
    Default:      "false"
    Description:  Flag indicating whether or not to coalesce states. 

    Key name:     compressed_infile
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  Flag indicating whether or not parallel restart wavefunction file 
                  uses compressed format. 

    Key name:     compressed_outfile
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  Flag indicating whether or not parallel output wavefunction file 
                  uses compressed format. 

    Key name:     description
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      ""
    Allowed:      
    Description:  Description of the run. 

    Key name:     energy_convergence_criterion
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    1.000000e-20
    Max value:    1.000000e-07
    Default:      1.000000e-10
    Description:  The RMS value of the estimated change in the total energy per step 
                  where we assume self consistency has been achieved. 

    Key name:     energy_output_units
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Hartrees"
    Allowed:      "Rydbergs" "Hartrees" 
    Description:  Units to be used when writing energy values to the output file. 
                  Hartrees or Rydbergs are available. 

    Key name:     exx_integrals_filepath
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "afqmc_rmg"
    Allowed:      
    Description:  File/path for exact exchange integrals. 

    Key name:     exx_mode
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Distributed fft"
    Allowed:      "Local fft" "Distributed fft" 
    Description:  FFT mode for exact exchange computations. 

    Key name:     exxdiv_treatment
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "gygi-baldereschi"
    Allowed:      "none" "gygi-baldereschi" 
    Description:  Exact exchange method for handling exx divergence at G=0. 

    Key name:     input_tddft_file
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Waves/wave_tddft.out"
    Allowed:      
    Description:  Input file/path to read wavefunctions and other binary data from 
                  on a restart. 

    Key name:     input_wave_function_file
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Waves/wave.out"
    Allowed:      
    Description:  Input file/path to read wavefunctions and other binary data from 
                  on a restart. 

    Key name:     interpolation_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "FFT"
    Allowed:      "FFT" "prolong" "Cubic Polynomial" 
    Description:  Interpolation method for transferring data between the potential 
                  grid and the wavefunction grid. Mostly for diagnostic purposes. 

    Key name:     max_exx_steps
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    2147483647
    Default:      100
    Description:  Maximum number of self consistent steps to perform with hybrid 
                  functionals. 

    Key name:     max_scf_steps
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    2147483647
    Default:      500
    Description:  Maximum number of self consistent steps to perform. Inner loop for 
                  hybrid functionals. 

    Key name:     noncollinear
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  if set true, calculate noncollinear 

    Key name:     nvme_orbitals
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Flag indicating whether or not orbitals should be mapped to disk. 

    Key name:     nvme_orbitals_filepath
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Orbitals/"
    Allowed:      
    Description:  File/path for runtime disk storage of orbitals. 

    Key name:     nvme_weights
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Flag indicating whether or not projector weights should be mapped 
                  to disk. 

    Key name:     nvme_weights_filepath
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Weights/"
    Allowed:      
    Description:  File/path for disk storage of projector weights. 

    Key name:     nvme_work
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Flag indicating whether or not work arrays should be mapped to 
                  disk. 

    Key name:     nvme_work_filepath
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Work/"
    Allowed:      
    Description:  File/path for disk storage of workspace. 

    Key name:     omp_threads_per_node
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    64
    Default:      0
    Description:  Number of Open MP threads each MPI process will use. A value of 0 
                  selects automatic setting. 

    Key name:     output_tddft_file
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Waves/wave_tddft.out"
    Allowed:      
    Description:  Output file/path to store wavefunctions and other binary data. 

    Key name:     output_wave_function_file
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Waves/wave.out"
    Allowed:      
    Description:  Output file/path to store wavefunctions and other binary data. 

    Key name:     pseudo_dir
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "."
    Allowed:      
    Description:  Directory where pseudopotentials are stored. 

    Key name:     qfunction_filepath
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Qfunctions/"
    Allowed:      
    Description:  File/path for runtime disk storage of qfunctions. 

    Key name:     read_serial_restart
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Directs RMG to read from serial restart files. Normally used when 
                  changing the sprocessor topology used during a restart run 

    Key name:     rms_convergence_criterion
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    1.000000e-03
    Default:      1.000000e-07
    Description:  The RMS value of the change in the total potential from step to 
                  step where we assume self consistency has been achieved. 

    Key name:     spinorbit
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  if set true, calculate with spinorbit coupling 

    Key name:     start_mode
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "LCAO Start"
    Allowed:      "Modified LCAO Start" "Restart TDDFT" "Start TDDFT" "Gaussian 
                  Start" "FIREBALL Start" "LCAO Start" "Restart From File" "Random 
                  Start" 
    Description:  Type of run. 

    Key name:     stress
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  flag to control stress cacluation 

    Key name:     stress_convergence_criterion
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    50.000000
    Default:      0.500000
    Description:  The stress criteria 

    Key name:     tddft_steps
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    2147483647
    Default:      2000
    Description:  Maximum number of tddft steps to perform. 

    Key name:     time_reversal
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  if false, no k -> -k symmetry 

    Key name:     vdw_corr
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "None"
    Allowed:      "DFT-D3" "DFT-D2" "Grimme-D2" "None" 
    Description:  Type of vdw correction 

    Key name:     vdwdf_kernel_filepath
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "vdW_kernel_table"
    Allowed:      
    Description:  File/path for vdW_kernel_table data. 

    Key name:     wannier90
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  set up informations for wannier90 interface 

    Key name:     wannier90_scdm_mu
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    -unlimited
    Max value:    unlimited
    Default:      0.000000e+00
    Description:  when wannier90 is used to build wannier functions, the energy 
                  window parameter 

    Key name:     write_data_period
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    5
    Max value:    500
    Default:      5
    Description:  How often to write checkpoint files during the initial quench in 
                  units of SCF steps. During structural relaxations of molecular 
                  dynamics checkpoints are written each ionic step. 

    Key name:     write_eigvals_period
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    100
    Default:      5
    Description:  How often to output eigenvalues in units of scf steps. 

    Key name:     write_pseudopotential_plots
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Flag to indicate whether or not to write pseudopotential plots. 

    Key name:     write_qmcpack_restart
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  If true then a QMCPACK restart file is written as well as a serial 
                  restart file. 

    Key name:     write_qmcpack_restart_localized
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  If true then a QMCPACK restart file for localized orbitals 

    Key name:     write_serial_restart
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  RMG normally writes parallel restart files. These require that 
                  restarts have the same processor topology. If write_serial_restart 
                  = "true" then RMG will also write a serial restart file that can 
                  be used with a different processor topology 

Cell parameter options

    Key name:     atomic_coordinate_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Absolute"
    Allowed:      "Absolute" "Cell Relative" 
    Description:  Flag indicated whether or not atomic coordinates are absolute or 
                  cell relative. 

    Key name:     bravais_lattice_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Orthorhombic Primitive"
    Allowed:      "Tetragonal Primitive" "Cubic Body Centered" "Orthorhombic 
                  Primitive" "Cubic Face Centered" "Hexagonal Primitive" "Cubic 
                  Primitive" "None" 
    Description:  Bravais Lattice Type. 

    Key name:     cell_movable
    Required:     no
    Key type:     integer array
    Expert:       No
    Experimental: No
    Default:      "0 0 0 0 0 0 0 0 0 "
    Description:  9 numbers to control cell relaxation 

    Key name:     crds_units
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Bohr"
    Allowed:      "Angstrom" "Bohr" 
    Description:  Units for the atomic coordinates. 

    Key name:     grid_spacing
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      0.350000
    Description:  Approximate grid spacing (bohr). 

    Key name:     kpoint_is_shift
    Required:     no
    Key type:     integer array
    Expert:       No
    Experimental: No
    Default:      "0 0 0 "
    Description:  Three-D layout of the kpoint shift. 

    Key name:     kpoint_mesh
    Required:     no
    Key type:     integer array
    Expert:       No
    Experimental: No
    Default:      "1 1 1 "
    Description:  Three-D layout of the kpoint mesh. 

    Key name:     lattice_units
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Bohr"
    Allowed:      "Angstrom" "Alat" "Bohr" 
    Description:  Units for the lattice vectors 

    Key name:     lattice_vector
    Required:     no
    Key type:     double array
    Expert:       No
    Experimental: No
    Default:      "Not done yet"
    Description:  Lattice vectors, a0, a1, a2 

    Key name:     potential_grid_refinement
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    4
    Default:      0
    Description:  Ratio of the potential grid density to the wavefunction grid 
                  density. For example if the wavefunction grid is (72,72,72) and 
                  potential_grid_refinement = "2" then the potential grid would be 
                  (144,144,144). The default value is 2 but it may sometimes be 
                  beneficial to adjust this. (For USPP the minimum value is also 2 
                  and it cannot be set lower. NCPP can be set to 1). 

    Key name:     processor_grid
    Required:     no
    Key type:     integer array
    Expert:       No
    Experimental: No
    Default:      "1 1 1 "
    Description:  Three-D (x,y,z) layout of the MPI processes. 

    Key name:     wavefunction_grid
    Required:     no
    Key type:     integer array
    Expert:       No
    Experimental: No
    Default:      "1 1 1 "
    Description:  Three-D (x,y,z) dimensions of the grid the wavefunctions are 
                  defined on. 

Pseudopotential related options

    Key name:     atomic_orbital_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "delocalized"
    Allowed:      "delocalized" "localized" 
    Description:  Atomic Orbital Type. Choices are localized and delocalized. 

    Key name:     energy_cutoff_parameter
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.600000
    Max value:    1.000000
    Default:      0.800000
    Description:  

    Key name:     filter_dpot
    Required:     no
    Key type:     boolean
    Expert:       Yes
    Experimental: No
    Default:      "false"
    Description:  Flag indicating whether or not to filter density dependent 
                  potentials. 

    Key name:     filter_factor
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.060000
    Max value:    1.000000
    Default:      0.250000
    Description:  Filtering factor. 

    Key name:     localize_localpp
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  The local potential associated with a particular ion also decays 
                  rapidly in real-space with increasing r. As with beta projectors 
                  truncating the real-space representation for large cells can lead 
                  to significant computational savings with a small loss of accuracy 
                  but it should be set to false for small cells. 

    Key name:     localize_projectors
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  The Beta function projectors for a particular ion decay rapidly in 
                  real-space with increasing r. For large cells truncating the 
                  real-space representation of the projector can lead to significant 
                  computational savings with a small loss of accuracy. For smaller 
                  cells the computational cost is the same for localized or 
                  delocalized projectors so it is better to set localize_projectors 
                  to false. 

    Key name:     max_nlradius
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    2.000000
    Max value:    10000.000000
    Default:      10000.000000
    Description:  maximum radius for non-local projectors 

    Key name:     max_qradius
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    2.000000
    Max value:    10000.000000
    Default:      10000.000000
    Description:  maximum radius for qfunc in ultra-pseudopotential 

    Key name:     min_nlradius
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    1.000000
    Max value:    10000.000000
    Default:      2.000000
    Description:  minimum radius for non-local projectors 

    Key name:     min_qradius
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    1.000000
    Max value:    10000.000000
    Default:      2.000000
    Description:  minimum radius for qfunc in ultra-pseudopotential 

    Key name:     projector_expansion_factor
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.500000
    Max value:    3.000000
    Default:      1.000000
    Description:  When using localized projectors the radius can be adjusted with 
                  this parameter. 

    Key name:     pseudopotential
    Required:     no
    Key type:     formatted
    Expert:       No
    Experimental: No
    Default:      ""
    Allowed:      
    Description:  External pseudopotentials may be specfied with this input key. The 
                  format uses the atomic symbol followed by the pseudopotential file 
                  name. pseudopotential = "Ni Ni.UPF O O.UPF" 

Kohn Sham solver options

    Key name:     davidson_max_steps
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    5
    Max value:    20
    Default:      8
    Description:  Maximum number of iterations for davidson diagonalization. 

    Key name:     davidson_multiplier
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    6
    Default:      0
    Description:  The davidson solver expands the eigenspace with the maximum 
                  expansion factor being set by the value of davidson_multiplier. 
                  Larger values often lead to faster convergence but because the 
                  computational cost of the davidson diagonalization step scales as 
                  the cube of the number of eigenvectors the optimal value based on 
                  the fastest time to solution depends on the number of orbitals. If 
                  not specified explicitly or set to 0 RMG uses the following 
                  algorithm to set the value. 
                  
                  Number of orbitals <= 600 davidson_multiplier= "4" 
                  600 < Number of orbitals <= 900 davidson_multiplier = "3" 
                  Number of orbitals > 900 davidson_multiplier = "2" 
                  
                  For very large problems the N^3 scaling makes even a factor of 2 
                  prohibitively costly and the multigrid solver is a better choice. 

    Key name:     kohn_sham_coarse_time_step
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.000000e+00
    Max value:    1.200000
    Default:      1.000000
    Description:  Time step to use in the kohn-sham multigrid solver on the coarse 
                  levels. 

    Key name:     kohn_sham_fd_order
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    6
    Max value:    10
    Default:      8
    Description:  RMG uses finite differencing to represent the kinetic energy 
                  operator and the accuracy of the representation is controllable by 
                  the kohn_sham_fd_order parameter. The default is 8 and is fine for 
                  most purposes but higher accuracy is obtainable with 10th order at 
                  the cost of some additional computational expense. 

    Key name:     kohn_sham_mg_levels
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    -1
    Max value:    6
    Default:      -1
    Description:  Number of multigrid levels to use in the kohn-sham multigrid 
                  preconditioner. 

    Key name:     kohn_sham_mg_timestep
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.000000e+00
    Max value:    2.000000
    Default:      0.666667
    Description:  timestep for multigrid correction. 

    Key name:     kohn_sham_mucycles
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    6
    Default:      2
    Description:  Number of mu (also known as W) cycles to use in the kohn-sham 
                  multigrid preconditioner. 

    Key name:     kohn_sham_post_smoothing
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    1
    Max value:    5
    Default:      2
    Description:  Number of global grid post-smoothing steps to perform after a 
                  multigrid preconditioner iteration. 

    Key name:     kohn_sham_pre_smoothing
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    1
    Max value:    5
    Default:      2
    Description:  Number of global grid pre-smoothing steps to perform before a 
                  multigrid preconditioner iteration. 

    Key name:     kohn_sham_solver
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "davidson"
    Allowed:      "davidson" "multigrid" 
    Description:  RMG supports a pure multigrid Kohn-Sham solver as well as a 
                  multigrid preconditioned davidson solver. The davidson solver is 
                  usually better for smaller problems with the pure multigrid solver 
                  often being a better choice for very large problems. 

    Key name:     kohn_sham_time_step
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.000000e+00
    Max value:    2.000000
    Default:      0.660000
    Description:  Smoothing timestep to use on the fine grid in the the kohn-sham 
                  multigrid preconditioner. 

    Key name:     unoccupied_tol_factor
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    1.000000
    Max value:    100000.000000
    Default:      1000.000000
    Description:  When using the Davidson Kohn-Sham solver unoccupied states are 
                  converged to a less stringent tolerance than occupied orbitals 
                  with the ratio set by this parameter. 

Exchange correlation options

    Key name:     exchange_correlation_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "AUTO_XC"
    Allowed:      "hartree-fock" "vdw-df-c09" "sla+pw+pbe+vdw1" "VDW-DF" "vdw-df" 
                  "gaupbe" "B3LYP" "hse" "mgga tb09" "AUTO_XC" "m06l" "VDW-DF-CX" 
                  "tpss" "ev93" "optbk88" "sogga" "wc" "HSE" "HCTH" "hcth" "Q2D" 
                  "q2d" "PBESOL" "tb09" "b86bpbe" "PW86PBE" "PBE0" "MGGA TB09" 
                  "pw86pbe" "REVPBE" "pbe" "revpbe" "GGA PBE" "BLYP" "pbe0" "pbesol" 
                  "blyp" "PBE" "GGA XP CP" "pw91" "GGA XB CP" "TB09" "optb86b" 
                  "olyp" "BP" "GGA BLYP" "bp" "b3lyp" "LDA" "vdw-df-cx" "PW91" "PZ" 
                  "pz" 
    Description:  Most pseudopotentials specify the exchange correlation type they 
                  were generated with and the default value of AUTO_XC means that 
                  the type specified in the pseudopotial is what RMG will use. That 
                  can be overridden by specifying a value here. 

    Key name:     exx_convergence_criterion
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    1.000000e-12
    Max value:    1.000000e-06
    Default:      1.000000e-09
    Description:  Convergence criterion for the EXX delta from step to step where we 
                  assume EXX consistency has been achieved. 

    Key name:     exx_fraction
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    -1.000000
    Max value:    1.000000
    Default:      -1.000000
    Description:  when hybrid functional is used, the fraction of Exx 

    Key name:     vexx_fft_threshold
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    1.000000e-14
    Max value:    0.100000
    Default:      1.000000e-14
    Description:  The value for the EXX delta where we switch from single to double 
                  precision ffts. Single precision ffts are generally accurate 
                  enough. 

    Key name:     x_gamma_extrapolation
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  if set true, use exx extrapolation to gamma 

Orbital occupation options

    Key name:     MP_order
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    5
    Default:      2
    Description:  order of Methefessel Paxton occupation. 

    Key name:     dos_broading
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    1.000000
    Default:      0.100000
    Description:  For DOS with Gaussian broading method 

    Key name:     dos_method
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "tetrahedra"
    Allowed:      "Gaussian" "tetrahedra" 
    Description:  tetrahedra or gauss smearing method for DOS calculation 

    Key name:     occupation_electron_temperature_eV
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    2.000000
    Default:      0.040000
    Description:  Target electron temperature when not using fixed occupations. 

    Key name:     occupation_number_mixing
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    1.000000
    Default:      1.000000
    Description:  Mixing parameter for orbital occupations when not using fixed 
                  occupations. 

    Key name:     occupations_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Fermi Dirac"
    Allowed:      "Error Function" "Gaussian" "Fermi Dirac" "MethfesselPaxton" "Cold 
                  Smearing" "Fixed" 
    Description:  RMG supports several different ways of specifying orbital 
                  occupations. For a spin polarized system one may specify the 
                  occupations for up and down separately. In the case of a non-zero 
                  electronic temperature these will be adjusted as the calculation 
                  proceeds based on this setting. 

    Key name:     states_count_and_occupation
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      ""
    Allowed:      
    Description:  Occupation string for states. Format for a system with 240 
                  electrons and 20 unoccupied states would be. "120 2.0 20 0.0" 

    Key name:     states_count_and_occupation_spin_down
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      ""
    Allowed:      
    Description:  Occupation string for spin down states. Format is the same as for 
                  states_count_and_occupation. Total number of states must match 
                  spin up occupation string. 

    Key name:     states_count_and_occupation_spin_up
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      ""
    Allowed:      
    Description:  Occupation string for spin up states. Format is the same as for 
                  states_count_and_occupation. Total number of states must match 
                  spin down occupation string. 

    Key name:     unoccupied_states_per_kpoint
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    2147483647
    Default:      10
    Description:  The number of unoccupied orbitals. A value that is 15-20% of the 
                  number of occupied orbitals generally works well. 

Charge density mixing options

    Key name:     charge_broyden_order
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    10
    Default:      5
    Description:  Number of previous steps to use when Broyden mixing is used to 
                  update the charge density. 

    Key name:     charge_broyden_scale
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    1.000000
    Default:      0.500000
    Description:  

    Key name:     charge_density_mixing
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    1.000000
    Default:      0.500000
    Description:  Proportion of the current charge density to replace with the new 
                  density after each scf step when linear mixing is used. 

    Key name:     charge_mixing_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Pulay"
    Allowed:      "Broyden" "Pulay" "Linear" 
    Description:  RMG supports Broyden, Pulay and Linear mixing When the davidson 
                  Kohn-Sham solver is selected Broyden or Pulay are preferred. For 
                  the multigrid solver Linear with potential acceleration is often 
                  (but not always) the best choice. 

    Key name:     charge_pulay_Gspace
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  if set true, charge density mixing the residual in G space 

    Key name:     charge_pulay_order
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    10
    Default:      5
    Description:  Number of previous steps to use when Pulay mixing is used to 
                  update the charge density. 

    Key name:     charge_pulay_refresh
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    2147483647
    Default:      100
    Description:  

    Key name:     charge_pulay_scale
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    1.000000
    Default:      0.500000
    Description:  

    Key name:     potential_acceleration_constant_step
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    4.000000
    Default:      0.000000e+00
    Description:  When set to a non-zero value this parameter causes RMG to perform 
                  a band by band update of the self-consistent potential during the 
                  course of an SCF step when the multigrid kohn_sham_solver is 
                  chosen. This means that updates to the lower energy orbitals are 
                  incorporated into the SCF potential seen by the higher energy 
                  orbitals as soon as they are computed. This can lead to faster 
                  convergence and better stability for many systems. The option 
                  should only be used with Linear mixing. Even when the davidson 
                  solver is chosen this parameter may be used since the first few 
                  steps with davidson usually uses the multigrid solver. 

Relaxation and Molecular dynamics options

    Key name:     dynamic_time_counter
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    0
    Default:      0
    Description:  

    Key name:     dynamic_time_delay
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    5
    Max value:    5
    Default:      5
    Description:  

    Key name:     force_grad_order
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    12
    Default:      8
    Description:  Atomic forces may be computed to varying degrees of accuracy 
                  depending on the requirements of a specific problem. A value of 0 
                  implies highest accuracy which is obtained by using FFTs in place 
                  of finite differencing. 

    Key name:     ionic_time_step
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      50.000000
    Description:  Ionic time step for use in molecular dynamics and structure 
                  optimizations. 

    Key name:     ionic_time_step_decrease
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    1.000000
    Default:      0.500000
    Description:  Factor by which ionic timestep is decreased when dynamic timesteps 
                  are enabled. 

    Key name:     ionic_time_step_increase
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    1.000000
    Max value:    3.000000
    Default:      1.100000
    Description:  Factor by which ionic timestep is increased when dynamic timesteps 
                  are enabled. 

    Key name:     max_ionic_time_step
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    150.000000
    Default:      150.000000
    Description:  Maximum ionic time step to use for molecular dynamics or 
                  structural optimizations. 

    Key name:     max_md_steps
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    2147483647
    Default:      100
    Description:  Maximum number of molecular dynamics steps to perform. 

    Key name:     md_integration_order
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "5th Beeman-Velocity Verlet"
    Allowed:      "5th Beeman-Velocity Verlet" "3rd Beeman-Velocity Verlet" "2nd 
                  Velocity Verlet" 
    Description:  Integration order for molecular dynamics. 

    Key name:     md_nose_oscillation_frequency_THz
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      15.590000
    Description:  

    Key name:     md_number_of_nose_thermostats
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    5
    Max value:    5
    Default:      5
    Description:  Number of Nose thermostats to use during Constant Volume and 
                  Temperature MD. 

    Key name:     md_randomize_velocity
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  The initial ionic velocities for a molecular dyanamics run are 
                  randomly initialized to the target temperature. 

    Key name:     md_temperature
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      300.000000
    Description:  Target MD Temperature. 

    Key name:     md_temperature_control
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Nose Hoover Chains"
    Allowed:      "Anderson Rescaling" "Nose Hoover Chains" 
    Description:  Type of temperature control method to use in molecular dynamics. 

    Key name:     relax_dynamic_timestep
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Flag indicating whether or not to use dynamic timesteps in 
                  relaxation mode. 

    Key name:     relax_mass
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Atomic"
    Allowed:      "Equal" "Atomic" 
    Description:  Mass to use for structural relaxation, either atomic masses, or 
                  the mass of carbon for all atoms. 

    Key name:     relax_max_force
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      2.500000e-03
    Description:  Force value at which an ionic relaxation is considered to be 
                  converged. 

    Key name:     relax_method
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Fast Relax"
    Allowed:      "LBFGS" "MD Min" "Quick Min" "FIRE" "Fast Relax" 
    Description:  Type of relaxation method to use for structural optimizations. 

    Key name:     renormalize_forces
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  Flag indicating whether or not to renormalize forces. 

    Key name:     tddft_time_step
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      0.200000
    Description:  TDDFT time step for use in TDDFT mode 

Diagonalization options

    Key name:     extra_random_lcao_states
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    2147483647
    Default:      0
    Description:  LCAO (Linear Combination of Atomic Orbitals) is the default 
                  startup method for RMG. The atomic orbitals are obtained from the 
                  pseudpotentials but in some cases better convergence may be 
                  obtained by adding extra random wavefunctions in addition to the 
                  atomic orbitals. 

    Key name:     folded_spectrum
    Required:     no
    Key type:     boolean
    Expert:       Yes
    Experimental: No
    Default:      "false"
    Description:  When the number of eigenvectors is large using folded_spectrum is 
                  substantially faster than standard diagonalization. It also tends 
                  to converge better for metallic systems. It works with the 
                  multigrid kohn_sham_solver but not the davidson solver. 

    Key name:     folded_spectrum_iterations
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    0
    Max value:    20
    Default:      2
    Description:  Number of folded spectrum iterations to perform. 

    Key name:     folded_spectrum_width
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.100000
    Max value:    1.000000
    Default:      0.300000
    Description:  Submatrix width to use as a fraction of the full spectrum. The 
                  folded spectrum width ranges from 0.10 to 1.0. For insulators and 
                  semiconductors a value of 0.3 is appropriate. For metals values 
                  between 0.15 to 0.2 tend to be better. The default value is 0.3 

    Key name:     initial_diagonalization
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  Perform initial subspace diagonalization. 

    Key name:     period_of_diagonalization
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    2147483647
    Default:      1
    Description:  Diagonalization period (per scf step). Mainly for debugging and 
                  should not be changed for production. 

    Key name:     scalapack_block_factor
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    4
    Max value:    512
    Default:      32
    Description:  Block size to use with scalapack. Optimal value is dependent on 
                  matrix size and system hardware. 

    Key name:     subdiag_driver
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "auto"
    Allowed:      "elpa" "cusolver" "auto" "scalapack" "magma" "lapack" 
    Description:  Driver type used for subspace diagonalization of the eigenvectors. 

Performance related options

    Key name:     mpi_queue_mode
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  Use mpi queue mode. 

    Key name:     non_local_block_size
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    64
    Max value:    40000
    Default:      512
    Description:  Block size to use when applying the non-local and S operators. 

    Key name:     preconditioner_threshold
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    1.000000e-09
    Max value:    0.100000
    Default:      0.100000
    Description:  The RMS value of the change in the total potential where we switch 
                  the preconditioner from single to double precision. 

    Key name:     require_huge_pages
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: Yes
    Default:      "false"
    Description:  If set RMG assumes that sufficient huge pages are available. Bad 
                  things may happen if this is not true. 

    Key name:     spin_manager_thread
    Required:     no
    Key type:     boolean
    Expert:       Yes
    Experimental: No
    Default:      "true"
    Description:  When mpi_queue_mode is enabled the manager thread spins instead of 
                  sleeping. 

    Key name:     spin_worker_threads
    Required:     no
    Key type:     boolean
    Expert:       Yes
    Experimental: No
    Default:      "true"
    Description:  When mpi_queue_mode is enabled the worker threads spin instead of 
                  sleeping. 

    Key name:     state_block_size
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    1
    Max value:    2147483647
    Default:      64
    Description:  state_block used in nlforce. 

    Key name:     use_alt_zgemm
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Flag indicating whether or not to use alternate zgemm 
                  implementation. 

    Key name:     use_async_allreduce
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  Use asynchronous allreduce if available. 

    Key name:     use_hwloc
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Use internal hwloc setup if available. If both this and use_numa 
                  are true hwloc takes precedence. 

    Key name:     use_numa
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  Use internal numa setup if available. 

LDAU options

    Key name:     Hubbard_U
    Required:     no
    Key type:     formatted
    Expert:       No
    Experimental: No
    Default:      ""
    Allowed:      
    Description:  Hubbard U parameter for each atomic species using the format 
                  Hubbard_U="Ni 6.5" 

    Key name:     ldaU_mode
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "None"
    Allowed:      "Simple" "None" 
    Description:  Type of lda+u implementation. 

    Key name:     ldaU_radius
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    1.000000
    Max value:    12.000000
    Default:      9.000000
    Description:  Max radius of atomic orbitals to be used in LDA+U projectors. 

Poisson solver options

    Key name:     hartree_max_sweeps
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    5
    Max value:    100
    Default:      10
    Description:  Maximum number of hartree iterations to perform per scf step. 

    Key name:     hartree_min_sweeps
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    5
    Default:      5
    Description:  Minimum number of hartree iterations to perform per scf step. 

    Key name:     hartree_rms_ratio
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    1000.000000
    Max value:    unlimited
    Default:      100000.000000
    Description:  Ratio between target RMS for get_vh and RMS total potential. 

    Key name:     poisson_coarse_time_step
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.400000
    Max value:    1.000000
    Default:      0.800000
    Description:  Time step to use in the poisson multigrid solver on the coarse 
                  levels. 

    Key name:     poisson_coarsest_steps
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    10
    Max value:    100
    Default:      25
    Description:  Number of smoothing steps to use on the coarsest level in the 
                  hartree multigrid solver. 

    Key name:     poisson_finest_time_step
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.400000
    Max value:    1.000000
    Default:      1.000000
    Description:  Time step to use in the poisson multigrid solver on the finest 
                  level. 

    Key name:     poisson_mucycles
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    1
    Max value:    4
    Default:      3
    Description:  Number of mu (also known as W) cycles to use in the hartree 
                  multigrid solver. 

    Key name:     poisson_post_smoothing
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    1
    Max value:    6
    Default:      1
    Description:  Number of global hartree grid post-smoothing steps to perform 
                  after a multigrid iteration. 

    Key name:     poisson_pre_smoothing
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    1
    Max value:    6
    Default:      2
    Description:  Number of global hartree grid pre-smoothing steps to perform 
                  before a multigrid iteration. 

    Key name:     poisson_solver
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "pfft"
    Allowed:      "pfft" "multigrid" 
    Description:  poisson solver. 

Testing options

    Key name:     test_energy
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    -1000000000.000000
    Max value:    1000000000.000000
    Default:      nan
    Description:  Expected final energy for testing. 

    Key name:     test_energy_tolerance
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    1.000000e-08
    Max value:    1.000000e-04
    Default:      1.000000e-07
    Description:  Test final energy tolerance. 

Miscellaneous options

    Key name:     E_POINTS
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    201
    Max value:    201
    Default:      201
    Description:  

    Key name:     Emax
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    -100.000000
    Max value:    100.000000
    Default:      0.000000e+00
    Description:  

    Key name:     Emin
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    -100.000000
    Max value:    100.000000
    Default:      -6.000000
    Description:  

    Key name:     ExxCholMax
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    64
    Default:      8
    Description:  maximum number of Exx integral cholesky vectors 

    Key name:     ExxIntCholosky
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  if set true, Exx integrals are Cholesky factorized to 3-index 

    Key name:     alt_laplacian
    Required:     no
    Key type:     boolean
    Expert:       Yes
    Experimental: No
    Default:      "true"
    Description:  Flag indicating whether or not to use alternate laplacian weights 
                  for some operators. 

    Key name:     boundary_condition_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Periodic"
    Allowed:      "Periodic" 
    Description:  Boundary condition type Only periodic is currently implemented. 

    Key name:     charge_analysis
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Voronoi"
    Allowed:      "Voronoi" "None" 
    Description:  Type of charge analysis to use. Only Voronoi deformation density 
                  is currently available. 

    Key name:     charge_analysis_period
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    500
    Default:      0
    Description:  How often to perform and write out charge analysis. 

    Key name:     cube_pot
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  if set true, total potential is printed out in cube format 

    Key name:     cube_rho
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  if set true, charge density is printed out in cube format 

    Key name:     cube_states_list
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      ""
    Allowed:      
    Description:  plot the states listed here 

    Key name:     cube_vh
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  if set true, hatree potential is printed out in cube format 

    Key name:     dftd3_version
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    2
    Max value:    6
    Default:      3
    Description:  Grimme's DFT-D3 versions, 

    Key name:     dipole_correction
    Required:     no
    Key type:     integer array
    Expert:       No
    Experimental: No
    Default:      "0 0 0 "
    Description:  (1,1,1) for molecule, dipole correction in all directions. 

    Key name:     dipole_moment
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Turns on calculation of dipole moment for the entire cell. 

    Key name:     ecutrho
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    10000.000000
    Default:      0.000000e+00
    Description:  ecut for rho in unit of Ry. 

    Key name:     ecutwfc
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    10000.000000
    Default:      0.000000e+00
    Description:  ecut for wavefunctions in unit of Ry. 

    Key name:     electric_field_magnitude
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      0.000000e+00
    Description:  Magnitude of external electric field. 

    Key name:     electric_field_vector
    Required:     no
    Key type:     double array
    Expert:       No
    Experimental: No
    Default:      "Not done yet"
    Description:  Components of the electric field. 

    Key name:     equal_initial_density
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Specifies whether to set initial up and down density to be equal. 

    Key name:     exx_int_flag
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  if set true, calculate the exact exchange integrals 

    Key name:     fast_density
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  Use a faster but less accurate method to generate the charge 
                  density from the electronic wavefunctions. As the cutoff 
                  (grid-density) increases this method improves in accuracy. This 
                  option should be set to false if you receive warnings about 
                  negative charge densities after interpolation. 

    Key name:     fd_allocation_limit
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1024
    Max value:    262144
    Default:      65536
    Description:  Allocation sizes in finite difference routines less than this 
                  value are stack rather than heap based. 

    Key name:     frac_symmetry
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  For supercell calculation, one can disable the fractional 
                  translation symmetry 

    Key name:     freeze_occupied
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: Yes
    Default:      "false"
    Description:  Flag indicating whether or not to freeze the density and occupied 
                  orbitals after a restart. 

    Key name:     gw_residual_convergence_criterion
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: Yes
    Min value:    1.000000e-14
    Max value:    4.000000e-04
    Default:      1.000000e-06
    Description:  The max value of the residual for unoccupied orbitals when 
                  performing a GW calculation. 

    Key name:     gw_residual_fraction
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: Yes
    Min value:    0.000000e+00
    Max value:    1.000000
    Default:      0.900000
    Description:  The residual value specified by gw_residual_convergence_criterion 
                  is applied to this fraction of the total spectrum. 

    Key name:     kohn_sham_ke_fft
    Required:     no
    Key type:     boolean
    Expert:       Yes
    Experimental: No
    Default:      "false"
    Description:  Special purpose flag which will force use of an FFT for the 
                  kinetic energy operator. 

    Key name:     kpoint_distribution
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    -2147483647
    Max value:    2147483647
    Default:      -1
    Description:  

    Key name:     laplacian_autocoeff
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  if set to true, we use LaplacianCoeff.cpp to generate coeff 

    Key name:     laplacian_offdiag
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  if set to true, we use LaplacianCoeff.cpp to generate coeff 

    Key name:     lcao_use_empty_orbitals
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Some pseudopotentials contain unbound atomic orbitals and this 
                  flag indicates whether or not they should be used for LCAO starts. 

    Key name:     md_steps_offset
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    0
    Default:      0
    Description:  

    Key name:     num_wanniers
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    2147483647
    Default:      0
    Description:  number of wannier functions to be used in wannier90 

    Key name:     output_rho_xsf
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Generate xsf format for electronic density. 

    Key name:     poisson_mg_levels
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    -1
    Max value:    6
    Default:      -1
    Description:  Number of multigrid levels to use in the hartree multigrid solver. 

    Key name:     restart_tddft
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  restart TDDFT 

    Key name:     rmg2bgw
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: Yes
    Default:      "false"
    Description:  Write wavefunction in G-space to BerkeleyGW WFN file. 

    Key name:     rmg_threads_per_node
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    64
    Default:      0
    Description:  Number of Multigrid/Davidson threads each MPI process will use. A 
                  value of 0 means set automatically. 

    Key name:     scf_steps_offset
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    0
    Default:      0
    Description:  

    Key name:     sqrt_interpolation
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Flag indicating whether or not to use square root technique for 
                  density interpolation. 

    Key name:     system_charge
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    -unlimited
    Max value:    unlimited
    Default:      0.000000e+00
    Description:  Number of excess holes in the system (useful for doped systems). 
                  Example, 2 means system is missing two electrons 

    Key name:     tddft_mode
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "electric field"
    Allowed:      "point charge" "electric field" 
    Description:  TDDFT mode 

    Key name:     tddft_qgau
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      1.000000
    Description:  Gaussian parameter for point charge to Gaussian charge 

    Key name:     tddft_qpos
    Required:     no
    Key type:     double array
    Expert:       No
    Experimental: No
    Default:      "Not done yet"
    Description:  cartesian coordinate of the point charge for tddft 

    Key name:     total_scf_steps_offset
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    0
    Default:      0
    Description:  

    Key name:     use_cpdgemr2d
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  if set to true, we use Cpdgemr2d to change matrix distribution 

    Key name:     use_symmetry
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  For non-gamma point, always true, for gamma point, optional 

    Key name:     vdwdf_grid_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Coarse"
    Allowed:      "Fine" "Coarse" 
    Description:  Type of grid to use when computing vdw-df correlation. 

    Key name:     verbose
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Flag for writing out extra information 

    Key name:     vxc_diag_nmax
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    10000
    Default:      1
    Description:  Maximum band index for diagonal Vxc matrix elements. 

    Key name:     vxc_diag_nmin
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    10000
    Default:      1
    Description:  Minimum band index for diagonal Vxc matrix elements. 

    Key name:     wannier90_scdm
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    -2147483647
    Max value:    2
    Default:      0
    Description:  use scdm method to set the trial wannier functions 

    Key name:     wannier90_scdm_sigma
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      1.000000
    Description:  when wannier90 is used to build wannier functions, the energy 
                  window parameter 

    Key name:     write_orbital_overlaps
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  If true the orbital overlap matrix from successive MD steps is 
                  written. 

    Key name:     write_pdos
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Flag to write partial density of states. 

    Key name:     z_average_output_mode
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "None"
    Allowed:      "potential and charge density" "wave functions" "None" 
    Description:  z_average_output_mode. 

''' 

undocumented_options = '''
Undocumented options

    Key name:     use_bessel_projectors
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"

    Key name:     kpoints
    Required:     no
    Key type:     formatted
    Expert:       No
    Experimental: No

    Key name:     kpoints_bandstructure
    Required:     no
    Key type:     formatted
    Expert:       No
    Experimental: No

    Key name:     atoms
    Required:     no
    Key type:     formatted
    Expert:       No
    Experimental: No
    
'''

deprecated_options = '''
Deprecated options

    Key name:     alpha
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No

    Key name:     beta
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No

    Key name:     gamma
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No

    Key name:     length_units
    Required:     no
    Key type:     string
    Allowed:      "Angstrom" "Bohr"
    Expert:       No
    Experimental: No

    Key name:     projector_mixing
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No

'''

raw_input_spec += undocumented_options
raw_input_spec += deprecated_options



def read_string(v):
    return v.strip().strip('"')
#end def read_string

truefalse_read = {'"true"':True,'"false"':False,'true':True,'false':False}
def read_boolean(v):
    return truefalse_read[v.lower()]
#end def read_boolean

def read_integer(v):
    return int(v.strip().strip('"'))
#end def read_integer

def read_double(v):
    return float(v.strip().strip('"'))
#end def read_double

def read_integer_array(v):
    return np.array(v.strip().strip('"').split(),dtype=int)
#end def read_integer_array

def read_double_array(v):
    return np.array(v.strip().strip('"').split(),dtype=float)
#end def read_double_array


def write_string(v):
    return '"'+v.strip()+'"'
#end def write_string

truefalse_write = {True:'"true"',False:'"false"'}
def write_boolean(v):
    return truefalse_write[v]
#end def write_boolean

def write_integer(v):
    return '"{}"'.format(v)
#end def write_integer

double_fmt = '{: 16.8f}'
double_fmt_exp = '{: 16.8e}'
def write_double(v):
    vs = double_fmt.format(v).strip()
    if abs((float(vs)-v))>1e-6*abs(v):
        vs = double_fmt_exp.format(v).strip()
    #end if
    return '"'+vs+'"'
#end def write_double

def write_integer_array(v):
    s = '"'
    if isinstance(v,np.ndarray):
        v = v.flatten()
    #end if
    for i in v:
        s+='{} '.format(i)
    #end for
    return s[:-1]+'"'
#end def write_integer_array

def write_double_array(v):
    s = '"'
    if isinstance(v,np.ndarray):
        v = v.flatten()
    #end if
    for d in v:
        s+=double_fmt.format(d).strip()+' '
    #end for
    return s[:-1]+'"'
#end def write_double_array


rmg_value_types = obj({
    'string'        : (str  , np.string_),
    'boolean'       : (bool , np.bool_  ),
    'integer'       : (int  , np.int_   ),
    'double'        : (float, np.float_ ),
    'integer array' : (tuple,list,np.ndarray),
    'double array'  : (tuple,list,np.ndarray),
    })

read_functions = obj({
    'string'        : read_string,
    'boolean'       : read_boolean,
    'integer'       : read_integer,
    'double'        : read_double,
    'integer array' : read_integer_array,
    'double array'  : read_double_array,
    })

write_functions = obj({
    'string'        : write_string,
    'boolean'       : write_boolean,
    'integer'       : write_integer,
    'double'        : write_double,
    'integer array' : write_integer_array,
    'double array'  : write_double_array,
    })

rmg_array_dtypes = obj({
    'integer array' : int,
    'double array'  : float,
    })


class RmgKeyword(DevBase):
    def __init__(self,key_spec,section_name):
        self.key_name     = None
        self.key_type     = None
        self.default      = None
        self.allowed      = None
        self.min_value    = None
        self.max_value    = None
        self.required     = None
        self.expert       = None
        self.experimental = None
        self.description  = None

        self.value_type   = None
        self.array_dtype  = None
        self.section_name = section_name

        spec = obj()
        name  = None
        value = None
        for line in key_spec.strip().splitlines():
            if ':' in line:
                if name is not None:
                    spec[name] = value
                #end if
                name,value = line.split(':',1)
                name  = name.strip().lower().replace(' ','_')
                value = value.strip()
            else:
                value += ' '+line.strip()
            #end if
        #end for
        if name is not None:
            spec[name] = value
        else:
            self.error('Invalid keyword specification text received.\nNo field names are present.\nInvalid spec: {}'.format(key_spec))
        #end if

        for name,value in spec.items():
            if name not in self:
                kname = 'unknown'
                if 'key_name' in spec:
                    kname = spec.key_name
                #end if
                self.error('Unrecognized keyword specification field.\nKeyword: {}\nField name: {}\nField value: {}'.format(kname,name,value))
            #end if
            self[name] = value
        #end for

        if self.key_name is None:
            self.error('Invalid keyword specification received.\nKey name must be defined.\nInvalid spec: {}'.format(key_spec))
        #end if
        if self.key_type is None:
            self.error('Invalid keyword specification received.\nKey type must be defined.\nInvalid spec: {}'.format(key_spec))
        #end if

        if self.key_type=='formatted':
            return
        #end if

        if self.key_type not in read_functions:
            self.error('Read function has not been implemented for key type "{}".'.format(self.key_type))
        #end if
        if self.key_type not in write_functions:
            self.error('Write function has not been implemented for key type "{}".'.format(self.key_type))
        #end if

        read_function = read_functions[self.key_type]

        yesno = dict(yes=True,no=False)
        if self.required is None:
            self.required = False
        else:
            self.required = yesno[self.required.lower()]
        #end if
        if self.expert is None:
            self.expert = False
        else:
            self.expert = yesno[self.expert.lower()]
        #end if
        if self.experimental is None:
            self.experimental = False
        else:
            self.experimental = yesno[self.experimental.lower()]
        #end if
        if self.min_value is not None:
            if 'unlimited' in self.min_value:
                self.min_value = None
            else:
                self.min_value = self.read(self.min_value)
            #end if
        #end if
        if self.max_value is not None:
            if 'unlimited' in self.max_value:
                self.max_value = None
            else:
                self.max_value = self.read(self.max_value)
            #end if
        #end if
        if self.default is not None and self.default!='"Not done yet"':
            self.default = self.read(self.default)
        #end if
        if self.allowed is not None:
            tokens = []
            i1 = -1
            i2 = -1
            for i,c in enumerate(self.allowed):
                if c=='"':
                    if i1==-1:
                        i1 = i
                    else:
                        i2 = i
                        tokens.append(self.read(self.allowed[i1:i2+1]))
                        i1 = -1
                        i2 = -1
                    #end if
                #end if
            #end for
            self.allowed = set(tokens)
            if len(self.allowed)==0:
                self.allowed = None
            #end if
        #end if

        self.value_type = rmg_value_types[self.key_type]
        if self.key_type in rmg_array_dtypes:
            self.array_dtype  = rmg_array_dtypes[self.key_type]
        #end if

    #end def __init__


    def read(self,value):
        return read_functions[self.key_type](value)
    #end def read


    def write(self,value):
        return write_functions[self.key_type](value)
    #end def write


    def assign(self,value):
        if not isinstance(value,self.value_type):
            self.error('cannot assign RMG keyword "{}".\nInvalid type encountered.\nType encoutered: {}\nType(s) expected: {}'.format(self.key_name,value.__class__.__name__,self.value_type))
        #end if
        if self.array_dtype is not None:
            return np.array(value,dtype=self.array_dtype)
        else:
            return value
        #end if
    #end def assign


    def valid(self,value,message=False):
        msg   = ''
        if not isinstance(value,self.value_type):
            msg += 'Keyword "{}" has the wrong type.\n  Type expected: {}\n  Type provided: {}\n'.format(self.key_name,self.key_type,value.__class__.__name__)
        else:
            if RmgInputSettings.enforce_min_value:
                if self.min_value is not None and value<self.min_value:
                    msg += 'Value for keyword "{}" is smaller than allowed.\n  Minimum value allowed: {}\n  Value provided: {}\n'.format(self.key_name,self.min_value,value)
                #end if
            #end if
            if RmgInputSettings.enforce_max_value:
                if self.max_value is not None and value>self.max_value:
                    self.warn('Value for keyword "{}" is larger than allowed.\n  Maximum value allowed: {}\n  Value provided: {}\n'.format(self.key_name,self.max_value,value))
                #end if
            #end if
            if RmgInputSettings.enforce_allowed:
                if self.allowed is not None and value not in self.allowed:
                    msg += 'Value for keyword "{}" is not allowed.\n  Value provided: {}\n  Allowed values: {}'.format(self.key_name,value,list(sorted(self.allowed)))
                #end if
            #end if
        #end if
        if not message:
            return len(msg)==0
        else:
            return len(msg)==0,msg
        #end if
    #end def valid
#end class RmgKeyword



class FormattedRmgKeyword(RmgKeyword):
    def read(self,value):
        self.not_implemented()
    #end def read

    def write(self,value):
        self.not_implemented()
    #end def write

    def assign(self,value):
        self.not_implemented()
    #end def assign

    def valid(self,value,message=False):
        valid = self.valid_no_msg(value)
        if not message:
            return valid
        else:
            return valid,'Data for keyword "{}" is invalid.\nInvalid value: {}'.format(self.key_name,value)
        #end if
    #end def valid

    def valid_no_msg(self,value):
        self.not_implemented()
    #end def valid_no_msg
#end class FormattedRmgKeyword



class FormattedTableRmgKeyword(FormattedRmgKeyword):
    array_options  = None
    array_types    = None
    exclude_fields = set()

    def assign(self,value):
        if isinstance(value,str):
            return value
        elif isinstance(value,(dict,obj)):
            for k,v in value.items():
                if isinstance(v,(tuple,list)):
                    value[k] = np.array(v)
                #end if
            #end for
            return value
        else:
            self.error('cannot assign RMG keyword "{}".\nInvalid type encountered.\nType encoutered: {}\nType(s) expected: str,dict,obj'.format(self.key_name,value.__class__.__name__))
        #end if
    #end def assign

    def valid_no_msg(self,value):
        cls = self.__class__
        if not isinstance(value,obj):
            return False
        #end if
        keys = set(value.keys())-set(cls.exclude_fields)
        match = False
        for key_set in cls.array_options:
            if keys==key_set:
                match = True
                break
            #end if
        #end if
        if not match:
            return False
        #end if
        lengths = [len(value[k]) for k in keys]
        if lengths[0]==0:
            return False
        #end if
        if len(set(lengths))!=1:
            return False
        #end if
        for k in keys:
            v = value[k]
            if not isinstance(v,np.ndarray):
                return False
            elif not isinstance(v.flatten()[0],cls.array_types[k]):
                return False
            #end if
        #end for
        return True
    #end def valid_no_msg
#end class FormattedTableRmgKeyword



class PseudopotentialKeyword(FormattedTableRmgKeyword):
    array_options = [
        set(('species','pseudos')),
        ]
    array_types   = obj(
        species = rmg_value_types.string,
        pseudos = rmg_value_types.string,
        )

    def read(self,value):
        d = np.array(value.split(),dtype=str)
        d.shape = len(d)//2,2
        species = d[:,0].flatten()
        pseudos = d[:,1].flatten()
        return obj(species=species,pseudos=pseudos)
    #end def read

    def write(self,value):
        v = value
        s = '"\n'
        for (sp,p) in zip(v.species,v.pseudos):
            s += '{:<4} {}\n'.format(sp,p)
        #end for
        s += '"'
        return s
    #end def write
#end class PseudopotentialKeyword



class KpointsKeyword(FormattedTableRmgKeyword):
    array_options = [
        set(('kpoints','weights')),
        ]
    array_types   = obj(
        kpoints = rmg_value_types.double,
        weights = rmg_value_types.double,
        )

    def read(self,value):
        d = np.array(value.split(),dtype=float)
        d.shape = len(d)//4,4
        kpoints = d[:,:3]
        weights = d[:,-1].flatten()
        return obj(kpoints=kpoints,weights=weights)
    #end def read

    def write(self,value):
        v = value
        s = '"\n'
        for (kp,w) in zip(v.kpoints,v.weights):
            s += '{: 16.12f} {: 16.12f} {: 16.12f} {: 16.12f}\n'.format(kp[0],kp[1],kp[2],w)
        #end for
        s += '"'
        return s
    #end def write
#end class KpointsKeyword



class KpointsBandstructureKeyword(FormattedTableRmgKeyword):
    array_options = [
        set(('kpoints','counts','labels')),
        ]
    array_types   = obj(
        kpoints = rmg_value_types.double,
        counts  = rmg_value_types.integer,
        labels  = rmg_value_types.string,
        )

    def read(self,value):
        d = np.array(value.split(),dtype=str)
        d.shape = len(d)//4,5
        kpoints = np.array(d[:,:3],dtype=float)
        counts  = np.array(d[:,3],dtype=int).flatten()
        labels  = d[:,-1].flatten()
        return obj(kpoints=kpoints,counts=counts,labels=labels)
    #end def read

    def write(self,value):
        v = value
        s = '"\n'
        for (kp,c,l) in zip(v.kpoints,v.counts,v.labels):
            s += '{: 16.12f} {: 16.12f} {: 16.12f}  {:>3}  {}\n'.format(kp[0],kp[1],kp[2],c,l)
        #end for
        s += '"'
        return s
    #end def write
#end class KpointsBandstructureKeyword



class AtomsKeyword(FormattedTableRmgKeyword):

    formats = ('basic','movable','movable_moment','moment','spin_ratio','full_spin')

    array_options = [
        set(('atoms','positions')),
        set(('atoms','positions','movable')),
        set(('atoms','positions','moments')),
        set(('atoms','positions','movable','moments')),
        set(('atoms','positions','movable','spin_ratio')),
        set(('atoms','positions','movable','spin_ratio','spin_theta','spin_phi')),
        ]
    array_types   = obj(
        atoms      = rmg_value_types.string,
        positions  = rmg_value_types.double,
        movable    = rmg_value_types.boolean,
        moments    = rmg_value_types.double,
        spin_ratio = rmg_value_types.double,
        spin_theta = rmg_value_types.double,
        spin_phi   = rmg_value_types.double,
        )
    exclude_fields = ['format']


    def read(self,value):
        # check if input data is empty
        value = value.strip()
        if len(value)==0:
            self.error('No data provided for "atoms".')
        #end if
        
        # determine the number of values per line
        if '\n' in value:
            first,rest = value.split('\n',1)
        else:
            first = value
        #end if
        nvals = len(first.split())

        # initial array parse of value table
        d = np.array(value.split(),dtype=str)
        d.shape = len(d)//nvals,nvals

        # extract universal atom labels and positions
        atom_labels = d[:,0]
        atom_labels = atom_labels.flatten()
        positions   = np.array(d[:,1:4],dtype=float)

        # extract remaining data and determine format
        v = obj(
            format    = None,
            atoms     = atom_labels,
            positions = positions,
            )

        boolset = set(['0','1'])
        invalid_format = False
        if nvals==4:
            v.format = 'basic'
        elif nvals==5:
            if len(set(d[:,4])-boolset)==0:
                v.movable = np.array(d[:,4],dtype=bool)
                v.format  = 'movable'
            else:
                try:
                    v.moments = np.array(d[:,4],dtype=float)
                    v.format  = 'moment'
                except:
                    invalid_format = True
                #end try
            #end if
        elif nvals==6:
            try:
                v.movable = np.array(d[:,4],dtype=bool)
                v.moments = np.array(d[:,5],dtype=float)
                v.format  = 'movable_moment'
            except:
                invalid_format = True
            #end try
        elif nvals==8:
            try:
                assert(len(set(d[:,4:7].flatten())-boolset)==0)
                v.movable    = np.array(d[:,4:7],dtype=bool)
                v.spin_ratio = np.array(d[:,7],dtype=float)
                v.format     = 'spin_ratio'
            except:
                invalid_format = True
            #end try
        elif nvals==10:
            try:
                assert(len(set(d[:,4:7].flatten())-boolset)==0)
                v.movable    = np.array(d[:,4:7],dtype=bool)
                v.spin_ratio = np.array(d[:,7],dtype=float)
                v.spin_theta = np.array(d[:,8],dtype=float)
                v.spin_phi   = np.array(d[:,9],dtype=float)
                v.format     = 'full_spin'
            except:
                invalid_format = True
            #end try
        #end if

        if invalid_format:
            self.error('Failed to read atoms data.\nPlease check the formatting:\n{}'.format(value))
        elif v.format is None or v.format not in AtomsKeyword.formats:
            self.error('Failed to read atoms data.\nThis is a developer error.\nPlease contact the developers.')
        #end if

        return v
    #end def read

    def write(self,value):
        v = value
        s = '"\n'
        if v.format=='basic':
            for (a,p) in zip(v.atoms,v.positions):
                s += '{:<4} {: 16.12f} {: 16.12f} {: 16.12f}\n'.format(a,p[0],p[1],p[2])
            #end for
        elif v.format=='movable':
            for (a,p,m) in zip(v.atoms,v.positions,v.movable):
                s += '{:<4} {: 16.12f} {: 16.12f} {: 16.12f}  {}\n'.format(a,p[0],p[1],p[2],int(m))
            #end for
        elif v.format=='moment':
            for (a,p,m) in zip(v.atoms,v.positions,v.moments):
                s += '{:<4} {: 16.12f} {: 16.12f} {: 16.12f}  {: 6.4f}\n'.format(a,p[0],p[1],p[2],m)
            #end for
        elif v.format=='movable_moment':
            for (a,p,mv,mo) in zip(v.atoms,v.positions,v.movable,v.moments):
                s += '{:<4} {: 16.12f} {: 16.12f} {: 16.12f}  {}  {: 6.4f}\n'.format(a,p[0],p[1],p[2],int(mv),mo)
            #end for
        elif v.format=='spin_ratio':
            for (a,p,m,s) in zip(v.atoms,v.positions,v.movable,v.spin_ratio):
                s += '{:<4} {: 16.12f} {: 16.12f} {: 16.12f} {} {} {} {: 6.4f}\n'.format(a,p[0],p[1],p[2],int(m[0]),int(m[1]),int(m[2]),s)
            #end for
        elif v.format=='full_spin':
            for (a,p,m,sr,st,sp) in zip(v.atoms,v.positions,v.movable,v.spin_ratio,v.spin_theta,v.spin_phi):
                s += '{:<4} {: 16.12f} {: 16.12f} {: 16.12f} {} {} {} {: 6.4f} {: 6.2f} {: 6.2f}\n'.format(a,p[0],p[1],p[2],int(m[0]),int(m[1]),int(m[2]),sr,st,sp)
            #end for
        else:
            self.error('Invalid atoms format encountered on write.\nInvalid format: {}\nValid options are: {}'.format(v.format,self.formats))
        #end if
        s += '"'
        return s
    #end def write
#end class AtomsKeyword



class HubbardUKeyword(RmgKeyword):
    def read(self,value):
        v = obj()
        tokens = read_string(value).split()
        for a,u in zip(tokens[::2],tokens[1::2]):
            v[a] = float(u)
        #end for
        return v
    #end def read

    def write(self,value):
        s = ''
        for a in sorted(value.keys()):
            s += ' {} {}'.format(a,value[a])
        #end for
        return write_string(s)
    #end def write

    def assign(self,value):
        if isinstance(value,str):
            return value
        elif isinstance(value,(dict,obj)):
            return obj(value)
        else:
            self.error('cannot assign RMG keyword "{}".\nInvalid type encountered.\nType encoutered: {}\nType(s) expected: str,dict,obj'.format(self.key_name,value.__class__.__name__))
        #end if
    #end def assign

    def valid(self,value,message=False):
        valid = True
        for k,v in value.items():
            if not isinstance(k,rmg_value_types.string):
                valid = False
                break
            elif not isinstance(v,rmg_value_types.double):
                valid = False
                break
            #end if
        #end for
        if not message:
            return valid
        else:
            return valid,'Data for keyword "{}" is invalid.\nInvalid value: {}'.format(self.key_name,value)
        #end if
    #end def valid
#end class HubbardUKeyword
    

formatted_keywords = obj(
    pseudopotential       = PseudopotentialKeyword,
    kpoints               = KpointsKeyword,
    kpoints_bandstructure = KpointsBandstructureKeyword,
    atoms                 = AtomsKeyword,
    Hubbard_U             = HubbardUKeyword,
    )



class RmgInputSpec(DevBase):
    def __init__(self):
        spec = raw_input_spec.strip()

        blocks = spec.split('\n\n')
    
        self.section_order    = []
        self.section_labels   = obj()
        self.section_contents = obj()
        self.keywords         = obj()

        section  = None
        sec_cont = None
        for b in blocks:
            if ':' not in b:
                b = b.strip().lower()
                if b.endswith('options'):
                    sec_cont = []
                    section  = b.replace(' options','').strip().replace('  ',' ').replace(' ','_')
                    self.section_order.append(section)
                    self.section_labels[section]   = b
                    self.section_contents[section] = sec_cont
                #end if
            else:
                k = RmgKeyword(b,section)
                if k.key_type=='formatted':
                    if k.key_name not in formatted_keywords:
                        self.error('unrecognized formatted keyword: "{}"'.format(k.key_name))
                    #end if
                    k = formatted_keywords[k.key_name](b,section)
                #end if
                self.keywords[k.key_name] = k
                sec_cont.append(k.key_name)
            #end if
        #end for
    #end def __init__
#end class RmgInputSpec

input_spec = RmgInputSpec()


class RmgCalcModes(DevBase):
    def __init__(self):
        self.full_calc = obj(
            scf         = 'Quench Electrons',
            exx         = 'Exx Only',
            neb         = 'NEB Relax', 
            band        = 'Band Structure Only',
            relax       = 'Relax Structure',
            dimer_relax = 'Dimer Relax',
            md_PE       = 'Constant Pressure And Energy',
            md_TE       = 'Constant Temperature And Energy',
            md_VE       = 'Constant Volume And Energy',
            tddft       = 'TDDFT',
            plot        = 'Plot',
            psi_plot    = 'Psi Plot',
            )

        self.short_calc = obj()
        for k,v in self.full_calc.items():
            self.short_calc[v] = k
        #end for
        self.full_calc_modes  = set(self.full_calc.values())
        self.short_calc_modes = set(self.short_calc.values())
    #end def __init__

    def is_full_mode(self,mode):
        return mode in self.full_calc_modes
    #end def is_full_mode

    def is_short_mode(self,mode):
        return mode in self.short_calc_modes
    #end def is_short_mode

    def full_mode(self,short_mode):
        mode = None
        if short_mode in self.full_calc:
            mode = self.full_calc[short_mode]
        #end if
        return mode
    #end def full_mode

    def short_mode(self,full_mode):
        mode = None
        if full_mode in self.short_calc:
            mode = self.short_calc[full_mode]
        #end if
        return mode
    #end def short_mode

    def mode_match(self,text,short=False):
        mode = None
        text = text.lower()
        for full_mode in self.full_calc_modes:
            if full_mode.lower() in text:
                if not short:
                    mode = full_mode
                else:
                    mode = self.short_mode(full_mode)
                #end if
                break
            #end if
        #end for
        return mode
    #end def mode_match
        
#end class RmgCalcModes

rmg_modes = RmgCalcModes()







class RmgInput(SimulationInput):
    def __init__(self,filepath=None):
        if filepath is not None:
            self.read(filepath)
        #end if
    #end def __init__


    def assign(self,**values):
        unrecognized = []
        for k,v in values.items():
            if k in input_spec.keywords:
                if isinstance(v,(str,np.string_)):
                    self[k] = input_spec.keywords[k].read(v)
                else:
                    self[k] = input_spec.keywords[k].assign(v)
                #end if
            else:
                unrecognized.append(k)
            #end if
        #end for
        if len(unrecognized)>0:
            unrec = obj(values).obj(unrecognized)
            self.error('Unrecognized keywords encountered during assignment.\nUnrecognized keywords: {}\nCorresponding values:\n{}'.format(list(sorted(unrecognized)),unrec))
        #end if
    #end def assign


    def read_text(self,contents,filepath=None):
        # remove comments and whitespace
        text = ''
        for line in contents.splitlines():
            i = line.find('#')
            if i!=-1:
                line = line[:i]
            #end if
            ls = line.strip()
            if len(ls)>0:
                text += ls+'\n'
            #end if
        #end for
        text = text.strip()

        # separate keywords and values
        values = obj()
        whitespace = ' \t\n'
        icur = 0
        k = 'some key'
        while k is not None:
            k  = None
            ie = text.find('=',icur)
            if ie!=-1:
                k   = text[icur:ie].strip()
                iv1 = text.find('"',ie)
                if iv1!=-1:
                    iv2 = text.find('"',iv1+1)
                    if iv2!=-1:
                        v         = text[iv1+1:iv2]
                        values[k] = v
                        icur      = iv2+1
                    #end if
                #end if
            #end if
        #end while

        # read the keyword values, checking for unrecognized ones
        self.assign(**values)
    #end def read_text


    def write_text(self,filepath=None):
        if RmgInputSettings.check_on_write:
            self.check_valid()
        #end if
        text = ''
        for section_name in input_spec.section_order:
            present = False
            for k in input_spec.section_contents[section_name]:
                if k in self:
                    if not present:
                        text += '\n\n# '+input_spec.section_labels[section_name]+'\n\n'
                        present = True
                    #end if
                    kw = input_spec.keywords[k]
                    text += '{:<22} = {}\n'.format(kw.key_name,kw.write(self[k]))
                #end if
            #end for
        #end for
        return text.lstrip()
    #end def write_text


    def check_valid(self,exit=True):
        msg = ''
        allowed = set(input_spec.keywords.keys())
        present = set(self.keys())
        unrecognized = present-allowed
        if len(unrecognized)>0:
            msg += 'Unrecognized keywords encountered.\n  Unrecognized keywords: {}\n  Valid keywords are: {}\n'.format(list(sorted(unrecognized)),list(sorted(allowed)))
        #end if
        recognized = present-unrecognized
        for k in sorted(recognized):
            kval,m = input_spec.keywords[k].valid(self[k],message=True)
            if not kval:
                msg += m+'\n'
            #end if
        #end if
        if len(msg)>0 and exit:
            self.log(msg)
            self.error('Input is invalid.\nPlease see messages above for specific issues.')
        #end if
        return len(msg)==0
    #end def check_valid


    def is_valid(self):
        return self.check_valid(exit=False)
    #end def is_valid


    def return_structure(self,units='B'):
        axes       = self.get('lattice_vector',None)
        axes_unit  = self.get('lattice_units','bohr')
        lattice    = self.get('bravais_lattice_type','orthorhombic primitive')
        a          = self.get('a_length',0.0)
        b          = self.get('b_length',0.0)
        c          = self.get('c_length',0.0)

        coord_type = self.get('atomic_coordinate_type','absolute')
        coord_unit = self.get('crds_units','bohr')
        atom_data  = self.get('atoms',obj())
        atoms      = atom_data.get('atoms',None)
        positions  = atom_data.get('positions',None)

        unit_dict = dict(angstrom='A',bohr='B')
        coord_unit = unit_dict[coord_unit.lower()]
        axes_unit  = unit_dict[axes_unit.lower()]

        if axes is not None:
            axes = np.array(axes,dtype=float)
        else:
            lattice_orig = lattice
            lattice = lattice.lower()
            if lattice=='cubic primitive':
                axes = np.diag((a,a,a))
            elif lattice=='tetragonal primitive':
                axes = np.diag((a,a,c))
            elif lattice=='orthorhombic primitive':
                axes = np.diag((a,b,c))
            elif lattice=='cubic body centered':
                axes = 0.5*a*np.array([[ 1, 1,-1],
                                       [-1, 1, 1],
                                       [ 1,-1, 1]],dtype=float)
            elif lattice=='cubic face centered':
                axes = 0.5*a*np.array([[ 1, 1, 0],
                                       [ 0, 1, 1],
                                       [ 1, 0, 1]],dtype=float)
            elif lattice=='hexagonal primitive':
                axes = np.array([[   a,              0, 0],
                                 [-a/2, np.sqrt(3)/2*a, 0],
                                 [   0,              0, c]],dtype=float)
            else:
                # cubic body centered, hexagonal primitive not yet supported
                self.error('Structure extraction failed.\nLattice type "{}" is currently unsupported.'.format(lattice_orig))
            #end if
        #end if
        axes = convert(axes,axes_unit,units)

        if atoms is None or positions is None:
            self.error('Structure extraction failed.\nEither atoms or positions could not be obtained.')
        #end if
        atoms     = np.array(atoms,dtype=object)
        positions = np.array(positions,dtype=float)
        if coord_type.lower()=='cell relative':
            positions = np.dot(positions,axes)
        else:
            positions = convert(positions,coord_unit,units)
        #end if
        
        s = generate_structure(
            units = units,
            axes  = axes,
            elem  = atoms,
            pos   = positions,
            )

        return s
    #end def return_structure
#end class RmgInput




def generate_rmg_input(**kwargs):
    selector = kwargs.pop('input_type','generic')
    if selector=='generic':
        return generate_any_rmg_input(**kwargs)
    else:
        RmgInput.class_error('Input type "{}" has not been implemented for RMG input generation.'.format(selector))
    #end if
#end def generate_rmg_input



generate_any_defaults = obj(
    none  = obj(),
    basic = obj(
        use_folded             = True,
        virtual_frac           = 0.20,
        ),
    #qmc   = obj(
    #    use_folded             = True,
    #    ),
    )

def generate_any_rmg_input(**kwargs):
    loc = 'generate_rmg_input'

    # set default values
    defaults = kwargs.pop('defaults','basic')
    kw = obj(**kwargs)
    kw.set_optional(generate_any_defaults[defaults])

    # extract keywords not appearing in RMG input file
    text            = kw.delete_optional('text'           , None   )
    wf_grid_spacing = kw.delete_optional('wf_grid_spacing', None   )
    pseudos         = kw.delete_optional('pseudos'        , None   )
    system          = kw.delete_optional('system'         , None   )
    copy_system     = kw.delete_optional('copy_system'    , True   )
    use_folded      = kw.delete_optional('use_folded'     , False  )
    virtual_frac    = kw.delete_optional('virtual_frac'   , None   )
    spin_polarized  = kw.delete_optional('spin_polarized' , None   )
    default_units   = kw.delete_optional('default_units'  , 'bohr' )

    default_units = dict(
        a        = 'angstrom',
        b        = 'bohr',
        angstrom = 'angstrom',
        bohr     = 'bohr'
        )[default_units.lower()]

    rmg_units_map = obj(
        angstrom = 'Angstrom',
        bohr     = 'Bohr',
        alat     = 'Alat',
        a        = 'Angstrom',
        b        = 'Bohr',
        )

    # generate RMG input
    ri = RmgInput()
    if text is not None:
        ri.read_text(text)
    #end if
    ri.assign(**kw)

    # incorporate pseudopotentials details provided via "pseudos"
    if pseudos is not None:
        species = []
        pps     = []
        for ppname in pseudos:
            label,element = pp_elem_label(ppname,guard=True)
            species.append(element)
            pps.append(ppname)
        #end for
        ri.pseudopotential = obj(
            species = np.array(species),
            pseudos = np.array(pps),
            )
    #end if
    
    # incorporate system details, if provided
    if system is not None:

        # add system details
        if copy_system:
            system = system.copy()
        #end if
        if use_folded:
            system = system.get_smallest()
        #end if
        system.check_folded_system()
        system.update_particles()

        # set atomic species, positions, magnetic moments and mobility
        if 'atomic_coordinate_type' not in ri:
            ri.atomic_coordinate_type = 'Absolute'
        #end if
        if 'crds_units' not in ri:
            cu = default_units
            
        else:
            cu = ri.crds_units.lower()
        #end if
        if cu=='angstrom':
            system.change_units('A')
        elif cu=='bohr':
            system.change_units('B')
        else:
            error('Invalid crds_units.\nExpected "Angstrom" or "Bohr".\nReceived: {}'.format(cu),loc)
        #end if
        rmg_length_units = rmg_units_map[cu]
        if 'crd_units' not in ri and 'atomic_coordinate_type' in ri and ri.atomic_coordinate_type=='Absolute':
            ri.crd_units = rmg_length_units
        #end if
        s = system.structure
        elem = np.array(s.elem)
        act = ri.atomic_coordinate_type.lower()
        if act=='absolute':
            pos = s.pos.copy()
        elif act=='cell relative':
            pos = s.pos_unit().copy()
        else:
            error('Invalid atomic_coordinate_type.\nExpected "Absolute" or "Cell Relative".\nReceived: {}'.format(cu),loc)
        #end if
        movable = None
        if s.frozen is not None:
            movable = ~s.is_frozen()
        #end if
        moments = None
        if s.mag is not None:
            moments = np.array(s.mag,dtype=float)
        #end if
        if movable is not None and moments is not None:
            ri.atoms = obj(
                format = 'movable_moment',
                atoms     = elem,
                positions = pos,
                movable   = movable,
                moments   = moments,
                )
        elif movable is not None:
            ri.atoms = obj(
                format = 'movable',
                atoms     = elem,
                positions = pos,
                movable   = movable,
                )
        else:
            ri.atoms = obj(
                format = 'basic',
                atoms     = elem,
                positions = pos,
                )
        #end if

        # set lattice vectors
        if 'a_length' not in ri:
            ri.lattice_units  = rmg_length_units
            ri.lattice_vector = s.axes.copy()
        #end if

        # set kpoints
        if len(s.kpoints)>0 and 'kpoint_mesh' not in ri:
            kpu = s.kpoints_unit()
            ri.kpoints = obj(
                kpoints = kpu.copy(),
                weights = s.kweights.copy(),
                )
            if 'kpoint_is_shift' in ri:
                del ri.kpoint_is_shift
            #end if
        #end if

        # set wavefunction grid
        if wf_grid_spacing is not None:
            wf_grid = []
            for a in system.structure.axes:
                g = int(np.ceil(np.linalg.norm(a)/wf_grid_spacing))
                wf_grid.append(g)
            #end for
            ri.assign(wavefunction_grid=wf_grid)
        #end if

        if spin_polarized is None and 'noncollinear' not in ri:
            spin_polarized = system.spin_polarized_orbitals()
        elif spin_polarized is None and 'noncollinear' in ri:
            if not ri.noncollinear:
                spin_polarized = system.spin_polarized_orbitals()
            #end if
        #end if

        # set occupations
        has_states = False
        has_states |= 'states_count_and_occupation' in ri
        has_states |= 'states_count_and_occupation_up' in ri and 'states_count_and_occupation_down' in ri
        if not has_states and virtual_frac is not None:
            states_keys = [
                'states_count_and_occupation',
                'states_count_and_occupation_up',
                'states_count_and_occupation_down',
                ]
            for k in states_keys:
                if k in ri:
                    del ri[k]
                #end if
            #end for
            nup,ndn = system.particles.electron_counts()
            nvirt = int(np.ceil(virtual_frac*max(nup,ndn)))
            nptot = max(nup,ndn) + nvirt
            nup_virt = nptot-nup
            ndn_virt = nptot-ndn
            if nup==ndn and not spin_polarized:
                occ_up = '{} 2.0 {} 0.0'.format(nup,nup_virt)
                ri.states_count_and_occupation = occ_up
            else:
                occ_up = '{} 1.0 {} 0.0'.format(nup,nup_virt)
                occ_dn = '{} 1.0 {} 0.0'.format(ndn,ndn_virt)
                ri.states_count_and_occupation_spin_up   = occ_up
                ri.states_count_and_occupation_spin_down = occ_dn
            #end if
        #end if

    #end if

    if spin_polarized is not None and spin_polarized:
        if 'states_count_and_occupation_spin_up' not in ri:
            error('System is spin polarized, but occupations not provided for up and down spins.',loc)
        #end if
    #end if

    return ri
#end def generate_any_rmg_input
