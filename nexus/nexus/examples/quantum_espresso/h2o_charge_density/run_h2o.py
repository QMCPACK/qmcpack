#! /usr/bin/env python

"""Example script using Quantum ESPRESSO to generate a charge density for water.

This script will do the following:
1. Create a vacuum box with a water molecule.
2. Center the molecule in the unit cell.
3. Create and run a geometry relaxation with Quantum ESPRESSO.
4. Create and run an SCF calculation with Quantum ESPRESSO.
5. Create and run a post-processing calculation to write the full charge
   density to an XSF file.

If you are a user using this script, you should check a few things first:
1. The paths to ``pw.x`` and ``pp.x``, if your QE ``build/bin``
   directory is not on ``$PATH``.
2. You also should check to make sure the core counts are correct, as
   this script was designed to use 6 cores on an 8-core laptop.
3. Check that you have enough memory for the calculation (~6 GB free).

   .. note::
      If you do not have enough memory, you can try to shrink the vacuum
      padding, though under 8 Å should be considered underconverged w.r.t.
      the final unit cell size.
"""

from __future__ import annotations

import numpy as np
from nexus import (
    PseudoSet,
    Structure,
    generate_physical_system,
    generate_pp,
    generate_pwscf,
    generate_structure,
    job,
    run_project,
    settings,
)


def calc_vac_box(structure: Structure, vacuum: float = 8.0) -> float:
    """Calculate the length of a unit cell for a provided vacuum padding."""
    coords = structure.pos

    min_x, min_y, min_z = np.min(coords, 0)

    coords[:,0] = coords[:,0] - min_x
    coords[:,1] = coords[:,1] - min_y
    coords[:,2] = coords[:,2] - min_z

    max_xyz = np.max(coords + (vacuum * 2), 0)

    return max(max_xyz)

# Nexus settings
settings(
    vdw_table     = "./",
    runs          = "runs",
    results       = "",
    generate_only = False,
    status_only   = False,
    sleep         = 0.5,
    machine       = "ws8",
    progress_tty  = True,
)

pseudos = PseudoSet.from_dir(
    pseudo_dir="../pseudopotentials",
)

#=====================#
#    User Settings    #
#=====================#

pwx_path = "pw.x"
ppx_path = "pp.x"
n_cores  = 6

#=====================#
#  Calculation Setup  #
#=====================#

structure = generate_structure(
    type       = "trimer",
    trimer     = ["O", "H", "H"],
    units      = "A",
    separation = [0.984, 0.984],
    angle      = 102.25,
)

# Distance to edge of vacuum box in Angstrom
# 2x equals distance from mol to mol.
vacuum_size = 8.00
cell_dim = calc_vac_box(structure, vacuum_size)
structure.set_axes(np.eye(3) * cell_dim)

# Move the molecule so it has equal padding on each side
structure.center_molecule()

# Create a `PhysicalSystem` with spin and charge information
mol = generate_physical_system(
    structure  = structure,
    net_charge = 0,
    **pseudos.get_Zeff(structure.elem),
)

job_name = "water"

relax = generate_pwscf(
    # Nexus Simulation Settings
    identifier = f"{job_name}_relax",
    path       = f"{job_name}_relax",
    job        = job(cores=n_cores, app=pwx_path),
    input_type = "generic",

    # QE Namelists
    # &CONTROL
    calculation   = "relax",
    max_seconds   = 1200, # Highly unlikely that this will take more than 20 min
    nstep         = 100,
    prefix        = "relax",
    outdir        = "relax_output",
    etot_conv_thr = 1.00e-6,
    forc_conv_thr = 1.00e-4,

    # &SYSTEM
    ecutwfc         = 80,
    ecutrho         = 320, # Default ecutrho in QE is 4x ecutwfc, but good to be explicit
    input_dft       = "pbe",
    assume_isolated = "mp",
    # Note: Makov-Payne correction is not always the best choice in production runs

    # &ELECTRONS
    conv_thr = 1e-7,

    # &IONS
    ion_dynamics = "bfgs",

    # QE Cards
    # ATOMIC_SPECIES
    pseudos = pseudos,

    # ATOMIC_SPECIES & ATOMIC_POSITIONS & CELL_PARAMETERS
    system = mol,

    # K_POINTS
    kgrid  = (1,1,1),
    kshift = (0,0,0),
)

scf = generate_pwscf(
    # Nexus Simulation Settings
    identifier   = f"{job_name}_scf",
    path         = f"{job_name}_scf",
    job          = job(cores=n_cores, app=pwx_path),
    input_type   = "generic",
    dependencies = (relax, "structure"),

    # QE Namelists
    # &CONTROL
    calculation = "scf",
    max_seconds = 1200,
    nstep       = 100,
    prefix      = "scf",
    outdir      = "scf_output",

    # &SYSTEM
    ecutwfc         = 100,
    ecutrho         = 400,
    input_dft       = "pbe",
    assume_isolated = "mp",
    # Note: Makov-Payne correction is not always the best choice in production runs

    # &ELECTRONS
    conv_thr = 1e-8, # Convergence threshold

    # QE Cards
    # ATOMIC_SPECIES
    pseudos = pseudos,

    # ATOMIC_SPECIES & ATOMIC_POSITIONS & CELL_PARAMETERS
    system = mol,

    # K_POINTS
    kgrid  = (1,1,1),
    kshift = (0,0,0),
)

post_proc = generate_pp(
    # Nexus Settings
    identifier    = f"{job_name}_chg_dens",
    path          = f"{job_name}_scf/chg_dens", # Store in SCF directory.
    job           = job(cores=n_cores, app=ppx_path),
    dependencies  = (scf, "other"), # 'other' means just wait to run until scf is done

    # QE PP Settings
    # &INPUTPP
    filplot       = f"{job_name}.pp",
    prefix        = "../scf_output/scf", # point to parent SCF directory.
    plot_num      = 0,

    # &PLOT
    fileout       = f"{job_name}_chg_dens.xsf",
    output_format = 5,
    iflag         = 3,
)

run_project()
