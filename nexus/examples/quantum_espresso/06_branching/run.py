#! /usr/bin/env python3
"""
Nexus Enhanced Workflow Example: Diamond SCF with Conditional Execution

This example demonstrates the new enhanced workflow features in Nexus:
1. Simulation dependencies (scf2 depends on scf1)
2. Error handlers with automatic resubmission
3. Conditional execution based on simulation results (energy threshold)
4. Branching workflow: scf3 or scf4 runs based on scf2's total energy

The workflow:
- scf1: Initial SCF calculation (no dependencies)
- scf2: SCF calculation that depends on scf1's output
- scf3: Runs if scf2's energy is below threshold (depends on scf2)
- scf4: Runs if scf2's energy is above threshold (depends on scf2)
- Only one of scf3 or scf4 will execute based on the energy condition
- All simulations have error handlers that retry on failure
"""

from nexus import settings, Job, run_project
from nexus import generate_physical_system
from nexus import generate_pwscf

# Add lib directory to path for enhanced workflow imports
import sys
import os
nexus_lib = os.path.join(os.path.dirname(__file__), '../../lib')
if nexus_lib not in sys.path:
    sys.path.insert(0, nexus_lib)

# Import enhanced workflow features
from enhanced_simulation import make_enhanced, create_branch
from error_handler import RetryErrorHandler

# Computer configuration
computer = 'baseline'

if computer == 'baseline':
    qe_modules = 'module purge; module load Core/25.05   gcc/12.4.0   openmpi/5.0.5   DefApps hdf5'
    qe_bin = '/ccsopen/home/ksu/SOFTWARE/qe/q-e-qe-7.4.1/build/bin'
    account = 'phy191'
else:
    print('Undefined computer')
    exit()

settings(
    pseudo_dir = '../pseudopotentials',
    results    = '',
    sleep      = 1,
    machine    = computer,
    account    = account,
    )

# Define physical system (diamond)
system = generate_physical_system(
    units    = 'A',
    axes     = '''1.785   1.785   0.000
                  0.000   1.785   1.785
                  1.785   0.000   1.785''',
    elem_pos = '''
               C  0.0000  0.0000  0.0000
               C  0.8925  0.8925  0.8925
               ''',
    C        = 4,
    )

# Simple example with two simulations and a dependency
# scf1 runs first, scf2 depends on scf1's output

# Create reusable job definition for baseline
qe_job = Job(nodes=1,
             queue='batch',
             hours=1,
             presub=qe_modules,
             app=qe_bin+'/pw.x')

# Simulation 1: Initial SCF calculation (no dependencies)
# Testing without make_enhanced wrapper - this will work but without error handlers
scf1 = generate_pwscf(
    identifier   = 'scf',
    path         = 'diamond/scf',
    job          = qe_job,
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

# Simulation 2: SCF calculation that depends on scf1
# This depends on scf1's charge_density output
scf2_base = generate_pwscf(
    identifier   = 'scf_dependent',
    path         = 'diamond/scf_dependent',
    job          = qe_job,
    input_type   = 'generic',
    calculation  = 'scf',
    input_dft    = 'lda', 
    ecutwfc      = 200,   
    conv_thr     = 1e-8, 
    system       = system,
    pseudos      = ['C.BFD.upf'],
    kgrid        = (4,4,4),
    kshift       = (0,0,0),
    dependencies = (scf1, 'charge_density'),
    )

# Convert to enhanced simulation with error handler
scf2 = make_enhanced(
    scf2_base,
    error_handlers=[RetryErrorHandler(max_retries=2)],
    # required_machine='andes',
    )

# Define energy threshold (in Ry, typical for Quantum ESPRESSO)
# Adjust this value based on your system's expected energy range
energy_threshold = -30.0

# Conditional function: scf3 runs if scf2 energy is below threshold
def energy_below_threshold(sim):
    """Check if scf2's energy is below threshold."""
    try:
        # Get scf2 from dependencies (scf3 depends on scf2)
        scf2_sim = None
        for dep in sim.dependencies.values():
            if dep.sim.identifier == 'scf_dependent':
                scf2_sim = dep.sim
                break
        
        if scf2_sim is None or not scf2_sim.finished:
            return False  # scf2 not finished yet
        
        # Load analyzer and get energy
        analyzer = scf2_sim.load_analyzer_image()
        if hasattr(analyzer, 'E') and analyzer.E != 0.0:
            energy = analyzer.E
            return energy < energy_threshold
        else:
            return False  # Energy not available
    except Exception:
        # If anything fails, don't run
        return False


# Simulation 3: Runs if scf2 energy is below threshold
scf3_base = generate_pwscf(
    identifier   = 'scf_low_energy',
    path         = 'diamond/scf_low_energy',
    job          = qe_job,
    input_type   = 'generic',
    calculation  = 'scf',
    input_dft    = 'lda', 
    ecutwfc      = 200,   
    conv_thr     = 1e-8, 
    system       = system,
    pseudos      = ['C.BFD.upf'],
    kgrid        = (4,4,4),
    kshift       = (0,0,0),
    # Note: dependencies will be set by create_branch
    )

# Simulation 4: Runs if scf2 energy is above threshold
scf4_base = generate_pwscf(
    identifier   = 'scf_high_energy',
    path         = 'diamond/scf_high_energy',
    job          = qe_job,
    input_type   = 'generic',
    calculation  = 'scf',
    input_dft    = 'lda', 
    ecutwfc      = 200,   
    conv_thr     = 1e-8, 
    system       = system,
    pseudos      = ['C.BFD.upf'],
    kgrid        = (4,4,4),
    kshift       = (0,0,0),
    # Note: dependencies will be set by create_branch
    )

# Create if/else branch using create_branch (more natural API)
scf3, scf4 = create_branch(
    parent=scf2,
    branches=[
        (scf3_base, lambda sim: not energy_below_threshold(sim), 'charge_density'),
        (scf4_base, energy_below_threshold, 'charge_density'),
    ],
    error_handlers=[RetryErrorHandler(max_retries=2)],
)

run_project()

