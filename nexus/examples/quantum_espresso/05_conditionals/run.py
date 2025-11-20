#! /usr/bin/env python3
"""
Using conditional utilities from conditionals.py

This example demonstrates the conditional utility functions provided in
nexus/lib/conditionals.py:

machine_conditional - Run simulation only on specific machines

The workflow:
- scf1: Initial SCF calculation (no dependencies)
- scf2: SCF that only runs on 'baseline' machine (machine_conditional)
- scf3: SCF that checks scf2's energy result (custom conditional function)
"""

from nexus import settings, Job, run_project
from nexus import generate_physical_system
from nexus import generate_pwscf

from enhanced_simulation import make_enhanced
from conditionals import (
    machine_conditional,
)

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

# Create reusable job definition for baseline
qe_job = Job(nodes=1,
             queue='batch',
             hours=1,
             presub=qe_modules,
             app=qe_bin+'/pw.x')

# Simulation 1: Initial SCF (no dependencies)
scf1 = generate_pwscf(
    identifier   = 'scf1',
    path         = 'diamond/scf1',
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

# Simulation 2: Only runs on 'baseline' machine (machine_conditional)
scf2_base = generate_pwscf(
    identifier   = 'scf2',
    path         = 'diamond/scf2',
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

# Use machine_conditional to only run on baseline
scf2 = make_enhanced(
    scf2_base,
    condition=machine_conditional('baseline', required=True),
    )

# Simulation 3: Checks scf2's energy result
scf3_base = generate_pwscf(
    identifier   = 'scf3',
    path         = 'diamond/scf3',
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
    dependencies = (scf2, 'charge_density'),
    )

# Use a custom conditional to check if scf2's energy is below threshold
def scf2_energy_below_threshold(sim):
    """Custom conditional that checks scf2's energy via analyzer."""
    try:
        scf2_sim = None
        for dep in sim.dependencies.values():
            if dep.sim.identifier == 'scf2':
                scf2_sim = dep.sim
                break
        if scf2_sim is None or not scf2_sim.finished:
            return False
        analyzer = scf2_sim.load_analyzer_image()
        if hasattr(analyzer, 'E') and analyzer.E != 0.0:
            return analyzer.E < -30.0
        return False
    except Exception:
        return False

scf3 = make_enhanced(
    scf3_base,
    condition=scf2_energy_below_threshold,
    )

run_project()