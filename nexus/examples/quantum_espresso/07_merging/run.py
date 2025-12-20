#! /usr/bin/env python3
"""
Branching and Merging Workflows

This example demonstrates:
1. Natural branching using create_branch() for if/else conditionals
2. Merging branches with different strategies (first, any, all, custom selector)

Branching:
- scf1 → scf2 (if condition A) OR scf3 (if condition B)
- Both scf2 and scf3 depend on scf1

Merging:
- scf4 depends on scf2 and scf3 with different merge strategies
- Demonstrates 'first', 'any', 'all', and custom selector options
"""

from nexus import settings, Job, run_project
from nexus import generate_physical_system
from nexus import generate_pwscf
import random

from enhanced_simulation import make_enhanced, create_branch

# Computer configuration
computer = 'ws16'

if computer == 'baseline':
    qe_modules = 'module purge; module load Core/25.05   gcc/12.4.0   openmpi/5.0.5   DefApps hdf5'
    qe_bin = '/ccsopen/home/ksu/SOFTWARE/qe/q-e-qe-7.4.1/build/bin'
    account = 'phy191'
elif computer == 'ws16':
    qe_modules = ''
    qe_bin = '/Users/ksu/Software/qe/q-e-qe-7.5/build/bin'
    account = 'ks'
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

# Parent simulation
scf1 = generate_pwscf(
    identifier   = 'scf1',
    path         = 'diamond/branch_merge/scf1',
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

# Create branches using create_branch (natural API)
scf2 = generate_pwscf(
    identifier   = 'scf2',
    path         = 'diamond/branch_merge/scf2',
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

scf3 = generate_pwscf(
    identifier   = 'scf3',
    path         = 'diamond/branch_merge/scf3',
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

# Now demonstrate merging strategies
# scf4 will merge scf2 and scf3 with different strategies

# Example 1: Merge with 'first' strategy (race condition)
scf4_first = generate_pwscf(
    identifier   = 'scf4_first',
    path         = 'diamond/branch_merge/scf4_first',
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

scf4_first_enhanced = make_enhanced(scf4_first)
# Run when FIRST of scf2 or scf3 completes
scf4_first_enhanced.depends((scf2, 'charge_density'), (scf3, 'charge_density'), strategy='first')


# Example 3: Merge with 'all' strategy (default, AND logic)
scf4_all = generate_pwscf(
    identifier   = 'scf4_all',
    path         = 'diamond/branch_merge/scf4_all',
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

scf4_all_enhanced = make_enhanced(scf4_all)
# Run when ALL of scf2 and scf3 complete (default behavior)
scf4_all_enhanced.depends((scf2, 'other'), (scf3, 'other'), strategy='all')

# Example 4: Merge with custom selector function
def select_lowest_energy(sims):
    """
    Custom selector: choose simulation with lowest energy.
    
    Args:
        sims: List of completed Simulation objects
        
    Returns:
        Simulation with lowest energy
    """
    if not sims:
        return None
    
    # Load analyzer for each sim and compare energies
    best_sim = None
    best_energy = float('inf')
    
    for sim in sims:
        try:
            if sim.finished:
                analyzer = sim.load_analyzer_image()
                if hasattr(analyzer, 'E') and analyzer.E != 0.0:
                    energy = analyzer.E
                    if energy < best_energy:
                        best_energy = energy
                        best_sim = sim
        except Exception:
            continue
    
    # Return best sim, or first sim if no energies found
    return best_sim if best_sim else sims[0]

scf4_custom = generate_pwscf(
    identifier   = 'scf4_custom',
    path         = 'diamond/branch_merge/scf4_custom',
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

scf4_custom_enhanced = make_enhanced(scf4_custom)
# Run when ALL complete, then select one with lowest energy
scf4_custom_enhanced.depends((scf2, 'charge_density'), (scf3, 'charge_density'), strategy='all', selector=select_lowest_energy)

# Run the workflow
run_project()
