#! /usr/bin/env python3
"""
Error Handler Demonstration with PwscfErrorHandler

This example demonstrates the error handler system, specifically error handling for Quantum ESPRESSO pw.x simulations.

Features demonstrated:
1. PwscfErrorHandler - QE-specific handler that modifies input on errors
2. RetryErrorHandler - Simple retry handler (for comparison)
3. Error strategies - Compact objects that combine diagnosis and handling
4. Modification strategies - Reusable strategies for input modification
5. Compatibility checking - Handlers check if they're compatible with simulations

The workflow:
- scf1: Default PwscfErrorHandler (auto-modifies input on convergence failures)
- scf2: RetryErrorHandler (simple retry, no diagnosis)
- scf3: Multiple handlers (compatibility checked automatically)
- scf4: Custom modification strategies (user-defined input modifications)
"""

from nexus import settings, Job, run_project
from nexus import generate_physical_system
from nexus import generate_pwscf

import os
import sys

nexus_lib = os.path.join(os.path.dirname(__file__), '../../lib')
if nexus_lib not in sys.path:
    sys.path.insert(0, nexus_lib)

from enhanced_simulation import make_enhanced
from error_handler import (
    RetryErrorHandler,
)
from pwscf_error_handler import (
    PwscfErrorHandler,
    RelaxConvergenceThreshold,
    IncreaseMaxSteps,
    AdjustMixingBeta,
    UseConjugateGradient,
    CombinedSimulationModifier,
)

# Computer configuration
computer = 'baseline'

if computer == 'baseline':
    qe_modules = 'module purge; module load Core/25.05   gcc/12.4.0   openmpi/5.0.5   DefApps hdf5'
    qe_bin = '/ccsopen/home/ksu/SOFTWARE/qe/q-e-qe-7.4.1/build/bin'
    account = 'phy191'
elif computer == 'ws4':
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

qe_job = Job(nodes=1,
             queue='batch',
             hours=0.25,
             presub=qe_modules,
             app=qe_bin+'/pw.x')

# Simulation 1: PwscfErrorHandler (default configuration)
# This handler uses error strategies to diagnose QE errors and automatically
# modifies input parameters to help the simulation succeed on retry.
# Default behavior: On convergence failures, it will:
# - Relax conv_thr by factor of 10
# - Increase electron_maxstep by 200
# - Reduce mixing_beta by 0.1
scf1_base = generate_pwscf(
    identifier   = 'scf1',
    path         = 'diamond/error_handlers/scf1',
    job          = qe_job,
    input_type   = 'generic',
    calculation  = 'scf',
    input_dft    = 'lda', 
    ecutwfc      = 200,   
    conv_thr     = 1e-8,
    electron_maxstep = 100,
    system       = system,
    pseudos      = ['C.BFD.upf'],
    kgrid        = (4,4,4),
    kshift       = (0,0,0),
    )

# Use PwscfErrorHandler with default modification strategies
# If convergence fails, input will be automatically modified before retry
scf1 = make_enhanced(
    scf1_base,
    error_handlers=[
        PwscfErrorHandler(
            max_retries=3,
            delay=1.0,  # Wait 1 second before retry
        )
    ],
)

# Simulation 2: RetryErrorHandler (simple)
# Simple retry handler - no error diagnosis, just retries up to max_retries
scf2_base = generate_pwscf(
    identifier   = 'scf2',
    path         = 'diamond/error_handlers/scf2',
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

scf2 = make_enhanced(
    scf2_base,
    error_handlers=[
        RetryErrorHandler(max_retries=2, delay=1.0)
    ],
)

# Simulation 3: Multiple handlers (compatibility checked)
# Demonstrates that handlers check compatibility automatically and are called in order.
# The FIRST handler that returns True will trigger resubmission and stop the chain.
scf3_base = generate_pwscf(
    identifier   = 'scf3',
    path         = 'diamond/error_handlers/scf3',
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

# Multiple handlers - compatibility is checked automatically
scf3 = make_enhanced(
    scf3_base,
    error_handlers=[
        PwscfErrorHandler(max_retries=2),  # Tried first: QE-specific with input modification
        RetryErrorHandler(max_retries=1),  # Last resort: simple retry
    ],
)

# Simulation 4: Custom modification strategies
# Demonstrates how to create custom modification strategies for specific error types
scf4_base = generate_pwscf(
    identifier   = 'scf4',
    path         = 'diamond/error_handlers/scf4',
    job          = qe_job,
    input_type   = 'generic',
    calculation  = 'scf',
    input_dft    = 'lda', 
    ecutwfc      = 200,   
    conv_thr     = 1e-8,  # Tight convergence - might fail
    electron_maxstep = 100,  # Limited steps - might need more
    system       = system,
    pseudos      = ['C.BFD.upf'],
    kgrid        = (4,4,4),
    kshift       = (0,0,0),
    dependencies = (scf3, 'charge_density'),
    )

# Use PwscfErrorHandler with custom modification strategies
# This demonstrates the compact strategy-based API:
# - Error strategies define what errors to detect
# - Modification strategies define how to fix them
# - Strategies can be combined and customized
scf4 = make_enhanced(
    scf4_base,
    error_handlers=[
        PwscfErrorHandler(
            max_retries=3,
            delay=1.0,
            modification_strategies={
                'convergence_failure': CombinedSimulationModifier([
                    RelaxConvergenceThreshold(factor=10.0),  # conv_thr *= 10
                    IncreaseMaxSteps(increase=200),          # electron_maxstep += 200
                    AdjustMixingBeta(reduction=0.1),        # mixing_beta -= 0.1 (min 0.1)
                    UseConjugateGradient(),                 # Use CG diagonalization
                ]),
            }
        )
    ],
)

run_project()
