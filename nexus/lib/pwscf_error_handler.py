"""
Pwscf-specific error handler for Quantum ESPRESSO pw.x simulations.

Provides error strategies and modification strategies specifically designed
for Pwscf simulations, including automatic input parameter modification.
"""

from typing import Dict, Any, Optional, List
from dataclasses import dataclass
import time

from error_handler import (
    SimulationError,
    SimulationModifier,
    CombinedSimulationModifier,
    read_simulation_files,
)


@dataclass
class PwscfError(SimulationError):
    """
    Pwscf-specific error class that inherits from SimulationError.
    
    This class could be extended in the future with Pwscf-specific
    error handling capabilities.
    - QE-specific error metadata
    - Pwscf input parameter extraction
    - Advanced error diagnosis based on QE output structure
    - Integration with Pwscf-specific recovery strategies
    
    """
    

# Pre-defined QE error strategies factory
def _create_default_qe_strategies() -> Dict[str, PwscfError]:
    """Create default QE error strategies using PwscfError."""
    return {
        'convergence_failure': PwscfError(
            pattern='convergence NOT achieved',
            error_type='convergence_failure',
            recoverable=True,
            should_retry=True,
        ),
        'time_limit': PwscfError(
            pattern='Maximum CPU time exceeded',
            error_type='time_limit',
            recoverable=True,
            should_retry=True,
        ),
        'user_stop': PwscfError(
            pattern='Program stopped by user request',
            error_type='user_stop',
            recoverable=False,
            should_retry=False,
        ),
        'error_in_routine': PwscfError(
            pattern='Error in routine',
            error_type='error_in_routine',
            recoverable=False,
            should_retry=False,
        )
    }


QE_ERROR_STRATEGIES = _create_default_qe_strategies()

# Pwscf-specific input modification strategies
@dataclass
class RelaxConvergenceThreshold:
    """Modification strategy: relax convergence threshold."""
    
    factor: float = 10.0
    
    def modify(self, sim: Any) -> None:
        """Relax conv_thr by multiplying by factor."""
        if not hasattr(sim, 'input') or not hasattr(sim.input, 'electrons'):
            return
        electrons = sim.input.electrons
        if hasattr(electrons, 'conv_thr'):
            current = getattr(electrons, 'conv_thr', 1e-6)
            new = current * self.factor
            electrons.conv_thr = new
            if hasattr(sim, 'log'):
                sim.log(f'Relaxed conv_thr: {current:.2e} -> {new:.2e}', n=2)


@dataclass
class IncreaseMaxSteps:
    """Modification strategy: increase maximum SCF steps."""
    
    increase: int = 200
    
    def modify(self, sim: Any) -> None:
        """Increase electron_maxstep by increase amount."""
        if not hasattr(sim, 'input') or not hasattr(sim.input, 'electrons'):
            return
        electrons = sim.input.electrons
        if hasattr(electrons, 'electron_maxstep'):
            current = getattr(electrons, 'electron_maxstep', 100)
            new = current + self.increase
            electrons.electron_maxstep = new
            if hasattr(sim, 'log'):
                sim.log(f'Increased electron_maxstep: {current} -> {new}', n=2)


@dataclass
class AdjustMixingBeta:
    """Modification strategy: adjust mixing beta for stability."""
    
    reduction: float = 0.1
    min_beta: float = 0.1
    
    def modify(self, sim: Any) -> None:
        """Reduce mixing_beta by reduction amount (with minimum)."""
        if not hasattr(sim, 'input') or not hasattr(sim.input, 'electrons'):
            return
        electrons = sim.input.electrons
        if hasattr(electrons, 'mixing_beta'):
            current = getattr(electrons, 'mixing_beta', 0.7)
            new = max(self.min_beta, current - self.reduction)
            electrons.mixing_beta = new
            if hasattr(sim, 'log'):
                sim.log(f'Reduced mixing_beta: {current:.3f} -> {new:.3f}', n=2)


@dataclass
class UseConjugateGradient:
    """Modification strategy: use conjugate gradient for convergence."""
    
    def modify(self, sim: Any) -> None:
        """Use conjugate gradient for convergence."""
        if not hasattr(sim, 'input') or not hasattr(sim.input, 'electrons'):
            return
        electrons = sim.input.electrons
        if hasattr(electrons, 'diagonalization'):
            electrons.diagonalization = 'cg'
            if hasattr(sim, 'log'):
                sim.log('Using conjugate gradient for convergence', n=2)



@dataclass
class PwscfErrorHandler:
    """
    Error handler for Pwscf that modifies input parameters based on error type.
    
    Uses error strategies to diagnose and handle errors. Modification strategies
    can be provided to customize how inputs are modified.
    
    Args:
        max_retries: Maximum number of retry attempts
        delay: Delay in seconds before retry
        error_strategies: List of PwscfError objects (defaults to QE_ERROR_STRATEGIES)
        modification_strategies: Dict mapping error_type to SimulationModifier
    """
    
    max_retries: int = 3
    delay: float = 0.0
    error_strategies: Optional[Dict[str, PwscfError]] = None
    modification_strategies: Optional[Dict[str, SimulationModifier]] = None
    
    def __post_init__(self):
        """Initialize with default QE error strategies if not provided."""
        if self.error_strategies is None:
            self.error_strategies = _create_default_qe_strategies()
        
        # Apply default modification strategies for convergence failures
        if self.modification_strategies is None:
            self.modification_strategies = {}
        
        # Default convergence failure modifications
        if 'convergence_failure' not in self.modification_strategies:
            self.modification_strategies['convergence_failure'] = CombinedSimulationModifier([
                RelaxConvergenceThreshold(factor=10.0),
                IncreaseMaxSteps(increase=200),
                AdjustMixingBeta(reduction=0.1),
            ])
    
    def diagnose_error(self, sim: Any) -> Optional[PwscfError]:
        """
        Diagnose error by checking all error strategies.
        
        Returns:
            Matching PwscfError or None if no match
        """
        # Read both output and error files using shared utility
        output, error = read_simulation_files(sim)
        
        # If no files found, can't diagnose
        if not output and not error:
            return None
        
        # Check strategies in order (job_done first, then others)
        for strategy in self.error_strategies.values():
            if strategy.matches(output or '', error or ''):
                return strategy
        
        return None
    
    def handle_error(self, sim: Any, error_info: Dict[str, Any]) -> bool:
        """
        Handle QE errors by modifying input and deciding whether to retry.
        
        Returns:
            True if should retry, False otherwise
        """
        attempt_count = getattr(sim, 'attempt_number', 0)
        
        # Check retry limit
        if attempt_count >= self.max_retries:
            return False
        
        # Diagnose error
        error_strategy = self.diagnose_error(sim)
        if error_strategy is None:
            return False  # Unknown error, don't retry
        
        # If job is done, no need to retry
        if error_strategy.error_type == 'completed':
            return False
        
        # Apply modification strategy if available
        if error_strategy.error_type in self.modification_strategies:
            mod_strategy = self.modification_strategies[error_strategy.error_type]
            if hasattr(mod_strategy, 'modify'):
                mod_strategy.modify(sim)
        
        # Also check if strategy has its own modification
        if error_strategy.modification_strategy is not None:
            error_strategy.modification_strategy.modify(sim)
        
        # Handle error using strategy
        should_retry = error_strategy.handle(sim, error_info)
        
        if should_retry and self.delay > 0:
            time.sleep(self.delay)
        
        return should_retry
    
    def is_compatible(self, sim: Any) -> bool:
        """Check if this handler is compatible with the simulation."""
        if hasattr(sim, '__class__'):
            class_name = sim.__class__.__name__
            if class_name == 'Pwscf' or (class_name.startswith('Enhanced') and 'Pwscf' in class_name):
                return True
            module_name = getattr(sim.__class__, '__module__', '')
            if 'pwscf' in module_name.lower():
                return True
        return False

