"""
Error handler infrastructure for DAG simulations.
"""

from typing import Protocol, Callable, Dict, Any, Optional, List
from dataclasses import dataclass
import time
import os

class ErrorHandler(Protocol):
    """Protocol for error handlers."""
    
    def diagnose_error(self, sim: Any) -> Optional[Dict[str, Any]]:
        """
        Diagnose an error for a simulation.
        
        Args:
            sim: Simulation instance that failed
        """
        ...
    
    def handle_error(self, sim: Any, error_info: Dict[str, Any]) -> bool:
        """
        Handle an error for a simulation.
        
        Args:
            sim: Simulation instance that failed
            error_info: Dictionary with error information (error_type, message, etc.)
            
        Returns:
            True if simulation should be resubmitted, False otherwise
        """
        ...
    
    def is_compatible(self, sim: Any) -> bool:
        """
        Check if this error handler is compatible with the given simulation.
        
        Args:
            sim: Simulation instance to check compatibility with
            
        Returns:
            True if compatible, False otherwise
        """
        # Default implementation: compatible with all simulations
        # Override in subclasses for simulation-specific handlers
        return True


@dataclass
class RetryErrorHandler:
    """Error handler that retries a simulation up to max_retries times."""
    
    max_retries: int = 3
    delay: float = 0.0  # Delay in seconds before retry
    
    def handle_error(self, sim: Any, error_info: Dict[str, Any]) -> bool:
        """Handle error by retrying if under max_retries."""
        attempt_count = getattr(sim, 'attempt_number', 0)
        if attempt_count < self.max_retries:
            if self.delay > 0:
                time.sleep(self.delay)
            return True  # Resubmit
        return False  # Don't resubmit
    
    def is_compatible(self, sim: Any) -> bool:
        """Generic retry handler is compatible with all simulations."""
        return True


@dataclass
class ModifyInputErrorHandler:
    """Error handler that modifies simulation input before resubmission."""
    
    modification_func: Callable[[Any, Dict[str, Any]], None]
    
    def handle_error(self, sim: Any, error_info: Dict[str, Any]) -> bool:
        """Handle error by modifying input and resubmitting."""
        self.modification_func(sim, error_info)
        return True  # Resubmit after modification
    
    def is_compatible(self, sim: Any) -> bool:
        """Generic input modifier is compatible with all simulations."""
        return True


@dataclass
class SkipDependentErrorHandler:
    """Error handler that skips dependent simulations when this one fails."""
    
    def handle_error(self, sim: Any, error_info: Dict[str, Any]) -> bool:
        """Handle error by blocking dependents."""
        if hasattr(sim, 'block_dependents'):
            sim.block_dependents(block_self=False)
        return False  # Don't resubmit
    
    def is_compatible(self, sim: Any) -> bool:
        """Generic handler is compatible with all simulations."""
        return True


@dataclass
class ConditionalErrorHandler:
    """Error handler that conditionally handles errors based on a function."""
    
    condition_func: Callable[[Any, Dict[str, Any]], bool]
    true_handler: ErrorHandler
    false_handler: Optional[ErrorHandler] = None
    
    def handle_error(self, sim: Any, error_info: Dict[str, Any]) -> bool:
        """Handle error conditionally."""
        if self.condition_func(sim, error_info):
            return self.true_handler.handle_error(sim, error_info)
        elif self.false_handler is not None:
            return self.false_handler.handle_error(sim, error_info)
        return False
    
    def is_compatible(self, sim: Any) -> bool:
        """Compatible if both handlers are compatible."""
        compatible = self.true_handler.is_compatible(sim) if hasattr(self.true_handler, 'is_compatible') else True
        if self.false_handler is not None and hasattr(self.false_handler, 'is_compatible'):
            compatible = compatible and self.false_handler.is_compatible(sim)
        return compatible


def read_file_safe(filepath: str, max_lines: Optional[int] = None) -> Optional[str]:
    """
    Safely read a file, returning None if it doesn't exist or can't be read.
    
    Args:
        filepath: Path to the file to read
        max_lines: If provided, only return the last max_lines lines
        
    Returns:
        File contents as string, or None if file doesn't exist or can't be read
    """
    if not os.path.exists(filepath):
        return None
    try:
        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            if max_lines is not None:
                lines = f.readlines()
                return ''.join(lines[-max_lines:])
            else:
                return f.read()
    except Exception:
        return None


def read_simulation_files(sim: Any, max_lines: Optional[int] = None) -> tuple[Optional[str], Optional[str]]:
    """
    Read output and error files from a simulation.
    
    Args:
        sim: Simulation instance
        max_lines: If provided, only return the last max_lines lines from each file
        
    Returns:
        Tuple of (output_content, error_content), either may be None
    """
    output_content = None
    error_content = None
    
    if hasattr(sim, 'outfile') and hasattr(sim, 'locdir'):
        outfile_path = os.path.join(sim.locdir, sim.outfile)
        output_content = read_file_safe(outfile_path, max_lines)
    
    if hasattr(sim, 'errfile') and hasattr(sim, 'locdir'):
        errfile_path = os.path.join(sim.locdir, sim.errfile)
        error_content = read_file_safe(errfile_path, max_lines)
    
    return output_content, error_content


class SimulationModifier(Protocol):
    """Protocol for input modification strategies."""
    
    def modify(self, sim: Any) -> None:
        """
        Modify simulation input to address an error.
        
        Args:
            sim: Simulation instance to modify
        """
        ...

@dataclass
class CombinedSimulationModifier:
    """Combine multiple modification strategies."""
    
    strategies: List[SimulationModifier]
    
    def modify(self, sim: Any) -> None:
        """Apply all modification strategies in order."""
        for strategy in self.strategies:
            if hasattr(strategy, 'modify'):
                strategy.modify(sim)


@dataclass
class SimulationError:
    """
    Base class for error strategies that combine diagnosis and handling.
    
    Each error strategy encapsulates:
    - Pattern to detect in output/error files
    - Whether it's recoverable
    - How to handle it (modify input, retry, etc.)
    """
    
    pattern: str  # Pattern to search for in output/error files
    error_type: str  # Name of error type
    recoverable: bool = True
    should_retry: bool = True
    modification_strategy: Optional[SimulationModifier] = None
    
    def matches(self, output: str, error: str) -> bool:
        """Check if this error strategy matches the output."""
        in_output = self.pattern in output if output else False
        in_error = self.pattern in error if error else False
        # Also check additional files if specified
        return in_output or in_error
    
    def handle(self, sim: Any, error_info: Dict[str, Any]) -> bool:
        """
        Handle the error by modifying input if needed and deciding whether to retry.
        
        Returns:
            True if should retry, False otherwise
        """
        # Apply modification strategy if provided
        if self.modification_strategy is not None:
            try:
                self.modification_strategy.modify(sim)
            except Exception as e:
                if hasattr(sim, 'warn'):
                    sim.warn(f'Error in modification strategy: {e}')
        
        return self.should_retry


